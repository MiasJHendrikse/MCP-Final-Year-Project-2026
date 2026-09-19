"""
`ScaledProblem`: the optimisation problem as SLSQP sees it (Phase 2, Stage A2).

What an optimiser is handed
----------------------------
    variable    u in [0, 1]^10,  d = bounds.to_physical(u)
    objective   fun(u) = J(u) / |J(u0)|,  J = -AEP [MWh/yr]
    constraint  the polar-cache Reynolds envelope, linear in u (below)
    gradient    jac_fd(u, h): central differences of `fun`
                jac_adjoint(u): the discrete adjoint (Phase 3), same units

The objective is normalised by |J| at the first point evaluated (recorded as
`J0`) because SLSQP's `ftol` is an *absolute* tolerance on the objective:
`ftol = 1e-8` means something for a quantity of order one and nothing for one
of order ten. `J0` is stored so every number can be undone back to MWh/yr.

The envelope constraint is mandatory, not optional
---------------------------------------------------
`CachedPolar` raises `PolarDomainError` outside the cache's Reynolds range
(`interpolant_for(polar_cache).re_values[0]` .. `[-1]`, 40 k .. 1 M for
SG6043). The box bounds let the optimiser leave that range: `chord_max_m =
0.30 m` at the tip gives Re ~ 1.04 M at 19.5 m/s under the 300 rpm ceiling
(and ~2 M with none). The solver's own Reynolds estimate
(`bem.rotor.solve_rotor`) is

    Re_{b,i} = W_i(V_b) c_i / nu,   W_i(V_b) = hypot(V_b, Omega_b r_i)

with `Omega_b = lambda_b V_b / R` the operating point's rotor speed on the
configured schedule (`objective.power.operating_points`). It depends on the
chord alone -- no state -- and `c = N_c d_c` is linear in the design vector,
so "stay inside the cache at every operating point" is a *linear*
constraint on `u` with a constant Jacobian:

    c_min,i (1 + mu) <= [N_c (lo + u * span)]_i <= c_max,i (1 - mu)

    c_max,i = Re_hi nu / max_b W_i(V_b),   c_min,i = Re_lo nu / min_b W_i(V_b)

The extremes are taken over the actual operating points rather than assumed
to sit at `V_min` and `V_max` (2026-09-19; until then every bin ran lambda =
6.5 and the two coincided). Under a ceiling `W_i` at the top bin is
`hypot(19.5, Omega_max r_i)`, well below the fixed-lambda value, so the
ceiling rows are looser than the conservative form the experiment reused;
`envelope_data()` reports which operating point set each row.

`mu = 0.05` is a stated margin, not a tuned one: SLSQP's intermediate
iterates may violate inequality constraints slightly, and an objective
evaluation that raises mid-line-search is an aborted run. If a
`PolarDomainError` still occurs the offending `u` is logged on the instance
and the error re-raised -- never a substitute value.

The angle-of-attack range of the cache (`xfoil_alpha_min/max`, -8..18 deg
for SG6043) is *state*-dependent (it needs the solved `phi`), so it is not a
constraint; `alpha_check` is a logged post-check on accepted iterates.

Bounds and the constraint set
------------------------------
`bounds` is a required argument with no default (`DesignBounds.from_config()`
since 2026-09-19; the 0.45 m placeholder that preceded it is retired).
`constraints()` is the full inequality set every optimisation run hands to
SLSQP -- the envelope and the configured solidity cap -- so no script can
assemble a different one; `active_set(u)` reports which rows and bounds are
active at a point.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

import numpy as np

from config import load_design_rotor, load_site
from gradients.finite_difference import central_difference
from objective.objective import objective
from objective.power import aerodynamic_power, operating_points
from polars.interpolant import PolarDomainError
from polars.polar import interpolant_for

#: Envelope margin, stated (see module docstring). Raised to 0.10 once, and
#: recorded, if a run still trips `PolarDomainError` -- never tuned beyond.
DEFAULT_ENVELOPE_MARGIN = 0.05

#: `|u|` or `|1 - u|` below this: the bound is active. `g(u)` below this: the
#: constraint row is active. One tolerance for every reporting script.
ACTIVE_TOL = 1e-6


class ScaledProblem:
    """
    The scaled, constrained optimisation problem over `u in [0, 1]^n`.

    Parameters
    ----------
    parameterisation : design.BladeParameterisation
    bounds : design.DesignBounds
        Required. The provisional bounds are constructed by the caller from
        `tests/test_parameterisation.py::PROVISIONAL_BOUNDS`.
    resource : objective.WeibullResource
        Required; `WeibullResource.from_config()` for the site.
    polar_cache : str
        The blade's polar cache, default `"sg6043"` -- the one default, and
        the same one `BladeParameterisation.to_geometry` carries.
    margin : float
        The envelope margin `mu`.

    Attributes
    ----------
    J0 : float or None
        `J` at the first point `fun` was evaluated at (MWh/yr, negative).
        `None` until then.
    n_fun_evals : int
        Objective evaluations made through this instance.
    domain_errors : list of list of float
        Every `u` at which the objective raised `PolarDomainError`.
    """

    def __init__(self, parameterisation, bounds, resource, polar_cache="sg6043",
                 margin=DEFAULT_ENVELOPE_MARGIN):
        if bounds.n_design_variables != parameterisation.n_design_variables:
            raise ValueError(
                f"bounds carry {bounds.n_design_variables} variables, the "
                f"parameterisation {parameterisation.n_design_variables}")

        self.parameterisation = parameterisation
        self.bounds = bounds
        self.resource = resource
        self.polar_cache = polar_cache
        self.margin = float(margin)

        self.design = load_design_rotor()
        self.site = load_site()

        interpolant = interpolant_for(polar_cache)
        self.reynolds_lo = float(interpolant.re_values[0])
        self.reynolds_hi = float(interpolant.re_values[-1])
        self.alpha_min_deg = float(interpolant.xfoil_alpha_min)
        self.alpha_max_deg = float(interpolant.xfoil_alpha_max)

        self.J0 = None
        self.n_fun_evals = 0
        self.n_adjoint_evals = 0
        self.domain_errors = []
        self._adjoint = None

    # -- variables ----------------------------------------------------------

    @property
    def n(self):
        return self.parameterisation.n_design_variables

    def physical(self, u):
        """`d = lo + u * span`."""

        return self.bounds.to_physical(u)

    def scaled(self, d):
        """`u = (d - lo) / span`."""

        return self.bounds.to_scaled(d)

    def operating_points(self):
        """
        Every operating point the objective solves at, `[(V, lambda)]`: the
        17 bin midpoints on the configured schedule. The rating is a
        configured number, so the rated speed is not one of them.
        """

        return operating_points(self.design)

    def operating_speeds(self):
        """The wind speeds of `operating_points()`."""

        return [v for v, _lam in self.operating_points()]

    # -- objective ----------------------------------------------------------

    def J(self, u):
        """`J(u) = -AEP(d(u))` in MWh/yr, unscaled."""

        d = self.physical(u)
        try:
            value = objective(d, resource=self.resource,
                              parameterisation=self.parameterisation,
                              polar_cache=self.polar_cache)
        except PolarDomainError:
            self.domain_errors.append([float(x) for x in np.asarray(u)])
            raise
        self.n_fun_evals += 1
        return float(value)

    def aep_mwh(self, u):
        """`-J(u)`: AEP in MWh/yr, for reporting."""

        return -self.J(u)

    def set_reference(self, u0):
        """Fix `J0 = J(u0)` explicitly (otherwise the first `fun` call does it)."""

        self.J0 = self.J(u0)
        return self.J0

    def fun(self, u):
        """
        `J(u) / |J0|` -- the objective SLSQP minimises.

        `J0` is recorded at the first call, so `fun(u0) == -1.0` exactly at
        the starting point: `x / |x|` for negative `x` is `-1.0` to the bit.
        """

        value = self.J(u)
        if self.J0 is None:
            self.J0 = value
        return value / abs(self.J0)

    def unscale(self, scaled_value):
        """Objective (or gradient) in `fun` units back to MWh/yr."""

        if self.J0 is None:
            raise RuntimeError("J0 not set: evaluate fun(u0) or set_reference(u0) first")
        return np.asarray(scaled_value) * abs(self.J0)

    def jac_fd(self, u, h):
        """
        Central-difference gradient of `fun` at `u`; `h` scalar or (n,).

        This is what `scipy.optimize.minimize(..., jac=...)` is given in the
        FD-driven run; the count of evaluations lands in `n_fun_evals`.
        """

        grad, _n_evals = central_difference(self.fun, u, h)
        return grad

    def adjoint_system(self):
        """
        The `adjoint.BEMSystem` over this problem's parameterisation, bounds,
        resource and polar cache, built on first use. Imported lazily so
        `gradients/` does not depend on `adjoint/` at import time (the
        dependency runs the other way for the FD path).
        """

        if self._adjoint is None:
            from adjoint.system import BEMSystem

            self._adjoint = BEMSystem(self.parameterisation, self.bounds,
                                      self.resource, self.polar_cache)
        return self._adjoint

    def jac_adjoint(self, u):
        """
        Discrete-adjoint gradient of `fun` at `u`:
        `gradient(d) * span / |J0|` (Phase 3, B3).

        One forward solve plus the partials -- the cost of about two
        objective evaluations, against `2n` for `jac_fd`. `J0` must be set
        (by `fun(u0)` or `set_reference`), as for `unscale`. A
        `PolarDomainError` is logged and re-raised exactly as in `J`.
        """

        if self.J0 is None:
            raise RuntimeError("J0 not set: evaluate fun(u0) or set_reference(u0) first")
        d = self.physical(u)
        try:
            result = self.adjoint_system().gradient(d)
        except PolarDomainError:
            self.domain_errors.append([float(x) for x in np.asarray(u)])
            raise
        self.n_adjoint_evals += 1
        return result.dJ_dd * self.bounds.span() / abs(self.J0)

    # -- the polar-cache envelope (linear, mandatory) ------------------------

    def _chord_affine(self):
        """`c(u) = A u + b` at every station, from `N_c` and the scaling."""

        n_c = self.parameterisation.dchord_dd()
        return n_c * self.bounds.span()[None, :], n_c @ self.bounds.lower()

    def envelope_data(self):
        """
        The rows of the envelope constraint, for reporting.

        Returns
        -------
        dict
            `radii`, `chord_min_m`, `chord_max_m` (the raw cache limits per
            station, before the margin), `margin`, `reynolds_lo/hi`, and per
            station the operating point (`v_at_min_w`, `v_at_max_w`, with
            `tsr_at_*`) whose `W` set each limit.
        """

        nu = float(self.site.kinematic_viscosity)
        R = self.parameterisation.radius_m
        radii = self.parameterisation.radii

        points = self.operating_points()
        # W[b, i] = hypot(V_b, Omega_b r_i), every operating point x station.
        w = np.array([np.hypot(v, lam * v * radii / R) for v, lam in points])
        b_min = np.argmin(w, axis=0)
        b_max = np.argmax(w, axis=0)
        stations = np.arange(len(radii))

        return {
            "radii": radii,
            "chord_min_m": self.reynolds_lo * nu / w[b_min, stations],
            "chord_max_m": self.reynolds_hi * nu / w[b_max, stations],
            "margin": self.margin,
            "reynolds_lo": self.reynolds_lo,
            "reynolds_hi": self.reynolds_hi,
            "v_at_min_w": np.array([points[b][0] for b in b_min]),
            "tsr_at_min_w": np.array([points[b][1] for b in b_min]),
            "v_at_max_w": np.array([points[b][0] for b in b_max]),
            "tsr_at_max_w": np.array([points[b][1] for b in b_max]),
        }

    def envelope_constraint(self):
        """
        SciPy inequality constraint `g(u) >= 0` keeping every station's
        Reynolds number inside the polar cache at every operating point.

        `g` has `2 * n_stations` rows: the first block is the floor
        (`c_i - c_min,i (1 + mu)`), the second the ceiling
        (`c_max,i (1 - mu) - c_i`). Constant Jacobian.
        """

        A, b = self._chord_affine()
        data = self.envelope_data()
        floor = data["chord_min_m"] * (1.0 + self.margin)
        ceiling = data["chord_max_m"] * (1.0 - self.margin)

        jac = np.vstack([A, -A])
        offset = np.concatenate([b - floor, ceiling - b])

        return {
            "type": "ineq",
            "fun": lambda u: jac @ np.asarray(u, dtype=float) + offset,
            "jac": lambda u: jac,
        }

    def envelope_row_labels(self):
        """`"floor r=..."` / `"ceiling r=..."` per row of `envelope_constraint`."""

        radii = self.parameterisation.radii
        return ([f"floor r={r:.4f}" for r in radii]
                + [f"ceiling r={r:.4f}" for r in radii])

    def solidity_row_labels(self):
        """`"solidity r=..."` per row of `solidity_constraint`."""

        return [f"solidity r={r:.4f}" for r in self.parameterisation.radii]

    def constraints(self):
        """
        The inequality set every optimisation run uses: the polar-cache
        envelope and the configured solidity cap, in that order. One place,
        so no script can quietly run a different problem.
        """

        return [self.envelope_constraint(), self.solidity_constraint()]

    def active_set(self, u, tol=ACTIVE_TOL):
        """
        Which bounds and constraint rows are active at `u`: a dict with
        `bounds_lower`, `bounds_upper` (variable indices), `envelope`,
        `solidity` (row labels), the tightest row of each constraint with
        its slack, and every station's solidity.
        """

        u = np.asarray(u, dtype=float)
        env = self.envelope_constraint()["fun"](u)
        sol = self.solidity_constraint()["fun"](u)
        env_labels = self.envelope_row_labels()
        sol_labels = self.solidity_row_labels()
        chord = self.parameterisation.chord(self.physical(u))
        sigma = self.design.n_blades * chord / (2.0 * math.pi * self.parameterisation.radii)
        return {
            "bounds_lower": [int(j) for j in range(self.n) if abs(u[j]) < tol],
            "bounds_upper": [int(j) for j in range(self.n) if abs(1.0 - u[j]) < tol],
            "envelope": [env_labels[k] for k in range(len(env)) if env[k] < tol],
            "solidity": [sol_labels[k] for k in range(len(sol)) if sol[k] < tol],
            "envelope_tightest": {"row": env_labels[int(np.argmin(env))],
                                  "slack_m": float(env.min())},
            "solidity_tightest": {"row": sol_labels[int(np.argmin(sol))],
                                  "slack": float(sol.min())},
            "sigma_max": float(sigma.max()),
            "sigma_cap": float(self.design.max_local_solidity),
        }

    # -- solidity -------------------------------------------------------------

    def solidity_constraint(self, cap=None):
        """
        SciPy inequality `cap - sigma_i >= 0`, `sigma_i = B c_i / (2 pi r_i)`.

        `cap` defaults to the configured `constraints.max_local_solidity`
        (0.5, resolved 2026-09-19 with `chord_max_m`). Linear, constant
        Jacobian. With the 0.30 m box on the control points and the spline's
        convex-hull property the cap cannot be reached (the first station
        would need 350 mm), so it reports inactive in every run; it is in
        the constraint set because it is the principled radius-aware form of
        the same limit, not because it shapes any result.
        """

        cap = float(self.design.max_local_solidity if cap is None else cap)
        A, b = self._chord_affine()
        factor = self.design.n_blades / (2.0 * math.pi * self.parameterisation.radii)
        jac = -A * factor[:, None]
        offset = cap - b * factor

        return {
            "type": "ineq",
            "fun": lambda u: jac @ np.asarray(u, dtype=float) + offset,
            "jac": lambda u: jac,
        }

    # -- post-checks on accepted iterates -----------------------------------

    def operating_state(self, u):
        """
        Per-operating-point station state at `u`: what the post-checks read.

        Solves the rotor at every point in `operating_points()` (17 solves,
        the cost of one objective evaluation). Returns a dict of lists, one
        entry per point: `speeds`, `tsr`, `alpha_deg` (25,), `reynolds`
        (25,), `phi` (25,), `a` (25,), `power_w`, `converged`.
        """

        d = self.physical(u)
        geometry = self.parameterisation.to_geometry(d, polar_cache=self.polar_cache)
        rho = float(self.site.air_density)
        nu = float(self.site.kinematic_viscosity)

        state = {"speeds": [], "tsr": [], "alpha_deg": [], "reynolds": [], "phi": [],
                 "a": [], "power_w": [], "converged": []}
        for v, tsr in self.operating_points():
            power, result = aerodynamic_power(geometry, v, tsr, rho, nu)
            stations = result["stations"]
            state["speeds"].append(v)
            state["tsr"].append(tsr)
            state["alpha_deg"].append([math.degrees(s["alpha"]) for s in stations])
            state["reynolds"].append([s["reynolds"] for s in stations])
            state["phi"].append([s["phi"] for s in stations])
            state["a"].append([s["a"] for s in stations])
            state["power_w"].append(float(power))
            state["converged"].append(bool(result["converged"]))
        return state

    def alpha_check(self, u):
        """
        Post-check: does every solved angle of attack lie inside the cache's
        XFOIL-converged band `[alpha_min_deg, alpha_max_deg]`?

        A logged check rather than a constraint because it needs the state.
        Returns a dict with the extreme angles, the limits, the margin to the
        nearest limit, the Reynolds extremes (the same pass answers whether
        the envelope ceiling is close, §4.5), and `within`.

        Reported twice: over every operating point, and over the *uncapped*
        points only. Under a rotor-speed ceiling the capped bins run
        lambda ~ 3-5 and their inboard stations sit in the Viterna
        extrapolation (alpha ~ 30 deg at 19.5 m/s); they contribute a
        constant to the objective and nothing to the gradient, so the
        uncapped figure is the one that says whether the *optimised* part of
        the objective is inside the validated polar band. Both are
        reported so neither can be mistaken for the other.
        """

        state = self.operating_state(u)
        alpha = np.array(state["alpha_deg"])
        reynolds = np.array(state["reynolds"])
        uncapped = np.array(state["power_w"]) <= float(self.design.rated_power_w)

        def band(rows):
            lo, hi = float(alpha[rows].min()), float(alpha[rows].max())
            return {
                "alpha_min_deg": lo,
                "alpha_max_deg": hi,
                "alpha_margin_deg": float(min(lo - self.alpha_min_deg, self.alpha_max_deg - hi)),
                "within": bool(lo >= self.alpha_min_deg and hi <= self.alpha_max_deg),
            }

        every = band(np.ones(len(alpha), dtype=bool))
        report = dict(every)
        report.update({
            "alpha_limits_deg": [self.alpha_min_deg, self.alpha_max_deg],
            "uncapped": band(uncapped) if uncapped.any() else None,
            "n_uncapped_points": int(uncapped.sum()),
            "reynolds_min": float(reynolds.min()),
            "reynolds_max": float(reynolds.max()),
            "reynolds_limits": [self.reynolds_lo, self.reynolds_hi],
            "all_converged": all(state["converged"]),
        })
        return report
