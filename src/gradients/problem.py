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
SG6043). The provisional box bounds let the optimiser leave that range:
`chord_max_m = 0.45 m` at the tip gives Re ~ 3 M at 19.5 m/s. The solver's
own Reynolds estimate (`bem.rotor.solve_rotor`) is

    Re_{b,i} = W_i(V_b) c_i / nu,   W_i(V) = hypot(V, lambda V r_i / R)

which depends on the chord alone -- no state -- and `c = N_c d_c` is linear
in the design vector, so "stay inside the cache at every operating point" is
a *linear* constraint on `u` with a constant Jacobian:

    c_min,i (1 + mu) <= [N_c (lo + u * span)]_i <= c_max,i (1 - mu)

    c_max,i = Re_hi nu / W_i(V_max),   c_min,i = Re_lo nu / W_i(V_min)

`mu = 0.05` is a stated margin, not a tuned one: SLSQP's intermediate
iterates may violate inequality constraints slightly, and an objective
evaluation that raises mid-line-search is an aborted run. If a
`PolarDomainError` still occurs the offending `u` is logged on the instance
and the error re-raised -- never a substitute value.

The angle-of-attack range of the cache (`xfoil_alpha_min/max`, -8..18 deg
for SG6043) is *state*-dependent (it needs the solved `phi`), so it is not a
constraint; `alpha_check` is a logged post-check on accepted iterates.

Provisional bounds
-------------------
`bounds` is a required argument with no default. The provisional numbers
live in one place, `tests/test_parameterisation.py::PROVISIONAL_BOUNDS`;
`chord_max_m = 0.45 m` there is a placeholder pending the hub-radius /
root-attachment decision, and every result produced through this class is
"under provisional bounds". Nothing here reads a bound from `config/`.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

import numpy as np

from config import load_design_rotor, load_site
from gradients.finite_difference import central_difference
from objective.objective import objective
from objective.power import aerodynamic_power, wind_speed_bins
from polars.interpolant import PolarDomainError
from polars.polar import interpolant_for

#: Envelope margin, stated (see module docstring). Raised to 0.10 once, and
#: recorded, if a run still trips `PolarDomainError` -- never tuned beyond.
DEFAULT_ENVELOPE_MARGIN = 0.05


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

    def operating_speeds(self):
        """
        Every wind speed the objective solves at: the 17 bin midpoints. The
        rating is a configured number, so the rated speed is not one of them.
        """

        _edges, midpoints, _width = wind_speed_bins()
        return [float(v) for v in midpoints]

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
            station, before the margin), `margin`, `reynolds_lo/hi`, and the
            speeds `v_min`, `v_max` the limits were computed at.
        """

        nu = float(self.site.kinematic_viscosity)
        tsr = float(self.design.design_tsr)
        R = self.parameterisation.radius_m
        radii = self.parameterisation.radii

        speeds = self.operating_speeds()
        v_min, v_max = min(speeds), max(speeds)
        w_at = lambda v: np.hypot(v, tsr * v * radii / R)  # noqa: E731

        return {
            "radii": radii,
            "chord_min_m": self.reynolds_lo * nu / w_at(v_min),
            "chord_max_m": self.reynolds_hi * nu / w_at(v_max),
            "margin": self.margin,
            "reynolds_lo": self.reynolds_lo,
            "reynolds_hi": self.reynolds_hi,
            "v_min": v_min,
            "v_max": v_max,
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

    # -- solidity (built, never run: no cap has been decided) ---------------

    def solidity_constraint(self, cap):
        """
        SciPy inequality `cap - sigma_i >= 0`, `sigma_i = B c_i / (2 pi r_i)`.

        `cap` is required and has no default: no solidity limit has been
        decided (it depends on the same undecided root-attachment concept as
        `chord_max_m`), so this constraint is not used in any run. It exists
        so that the day a value lands, it is one argument away.
        """

        cap = float(cap)
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

        Solves the rotor at every speed in `operating_speeds()` (17 solves,
        the cost of one objective evaluation). Returns a dict of lists, one
        entry per speed: `speeds`, `alpha_deg` (25,), `reynolds` (25,),
        `phi` (25,), `a` (25,), `power_w`, `converged`.
        """

        d = self.physical(u)
        geometry = self.parameterisation.to_geometry(d, polar_cache=self.polar_cache)
        tsr = float(self.design.design_tsr)
        rho = float(self.site.air_density)
        nu = float(self.site.kinematic_viscosity)

        state = {"speeds": [], "alpha_deg": [], "reynolds": [], "phi": [], "a": [],
                 "power_w": [], "converged": []}
        for v in self.operating_speeds():
            power, result = aerodynamic_power(geometry, v, tsr, rho, nu)
            stations = result["stations"]
            state["speeds"].append(v)
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
        """

        state = self.operating_state(u)
        alpha = np.array(state["alpha_deg"])
        reynolds = np.array(state["reynolds"])

        alpha_min, alpha_max = float(alpha.min()), float(alpha.max())
        return {
            "alpha_min_deg": alpha_min,
            "alpha_max_deg": alpha_max,
            "alpha_limits_deg": [self.alpha_min_deg, self.alpha_max_deg],
            "alpha_margin_deg": float(min(alpha_min - self.alpha_min_deg,
                                          self.alpha_max_deg - alpha_max)),
            "within": bool(alpha_min >= self.alpha_min_deg and alpha_max <= self.alpha_max_deg),
            "reynolds_min": float(reynolds.min()),
            "reynolds_max": float(reynolds.max()),
            "reynolds_limits": [self.reynolds_lo, self.reynolds_hi],
            "all_converged": all(state["converged"]),
        }
