"""
`ScaledProblem`: the optimisation problem as SLSQP sees it.

What an optimiser is handed
----------------------------
    variable    u in [0, 1]^10,  d = bounds.to_physical(u)
    objective   fun(u) = J(u) / |J(u0)|,  J = -AEP [MWh/yr]
    constraint  the polar-cache Reynolds envelope, linear in u (below)
    gradient    jac_fd(u, h): central differences of `fun`
                jac_adjoint(u): the discrete adjoint, same units

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

The mass problem (2026-09-20)
------------------------------
`docs/DESIGN-BASIS.md` §6. The energy objective above becomes
a constraint and the objective becomes the material proxy:

    minimise    mass(u)  = A_shell(d) / A_shell(x0)          `objective.mass`, no BEM
    subject to  AEP floor      -fun(u) - (1 - delta)  >= 0   `aep_floor_constraint(delta)`
                moment cap     KS0 - KS(u)            >= 0   `moment_constraint(0)`
                stress proxy   KS0/c00^2 - KS(u)/c0^2 >= 0   `stress_constraint()`
                deflection     D0 - D(u)              >= 0   `deflection_constraint()`
                monotone chord and twist control points     `manufacturing_constraints()`
                envelope, solidity, box                      unchanged

with `x0`, the committed Schmitz baseline, the reference of every relative
row (it has zero slack on each of them by construction) and
`constraints_for_mass_problem(delta)` the one assembly. The AEP row's
Jacobian is `-jac_adjoint`; the moment and stress rows read the one cached
`RootMomentSystem.gradient` at `u`; the deflection row reads a
`DeflectionSystem.gradient` at `u` built on the *same* forward state (one
9-point solve for three rows). The objective is geometry and cannot fail; the
four state rows can (`PolarDomainError`, or an unconverged solve), and inside
the mass problem their `fun` reports a large violation (`FAILED_SLACK`) and
logs the point in `evaluation_failures` so SLSQP's line search backs off,
while their `jac` raises -- a stale Jacobian is never substituted, and an
accepted iterate that cannot be evaluated is a failed start.

The absolute rows (2026-09-20, evening)
----------------------------------------
With the laminate resolved (`config/rotor_design.yaml::structure`,
`docs/OUTSTANDING-INPUTS.md` section 11) the two proxies have absolute
counterparts on the same adjoints, with the constants carried instead of
cancelled:

    root stress     sigma(u)     = KS(u) M_ref / (k_Z c0(u)^2 t)          [Pa]
    tip deflection  delta_tip(u) = delta_ref D(u) / (E k_I t)             [m]

    stress_absolute      1 - sigma(u) / (sigma_allow / SF)      >= 0
    deflection_absolute  1 - delta_tip(u) / tip_clearance_m     >= 0

Each is the natural inequality (`sigma_allow/SF - sigma >= 0`,
`delta_allow - delta_tip >= 0`) divided by its allowable, so that the row
SLSQP sees is of order one like every other row rather than of order
1e8 Pa; the feasible set is the same and `mass_problem_slacks` quotes the
MPa and mm. The stress Jacobian is the relative row's chain rule
(`dKS/dd / c0^2 - 2 KS / c0^3 e_0`) times `M_ref / (k_Z t sigma_design)`;
the deflection Jacobian is `dD/dd` times `delta_ref / (E k_I t delta_allow)`
-- the same cached `RootMomentSystem` and `DeflectionSystem` results, no new
derivation, so Tiers 1-2 are inherited and Tier 3 is re-run at the new
constants (`tests/test_mass_problem.py`, `verification/absolute_material/`).

WHICH ROWS ARE THE DESIGN ROWS, AND WHY. The relative rows stay the
production set (`MASS_PROBLEM_ROWS`, the committed `x_m`), and the absolute
rows are ADDITIONAL (`ABSOLUTE_ROWS`, assembled by name), for a reason that
is mechanics, not caution: the load set `L` is the *operating* set (the
uncapped bins and the rated point), and at those loads the thin-shell root
stress of `x0` is 22.6 MPa against a 196.5 MPa design allowable -- the
absolute row is slack by a factor of nine and would not size anything. The
loads that size a small blade's root are the ones IEC 61400-2 prescribes and
this model does not compute (parked at the 50-year gust, fatigue), and the
relative row `KS/c0^2 <= KS0/c00^2` is the statement that survives that
gap: whatever the extreme load is, the optimum's root carries no more stress
per unit of it than the reference's. Replacing the relative stress row by
the absolute one therefore removes the stress constraint in effect (the
optimum goes to the "no stress row" ablation of
`verification/mass_optimisation/`), which is measured and recorded in
`verification/absolute_material/`, not assumed. The moment cap `KS <= KS0`
is likewise kept: it is a LOAD cap -- the hub, shaft and tower see no more
flapwise moment than with the reference -- and not a strength check, which
is what the stress rows are. `deflection_absolute` needs `tip_clearance_m`,
which is machine geometry and still TODO; the row is built and tested with
an injected clearance and `absolute_rows_available()` says why it is not
assembled until then.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math

import numpy as np

from config import is_resolved, load_design_rotor, load_site
from gradients.finite_difference import central_difference
from objective.objective import objective
from objective.power import aerodynamic_power, operating_points, tsr_schedule
from polars.interpolant import PolarDomainError
from polars.polar import interpolant_for

#: Envelope margin, stated (see module docstring). Raised to 0.10 once, and
#: recorded, if a run still trips `PolarDomainError` -- never tuned beyond.
DEFAULT_ENVELOPE_MARGIN = 0.05

#: KS stiffness for the root-moment aggregate. The default the load
#: constraint is built at; `moment_report` also quotes 30 and 300.
DEFAULT_KS_RHO = 100.0

#: `|u|` or `|1 - u|` below this: the bound is active. `g(u)` below this: the
#: constraint row is active. One tolerance for every reporting script.
ACTIVE_TOL = 1e-6

#: The slack a state-dependent row of the mass problem reports at a trial
#: point it cannot evaluate: a violation large enough that SLSQP's merit
#: function rejects the step. Reported, never silently absorbed.
FAILED_SLACK = 1.0e3

#: The rows of the committed mass problem, in assembly order -- the relative
#: set every `verification/mass_optimisation/` artefact was run with.
MASS_PROBLEM_ROWS = ("envelope", "solidity", "manufacturing", "aep_floor", "moment",
                     "stress", "deflection")

#: The absolute rows (module docstring): additional, assembled by name, never
#: part of `MASS_PROBLEM_ROWS` so the committed runs reproduce as run.
ABSOLUTE_ROWS = ("stress_absolute", "deflection_absolute")

#: Every row `mass_problem_rows` knows, in assembly order.
ALL_MASS_PROBLEM_ROWS = MASS_PROBLEM_ROWS + ABSOLUTE_ROWS


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

        # Load constraint: built on first use from the
        # committed baseline. KS0 is `KS_rho(u0)` alongside J0; n_moment_solves
        # counts the 9-point forward solves the constraint's fun/jac need.
        self.KS0 = None
        self.m_ref = None
        self.u0 = None
        self.n_moment_solves = 0
        self._systems = {}
        self._ks0 = {}
        self._moment_cache = None
        self._moment_cache_key = None

        # The mass problem (2026-09-20): the material model and its reference,
        # the deflection system with its reference, its own cache, and the
        # log of trial points a state row could not evaluate.
        self._material = None
        self.material_ref = None
        self._deflection = None
        self.delta_ref = None
        self.D0 = None
        self.n_deflection_solves = 0
        self._deflection_cache = None
        self._deflection_cache_key = None
        self.evaluation_failures = []
        # `(u bytes, error)` of the last trial point whose forward load solve
        # failed: the moment, stress and deflection rows share that solve, so
        # the second and third row re-raise instead of repeating it.
        self._load_failure = None

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
        `gradient(d) * span / |J0|`.

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

    # -- the root-moment functional -----------------------

    def _system_for(self, rho):
        """
        The `adjoint.loads.RootMomentSystem` at KS stiffness `rho`, built on
        first use and cached. The default (`rho = 100`) is the one
        `moment_constraint` uses; the others exist so `moment_report` can quote
        the conservatism at 30 and 300 from one set of moments.

        `M_ref = M(x0)` and the starting point are read once from the committed
        Schmitz baseline (`design.baseline.build_schmitz_baseline`, which
        reproduces `verification/baseline/x0.json` bit for bit), so the
        normalisation is one number in one place. `KS0 = KS_rho(u0)` is measured
        in the same pass and cached; that is what makes the `eps = 0`
        constraint exactly zero-slack at `u0`.
        """

        if rho in self._systems:
            return self._systems[rho]

        from adjoint.loads import RootMomentSystem
        from objective.loads import root_moment

        if self.m_ref is None:
            x0 = self.reference_design()
            geometry = self.parameterisation.to_geometry(x0, polar_cache=self.polar_cache)
            v_rated = float(self.design.rated_wind_speed_ms)
            tsr_rated = tsr_schedule(v_rated, self.design.design_tsr,
                                     self.design.max_tip_speed_ms)
            self.m_ref = float(root_moment(
                geometry, v_rated, tsr_rated, float(self.site.air_density),
                float(self.site.kinematic_viscosity))[0])

        system = RootMomentSystem(self.parameterisation, self.bounds, self.resource,
                                  polar_cache=self.polar_cache, rho=rho,
                                  m_ref_nm=self.m_ref)
        self._systems[rho] = system

        # The forward state at u0 does not depend on rho: a second stiffness
        # reuses the solve the first one made rather than repeating it.
        u0_key = np.asarray(self.u0, dtype=float).tobytes()
        state = None
        if self._moment_cache is not None and self._moment_cache_key[1] == u0_key:
            state = self._moment_cache.state
        else:
            self.n_moment_solves += 1
        result = system.gradient(self.physical(self.u0), state=state)
        self._ks0[rho] = float(result.KS)
        self._moment_cache_key = (float(rho), u0_key)
        self._moment_cache = result
        if rho == DEFAULT_KS_RHO:
            self.KS0 = float(result.KS)
        return system

    def reference_design(self):
        """
        `x0`, the committed Schmitz baseline (`design.baseline.build_schmitz_baseline`,
        which reproduces `verification/baseline/x0.json` bit for bit), and the
        one place `self.u0` is set. Every relative row normalises on it.
        """

        from design.baseline import build_schmitz_baseline

        x0 = np.asarray(build_schmitz_baseline().design_vector, dtype=float)
        if self.u0 is None:
            self.u0 = self.scaled(x0)
        return x0

    def load_system(self):
        """The default (`rho = 100`) `RootMomentSystem` over the load set `L`."""

        return self._system_for(DEFAULT_KS_RHO)

    def _shared_load_state(self, u_key):
        """
        A forward state over `L` already solved at this exact `u`, from either
        the moment or the deflection cache, or `None`. The two systems share
        the points, so a state one of them solved is the other's.
        """

        if self._moment_cache is not None and self._moment_cache_key[1] == u_key:
            return self._moment_cache.state
        if self._deflection_cache is not None and self._deflection_cache_key == u_key:
            return self._deflection_cache.state
        return None

    def _moment_at(self, u, rho):
        """
        The cached `RootMomentSystem.gradient(d(u))` at stiffness `rho`.

        SLSQP evaluates a constraint's `fun` and `jac` separately at the same
        `u` and `fun` again in the line search; one cache entry keyed on the
        exact `u` bytes makes that cost one 9-point solve, not three.
        """

        u = np.asarray(u, dtype=float)
        key = (float(rho), u.tobytes())
        if key == self._moment_cache_key:
            return self._moment_cache
        self._raise_if_failed(key[1])
        system = self._system_for(rho)
        state = self._shared_load_state(key[1])
        if state is None:
            self.n_moment_solves += 1
        try:
            result = system.gradient(self.physical(u), state=state)
        except (PolarDomainError, RuntimeError) as error:
            self._load_failure = (key[1], error)
            raise
        self._moment_cache_key = key
        self._moment_cache = result
        return result

    def _raise_if_failed(self, u_key):
        """Re-raise the remembered failure of the forward load solve at this `u`."""

        if self._load_failure is not None and self._load_failure[0] == u_key:
            raise self._load_failure[1]

    def moment_state(self, u):
        """The cached moment-adjoint result at `u` (default stiffness)."""

        return self._moment_at(u, DEFAULT_KS_RHO)

    def moment_ks(self, u, rho=None):
        """`KS_rho(u)`, as a float, at the default stiffness (or `rho`)."""

        rho = DEFAULT_KS_RHO if rho is None else float(rho)
        return float(self._moment_at(u, rho).KS)

    def _moment_reference(self, rho):
        """`KS0 = KS_rho(u0)` at the committed baseline, measured once."""

        if rho not in self._ks0:
            self._system_for(rho)
        return self._ks0[rho]

    def moment_constraint(self, eps, rho=None):
        """
        SciPy inequality for the root-moment load cap:

            g_eps(u) = (1 - eps) KS0 - KS(u)  =  (KS0 - KS(u)) - eps KS0  >= 0

        written in the regrouped form so the `eps = 0` slack at `u0` and the
        `eps`-linear slack are exact to the bit (the two groupings differ in
        the last bit). `jac` is the adjoint Jacobian `-dKS/dd * span`.
        """

        rho = DEFAULT_KS_RHO if rho is None else float(rho)
        ks0 = self._moment_reference(rho)
        eps = float(eps)
        span = self.bounds.span()

        def fun(u):
            return np.array([(ks0 - self._moment_at(u, rho).KS) - eps * ks0])

        def jac(u):
            return (-self._moment_at(u, rho).dKS_dd * span)[None, :]

        return {"type": "ineq", "fun": fun, "jac": jac}

    def constraints_with_moment(self, eps):
        """The full inequality set: envelope, solidity, moment."""

        return self.constraints() + [self.moment_constraint(eps)]

    def moment_active(self, u, eps, rho=None, tol=ACTIVE_TOL):
        """The moment row's slack at `u` for reduction fraction `eps`."""

        rho = DEFAULT_KS_RHO if rho is None else float(rho)
        ks0 = self._moment_reference(rho)
        ks = float(self._moment_at(u, rho).KS)
        slack = (ks0 - ks) - float(eps) * ks0
        return {"row": "moment", "eps": float(eps), "rho": float(rho),
                "KS": ks, "KS0": ks0, "slack": float(slack),
                "active": bool(slack < tol)}

    def moment_report(self, u):
        """
        Everything the load-constraint artefacts record at `u`, forward path
        only: the per-point moments and softmax weights, `KS` and its
        conservatism `KS - max` over `rho in {30, 100, 300}`, the
        design-condition thrust and `Ct`, and the cut-out
        post-check at 20 m/s on the ceiling (reported, never constrained).
        """

        from objective.loads import ks, ks_weights, root_moment

        u = np.asarray(u, dtype=float)
        system = self.load_system()
        d = self.physical(u)
        state = system.solve(d)
        parts = system.partials(state.phi, state.d)
        moments = system.moments_from_m(parts["m"])
        m_ref = float(system.m_ref)
        normalised = moments / m_ref
        weights = ks_weights(normalised, system.rho)
        peak = float(np.max(normalised))

        points = []
        for b, (v, lam) in enumerate(system.points):
            points.append({
                "v_ms": float(v),
                "tsr": float(lam),
                "rpm": float(system.omega[b] * 60.0 / (2.0 * math.pi)),
                "moment_nm": float(moments[b]),
                "moment_normalised": float(normalised[b]),
                "softmax_weight": float(weights[b]),
            })

        rhos = (30.0, 100.0, 300.0)
        ks_by_rho = {f"{rho:g}": float(ks(normalised, rho)) for rho in rhos}
        conservatism = {f"{rho:g}": float(ks(normalised, rho) - peak) for rho in rhos}

        geometry = self.parameterisation.to_geometry(d, polar_cache=self.polar_cache)
        rho_air = float(self.site.air_density)
        nu = float(self.site.kinematic_viscosity)
        area = math.pi * self.parameterisation.radius_m ** 2

        v_rated, tsr_rated = system.points[-1]
        rated_moment, rated_result = root_moment(geometry, v_rated, tsr_rated,
                                                 rho_air, nu)
        Ct_rated = float(rated_result["Ct"])

        v_cut = float(self.design.cut_out_wind_speed_ms)
        tsr_cut = tsr_schedule(v_cut, self.design.design_tsr,
                               self.design.max_tip_speed_ms)
        cut_moment, cut_result = root_moment(geometry, v_cut, tsr_cut, rho_air, nu)
        cut_alpha = [math.degrees(s["alpha"]) for s in cut_result["stations"]]
        Ct_cut = float(cut_result["Ct"])
        alpha_max_deg = float(max(cut_alpha))

        return {
            "u": [float(x) for x in u],
            "x": [float(x) for x in d],
            "m_ref_nm": m_ref,
            "m_ref_label": ("M(x0) at 11 m/s, lambda = 5.711986642890533 "
                            "(300 rpm), forward path"),
            "KS0": None if self.KS0 is None else float(self.KS0),
            "rho": float(system.rho),
            "points": points,
            "peak_normalised": peak,
            "KS_by_rho": ks_by_rho,
            "conservatism_by_rho": conservatism,
            "softmax_weight_on_rated": float(weights[-1]),
            "design_condition": {
                "v_ms": float(v_rated),
                "tsr": float(tsr_rated),
                "rpm": float(system.omega[-1] * 60.0 / (2.0 * math.pi)),
                "Ct": Ct_rated,
                "thrust_n": float(Ct_rated * 0.5 * rho_air * float(v_rated) ** 2 * area),
                "moment_nm": float(rated_moment),
                "moment_over_ref": float(rated_moment / m_ref),
            },
            "cut_out": {
                "label": ("B3-dependent: the model holds P = P_rated with no "
                          "mechanism, so the state here is not the machine's; "
                          f"alpha up to {alpha_max_deg:.1f} deg on Viterna"),
                "B3_dependent": True,
                "v_ms": v_cut,
                "tsr": float(tsr_cut),
                "rpm": float(tsr_cut * v_cut / self.parameterisation.radius_m
                             * 60.0 / (2.0 * math.pi)),
                "Ct": Ct_cut,
                "thrust_n": float(Ct_cut * 0.5 * rho_air * v_cut ** 2 * area),
                "moment_nm": float(cut_moment),
                "alpha_min_deg": float(min(cut_alpha)),
                "alpha_max_deg": alpha_max_deg,
                "converged": bool(cut_result["converged"]),
            },
        }

    # -- the mass problem (2026-09-20) ------------------------------------------

    def material_model(self):
        """
        The `objective.mass.MaterialModel` over this parameterisation (the
        configured `objective.mass_model`), built on first use; the
        normaliser `material_ref = value(x0)` is measured in the same pass.
        """

        if self._material is None:
            from objective.mass import MaterialModel

            self._material = MaterialModel(self.parameterisation, polar_cache=self.polar_cache,
                                           design=self.design)
            # On `physical(u0)`, not `x0`: the scaled round trip differs from
            # `x0` in the last bit, and the objective must be 1.0 at `u0`
            # exactly as `fun(u0)` is -1.0 (KS0 and D0 are measured the same way).
            self.reference_design()
            self.material_ref = float(self._material.value(self.physical(self.u0)))
        return self._material

    def mass(self, u):
        """
        The mass-problem objective: the material proxy at `u` as a fraction
        of the reference's, `value(d) / value(x0)`. Exactly `1.0` at `u0`.
        Geometry only -- no BEM, cannot fail.
        """

        model = self.material_model()
        return float(model.value(self.physical(u))) / self.material_ref

    def mass_jac(self, u):
        """`d mass / du = gradient(d) * span / value(x0)` -- exact, constant for the shell."""

        model = self.material_model()
        return model.gradient(self.physical(u)) * self.bounds.span() / self.material_ref

    def material_report(self, u):
        """Both proxies at `u`, and as fractions of `x0`'s (`MaterialModel.report`)."""

        model = self.material_model()
        return model.report(self.physical(u), reference=self.physical(self.u0))

    def aep_floor_constraint(self, delta):
        """
        SciPy inequality for the energy floor `AEP(u) >= (1 - delta) AEP(x0)`:

            g_delta(u) = -fun(u) - (1 - delta)  >=  0

        `fun = J / |J0|` with `J0 = J(u0)` (`set_reference(u0)` must have been
        called with the reference), so `-fun(u0) = 1.0` to the bit and the
        `delta = 0` slack at `u0` is exactly zero. `jac` is `-jac_adjoint`,
        the verified objective adjoint as a constraint Jacobian.
        """

        if self.J0 is None:
            raise RuntimeError("J0 not set: call set_reference(u0) with the reference blade first")
        delta = float(delta)
        if not 0.0 <= delta < 1.0:
            raise ValueError(f"delta must be in [0, 1), got {delta}")
        floor = 1.0 - delta

        def fun(u):
            return np.array([-self.fun(u) - floor])

        def jac(u):
            return (-self.jac_adjoint(u))[None, :]

        return {"type": "ineq", "fun": fun, "jac": jac}

    def _root_chord(self, u):
        """`c0 = d_0`: the chord at the clamp `r_hub` (`s = 0`), which the clamped
        spline interpolates, so its basis row is `e_0`."""

        return float(self.physical(u)[0])

    def _stress_reference(self):
        """`KS0 / c00^2` at the committed reference, measured once."""

        ks0 = self._moment_reference(DEFAULT_KS_RHO)
        return ks0, self._root_chord(self.u0)

    def stress_ratio(self, u):
        """`KS(u) / c0(u)^2` over `KS0 / c00^2`: the root-stress proxy relative to `x0`."""

        ks0, c00 = self._stress_reference()
        return (self.moment_ks(u) / self._root_chord(u) ** 2) / (ks0 / c00 ** 2)

    def stress_constraint(self):
        """
        SciPy inequality for the root-stress proxy, `sigma ~ KS_rho(M) / Z`
        with `Z ~ c0^2 t` for a thin shell of constant laminate thickness:

            g(u) = KS0 / c00^2 - KS(u) / c0(u)^2  >=  0

        zero at `u0` to the bit. `jac = -(dKS/dd / c0^2 - 2 KS / c0^3 e_0) * span`
        with `dKS/dd`, `KS` from the cached moment state at `u`.
        """

        ks0, c00 = self._stress_reference()
        reference = ks0 / c00 ** 2
        span = self.bounds.span()
        e0 = np.zeros(self.n)
        e0[0] = 1.0

        def fun(u):
            c0 = self._root_chord(u)
            return np.array([reference - self.moment_state(u).KS / c0 ** 2])

        def jac(u):
            c0 = self._root_chord(u)
            state = self.moment_state(u)
            g = state.dKS_dd / c0 ** 2 - 2.0 * state.KS / c0 ** 3 * e0
            return (-g * span)[None, :]

        return {"type": "ineq", "fun": fun, "jac": jac}

    def _deflection_system(self):
        """
        The `adjoint.deflection.DeflectionSystem` over `L` at the default
        stiffness, built on first use. `delta_ref` is the rated-point
        deflection of the committed `x0` through the forward path
        (`objective.loads.tip_deflection_at`, as `m_ref` is `root_moment`);
        `D0 = D(u0)` is measured in the same pass on the moment row's state.
        """

        if self._deflection is not None:
            return self._deflection

        from adjoint.deflection import DeflectionSystem
        from objective.loads import tip_deflection_at

        x0 = self.reference_design()
        self._system_for(DEFAULT_KS_RHO)  # m_ref, KS0 and the u0 state
        geometry = self.parameterisation.to_geometry(x0, polar_cache=self.polar_cache)
        v_rated = float(self.design.rated_wind_speed_ms)
        tsr_rated = tsr_schedule(v_rated, self.design.design_tsr, self.design.max_tip_speed_ms)
        self.delta_ref = float(tip_deflection_at(
            geometry, v_rated, tsr_rated, float(self.site.air_density),
            float(self.site.kinematic_viscosity))[0])

        system = DeflectionSystem(self.parameterisation, self.bounds, self.resource,
                                  polar_cache=self.polar_cache, rho=DEFAULT_KS_RHO,
                                  m_ref_nm=self.m_ref, delta_ref=self.delta_ref)
        self._deflection = system
        u0_key = np.asarray(self.u0, dtype=float).tobytes()
        state = self._shared_load_state(u0_key)
        if state is None:
            self.n_deflection_solves += 1
        result = system.gradient(self.physical(self.u0), state=state)
        self._deflection_cache_key = u0_key
        self._deflection_cache = result
        self.D0 = float(result.D)
        return system

    def _deflection_at(self, u):
        """The cached `DeflectionSystem.gradient(d(u))`, on the moment row's state if it has one."""

        u = np.asarray(u, dtype=float)
        key = u.tobytes()
        system = self._deflection_system()
        if key == self._deflection_cache_key:
            return self._deflection_cache
        self._raise_if_failed(key)
        state = self._shared_load_state(key)
        if state is None:
            self.n_deflection_solves += 1
        try:
            result = system.gradient(self.physical(u), state=state)
        except (PolarDomainError, RuntimeError) as error:
            self._load_failure = (key, error)
            raise
        self._deflection_cache_key = key
        self._deflection_cache = result
        return result

    def deflection_state(self, u):
        """The cached deflection-adjoint result at `u`."""

        return self._deflection_at(u)

    def deflection_ks(self, u):
        """`D(u) = KS_rho(delta / delta_ref)` as a float."""

        return float(self._deflection_at(u).D)

    def deflection_constraint(self):
        """
        SciPy inequality for the tip-deflection proxy, `g(u) = D0 - D(u) >= 0`,
        zero at `u0`; `jac = -dD/dd * span` by the deflection adjoint.
        """

        self._deflection_system()
        d0 = self.D0
        span = self.bounds.span()

        def fun(u):
            return np.array([d0 - self._deflection_at(u).D])

        def jac(u):
            return (-self._deflection_at(u).dD_dd * span)[None, :]

        return {"type": "ineq", "fun": fun, "jac": jac}

    # -- the absolute rows (2026-09-20, evening) -----------------------------------

    def mass_kg(self, u):
        """The configured proxy at `u` as kg per blade (`MaterialModel.mass_kg`), or `None`."""

        return self.material_model().mass_kg(self.physical(u))

    def design_allowable_stress_pa(self):
        """`allowable_stress_pa / safety_factor`; raises while either is TODO."""

        return float(self.design.design_allowable_stress_pa)

    def _stress_scale_pa(self):
        """`M_ref / (k_Z t)`: what multiplies `KS / c0^2` to make pascals."""

        self._system_for(DEFAULT_KS_RHO)  # m_ref
        factor = self.material_model().section_modulus_factor_m()
        if factor is None:
            float(self.design.shell_thickness_m)  # raises, naming the field
        return float(self.m_ref) / factor

    def stress_absolute_pa(self, u):
        """
        The thin-shell root stress at `u`, Pa: `KS(u) M_ref / (k_Z c0(u)^2 t)`
        -- the KS aggregate of the load set's root moments over the section
        modulus of the root chord (module docstring).
        """

        return self._stress_scale_pa() * self.moment_ks(u) / self._root_chord(u) ** 2

    def stress_constraint_absolute(self):
        """
        SciPy inequality for the absolute root stress,

            g(u) = 1 - sigma(u) / sigma_design,  sigma_design = sigma_allow / SF

        i.e. `sigma_allow/SF - sigma(u) >= 0` divided by the allowable.
        `jac = -(M_ref / (k_Z t sigma_design)) (dKS/dd / c0^2 - 2 KS / c0^3 e_0) span`
        on the cached moment state -- the relative row's chain rule with the
        constants carried.
        """

        scale = self._stress_scale_pa() / self.design_allowable_stress_pa()
        span = self.bounds.span()
        e0 = np.zeros(self.n)
        e0[0] = 1.0

        def fun(u):
            c0 = self._root_chord(u)
            return np.array([1.0 - scale * self.moment_state(u).KS / c0 ** 2])

        def jac(u):
            c0 = self._root_chord(u)
            state = self.moment_state(u)
            g = state.dKS_dd / c0 ** 2 - 2.0 * state.KS / c0 ** 3 * e0
            return (-scale * g * span)[None, :]

        return {"type": "ineq", "fun": fun, "jac": jac}

    def _deflection_scale_m(self):
        """`delta_ref / (E k_I t)`: what multiplies `D` to make metres."""

        self._deflection_system()  # delta_ref
        factor = self.material_model().stiffness_factor_pa_m()
        if factor is None:
            float(self.design.youngs_modulus_pa)
            float(self.design.shell_thickness_m)
        return float(self.delta_ref) / factor

    def deflection_absolute_m(self, u):
        """
        The thin-shell tip deflection at `u`, m: `delta_ref D(u) / (E k_I t)`,
        the KS aggregate over the load set of the rated-point-normalised
        deflections, made absolute (module docstring).
        """

        return self._deflection_scale_m() * self.deflection_ks(u)

    def deflection_constraint_absolute(self):
        """
        SciPy inequality for the absolute tip deflection,

            g(u) = 1 - delta_tip(u) / tip_clearance_m

        i.e. `tip_clearance_m - delta_tip(u) >= 0` divided by the clearance;
        `jac = -(delta_ref / (E k_I t delta_allow)) dD/dd span` on the cached
        deflection state. Raises `UnresolvedConfigError` naming
        `tip_clearance_m` while the clearance is TODO.
        """

        clearance = float(self.design.tip_clearance_m)
        scale = self._deflection_scale_m() / clearance
        span = self.bounds.span()

        def fun(u):
            return np.array([1.0 - scale * self._deflection_at(u).D])

        def jac(u):
            return (-scale * self._deflection_at(u).dD_dd * span)[None, :]

        return {"type": "ineq", "fun": fun, "jac": jac}

    def absolute_rows_available(self):
        """
        Which of `ABSOLUTE_ROWS` can be assembled from the config as it
        stands: `{row: True}` or `{row: "TODO: <the field's note>"}`. What an
        artefact records instead of silently omitting a row.
        """

        needs = {
            "stress_absolute": ("shell_thickness_m", "allowable_stress_pa", "safety_factor"),
            "deflection_absolute": ("shell_thickness_m", "youngs_modulus_pa", "tip_clearance_m"),
        }
        out = {}
        for row, fields in needs.items():
            missing = [f for f in fields if not is_resolved(getattr(self.design, f))]
            out[row] = True if not missing else "TODO: " + "; ".join(
                f"{f}: {getattr(self.design, f).note}" for f in missing)
        return out

    def manufacturing_rows(self):
        """
        The difference matrix `D` (rows x n) of the enabled monotone blocks:
        row `c_i - c_{i+1}` for chord, `theta_i - theta_{i+1}` for twist, per
        `manufacturing.monotone_chord` / `monotone_twist` in config. Non-
        increasing control points give a non-increasing B-spline (variation
        diminishing), so with the min-chord floor (`min_chord_rows`) these
        linear rows are the whole manufacturability set. Empty (0 x n) when
        both are off.
        """

        n_c = self.parameterisation.n_chord
        n_t = self.parameterisation.n_twist
        rows, labels = [], []
        blocks = []
        if self.design.monotone_chord:
            blocks.append(("chord", 0, n_c))
        if self.design.monotone_twist:
            blocks.append(("twist", n_c, n_t))
        for name, offset, count in blocks:
            for i in range(count - 1):
                row = np.zeros(self.n)
                row[offset + i] = 1.0
                row[offset + i + 1] = -1.0
                rows.append(row)
                labels.append(f"monotone {name}_{i}-{name}_{i + 1}")
        matrix = np.vstack(rows) if rows else np.zeros((0, self.n))
        return matrix, labels

    def min_chord_rows(self):
        """
        The buildable-tip floor `c_i - min_chord_m >= 0` on every chord
        control point (`manufacturing.min_chord_m`, 2026-09-20): `(matrix,
        rhs, labels)` with `matrix = [I 0]` and `rhs = min_chord_m`. A row
        rather than the box bound so the scaling of `u` stays as
        committed. Empty when the floor is at or below the box bound.
        """

        n_c = self.parameterisation.n_chord
        floor = float(self.design.min_chord_m)
        if floor <= self.bounds.chord_min_m:
            return np.zeros((0, self.n)), np.zeros(0), []
        matrix = np.zeros((n_c, self.n))
        matrix[np.arange(n_c), np.arange(n_c)] = 1.0
        return matrix, np.full(n_c, floor), [f"min chord chord_{i}" for i in range(n_c)]

    def manufacturing_row_labels(self):
        return self.manufacturing_rows()[1] + self.min_chord_rows()[2]

    def manufacturing_constraints(self):
        """
        SciPy inequality `A (lo + u span) - b >= 0` over the monotone rows
        (`b = 0`) and the min-chord rows (`b = min_chord_m`), in that order;
        constant Jacobian `A span`.
        """

        monotone, _labels = self.manufacturing_rows()
        floor, rhs, _floor_labels = self.min_chord_rows()
        matrix = np.vstack([monotone, floor])
        b = np.concatenate([np.zeros(monotone.shape[0]), rhs])
        jac = matrix * self.bounds.span()[None, :]
        offset = matrix @ self.bounds.lower() - b

        return {
            "type": "ineq",
            "fun": lambda u: jac @ np.asarray(u, dtype=float) + offset,
            "jac": lambda u: jac,
        }

    def _guarded(self, label, constraint):
        """
        The mass-problem wrapper for a state-dependent row: `fun` returns
        `-FAILED_SLACK` and logs the point when the row cannot be evaluated
        (outside the polar cache, or an unconverged solve); `jac` is left to
        raise. See the module docstring.
        """

        fun, jac = constraint["fun"], constraint["jac"]

        def guarded_fun(u):
            try:
                return fun(u)
            except (PolarDomainError, RuntimeError) as error:
                self.evaluation_failures.append({
                    "row": label,
                    "u": [float(x) for x in np.asarray(u, dtype=float)],
                    "error": f"{type(error).__name__}: {error}",
                })
                return np.array([-FAILED_SLACK])

        return {"type": "ineq", "fun": guarded_fun, "jac": jac}

    def mass_problem_rows(self, delta, include=MASS_PROBLEM_ROWS):
        """
        The rows of the mass problem as `(label, constraint)` pairs, in
        `MASS_PROBLEM_ROWS` order; `include` drops rows for an ablation, never
        reorders them. The state rows are guarded.
        """

        builders = {
            "envelope": lambda: self.envelope_constraint(),
            "solidity": lambda: self.solidity_constraint(),
            "manufacturing": lambda: self.manufacturing_constraints(),
            "aep_floor": lambda: self._guarded("aep_floor", self.aep_floor_constraint(delta)),
            "moment": lambda: self._guarded("moment", self.moment_constraint(0.0)),
            "stress": lambda: self._guarded("stress", self.stress_constraint()),
            "deflection": lambda: self._guarded("deflection", self.deflection_constraint()),
            "stress_absolute": lambda: self._guarded("stress_absolute",
                                                     self.stress_constraint_absolute()),
            "deflection_absolute": lambda: self._guarded("deflection_absolute",
                                                         self.deflection_constraint_absolute()),
        }
        unknown = set(include) - set(ALL_MASS_PROBLEM_ROWS)
        if unknown:
            raise ValueError(f"unknown mass-problem rows {sorted(unknown)}")
        return [(name, builders[name]()) for name in ALL_MASS_PROBLEM_ROWS if name in include]

    def constraints_for_mass_problem(self, delta, include=MASS_PROBLEM_ROWS):
        """The full inequality set of the mass problem at energy floor `delta`."""

        return [constraint for _label, constraint in self.mass_problem_rows(delta, include)]

    def mass_problem_slacks(self, u, delta, tol=ACTIVE_TOL):
        """
        Every scalar row of the mass problem at `u`: its slack, whether it is
        active, and the ratios the artefacts report (`aep_over_ref`,
        `ks_over_ks0`, `stress_ratio`, `deflection_ratio`), plus the
        absolute rows with their MPa and mm (`stress_absolute`,
        `deflection_absolute`; each row's `available` says whether it could
        be assembled, and a row that could not carries its numbers where the
        inputs allow and nothing where they do not). Evaluated unguarded --
        a point that cannot be evaluated raises here.
        """

        aep_slack = float(self.aep_floor_constraint(delta)["fun"](u)[0])
        moment = self.moment_active(u, 0.0, tol=tol)
        stress_slack = float(self.stress_constraint()["fun"](u)[0])
        deflection_slack = float(self.deflection_constraint()["fun"](u)[0])
        labels = self.manufacturing_row_labels()
        mfg = self.manufacturing_constraints()["fun"](u)
        available = self.absolute_rows_available()

        material = self.material_model()
        stress_abs = {"available": available["stress_absolute"]}
        if material.section_modulus_factor_m() is not None:
            stress_abs["stress_mpa"] = float(self.stress_absolute_pa(u)) / 1e6
        if available["stress_absolute"] is True:
            allowable = self.design_allowable_stress_pa()
            slack = float(self.stress_constraint_absolute()["fun"](u)[0])
            stress_abs.update({"allowable_mpa": allowable / 1e6,
                               "ultimate_mpa": float(self.design.allowable_stress_pa) / 1e6,
                               "safety_factor": float(self.design.safety_factor),
                               "slack": slack, "slack_mpa": slack * allowable / 1e6,
                               "active": bool(slack < tol)})

        deflection_abs = {"available": available["deflection_absolute"]}
        if material.stiffness_factor_pa_m() is not None:
            deflection_abs["tip_deflection_mm"] = float(self.deflection_absolute_m(u)) * 1e3
        if available["deflection_absolute"] is True:
            clearance = float(self.design.tip_clearance_m)
            slack = float(self.deflection_constraint_absolute()["fun"](u)[0])
            deflection_abs.update({"clearance_mm": clearance * 1e3, "slack": slack,
                                   "slack_mm": slack * clearance * 1e3,
                                   "active": bool(slack < tol)})
        return {
            "aep_floor": {"delta": float(delta), "slack": aep_slack,
                          "aep_over_ref": float(-self.fun(u)),
                          "active": bool(aep_slack < tol)},
            "moment": {"slack": moment["slack"], "ks_over_ks0": moment["KS"] / moment["KS0"],
                       "active": moment["active"]},
            "stress": {"slack": stress_slack, "stress_ratio": float(self.stress_ratio(u)),
                       "active": bool(stress_slack < tol)},
            "deflection": {"slack": deflection_slack,
                           "deflection_ratio": float(self.deflection_ks(u) / self.D0),
                           "active": bool(deflection_slack < tol)},
            "manufacturing": {"rows": [labels[k] for k in range(len(labels)) if mfg[k] < tol],
                              "slack_min": float(mfg.min()) if len(mfg) else None},
            "stress_absolute": stress_abs,
            "deflection_absolute": deflection_abs,
        }

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

    def active_set(self, u, tol=ACTIVE_TOL, moment_eps=None, mass_delta=None):
        """
        Which bounds and constraint rows are active at `u`: a dict with
        `bounds_lower`, `bounds_upper` (variable indices), `envelope`,
        `solidity` (row labels), the tightest row of each constraint with
        its slack, and every station's solidity.

        With `moment_eps` given, the moment row is attached too (its
        slack is `(KS0 - KS(u)) - eps KS0`, from `moment_active`). With
        `mass_delta` given, every scalar row of the mass problem is attached
        under `mass_problem` (`mass_problem_slacks`).
        """

        u = np.asarray(u, dtype=float)
        env = self.envelope_constraint()["fun"](u)
        sol = self.solidity_constraint()["fun"](u)
        env_labels = self.envelope_row_labels()
        sol_labels = self.solidity_row_labels()
        chord = self.parameterisation.chord(self.physical(u))
        sigma = self.design.n_blades * chord / (2.0 * math.pi * self.parameterisation.radii)
        report = {
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
        if moment_eps is not None:
            moment = self.moment_active(u, moment_eps, tol=tol)
            report["moment"] = moment
            report["moment_tightest"] = {"row": "moment", "slack": moment["slack"]}
        if mass_delta is not None:
            report["mass_problem"] = self.mass_problem_slacks(u, mass_delta, tol=tol)
        return report

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
