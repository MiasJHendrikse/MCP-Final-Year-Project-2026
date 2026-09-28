"""
EXPERIMENT (2026-09-19, not a project result): the AEP objective under a
rotor-speed ceiling, and its discrete adjoint.

The project's objective (`src/objective`) is already `J = -AEP`, but with
`lambda = 6.5` at every wind-speed bin, which makes AEP a fixed convex
combination of `Cp(6.5, Re_b; d)` and the optimum a Cp-at-one-TSR optimum
(`docs/DESIGN-BASIS.md` section 2). This module puts the one thing
that gives AEP a TSR dimension -- a maximum rotor speed -- into the
operating law,

    lambda_b = min(lambda_opt, V_tip,max / V_b),    lambda_opt = 6.5,

and extends the discrete adjoint to it. Nothing in `src/` is touched: the
two classes below are subclasses that override exactly what the ceiling
changes and inherit everything else.

What the ceiling changes in the adjoint, and what it does not
--------------------------------------------------------------
`adjoint.system.BEMSystem` carries one `Omega_b = lambda V_b / R` per
operating point in `self.omega`; every place lambda enters -- the local
speed ratio `lam_r = Omega_b r_i / V_b`, the zero-induction Reynolds
estimate `W = hypot(V_b, Omega_b r_i)`, and the power assembly
`P_b = Omega_b sum_i t_i q_{b,i}` -- reads that list. With `lambda_b` per
bin, `Omega_b = lambda_b V_b / R` and the same code runs. The state
Jacobian stays diagonal, `dR/dd` keeps its form, and

    dJ/dd = sum_b omega_b dP_b/dd,    omega_b = -(T/1e6) m_b  (0 if capped)

is the weighted sum of per-bin adjoint sensitivities the experiment
needs -- the parent's `gradient()` already assembles it that way. The
forward solve is the one method that has to be re-stated, because the
parent passes the scalar `self.tsr` to `aerodynamic_power`.

Independence of the two AEP evaluations
----------------------------------------
`evaluate_aep` (the optimiser's objective, and the report's numbers) goes
`aerodynamic_power -> solve_rotor -> Cp`, the project's forward chain.
`CeilingBEMSystem.J(phi, d)` re-assembles the same energy from the station
integrand `q` through `adjoint.kernels.station_partials`, the trapezoid
weights and `Omega_b` -- a second code path with no shared arithmetic
above the residual. The two agreeing at every optimum (to ~1e-12) is the
check that the adjoint's own picture of J is the thing being optimised.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math

import numpy as np

from adjoint.system import BEMSystem, ForwardState
from gradients import ScaledProblem
from objective.objective import HOURS_PER_YEAR
from objective.power import aerodynamic_power, wind_speed_bins
from polars.interpolant import PolarDomainError


def tsr_schedule(v, lambda_opt, vtip_max):
    """`lambda(V) = min(lambda_opt, V_tip,max / V)`; `vtip_max = None` is no ceiling."""

    if vtip_max is None:
        return float(lambda_opt)
    return float(min(lambda_opt, vtip_max / v))


class CeilingBEMSystem(BEMSystem):
    """
    `BEMSystem` with `lambda_b = min(lambda_opt, V_tip,max / V_b)` per bin.

    `vtip_max=None` reproduces the parent exactly (same `omega`, same
    forward solve), which is the control case.
    """

    def __init__(self, parameterisation, bounds, resource, vtip_max=None,
                 polar_cache="sg6043"):
        super().__init__(parameterisation, bounds, resource, polar_cache)
        self.vtip_max = None if vtip_max is None else float(vtip_max)
        self.lambda_opt = self.tsr
        self.tsr_per_bin = [tsr_schedule(v, self.lambda_opt, self.vtip_max)
                            for v in self.speeds]
        # The one list every lambda-dependent quantity in the parent reads.
        self.omega = [lam * v / self.R for lam, v in zip(self.tsr_per_bin, self.speeds)]

    def solve(self, d):
        """The parent's forward solve with the per-bin TSR in place of `self.tsr`."""

        d = np.asarray(d, dtype=float)
        geometry = self.parameterisation.to_geometry(d, polar_cache=self.polar_cache)

        phi, power, a, reynolds, converged = [], [], [], [], True
        for v, lam in zip(self.speeds, self.tsr_per_bin):
            p, result = aerodynamic_power(geometry, v, lam, self.air_density,
                                          self.kinematic_viscosity)
            stations = result["stations"]
            phi.append([s["phi"] for s in stations])
            a.append([s["a"] for s in stations])
            reynolds.append([s["reynolds"] for s in stations])
            power.append(float(p))
            converged = converged and bool(result["converged"])

        if not converged:
            raise RuntimeError(f"forward solve did not converge at d = {d.tolist()}")

        power = np.array(power)
        limited = power > self.p_rated
        return ForwardState(d=d, phi=np.array(phi), power_w=power, limited=limited,
                            a=np.array(a), reynolds=np.array(reynolds), converged=True)


def evaluate_aep(x, parameterisation, resource, design, site, vtip_max,
                 polar_cache="sg6043"):
    """
    AEP [MWh/yr] of design vector `x` under the ceiling, through the
    project's forward chain, with the per-bin record the report needs.

    Bins, midpoint rule, CDF bin masses, the fixed rating and `T = 8766 h`
    are the project's (`objective.power`, `objective.objective`); only the
    per-bin TSR differs from `annual_energy_mwh`.
    """

    edges, mids, _width = wind_speed_bins()
    mass = np.asarray(resource.probability_between(edges[:-1], edges[1:]), dtype=float)
    lam_opt = float(design.design_tsr)
    p_rated = float(design.rated_power_w)
    air = (float(site.air_density), float(site.kinematic_viscosity))
    radius = float(design.radius_m)

    geometry = parameterisation.to_geometry(np.asarray(x, dtype=float), polar_cache=polar_cache)
    rows, total, converged = [], 0.0, True
    for v, m in zip(mids, mass):
        v = float(v)
        lam = tsr_schedule(v, lam_opt, vtip_max)
        p, res = aerodynamic_power(geometry, v, lam, *air)
        converged = converged and bool(res["converged"])
        limited = p > p_rated
        energy = min(p, p_rated) * m * HOURS_PER_YEAR / 1e6
        total += energy
        rows.append({"v_ms": v, "tsr": lam, "rpm": lam * v / radius * 60.0 / (2.0 * math.pi),
                     "mass": float(m), "Cp": float(res["Cp"]), "Ct": float(res["Ct"]),
                     "power_w": float(p), "limited": bool(limited),
                     "energy_mwh": float(energy), "converged": bool(res["converged"])})
    return float(total), rows, converged


class CeilingProblem(ScaledProblem):
    """
    `ScaledProblem` with the ceiling objective and its adjoint.

    Inherited unchanged: the scaling `u <-> d`, the normalisation by
    `|J0|`, the polar-cache envelope constraint (which assumes
    `lambda = 6.5` at every speed and therefore over-estimates `W` under a
    ceiling -- conservative, and the same choice `reoptimise.py` made),
    the counters and the domain-error log.
    """

    def __init__(self, parameterisation, bounds, resource, vtip_max=None,
                 polar_cache="sg6043", margin=0.05):
        super().__init__(parameterisation, bounds, resource, polar_cache, margin)
        self.vtip_max = None if vtip_max is None else float(vtip_max)

    def J(self, u):
        d = self.physical(u)
        try:
            aep, _rows, converged = evaluate_aep(d, self.parameterisation, self.resource,
                                                 self.design, self.site, self.vtip_max,
                                                 self.polar_cache)
        except PolarDomainError:
            self.domain_errors.append([float(x) for x in np.asarray(u)])
            raise
        if not converged:
            raise RuntimeError(f"BEM did not converge at u = {np.asarray(u).tolist()}")
        self.n_fun_evals += 1
        return -aep

    def adjoint_system(self):
        if self._adjoint is None:
            self._adjoint = CeilingBEMSystem(self.parameterisation, self.bounds,
                                             self.resource, self.vtip_max, self.polar_cache)
        return self._adjoint

    def operating_state(self, u):
        """The parent's post-check state, solved on the ceiling schedule."""

        d = self.physical(u)
        geometry = self.parameterisation.to_geometry(d, polar_cache=self.polar_cache)
        lam_opt = float(self.design.design_tsr)
        rho = float(self.site.air_density)
        nu = float(self.site.kinematic_viscosity)

        state = {"speeds": [], "tsr": [], "alpha_deg": [], "reynolds": [], "phi": [],
                 "a": [], "power_w": [], "converged": []}
        for v in self.operating_speeds():
            lam = tsr_schedule(v, lam_opt, self.vtip_max)
            power, result = aerodynamic_power(geometry, v, lam, rho, nu)
            stations = result["stations"]
            state["speeds"].append(v)
            state["tsr"].append(lam)
            state["alpha_deg"].append([math.degrees(s["alpha"]) for s in stations])
            state["reynolds"].append([s["reynolds"] for s in stations])
            state["phi"].append([s["phi"] for s in stations])
            state["a"].append([s["a"] for s in stations])
            state["power_w"].append(float(power))
            state["converged"].append(bool(result["converged"]))
        return state
