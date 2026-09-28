"""
SLSQP re-optimisation under a modified operating strategy, for
docs/DESIGN-BASIS.md section 2 (an earlier investigation, kept as a record).

A4's optimiser (`verification/fd_optimisation/`) -- central FD at h = 3e-6 in
u, SLSQP with ftol = 1e-8, the polar-cache envelope constraint at margin
0.05, Bounds(0, 1), start at x0 -- on an objective that differs from the
project's in two stated ways:

  * the generator rating may be FIXED at x0's aerodynamic rated power
    (3822.2 W) instead of floating with the design as P_aero(V_rated; d);
  * an optional tip-speed ceiling: below V_c = V_tip,max / lambda the strategy
    is unchanged (lambda = 6.5); above it the rotor holds Omega_max and
    lambda = V_tip,max / V (Region 2.5), until the power cap.

The envelope constraint is reused as-is: it assumes lambda = 6.5 at every
speed, which over-estimates Re under a ceiling, so it is conservative.

    python verification/aep_gain_audit/reoptimise.py <vtip_max_or_none> <fixed|float> <out.json>

e.g.  reoptimise.py none fixed opt_none_fixed.json      (~2.5 min)
      reoptimise.py 55 fixed opt_55_fixed.json          (~3-5 min)

Under provisional bounds (chord_max_m = 0.45 m provisional). Not a project
result: an audit measurement of what the objective would do under a
machine model that hadn't been specified.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import Bounds, minimize

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))

from config import load_design_rotor, load_site  # noqa: E402
from design import BladeParameterisation, DesignBounds  # noqa: E402
from gradients import ScaledProblem  # noqa: E402
from objective import WeibullResource  # noqa: E402
from objective.objective import HOURS_PER_YEAR  # noqa: E402
from objective.power import aerodynamic_power, wind_speed_bins  # noqa: E402

# The bounds this 2026-09-13 audit ran under, kept here verbatim so the script
# stays a reproducible record. They were `tests/test_parameterisation.py::
# PROVISIONAL_BOUNDS` at the time; that set was retired on 2026-09-19 when
# `chord_max_m = 0.30 m` went into config/rotor_design.yaml. Do not read
# 0.45 m as a current bound.
PROVISIONAL_BOUNDS = {
    "chord_min_m": 0.045, "chord_max_m": 0.45,
    "twist_min_rad": math.radians(-2.0), "twist_max_rad": math.radians(35.0),
}

#: x0's aerodynamic rated power under the project's strategy, rounded to
#: 0.1 W (verification/baseline/baseline_reference.json: 3822.19 W).
P_GEN_W = 3822.2
FD_STEP = 3e-6


def main(argv):
    vtip = None if argv[0].lower() == "none" else float(argv[0])
    rating = argv[1]
    if rating not in ("fixed", "float"):
        raise SystemExit("rating must be 'fixed' or 'float'")
    out_path = argv[2]

    design, site = load_design_rotor(), load_site()
    air = (site.air_density, site.kinematic_viscosity)
    resource = WeibullResource.from_config()
    par = BladeParameterisation()
    bounds = DesignBounds(n_chord=par.n_chord, n_twist=par.n_twist, **PROVISIONAL_BOUNDS)
    problem = ScaledProblem(par, bounds, resource, margin=0.05)

    with open(os.path.join(REPO, "verification", "baseline", "x0.json"), encoding="utf-8") as handle:
        d = json.load(handle)
    x0 = np.array(d["chord_control_points_m"] + d["twist_control_points_rad"], dtype=float)
    u0 = bounds.to_scaled(x0)

    edges, mids, _ = wind_speed_bins()
    mass = resource.probability_between(edges[:-1], edges[1:])
    lam, v_rated = design.design_tsr, design.rated_wind_speed_ms

    def tsr(v):
        return lam if vtip is None else min(lam, vtip / v)

    def aep(u):
        g = par.to_geometry(bounds.to_physical(u))
        if rating == "fixed":
            p_rated = P_GEN_W
        else:
            p_rated = aerodynamic_power(g, v_rated, tsr(v_rated), *air)[0]
        total = 0.0
        for v, m in zip(mids, mass):
            p, _ = aerodynamic_power(g, float(v), tsr(float(v)), *air)
            total += min(p, p_rated) * m
        return total * HOURS_PER_YEAR / 1e6

    aep0 = aep(u0)
    n_eval = [1]

    def fun(u):
        n_eval[0] += 1
        return -aep(u) / aep0

    def jac(u):
        grad = np.zeros(len(u))
        for j in range(len(u)):
            up, um = u.copy(), u.copy()
            up[j] += FD_STEP
            um[j] -= FD_STEP
            grad[j] = (fun(up) - fun(um)) / (2 * FD_STEP)
        return grad

    started = time.perf_counter()
    result = minimize(fun, u0, jac=jac, method="SLSQP", bounds=Bounds(0, 1),
                      constraints=[problem.envelope_constraint()],
                      options={"ftol": 1e-8, "maxiter": 200})
    aep_opt = aep(result.x)

    record = {
        "provisional_bounds": "under provisional bounds (chord_max_m = 0.45 m provisional)",
        "vtip_max_ms": vtip,
        "rating": rating,
        "p_gen_w": P_GEN_W if rating == "fixed" else None,
        "aep_x0_mwh": aep0,
        "aep_opt_mwh": aep_opt,
        "gain_pct": 100 * (aep_opt / aep0 - 1),
        "nit": int(result.nit),
        "nfev_total": n_eval[0],
        "status": int(result.status),
        "message": result.message,
        "wall_time_s": time.perf_counter() - started,
        "x_opt": [float(v) for v in bounds.to_physical(result.x)],
        "u_opt": [float(v) for v in result.x],
        "x0": [float(v) for v in x0],
    }
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(json.dumps({k: v for k, v in record.items() if k not in ("x_opt", "u_opt", "x0")}))


if __name__ == "__main__":
    main(sys.argv[1:])
