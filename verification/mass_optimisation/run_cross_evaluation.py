"""
Cross-evaluation of `x0`, `x_c` and `x_m` under other
resources and other rotor-speed ceilings. No optimisation: three fixed
blades, evaluated.

The mass optimum `x_m` holds AEP(x0) to zero slack under ONE resource
(`config/site.yaml`, `k = 1.709, c = 7.274 m/s`) and ONE operating law
(300 rpm ceiling). This script asks how much of that energy parity survives
when either is changed, which is the honest caveat on a "-4.2 % material for
0 % energy" headline:

  * **resource corners** -- `WeibullResource(k, c)` at the central case and
    at the four corners of `c in {7.049, 7.633} m/s` (the log-law cross-check
    band of `verification/wind_resource/README.md`, `z0 = 0.05 .. 0.5 m`) x
    `k in {1.5, 2.0}` (about +/- 0.25 around the extrapolated 1.709; the GWA
    source value at 50 m is 1.87). THIS BAND IS AN ASSUMPTION OF THIS
    ARTEFACT: no sensitivity band has been formally chosen (the README says
    so) and nothing here chooses it. The loads do not depend on the
    resource, so only AEP is evaluated here.
  * **rotor-speed ceilings** -- `V_tip,max in {62.83 (300 rpm, config), 60,
    55, 50} m/s`, the operating law `lambda(V) = min(6.5, V_tip,max / V)` and
    the fixed rating, exactly `verification/aep_gain_audit/cross_evaluate_xstar.py`.
    The rated point moves with the ceiling, so the rated-point root moment and
    tip deflection are evaluated too.

Outputs, next to this script: `cross_evaluation.json`.

Run from the repo root (about one minute):

    python verification/mass_optimisation/run_cross_evaluation.py

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import datetime
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import _common as C  # noqa: E402
from objective import WeibullResource  # noqa: E402
from objective.loads import root_moment, tip_deflection_at  # noqa: E402
from objective.objective import HOURS_PER_YEAR  # noqa: E402
from objective.power import aerodynamic_power, wind_speed_bins  # noqa: E402

RESULT_PATH = os.path.join(_HERE, "result_delta0.json")
OUT_PATH = os.path.join(_HERE, "cross_evaluation.json")

C_BAND_MS = (7.049, 7.633)
K_BAND = (1.5, 2.0)
CEILINGS_MS = (60.0, 55.0, 50.0)


def aep_under(problem, geometry, resource, vtip):
    """AEP [MWh/yr] of one blade under `resource` and ceiling `vtip` (None: no ceiling)."""

    lam = float(problem.design.design_tsr)
    p_rated = float(problem.design.rated_power_w)
    air = (float(problem.site.air_density), float(problem.site.kinematic_viscosity))
    edges, mids, _ = wind_speed_bins()
    mass = resource.probability_between(edges[:-1], edges[1:])
    total, converged = 0.0, True
    for v, m in zip(mids, mass):
        v = float(v)
        tsr = lam if vtip is None else min(lam, vtip / v)
        p, res = aerodynamic_power(geometry, v, tsr, *air)
        converged = converged and bool(res["converged"])
        total += min(p, p_rated) * m * HOURS_PER_YEAR / 1e6
    return float(total), converged


def loads_under(problem, geometry, vtip):
    """Rated-point root moment and tip deflection under ceiling `vtip`."""

    v = float(problem.design.rated_wind_speed_ms)
    lam = float(problem.design.design_tsr)
    tsr = lam if vtip is None else min(lam, vtip / v)
    air = (float(problem.site.air_density), float(problem.site.kinematic_viscosity))
    moment = float(root_moment(geometry, v, tsr, *air)[0])
    deflection = float(tip_deflection_at(geometry, v, tsr, *air)[0])
    return {"tsr_rated": float(tsr), "root_bending_moment_nm": moment,
            "tip_deflection_per_unit_stiffness": deflection}


def main(argv=None):
    problem, u0 = C.prepared_problem()
    artefact = C.load_json(RESULT_PATH)
    blades = {"x0": C.load_x0(), "x_c": C.load_xc(), "x_m": np.array(artefact["x_m"], dtype=float)}
    geometries = {name: problem.parameterisation.to_geometry(x, polar_cache=problem.polar_cache)
                  for name, x in blades.items()}
    site = problem.resource
    vtip_config = float(problem.design.max_tip_speed_ms)

    resource_cases = [("central", site.k, site.c)]
    for k in K_BAND:
        for c in C_BAND_MS:
            resource_cases.append((f"k={k:g}, c={c:g}", k, c))
    resources = {}
    for label, k, c in resource_cases:
        resource = WeibullResource(k=k, c=c)
        case = {"k": float(k), "c_ms": float(c), "mean_wind_speed_ms": float(c * math.gamma(1.0 + 1.0 / k))}
        for name, geometry in geometries.items():
            aep, ok = aep_under(problem, geometry, resource, vtip_config)
            case[name] = {"aep_mwh": aep, "all_converged": ok}
        a0 = case["x0"]["aep_mwh"]
        case["x_c_pct_vs_x0"] = 100.0 * (case["x_c"]["aep_mwh"] / a0 - 1.0)
        case["x_m_pct_vs_x0"] = 100.0 * (case["x_m"]["aep_mwh"] / a0 - 1.0)
        resources[label] = case
        print(f"resource {label:>16}: AEP x0 {a0:.4f}  x_c {case['x_c_pct_vs_x0']:+.4f} %  "
              f"x_m {case['x_m_pct_vs_x0']:+.4f} %", flush=True)

    ceiling_cases = [("none", None), (f"{problem.design.max_rotor_speed_rpm:g} rpm", vtip_config),
                     *((f"{c:g} m/s", c) for c in CEILINGS_MS)]
    ceilings = {}
    for label, vtip in ceiling_cases:
        case = {"vtip_max_ms": vtip, "v_c_ms": None if vtip is None else vtip / float(problem.design.design_tsr)}
        for name, geometry in geometries.items():
            aep, ok = aep_under(problem, geometry, site, vtip)
            case[name] = {"aep_mwh": aep, "all_converged": ok, "rated_point": loads_under(problem, geometry, vtip)}
        a0 = case["x0"]["aep_mwh"]
        m0 = case["x0"]["rated_point"]["root_bending_moment_nm"]
        d0 = case["x0"]["rated_point"]["tip_deflection_per_unit_stiffness"]
        for name in ("x_c", "x_m"):
            case[f"{name}_pct_vs_x0"] = 100.0 * (case[name]["aep_mwh"] / a0 - 1.0)
            case[f"{name}_moment_pct_vs_x0"] = 100.0 * (case[name]["rated_point"]["root_bending_moment_nm"] / m0 - 1.0)
            case[f"{name}_deflection_pct_vs_x0"] = 100.0 * (
                case[name]["rated_point"]["tip_deflection_per_unit_stiffness"] / d0 - 1.0)
        ceilings[label] = case
        print(f"ceiling {label:>8}: AEP x0 {a0:.4f}  x_c {case['x_c_pct_vs_x0']:+.4f} %  "
              f"x_m {case['x_m_pct_vs_x0']:+.4f} %  | rated moment x_m {case['x_m_moment_pct_vs_x0']:+.2f} %  "
              f"deflection x_m {case['x_m_deflection_pct_vs_x0']:+.2f} %", flush=True)

    xm_resource = [case["x_m_pct_vs_x0"] for case in resources.values()]
    xm_ceiling = [case["x_m_pct_vs_x0"] for case in ceilings.values()]
    record = {
        "problem": C.PROBLEM_LABEL,
        "description": __doc__.split("\n\n")[0].strip(),
        "command": "python verification/mass_optimisation/run_cross_evaluation.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "bounds_label": C.BOUNDS_LABEL, "law_label": C.LAW_LABEL,
        "resource_band_is_an_assumption": ("c from the log-law cross-check band of "
                                           "verification/wind_resource/README.md; k +/- ~0.25 "
                                           "around the extrapolated 1.709. Plan 1.3's sensitivity "
                                           "band is not chosen; this does not choose it."),
        "c_band_ms": list(C_BAND_MS), "k_band": list(K_BAND),
        "ceilings_ms": [vtip for _, vtip in ceiling_cases],
        "blades": {k: [float(v) for v in x] for k, x in blades.items()},
        "sources": {"x0": "verification/baseline/x0.json",
                    "x_c": "verification/load_constraint/result_eps0.json (x_c)",
                    "x_m": "verification/mass_optimisation/result_delta0.json (x_m)"},
        "resources": resources,
        "ceilings": ceilings,
        "x_m_aep_pct_vs_x0_range_over_resources": [min(xm_resource), max(xm_resource)],
        "x_m_aep_pct_vs_x0_range_over_ceilings": [min(xm_ceiling), max(xm_ceiling)],
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"\nx_m AEP vs x0: over resources [{min(xm_resource):+.4f}, {max(xm_resource):+.4f}] %; "
          f"over ceilings [{min(xm_ceiling):+.4f}, {max(xm_ceiling):+.4f}] %")
    print(f"wrote {OUT_PATH}")
    return record


if __name__ == "__main__":
    main()
