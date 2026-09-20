"""
Cross-evaluation of the optimum under a rotor-speed ceiling
(2026-09-19; docs/journal/CROSS_EVALUATION_ROTOR_SPEED_CEILING.md).

Originally: does the blade optimised WITHOUT a rotor-speed ceiling (the
adjoint-driven SLSQP optimum x*, verification/adjoint_optimisation/) still
beat the Schmitz baseline x0 once a ceiling is imposed? No re-optimisation;
pure evaluation of two fixed blades under several operating strategies.

RE-RUN 2026-09-19 under the 300 rpm operating law (B1) and the configured
bounds. The machine ceiling is now a config fact, so x* itself is the
CEILING-CONSTRAINED optimum and the first case below is its own law rather
than a counterfactual. The `none` case survives as the diagnostic it always
was: what that same blade would do with the ceiling removed.

Operating strategy per case, exactly as reoptimise.py (audit section 3.1):

    lambda(V) = 6.5                          no ceiling ("none")
    lambda(V) = min(6.5, V_tip,max / V)      ceiling at V_tip,max
                                             in {62.83 (300 rpm, config),
                                                 60, 55, 50} m/s

Above rated, power is held at the FIXED generator rating from config
(operating.rated_power_w, the same number for every blade and every case).
Bins, masses, air properties, polar cache: the project's. The audit's
ceiling-optimised blades (opt_{60,55,50}_fixed.json) are evaluated too, as
context for how far x* is from the blade a ceiling actually wants.

    python verification/aep_gain_audit/cross_evaluate_xstar.py      # ~1 min

The 62.83 m/s case uses `design.max_tip_speed_ms` from config, never a
hard-coded number, so this file follows the law rather than restating it.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import datetime
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))

from config import load_design_rotor, load_site  # noqa: E402
from design import BladeParameterisation  # noqa: E402
from objective import WeibullResource  # noqa: E402
from objective.objective import HOURS_PER_YEAR  # noqa: E402
from objective.power import aerodynamic_power, wind_speed_bins  # noqa: E402

CEILINGS_MS = (60.0, 55.0, 50.0)
OUT_PATH = os.path.join(_HERE, "cross_evaluate_xstar.json")


def load_vector(path, chord_key, twist_key=None):
    with open(path, encoding="utf-8") as handle:
        d = json.load(handle)
    if twist_key is None:
        return np.array(d[chord_key], dtype=float)
    return np.array(d[chord_key] + d[twist_key], dtype=float)


def main():
    design, site = load_design_rotor(), load_site()
    air = (site.air_density, site.kinematic_viscosity)
    resource = WeibullResource.from_config()
    par = BladeParameterisation()
    lam = float(design.design_tsr)
    p_rated = float(design.rated_power_w)
    radius = float(design.radius_m)

    edges, mids, _ = wind_speed_bins()
    mass = resource.probability_between(edges[:-1], edges[1:])

    blades = {
        "x0": load_vector(os.path.join(REPO, "verification", "baseline", "x0.json"),
                          "chord_control_points_m", "twist_control_points_rad"),
        "x_star": load_vector(os.path.join(REPO, "verification", "adjoint_optimisation", "result.json"),
                              "x_star"),
    }
    context = {}
    for v in (60, 55, 50):
        path = os.path.join(_HERE, f"opt_{v}_fixed.json")
        if os.path.exists(path):
            context[f"opt_{v}"] = load_vector(path, "x_opt")

    def evaluate(x, vtip):
        g = par.to_geometry(x)
        rows, total, converged = [], 0.0, True
        for v, m in zip(mids, mass):
            v = float(v)
            tsr = lam if vtip is None else min(lam, vtip / v)
            p, res = aerodynamic_power(g, v, tsr, *air)
            converged = converged and bool(res["converged"])
            limited = p > p_rated
            energy = min(p, p_rated) * m * HOURS_PER_YEAR / 1e6
            total += energy
            rows.append({"v": v, "tsr": tsr, "rpm": tsr * v / radius * 60 / (2 * math.pi),
                         "Cp": float(res["Cp"]), "Ct": float(res["Ct"]),
                         "power_w": float(p), "limited": bool(limited), "energy_mwh": float(energy)})
        return total, rows, converged

    results = {}
    cases = (("none", None),
             (f"{design.max_rotor_speed_rpm:g} rpm", float(design.max_tip_speed_ms)),
             *((f"{c:g} m/s", c) for c in CEILINGS_MS))
    for label, vtip in cases:
        case = {"vtip_max_ms": vtip}
        for name, x in list(blades.items()) + list(context.items()):
            aep, rows, ok = evaluate(x, vtip)
            case[name] = {"aep_mwh": aep, "all_converged": ok, "bins": rows}
        a0, a1 = case["x0"]["aep_mwh"], case["x_star"]["aep_mwh"]
        case["gain_xstar_over_x0_pct"] = 100 * (a1 / a0 - 1)
        case["delta_xstar_minus_x0_mwh"] = a1 - a0
        for name in context:
            case[f"gain_{name}_over_x0_pct"] = 100 * (case[name]["aep_mwh"] / a0 - 1)
        # where the ceiling bites: the first bin whose lambda is below 6.5
        first = next((r["v"] for r in case["x0"]["bins"] if r["tsr"] < lam - 1e-12), None)
        case["first_bin_below_design_tsr_ms"] = first
        case["v_c_ms"] = None if vtip is None else vtip / lam
        results[label] = case

        print(f"ceiling {label:>4}: AEP x0 {a0:.6f}  x* {a1:.6f}  "
              f"x* vs x0 {case['gain_xstar_over_x0_pct']:+.3f} %"
              + "".join(f"  | {n} {case[f'gain_{n}_over_x0_pct']:+.3f} %" for n in context))
        for r0, r1 in zip(case["x0"]["bins"], case["x_star"]["bins"]):
            print(f"    V={r0['v']:4.1f}  lam={r0['tsr']:.3f}  rpm={r0['rpm']:6.1f}  "
                  f"Cp x0 {r0['Cp']:.4f}  x* {r1['Cp']:.4f}  dCp {100*(r1['Cp']/r0['Cp']-1):+.2f} %  "
                  f"{'capped' if r0['limited'] and r1['limited'] else ''}")

    record = {
        "description": __doc__.split("\n\n")[0].strip(),
        "operating_law": ("lambda(V) = min(design_tsr, V_tip,max / V); the "
                          "'300 rpm' case is the machine ceiling from config "
                          "(max_tip_speed_ms)"),
        "bounds": ("configured bounds (chord_max_m = 0.30 m, "
                   "max_local_solidity = 0.5). x* and x0 satisfy them; the "
                   "audit's ceiling-optimised context blades opt_60/55/50 do "
                   "not -- they were optimised under the retired 0.45 m "
                   "provisional cap, so they are context, not candidates"),
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": "python verification/aep_gain_audit/cross_evaluate_xstar.py",
        "rated_power_w": p_rated,
        "design_tsr": lam,
        "ceilings_ms": [vtip for _, vtip in cases],
        "blades": {k: [float(v) for v in x] for k, x in list(blades.items()) + list(context.items())},
        "sources": {
            "x0": "verification/baseline/x0.json",
            "x_star": "verification/adjoint_optimisation/result.json (x_star)",
            **{k: f"verification/aep_gain_audit/{k}_fixed.json (x_opt)" for k in context},
        },
        "cases": results,
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print("wrote", OUT_PATH)


if __name__ == "__main__":
    main()
