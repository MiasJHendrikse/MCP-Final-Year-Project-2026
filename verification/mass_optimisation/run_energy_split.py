"""
How much of the annual energy can blade shape move at all?

Review roadmap item 6 (2026-09-26). Under the fixed rating, a bin whose power
reaches `P_rated` contributes `P_rated` whatever the blade does, so part of
AEP is fixed by the rating and not by the shape. This script splits the AEP
of `x0`, `x_c` and `x_m` into the capped (design-invariant) part and the
uncapped (design-dependent) part, and restates the report's two headline
energy figures on the design-dependent basis:

  * the energy optimum's gain over `x0` (+0.1465 % of total AEP);
  * the marginal exchange rate at the reference energy (the AEP-floor KKT
    multiplier, 11.0 % shell material per 1 % of total AEP).

It also restates every KKT multiplier of `x_m` in comparable units (review
roadmap, minor item on multiplier scaling). The rows are written in different
units -- the AEP floor as a fraction of AEP(x0), the stress row in
normalised-KS per m^2, the deflection row as a fraction of D0, the floor rows
in metres of chord -- so the raw multipliers cannot be compared. Each is
multiplied by its row's reference scale, giving percent of shell material per
percent of that row's reference value (per millimetre for the floor rows).

Evaluation only. Inputs: `verification/baseline/x0.json`,
`load_constraint/result_eps0.json` (`x_c`), `result_delta0.json` (`x_m` and
its KKT record). Output: `energy_split.json`, next to this script.

    python verification/mass_optimisation/run_energy_split.py     # ~10 s

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import datetime
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import _common as C  # noqa: E402
from objective.objective import HOURS_PER_YEAR  # noqa: E402
from objective.power import aerodynamic_power, wind_speed_bins  # noqa: E402

RESULT_PATH = os.path.join(_HERE, "result_delta0.json")
OUT_PATH = os.path.join(_HERE, "energy_split.json")


def split(problem, d):
    """(capped MWh/yr, uncapped MWh/yr, per-bin record) for one blade."""

    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    lam = float(problem.design.design_tsr)
    vtip = float(problem.design.max_tip_speed_ms)
    p_rated = float(problem.design.rated_power_w)
    air = (float(problem.site.air_density), float(problem.site.kinematic_viscosity))
    edges, mids, _ = wind_speed_bins()
    mass = problem.resource.probability_between(edges[:-1], edges[1:])
    capped = uncapped = 0.0
    bins = []
    for v, m in zip(mids, mass):
        v = float(v)
        p, res = aerodynamic_power(geometry, v, min(lam, vtip / v), *air)
        is_capped = p >= p_rated
        e = min(p, p_rated) * float(m) * HOURS_PER_YEAR / 1e6
        if is_capped:
            capped += e
        else:
            uncapped += e
        bins.append({"v_ms": v, "mass": float(m), "power_w": float(p), "capped": bool(is_capped),
                     "energy_mwh": float(e), "converged": bool(res["converged"])})
    return capped, uncapped, bins


def rated_speed(problem, d):
    """The wind speed at which a blade first reaches the fixed rating under the operating law."""

    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    lam = float(problem.design.design_tsr)
    vtip = float(problem.design.max_tip_speed_ms)
    p_rated = float(problem.design.rated_power_w)
    air = (float(problem.site.air_density), float(problem.site.kinematic_viscosity))

    def excess(v):
        return aerodynamic_power(geometry, v, min(lam, vtip / v), *air)[0] - p_rated

    v = brentq(excess, 10.5, 12.5, xtol=1e-4)
    return float(v), float(vtip / v)


def main():
    problem, _ = C.prepared_problem()
    artefact = C.load_json(RESULT_PATH)
    blades = {"x0": C.load_x0(), "x_c": C.load_xc(), "x_m": np.array(artefact["x_m"], dtype=float)}
    records = {}
    for name, d in blades.items():
        capped, uncapped, bins = split(problem, d)
        records[name] = {"aep_mwh": capped + uncapped, "capped_mwh": capped,
                         "uncapped_mwh": uncapped, "capped_share": capped / (capped + uncapped),
                         "n_capped_bins": sum(b["capped"] for b in bins), "bins": bins}
        v_rated, tsr_rated = rated_speed(problem, d)
        records[name]["rating_first_reached_at_ms"] = v_rated
        records[name]["tsr_there"] = tsr_rated
        print(f"{name:3s}: AEP {capped + uncapped:.4f} MWh/yr = capped {capped:.4f} "
              f"({100 * capped / (capped + uncapped):.1f} %) + uncapped {uncapped:.4f}; "
              f"{records[name]['n_capped_bins']} capped bins; rating reached at "
              f"{records[name]['rating_first_reached_at_ms']:.3f} m/s")

    x0 = records["x0"]
    share_dep = x0["uncapped_mwh"] / x0["aep_mwh"]
    gain_total = 100.0 * (records["x_c"]["aep_mwh"] / x0["aep_mwh"] - 1.0)
    gain_dep = 100.0 * (records["x_c"]["uncapped_mwh"] - x0["uncapped_mwh"]) / x0["uncapped_mwh"]
    multipliers = artefact["kkt"]["multipliers"]
    aep_multiplier = float(next(v for k, v in multipliers.items() if "aep" in k.lower()))
    rate_dep = aep_multiplier * share_dep
    print(f"design-dependent share of x0's AEP: {100 * share_dep:.1f} %")
    print(f"x_c gain: {gain_total:+.4f} % of total AEP, {gain_dep:+.4f} % of design-dependent energy")
    print(f"exchange rate: {aep_multiplier:.2f} % per 1 % total AEP = {rate_dep:.2f} % per 1 % "
          "of design-dependent energy")
    ks0, c00 = problem._stress_reference()
    stress_scale = ks0 / c00 ** 2
    deflection_scale = float(problem.D0)
    normalised = {
        "aep_floor_pct_per_pct": aep_multiplier * 1.0,
        "stress_pct_per_pct": float(multipliers["stress"]) * stress_scale,
        "deflection_pct_per_pct": float(multipliers["deflection"]) * deflection_scale,
        "min_chord_chord_3_pct_per_mm": 100.0 * float(multipliers["min chord chord_3"]) * 1e-3,
        "min_chord_chord_4_pct_per_mm": 100.0 * float(multipliers["min chord chord_4"]) * 1e-3,
        "monotone_chord_3_chord_4": float(multipliers["monotone chord_3-chord_4"]),
    }
    print("normalised multipliers:", {k: round(v, 4) for k, v in normalised.items()})
    record = {
        "description": __doc__.split("\n\n")[0].strip(),
        "command": "python verification/mass_optimisation/run_energy_split.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "rated_power_w": float(problem.design.rated_power_w),
        "blades": records,
        "x0_design_dependent_share": share_dep,
        "x0_capped_share": 1.0 - share_dep,
        "x_c_gain_pct_total_aep": gain_total,
        "x_c_gain_pct_design_dependent": gain_dep,
        "exchange_rate_pct_per_pct_total_aep": aep_multiplier,
        "exchange_rate_pct_per_pct_design_dependent": rate_dep,
        "kkt_multipliers_raw": {k: float(v) for k, v in multipliers.items()},
        "kkt_row_scales": {"stress_KS0_over_c00_sq_per_m2": float(stress_scale),
                           "deflection_D0": deflection_scale},
        "kkt_multipliers_normalised": normalised,
        "kkt_normalisation_note": ("percent of x0's shell material per percent of the row's "
                                   "reference value (AEP(x0), KS0/c00^2, D0), or per mm of "
                                   "chord for the floor rows; each is the local sensitivity "
                                   "of the optimum to relaxing that row"),
        "note": ("A bin is capped when its aerodynamic power reaches the fixed rating; it then "
                 "contributes P_rated whatever the shape. The design-dependent rate is the KKT "
                 "multiplier times the uncapped share of x0's AEP: 1 % of the design-dependent "
                 "energy is that share of 1 % of total AEP."),
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
