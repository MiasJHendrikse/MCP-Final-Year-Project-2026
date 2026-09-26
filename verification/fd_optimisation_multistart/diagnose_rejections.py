"""
Why the rejected random draws of `starts.json` do not converge.

Review roadmap item 3 (2026-09-26). `run_multistart.py` rejected 34 of its 58
draws with the single reason "a station did not converge". This script
re-solves each of those draws at all seventeen operating points and records,
per failed station, the solver's own `failure` string, the station, the
operating point and the local geometry, so the report can name the failure
mode instead of asserting that the solver converges everywhere.

No optimisation and nothing written anywhere but `rejection_diagnosis.json`,
next to this script. Run from the repo root (about a minute):

    python verification/fd_optimisation_multistart/diagnose_rejections.py

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import collections
import datetime
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "verification", "mass_optimisation"))

import _common as C  # noqa: E402
from bem.rotor import solve_rotor  # noqa: E402
from objective.power import wind_speed_bins  # noqa: E402

STARTS_PATH = os.path.join(_HERE, "starts.json")
OUT_PATH = os.path.join(_HERE, "rejection_diagnosis.json")


def failure_kind(text):
    """The solver's failure string, reduced to its cause."""

    if text.startswith("no sign change in the momentum region"):
        return "no bracket in the momentum region (propeller-brake state)"
    if text.startswith("brentq did not converge"):
        return "brentq iteration limit"
    if text.startswith("residual not reduced"):
        return "residual not reduced"
    return text[:80]


def main():
    problem, _ = C.prepared_problem()
    starts = C.load_json(STARTS_PATH)
    rejected = [r for r in starts["rejected"] if r["reason"] == "a station did not converge"]
    lam = float(problem.design.design_tsr)
    vtip = float(problem.design.max_tip_speed_ms)
    rho = float(problem.site.air_density)
    nu = float(problem.site.kinematic_viscosity)
    _, mids, _ = wind_speed_bins()

    kinds = collections.Counter()
    draws = []
    for entry in rejected:
        u = np.array(entry["u"], dtype=float)
        d = problem.physical(u)
        geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
        failures = []
        for v in mids:
            v = float(v)
            tsr = min(lam, vtip / v)
            result = solve_rotor(geometry, tsr=tsr, v_inf=v, air_density=rho,
                                 kinematic_viscosity=nu)
            for i, station in enumerate(result["stations"]):
                if station["converged"]:
                    continue
                kind = failure_kind(station["failure"])
                kinds[kind] += 1
                failures.append({
                    "v_ms": v, "tsr": float(tsr), "station": i,
                    "r_m": float(geometry.r[i]),
                    "r_over_R": float(geometry.r[i] / geometry.R),
                    "chord_mm": float(1e3 * geometry.chord[i]),
                    "twist_deg": float(math.degrees(geometry.twist[i])),
                    # the inflow angle with no induction, atan(1 / lambda_r):
                    # a section twisted past it meets the wind at a negative
                    # angle of attack before any induction is added
                    "twist_minus_zero_induction_inflow_deg": float(
                        math.degrees(geometry.twist[i])
                        - math.degrees(math.atan(1.0 / (tsr * geometry.r[i] / geometry.R)))),
                    "kind": kind,
                })
        stations = sorted({f["station"] for f in failures})
        draws.append({
            "draw": entry["draw"],
            "control_points": C.control_points_record(problem, d),
            "n_failed_station_points": len(failures),
            "failed_stations": stations,
            "failed_r_over_R": sorted({round(f["r_over_R"], 3) for f in failures}),
            "failed_wind_speeds_ms": sorted({f["v_ms"] for f in failures}),
            "failed_twist_deg_range": ([min(f["twist_deg"] for f in failures),
                                        max(f["twist_deg"] for f in failures)]
                                       if failures else None),
            "kinds": dict(collections.Counter(f["kind"] for f in failures)),
            "failures": failures,
        })
        print(f"draw {entry['draw']:>2}: {len(failures):>3} failed station-points, "
              f"r/R {draws[-1]['failed_r_over_R'][:3]}..., kinds {draws[-1]['kinds']}", flush=True)

    tip_twists = [dr["control_points"]["twist_deg"][-1] for dr in draws]
    accepted_tip_twists = [s["x"][-1] * 180.0 / math.pi for s in starts["starts"]]
    excess = [f["twist_minus_zero_induction_inflow_deg"] for dr in draws for f in dr["failures"]]
    record = {
        "description": __doc__.split("\n\n")[0].strip(),
        "command": "python verification/fd_optimisation_multistart/diagnose_rejections.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "n_rejected_for_convergence": len(rejected),
        "failure_kinds_station_points": dict(kinds),
        "all_failed_stations_outboard_of_r_over_R": min(
            min(dr["failed_r_over_R"]) for dr in draws if dr["failed_r_over_R"]),
        "tip_twist_deg_rejected": [min(tip_twists), max(tip_twists)],
        "tip_twist_deg_accepted": [min(accepted_tip_twists), max(accepted_tip_twists)],
        "twist_minus_zero_induction_inflow_deg_at_failures": {
            "min": min(excess), "median": float(np.median(excess)), "max": max(excess),
            "fraction_positive": float(np.mean(np.array(excess) > 0.0))},
        "draws": draws,
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"\nfailure kinds over all station-points: {dict(kinds)}")
    print(f"innermost failed station r/R = {record['all_failed_stations_outboard_of_r_over_R']:.3f}")
    print(f"tip twist: rejected {record['tip_twist_deg_rejected']}, accepted {record['tip_twist_deg_accepted']}")
    print(f"twist minus zero-induction inflow at failures: {record['twist_minus_zero_induction_inflow_deg_at_failures']}")
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
