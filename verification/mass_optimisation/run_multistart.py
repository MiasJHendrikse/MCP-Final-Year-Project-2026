"""
Phase 5 (2026-09-20) -- multi-start of the mass problem at delta = 0.

Starts: `x0` (the production start), `x_c` (the Phase 4 energy optimum,
AEP slack +0.146 %, below the 60 mm tip floor it predates), and the eight
committed random starts of `verification/fd_optimisation_multistart/starts.json`
(drawn for the energy problem, so most are infeasible for the mass problem
and SLSQP has to find the feasible set first). A start whose first
evaluation fails, or whose run raises from an accepted iterate, is recorded
as failed, never patched.

Agreement criterion: the committed multi-start spread of the energy problem
in `u`, `spread_of_optima_u_inf = 0.0166` (results.json). Every successful
start must land within it of the best (lowest-mass, feasible) optimum, or
the README says which did not and why.

Outputs, next to this script: `multistart_delta0.json`, `multistart_delta0.png`.

Run from the repo root (about five minutes):

    python verification/mass_optimisation/run_multistart.py

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import datetime
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import _common as C  # noqa: E402
from gradients.problem import ACTIVE_TOL, MASS_PROBLEM_ROWS  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402

DELTA = 0.0
FEASIBILITY_TOL = 1e-6
OUT_PATH = os.path.join(_HERE, "multistart_delta0.json")
FIGURE_PATH = os.path.join(_HERE, "multistart_delta0.png")


def starts(problem):
    yield "x0", np.asarray(problem.u0, dtype=float)
    yield "x_c", problem.scaled(C.load_xc())
    for entry in C.load_json(C.STARTS_PATH)["starts"]:
        yield f"start_{entry['start']}", np.array(entry["u"], dtype=float)


def feasible(record, tol=FEASIBILITY_TOL):
    worst = C.worst_slack(record["slacks"])
    return bool(worst >= -tol), worst


def run_start(label, u_start, names):
    problem, u0 = C.prepared_problem()
    started = time.perf_counter()
    try:
        start_record = {"mass": float(problem.mass(u_start)),
                        "aep_over_x0": float(-problem.fun(u_start))}
    except (PolarDomainError, RuntimeError) as error:
        return {"label": label, "u_start": [float(x) for x in u_start], "status": "unevaluable",
                "error": f"{type(error).__name__}: {error}"}
    try:
        result, recorder, wall = C.run_mass_slsqp(problem, u_start, DELTA, MASS_PROBLEM_ROWS,
                                                  verbose=False)
    except (PolarDomainError, RuntimeError) as error:
        return {"label": label, "u_start": [float(x) for x in u_start], "status": "failed",
                "start": start_record, "error": f"{type(error).__name__}: {error}",
                "evaluation_failures": problem.evaluation_failures,
                "wall_time_s": time.perf_counter() - started}
    u_m = np.array(result.x, dtype=float)
    try:
        record = C.blade_record(problem, u_m, DELTA, label, u0=u0)
    except (PolarDomainError, RuntimeError) as error:
        return {"label": label, "u_start": [float(x) for x in u_start], "status": "failed",
                "start": start_record, "error": f"optimum unevaluable: {type(error).__name__}: {error}",
                "u_m": [float(x) for x in u_m], "slsqp_status": int(result.status),
                "message": str(result.message), "wall_time_s": time.perf_counter() - started}
    is_feasible, worst = feasible(record)
    kkt = C.kkt_report(problem, u_m, DELTA, MASS_PROBLEM_ROWS, names)
    print(f"  [{label:>8}] exit {result.status} nit {result.nit:3d}  mass {record['mass_objective']:.6f} "
          f"(shell {record['shell_pct_vs_x0']:+.3f} %)  AEP {record['aep_pct_vs_x0']:+.5f} %  "
          f"feasible {is_feasible} (worst slack {worst:.1e})  failed trials "
          f"{len(recorder.failures)}  {wall:.0f}s", flush=True)
    return {
        "label": label, "u_start": [float(x) for x in u_start], "status": "ok",
        "start": start_record,
        "success": bool(result.success), "slsqp_status": int(result.status),
        "message": str(result.message), "nit": int(result.nit), "nfev": int(result.nfev),
        "u_m": record["u"], "x_m": record["x"], "control_points": record["control_points"],
        "mass": record["mass_objective"],
        "shell_pct_vs_x0": record["shell_pct_vs_x0"], "solid_pct_vs_x0": record["solid_pct_vs_x0"],
        "aep_pct_vs_x0": record["aep_pct_vs_x0"],
        "slacks": record["slacks"], "feasible": is_feasible, "worst_slack": worst,
        "active_rows": kkt["active_rows"], "multipliers": kkt["multipliers"],
        "reynolds_min": record["reynolds_min"],
        "n_evaluation_failures": len(recorder.failures),
        "wall_time_s": float(wall),
    }


def plot(runs, spread, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok = [r for r in runs if r["status"] == "ok"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    labels = [r["label"] for r in ok]
    masses = [r["shell_pct_vs_x0"] for r in ok]
    colours = ["#1f5fbf" if r["feasible"] else "#e08a1e" for r in ok]
    axes[0].bar(range(len(ok)), masses, color=colours)
    axes[0].set_xticks(range(len(ok)))
    axes[0].set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    axes[0].set_ylabel("shell material vs x0 [%]")
    axes[0].set_title("optimum reached from each start (orange: infeasible)", fontsize=11)
    axes[0].grid(True, axis="y", color="#dddddd", lw=0.6)

    du = [r["du_inf_from_best"] for r in ok]
    axes[1].bar(range(len(ok)), du, color=colours)
    axes[1].axhline(spread, color="#d1495b", lw=1.2, ls="--",
                    label=f"energy-problem multi-start spread {spread:.4f}")
    axes[1].set_xticks(range(len(ok)))
    axes[1].set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    axes[1].set_ylabel("|u - u_best|_inf")
    axes[1].set_yscale("log")
    axes[1].set_title("distance in u from the best feasible optimum", fontsize=11)
    axes[1].grid(True, axis="y", color="#dddddd", lw=0.6)
    axes[1].legend(fontsize=8, frameon=False)
    fig.suptitle("Mass problem multi-start, delta = 0 -- " + C.BOUNDS_LABEL, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    spread = float(C.load_json(C.MULTISTART_PATH)["spread_of_optima_u_inf"])
    probe, _u0 = C.prepared_problem()
    names = C.variable_names(probe)

    started = time.perf_counter()
    runs = [run_start(label, u, names) for label, u in starts(probe)]
    wall = time.perf_counter() - started

    ok = [r for r in runs if r["status"] == "ok"]
    feasible_runs = [r for r in ok if r["feasible"]]
    best = min(feasible_runs, key=lambda r: r["mass"]) if feasible_runs else None
    for r in ok:
        r["du_inf_from_best"] = (None if best is None else
                                 float(np.max(np.abs(np.array(r["u_m"]) - np.array(best["u_m"])))))
        r["agrees_to_spread"] = (None if best is None else bool(r["du_inf_from_best"] <= spread))
        r["mass_minus_best"] = None if best is None else float(r["mass"] - best["mass"])
    agreeing = [r["label"] for r in feasible_runs if r["agrees_to_spread"]]
    disagreeing = [r["label"] for r in feasible_runs if not r["agrees_to_spread"]]
    mass_spread = (float(max(r["mass"] for r in feasible_runs) - best["mass"])
                   if feasible_runs else None)

    plot(runs, spread, FIGURE_PATH)
    summary = {
        "problem": C.PROBLEM_LABEL,
        "bounds_label": C.BOUNDS_LABEL,
        "law_label": C.LAW_LABEL,
        "command": "python verification/mass_optimisation/run_multistart.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "delta": DELTA,
        "rows": list(MASS_PROBLEM_ROWS),
        "slsqp_options": {"ftol": C.SLSQP_FTOL, "maxiter": C.SLSQP_MAXITER},
        "feasibility_tol": FEASIBILITY_TOL,
        "agreement_criterion": f"|u - u_best|_inf <= {spread} (energy-problem multi-start spread)",
        "multistart_u_spread_inf": spread,
        "n_starts": len(runs),
        "n_ok": len(ok),
        "n_feasible": len(feasible_runs),
        "n_failed": len([r for r in runs if r["status"] != "ok"]),
        "best_start": None if best is None else best["label"],
        "best_mass": None if best is None else best["mass"],
        "best_shell_pct_vs_x0": None if best is None else best["shell_pct_vs_x0"],
        "spread_of_feasible_optima_mass": mass_spread,
        "spread_of_feasible_optima_u_inf": (float(max(r["du_inf_from_best"] for r in feasible_runs))
                                            if feasible_runs else None),
        "starts_agreeing_to_spread": agreeing,
        "starts_not_agreeing": disagreeing,
        "all_feasible_agree": bool(feasible_runs) and not disagreeing,
        "runs": runs,
        "wall_time_s": wall,
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    print(f"\nbest start {summary['best_start']}: mass {summary['best_mass']} "
          f"({summary['best_shell_pct_vs_x0']:+.3f} %); feasible optima spread in u "
          f"{summary['spread_of_feasible_optima_u_inf']} vs {spread:.4f}; "
          f"agree {agreeing}; disagree {disagreeing}; failed {summary['n_failed']}  ({wall:.0f} s)")
    print(f"wrote {OUT_PATH}\n      {FIGURE_PATH}")
    return summary


if __name__ == "__main__":
    main()
