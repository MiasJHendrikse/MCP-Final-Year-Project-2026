"""
Multi-start of the mass problem at delta = 0.

Starts: `x0` (the production start), `x_c` (the energy optimum,
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
    python verification/mass_optimisation/run_multistart.py --replot

`--replot` redraws `multistart_delta0.png` from the committed
`multistart_delta0.json`; no solve is run and no JSON is written.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import argparse
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


def _start_label(label):
    """A run's tick label: the blade symbols for x0 / x_c, 'start k' otherwise."""

    if label == "x0":
        return r"$\mathbf{x}_0$"
    if label == "x_c":
        return r"$\mathbf{x}_c$"
    return label.replace("_", " ")


def plot_multistart(runs, spread, path):
    """
    The optima the ten starts reached, one panel, written through
    `figstyle.save`: `|u - u_best|_inf` for each start on a log axis, with
    the committed agreement criterion as a dashed line. The best run itself
    is at distance zero, which a log axis cannot draw, so its slot is
    labelled "best" instead of carrying a bar.

    The objective values are not plotted: the ten runs agree to about
    1e-11 in the normalised objective, so bars of the objective are ten
    identical bars and carry no information.

    `runs` is `multistart_delta0.json`'s `runs`; the best run is the
    feasible one with the lowest objective, and the distance is taken from
    its `u_m`. Returns the path.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from plotting import figstyle

    figstyle.apply()

    ok = [r for r in runs if r["status"] == "ok" and r.get("u_m")]
    feasible = [r for r in ok if r.get("feasible")]
    best = min(feasible or ok, key=lambda r: r["mass"]) if ok else None
    distance = []
    for r in ok:
        if best is None:
            distance.append(float("nan"))
        else:
            u = np.array(r["u_m"], dtype=float)
            u_best = np.array(best["u_m"], dtype=float)
            distance.append(float(np.max(np.abs(u - u_best))))
    x = np.arange(len(ok))
    blue, red = figstyle.PALETTE[0], figstyle.PALETTE[1]
    orange = figstyle.PALETTE[3]
    colours = [blue if r.get("feasible") else orange for r in ok]
    labels = [_start_label(r["label"]) for r in ok]

    fig, ax = plt.subplots(figsize=(4.8, 2.6))

    ax.bar(x, distance, color=colours)
    ax.axhline(spread, color=red, lw=1.2, ls="--")
    # The criterion is named on its line, below it at the right-hand end,
    # where the bars (four decades lower) leave the panel empty.
    ax.annotate("agreement criterion", (1.0, spread),
                xycoords=ax.get_yaxis_transform(), xytext=(-6, -4),
                textcoords="offset points", ha="right", va="top",
                fontsize=8, color=red)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel(figstyle.label(
        r"$\|\mathbf{u} - \mathbf{u}_{\mathrm{best}}\|_\infty$", None, "--"))
    ax.set_yscale("log")
    # Fix the lower limit before labelling the best slot, so the label sits
    # just above the axis whatever the bars' range.
    positive = [d for d in distance if d > 0]
    if positive:
        ax.set_ylim(min(positive) / 3.0, spread * 3.0)
    for i, d in enumerate(distance):
        if d == 0.0:
            ax.annotate("best", (i, 0.0), xycoords=("data", "axes fraction"),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7, color=blue)

    if any(r.get("feasible") for r in ok) and any(not r.get("feasible") for r in ok):
        handles = [Line2D([], [], color=blue, lw=6, label="feasible optimum"),
                   Line2D([], [], color=orange, lw=6, label="infeasible optimum")]
        ax.legend(handles=handles, loc="center right", fontsize=7)

    fig.tight_layout()
    figstyle.save(fig, path)
    plt.close(fig)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--replot", action="store_true",
                        help="redraw multistart_delta0.png from the committed "
                             "multistart_delta0.json; no solve is run and no "
                             "JSON is written")
    args = parser.parse_args(argv)

    if args.replot:
        artefact = C.load_json(OUT_PATH)
        plot_multistart(artefact["runs"], artefact["multistart_u_spread_inf"],
                        FIGURE_PATH)
        print(f"redrew {FIGURE_PATH} from {OUT_PATH}; no solve", flush=True)
        return artefact

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

    plot_multistart(runs, spread, FIGURE_PATH)
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
