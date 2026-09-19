"""
Multi-start FD-driven SLSQP: is A4's single-start optimum the optimum?

Re-runs `verification/fd_optimisation/run_fd_slsqp.py`'s optimisation --
same objective scaling (`fun = J / |J(x0)|`), same envelope constraint
(margin 0.05), same `h*`, same SLSQP options (`ftol = 1e-8`, `maxiter = 200`)
-- from starting points sampled uniformly across the provisional bounds
instead of from `x0`.

Starting points
----------------
`u ~ U[0, 1]^10` with a fixed seed, accepted only if the objective can be
evaluated there: the envelope constraint holds (otherwise `CachedPolar`
raises before SLSQP takes a step), every station converges at every
operating point, and `J` is finite. Rejected candidates are recorded with the
reason. The envelope ceiling binds the outboard chord control points hard
(`c_4 <= 0.148 m` at the tip against a 0.45 m box), so "across the bounds"
means across the *feasible* part of the box; that is stated in the README.

Usage (from the repo root)
--------------------------
    python verification/fd_optimisation_multistart/run_multistart.py --sample
        print the candidate list (fast) and write starts.json
    python verification/fd_optimisation_multistart/run_multistart.py --start K
        run start K (one SLSQP run, ~5-10 min) -> start_K.json
    python verification/fd_optimisation_multistart/run_multistart.py --collect
        aggregate start_*.json -> results.json, multistart.png

The starts are independent, so `--start K` for every K can run in parallel
processes.

Provisional bounds: `chord_max_m = 0.45 m` is a placeholder pending the
hub-radius / root-attachment decision; every optimum here is "under
provisional bounds".

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import datetime
import glob
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))
sys.path.insert(0, os.path.join(REPO_ROOT, "verification", "fd_optimisation"))

from gradients.problem import DEFAULT_ENVELOPE_MARGIN  # noqa: E402
from run_fd_slsqp import (  # noqa: E402
    ACTIVE_TOL, ENVELOPE_ACTIVE_TOL, BOUNDS_LABEL,
    build_problem, load_h_star, load_x0, run,
)

A4_RESULT_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation", "result.json")
STARTS_PATH = os.path.join(_HERE, "starts.json")
RESULTS_PATH = os.path.join(_HERE, "results.json")
FIGURE_PATH = os.path.join(_HERE, "multistart.png")

SEED = 20260913
N_STARTS = 8
MAX_CANDIDATES = 20000


def _names(problem):
    p = problem.parameterisation
    return [f"chord_{j}" for j in range(p.n_chord)] + [f"twist_{j}" for j in range(p.n_twist)]


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------

def sample_starts(problem, n_starts=N_STARTS, seed=SEED):
    """
    The first `n_starts` feasible draws of `U[0,1]^10` under `seed`, and every
    rejected draw with its reason. Deterministic, so `--start K` in a separate
    process gets the same point K.
    """

    rng = np.random.default_rng(seed)
    envelope = problem.envelope_constraint()
    accepted, rejected = [], []

    for draw in range(MAX_CANDIDATES):
        u = rng.uniform(0.0, 1.0, problem.n)
        g = envelope["fun"](u)
        if np.any(g < 0.0):
            rejected.append({"draw": draw, "u": u.tolist(),
                             "reason": f"envelope violated, min row {float(g.min()):.4f} m"})
            continue
        try:
            check = problem.alpha_check(u)
        except Exception as error:  # the objective is unevaluable here
            rejected.append({"draw": draw, "u": u.tolist(), "reason": f"{type(error).__name__}: {error}"})
            continue
        if not check["all_converged"]:
            rejected.append({"draw": draw, "u": u.tolist(), "reason": "a station did not converge"})
            continue
        J = problem.J(u)
        if not np.isfinite(J):
            rejected.append({"draw": draw, "u": u.tolist(), "reason": "J not finite"})
            continue
        accepted.append({"start": len(accepted), "draw": draw, "u": u.tolist(),
                         "x": problem.physical(u).tolist(), "J_mwh_per_yr": J,
                         "aep_mwh_per_yr": -J, "alpha_check": check})
        if len(accepted) == n_starts:
            break

    return accepted, rejected


def do_sample(args):
    problem = build_problem(DEFAULT_ENVELOPE_MARGIN)
    accepted, rejected = sample_starts(problem, args.n_starts, args.seed)
    payload = {
        "bounds_label": BOUNDS_LABEL,
        "seed": args.seed,
        "n_starts": len(accepted),
        "acceptance": "envelope satisfied, all stations converged at all 18 operating "
                      "points, J finite",
        "variables": _names(problem),
        "starts": accepted,
        "rejected": rejected,
    }
    with open(STARTS_PATH, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)
    print(f"{len(accepted)} starts accepted from {len(accepted) + len(rejected)} draws "
          f"({len(rejected)} rejected)")
    for s in accepted:
        print(f"  start {s['start']}: AEP {s['aep_mwh_per_yr']:.4f} MWh/yr  "
              f"u = {np.round(s['u'], 3).tolist()}")


# ---------------------------------------------------------------------------
# one start
# ---------------------------------------------------------------------------

def do_start(args):
    with open(STARTS_PATH, encoding="utf-8") as handle:
        starts = json.load(handle)["starts"]
    start = starts[args.start]
    u_start = np.array(start["u"])

    h = load_h_star()
    problem = build_problem(DEFAULT_ENVELOPE_MARGIN)
    u0 = problem.scaled(load_x0())
    problem.set_reference(u0)          # fun = J / |J(x0)|, exactly as in A4
    names = _names(problem)

    result, recorder, wall = run(problem, u_start, h, args.maxiter, args.ftol)

    u_star = np.array(result.x, dtype=float)
    aep0 = -problem.J0
    aep_start = -start["J_mwh_per_yr"]
    aep_star = -problem.J(u_star)
    g_env = problem.envelope_constraint()["fun"](u_star)
    labels = problem.envelope_row_labels()

    with open(A4_RESULT_PATH, encoding="utf-8") as handle:
        a4 = json.load(handle)
    u_a4 = np.array(a4["u_star"])

    record = {
        "bounds_label": BOUNDS_LABEL,
        "start": args.start,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "fd_step": h,
        "envelope_margin": problem.margin,
        "slsqp_options": {"ftol": args.ftol, "maxiter": args.maxiter},
        "u_start": u_start.tolist(),
        "x_start": problem.physical(u_start).tolist(),
        "aep_start_mwh_per_yr": aep_start,
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nit": int(result.nit),
        "nfev": int(result.nfev),
        "njev": int(result.njev),
        "objective_evaluations_total": int(problem.n_fun_evals),
        "wall_time_s": wall,
        "domain_errors": problem.domain_errors,
        "u_star": u_star.tolist(),
        "x_star": problem.physical(u_star).tolist(),
        "aep_optimum_mwh_per_yr": aep_star,
        "delta_aep_pct_vs_x0": 100.0 * (aep_star - aep0) / aep0,
        "delta_aep_pct_vs_start": 100.0 * (aep_star - aep_start) / aep_start,
        "gradient_star_mwh_per_u": problem.unscale(result.jac).tolist(),
        "active_bounds": {"lower": [names[i] for i in range(problem.n) if abs(u_star[i]) < ACTIVE_TOL],
                          "upper": [names[i] for i in range(problem.n) if abs(1 - u_star[i]) < ACTIVE_TOL]},
        "active_envelope_rows": [labels[i] for i in np.flatnonzero(g_env < ENVELOPE_ACTIVE_TOL)],
        "envelope_min_row_value_m": float(g_env.min()),
        "distance_to_a4_u_inf": float(np.max(np.abs(u_star - u_a4))),
        "distance_to_a4_x": {names[i]: float(v) for i, v in
                             enumerate(problem.physical(u_star) - np.array(a4["x_star"]))},
        "post_check_optimum": problem.alpha_check(u_star),
        "iterates": recorder.iterates,
    }
    path = os.path.join(_HERE, f"start_{args.start}.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"start {args.start}: {result.message}; nit={result.nit} AEP {aep_start:.4f} -> "
          f"{aep_star:.6f} ({record['delta_aep_pct_vs_x0']:+.3f} % vs x0); "
          f"|u*-u*_A4|_inf = {record['distance_to_a4_u_inf']:.4f}; wrote {path}")


# ---------------------------------------------------------------------------
# collect
# ---------------------------------------------------------------------------

def do_collect(args):
    with open(A4_RESULT_PATH, encoding="utf-8") as handle:
        a4 = json.load(handle)
    with open(STARTS_PATH, encoding="utf-8") as handle:
        starts_payload = json.load(handle)

    runs = []
    for path in sorted(glob.glob(os.path.join(_HERE, "start_*.json")),
                       key=lambda p: int(os.path.basename(p)[6:-5])):
        with open(path, encoding="utf-8") as handle:
            record = json.load(handle)
        record.pop("iterates")
        record.pop("post_check_optimum")
        runs.append(record)
    if not runs:
        raise SystemExit("no start_*.json found; run --start K first")

    aep_a4 = a4["aep_optimum_mwh_per_yr"]
    best = max(runs, key=lambda r: r["aep_optimum_mwh_per_yr"])
    u_stars = np.array([r["u_star"] for r in runs])
    x_stars = np.array([r["x_star"] for r in runs])

    summary = {
        "bounds_label": BOUNDS_LABEL,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": "python verification/fd_optimisation_multistart/run_multistart.py "
                   "--sample | --start K | --collect",
        "seed": starts_payload["seed"],
        "n_starts": len(runs),
        "n_rejected_draws": len(starts_payload["rejected"]),
        "a4_single_start": {"aep_optimum_mwh_per_yr": aep_a4,
                            "delta_aep_pct_vs_x0": a4["delta_aep_pct_vs_x0"],
                            "nit": a4["nit"], "u_star": a4["u_star"], "x_star": a4["x_star"]},
        "best_start": best["start"],
        "best_aep_mwh_per_yr": best["aep_optimum_mwh_per_yr"],
        "best_minus_a4_mwh_per_yr": best["aep_optimum_mwh_per_yr"] - aep_a4,
        "best_minus_a4_pct_of_a4": 100.0 * (best["aep_optimum_mwh_per_yr"] - aep_a4) / aep_a4,
        "spread_of_optima_aep_mwh_per_yr": float(np.ptp([r["aep_optimum_mwh_per_yr"] for r in runs])),
        "spread_of_optima_u_inf": float(np.max(np.ptp(u_stars, axis=0))),
        "spread_of_optima_u_per_variable": np.ptp(u_stars, axis=0).tolist(),
        "spread_of_optima_x_per_variable": np.ptp(x_stars, axis=0).tolist(),
        "all_converged": all(r["success"] for r in runs),
        "runs": runs,
    }
    with open(RESULTS_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)

    plot(runs, a4, starts_payload, FIGURE_PATH)

    print(f"{'start':>5} {'AEP start':>10} {'nit':>4} {'nfev':>5} {'AEP*':>10} {'dAEP% x0':>9} "
          f"{'|u*-uA4|':>9} {'wall s':>7}  status")
    for r in runs:
        print(f"{r['start']:>5} {r['aep_start_mwh_per_yr']:>10.4f} {r['nit']:>4} {r['nfev']:>5} "
              f"{r['aep_optimum_mwh_per_yr']:>10.6f} {r['delta_aep_pct_vs_x0']:>+9.3f} "
              f"{r['distance_to_a4_u_inf']:>9.4f} {r['wall_time_s']:>7.0f}  {r['message']}")
    print(f"A4 single start: AEP* {aep_a4:.6f} ({a4['delta_aep_pct_vs_x0']:+.3f} %), nit {a4['nit']}")
    print(f"best start {best['start']}: {best['aep_optimum_mwh_per_yr']:.6f} MWh/yr, "
          f"{summary['best_minus_a4_mwh_per_yr']:+.6f} vs A4 "
          f"({summary['best_minus_a4_pct_of_a4']:+.4f} %)")
    print(f"wrote {RESULTS_PATH}\n      {FIGURE_PATH}")


def plot(runs, a4, starts_payload, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    problem = build_problem(DEFAULT_ENVELOPE_MARGIN)
    p = problem.parameterisation
    r = p.radii
    x0 = load_x0()
    x_a4 = np.array(a4["x_star"])

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for ax, (label, unit, get) in zip(axes[:2], [
        ("chord", "m", p.chord),
        ("twist", "deg", lambda d: np.degrees(p.twist(d))),
    ]):
        for s in starts_payload["starts"]:
            ax.plot(r, get(np.array(s["x"])), "-", color="#e9a03b", lw=0.7, alpha=0.5,
                    label="starting points" if s["start"] == 0 else None)
        for run_ in runs:
            ax.plot(r, get(np.array(run_["x_star"])), "-", color="#2a9d8f", lw=1.0, alpha=0.8,
                    label="multi-start optima" if run_ is runs[0] else None)
        ax.plot(r, get(x0), "--", color="#888888", lw=1.4, label="x0 (fitted Schmitz)")
        ax.plot(r, get(x_a4), "-", color="#1f5fbf", lw=1.8, label="A4 single-start optimum")
        ax.set_xlabel("radius r [m]")
        ax.set_ylabel(f"{label} [{unit}]")
        ax.set_title(f"{label} distribution", fontsize=11)
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    ax = axes[2]
    for run_ in runs:
        path_json = os.path.join(_HERE, f"start_{run_['start']}.json")
        with open(path_json, encoding="utf-8") as handle:
            iterates = json.load(handle)["iterates"]
        ax.plot([it["k"] for it in iterates], [-it["J_mwh_per_yr"] for it in iterates],
                "-", color="#2a9d8f", lw=1.0, alpha=0.8)
    ax.axhline(a4["aep_optimum_mwh_per_yr"], color="#1f5fbf", lw=1.4, label="A4 optimum")
    ax.axhline(a4["aep_x0_mwh_per_yr"], color="#888888", lw=1.2, ls="--", label="x0")
    ax.set_xlabel("SLSQP iteration k")
    ax.set_ylabel("AEP [MWh/yr]")
    ax.set_title("convergence of every start", fontsize=11)
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Multi-start FD-driven SLSQP -- " + BOUNDS_LABEL, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample", action="store_true")
    group.add_argument("--start", type=int, metavar="K")
    group.add_argument("--collect", action="store_true")
    parser.add_argument("--n-starts", type=int, default=N_STARTS)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--ftol", type=float, default=1e-8)
    args = parser.parse_args(argv)

    if args.sample:
        do_sample(args)
    elif args.start is not None:
        do_start(args)
    else:
        do_collect(args)


if __name__ == "__main__":
    main()
