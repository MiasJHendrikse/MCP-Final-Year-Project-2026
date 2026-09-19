"""
Phase 4, Step 2d -- the eps-sweep Pareto front for the root-moment constraint.

For `eps in {0, 0.02, 0.05, 0.10}` the constrained problem of
`run_constrained_slsqp.py` is solved twice: **cold** from `x0` and **warm**
from the previous step's optimum. The two must agree to the multi-start
study's spread in `u` (`verification/fd_optimisation_multistart/results.json`,
0.0166); where they do not, the better optimum (the higher AEP) is chosen and
the README says why.

The Pareto figure plots AEP cost [%] against rated-point moment reduction [%],
with `x0` and the unconstrained `u*` marked (`u*` sits on the wrong side of the
origin: -0.147 % cost for +0.31 % moment). A second and third panel show where
the load came off the blade -- chord and twist for `x0`, `u*` and each `eps`.

Outputs, next to this script: `pareto.json`, `pareto.png`.

Run from the repo root (about ten minutes: eight SLSQP runs):

    python verification/load_constraint/run_pareto.py

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import datetime
import importlib.util
import json
import os
import subprocess
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

RUNNER_PATH = os.path.join(_HERE, "run_constrained_slsqp.py")
MULTISTART_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation_multistart",
                               "results.json")
ADJOINT_RESULT_PATH = os.path.join(REPO_ROOT, "verification", "adjoint_optimisation",
                                   "result.json")
PARETO_PATH = os.path.join(_HERE, "pareto.json")
FIGURE_PATH = os.path.join(_HERE, "pareto.png")

EPS_STEPS = (0.0, 0.02, 0.05, 0.10)


def load_runner():
    spec = importlib.util.spec_from_file_location("run_constrained_slsqp", RUNNER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def src_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       cwd=REPO_ROOT, text=True).strip()
    except Exception:  # pragma: no cover - reporting only
        return None


def run_point(rc, eps, x_start, label):
    """One constrained solve, returning its record (no files written)."""

    x0 = rc.load_x0()
    problem = rc.build_problem()
    u0 = problem.scaled(x0)
    problem.set_reference(u0)
    problem.load_system()
    started = time.perf_counter()
    result, recorder, wall = rc.run_constrained(problem, u0, eps,
                                                maxiter=200, ftol=1e-8,
                                                x_start=x_start)
    u_c = np.array(result.x, dtype=float)
    x_c = problem.physical(u_c)
    aep_c = -problem.J(u_c)
    aep_star = float(load_json(ADJOINT_RESULT_PATH)["aep_optimum_mwh_per_yr"])
    ks_c = problem.moment_ks(u_c)
    moment_rated = rc.loads_at(problem, x_c)["root_bending_moment_nm"]
    record = {
        "label": label,
        "eps": float(eps),
        "warm_started": x_start is not None,
        "u_c": [float(x) for x in u_c],
        "x_c": [float(x) for x in x_c],
        "aep_mwh_per_yr": float(aep_c),
        "change_pct_vs_x0": float(100.0 * (aep_c + problem.J0) / -problem.J0),
        "aep_cost_pct_vs_unconstrained": float(100.0 * (aep_star - aep_c) / aep_star),
        "KS_c": float(ks_c),
        "KS0": float(problem.KS0),
        "KS_reduction_pct": float(100.0 * (problem.KS0 - ks_c) / problem.KS0),
        "moment_rated_nm": float(moment_rated),
        "moment_reduction_pct": float(100.0 * (problem.m_ref - moment_rated)
                                      / problem.m_ref),
        "success": bool(result.success), "status": int(result.status),
        "nit": int(result.nit), "message": str(result.message),
        "moment_solves": int(problem.n_moment_solves),
        "wall_time_s": float(wall),
    }
    print(f"  [{label}] eps={eps:g} AEP={aep_c:.6f} "
          f"(vs x0 {record['change_pct_vs_x0']:+.3f} %, cost vs u* "
          f"{record['aep_cost_pct_vs_unconstrained']:+.3f} %) KS red "
          f"{record['KS_reduction_pct']:+.3f} % moment red "
          f"{record['moment_reduction_pct']:+.3f} %", flush=True)
    return record


def plot(rc, probe, steps, x0, x_star, star_reduction_pct, star_cost_pct, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

    # (a) Pareto: AEP cost against rated-point moment reduction.
    axes[0].plot([0.0], [0.0], "o", color="#888888", ms=8, label="x0")
    axes[0].plot([star_reduction_pct], [star_cost_pct], "s", color="#d1495b", ms=8,
                 label="unconstrained u*")
    for step in steps:
        chosen = step["chosen"]
        axes[0].plot([chosen["moment_reduction_pct"]], [chosen["change_pct_vs_x0"]],
                     "o-", color="#1f5fbf", ms=6)
        axes[0].annotate(f"eps={step['eps']:g}", (chosen["moment_reduction_pct"],
                                                  chosen["change_pct_vs_x0"]),
                         textcoords="offset points", xytext=(6, 5), fontsize=8)
        if not step["agrees_to_multistart_spread"]:
            axes[0].plot([step["cold"]["moment_reduction_pct"]],
                         [step["cold"]["change_pct_vs_x0"]], "x", color="#e08a1e", ms=7)
            axes[0].plot([step["warm"]["moment_reduction_pct"]],
                         [step["warm"]["change_pct_vs_x0"]], "x", color="#e08a1e", ms=7)
    axes[0].axhline(0.0, color="#bbbbbb", lw=0.8)
    axes[0].axvline(0.0, color="#bbbbbb", lw=0.8)
    axes[0].set_xlabel("rated-point moment reduction [%]  (positive = less load)")
    axes[0].set_ylabel("AEP change vs x0 [%]  (negative = gain; gap to u* = cost)")
    axes[0].set_title("Pareto front", fontsize=11)
    axes[0].grid(True, color="#dddddd", lw=0.6)
    axes[0].legend(fontsize=8, frameon=False)

    # (b) and (c): where the load came off -- chord and twist.
    p = probe.parameterisation
    radii = p.radii
    for ax, (field, unit, get, block) in zip(axes[1:], [
        ("chord", "m", p.chord, slice(0, p.n_chord)),
        ("twist", "deg", lambda d: np.degrees(p.twist(d)), slice(p.n_chord, None)),
    ]):
        ax.plot(radii, get(x0), "-", color="#888888", lw=1.6, label="x0")
        ax.plot(radii, get(x_star), "--", color="#d1495b", lw=1.6,
                label="unconstrained u*")
        for step in steps:
            ax.plot(radii, get(np.array(step["chosen"]["x_c"], dtype=float)), "-",
                    lw=1.6, label=f"eps={step['eps']:g}")
        ax.set_xlabel("radius r [m]")
        ax.set_ylabel(f"{field} [{unit}]")
        ax.set_title(f"{field} distribution", fontsize=11)
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Phase 4 root-moment constraint: AEP cost vs load reduction -- "
                 + rc.BOUNDS_LABEL, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    rc = load_runner()
    x0 = rc.load_x0()
    u_star = np.array(load_json(ADJOINT_RESULT_PATH)["u_star"], dtype=float)
    spread = float(load_json(MULTISTART_PATH)["spread_of_optima_u_inf"])
    unconstrained = load_json(ADJOINT_RESULT_PATH)
    probe = rc.build_problem()
    x_star = probe.physical(u_star)
    probe_moment0 = rc.loads_at(probe, x0)["root_bending_moment_nm"]
    probe_moment_star = rc.loads_at(probe, x_star)["root_bending_moment_nm"]
    star_reduction_pct = 100.0 * (probe_moment0 - probe_moment_star) / probe_moment0
    star_cost_pct = -float(unconstrained["delta_aep_pct_vs_x0"])

    steps = []
    previous = None
    for eps in EPS_STEPS:
        print(f"\n== eps = {eps:g} ==", flush=True)
        cold = run_point(rc, eps, None, "cold")
        warm = run_point(rc, eps, previous, "warm")
        du_inf = float(np.max(np.abs(np.array(cold["u_c"]) - np.array(warm["u_c"]))))
        agrees = bool(du_inf <= spread)
        chosen = cold if cold["aep_mwh_per_yr"] >= warm["aep_mwh_per_yr"] else warm
        reason = None
        if not agrees:
            reason = (f"cold and warm optima differ by {du_inf:.3e} in u, above the "
                      f"multi-start spread {spread:.4f}; the higher-AEP optimum "
                      f"({chosen['label']}) is kept")
        steps.append({
            "eps": float(eps),
            "cold": cold, "warm": warm,
            "du_inf": du_inf, "agrees_to_multistart_spread": agrees,
            "chosen": chosen, "chosen_label": chosen["label"], "reason": reason,
        })
        previous = np.array(chosen["u_c"], dtype=float)

    # Monotone check: the AEP cost relative to the unconstrained optimum must
    # not fall as eps rises (section 3.1). The change vs x0 turns negative at
    # eps = 0.10, so it is not the right axis for the monotonicity statement.
    costs = [s["chosen"]["aep_cost_pct_vs_unconstrained"] for s in steps]
    monotone = all(costs[k] <= costs[k + 1] + 1e-9 for k in range(len(costs) - 1))

    plot(rc, probe, steps, x0, x_star, star_reduction_pct, star_cost_pct, FIGURE_PATH)

    summary = {
        "description": "Phase 4 eps-sweep: cold and warm constrained SLSQP at "
                       "eps in {0, 0.02, 0.05, 0.10}.",
        "bounds_label": rc.BOUNDS_LABEL,
        "law_label": rc.LAW_LABEL,
        "command": "python verification/load_constraint/run_pareto.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": src_commit(),
        "multistart_u_spread_inf": spread,
        "eps_steps": steps,
        "cost_pct_monotone": bool(monotone),
        "all_agree_to_multistart_spread": bool(all(s["agrees_to_multistart_spread"]
                                                   for s in steps)),
        "x0": [float(x) for x in x0],
        "u_star_unconstrained": [float(x) for x in u_star],
    }
    with open(PARETO_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)

    print(f"\ncost monotone in eps: {monotone}; all agree to spread "
          f"{spread:.4f}: {summary['all_agree_to_multistart_spread']}")
    print(f"wrote {PARETO_PATH}\n      {FIGURE_PATH}")
    return summary


if __name__ == "__main__":
    main()
