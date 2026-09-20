"""
Phase 5 (2026-09-20) -- the energy-floor sweep and the ablation of the mass
problem.

**Sweep.** For `delta in {0, 0.0025, 0.005, 0.01, 0.02}` the problem of
`run_mass_slsqp.py` is solved twice: **cold** from `x0` and **warm** from the
previous floor's optimum. The two must agree to the multi-start spread in
`u` (`verification/fd_optimisation_multistart/results.json`, 0.0166); where
they do not, the lighter feasible optimum is chosen and the README says why.
The KKT multiplier of the AEP row at each optimum is the exchange rate at
that floor (% material per % energy); it is drawn as the tangent of the
front in `pareto.png`, so the front's slope and the adjoint's multiplier
check each other.

**Ablation at delta = 0.** The full row set, then without the deflection
row, without the stress row, without both, without the manufacturing block
(the monotone rows and the 60 mm min-chord rows together -- the 45 mm box
result reappears there), and without the moment cap.
Removing rows can only enlarge the feasible set, so the saving must grow
(or stay) as rows are removed -- the check that says which row is doing the
work, and the table the report's "what makes the number defensible"
paragraph is built from.

Outputs, next to this script: `pareto.json`, `pareto.png`, `ablation.json`.

Run from the repo root (about five minutes):

    python verification/mass_optimisation/run_sweep.py

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
from gradients.problem import MASS_PROBLEM_ROWS  # noqa: E402

PARETO_PATH = os.path.join(_HERE, "pareto.json")
FIGURE_PATH = os.path.join(_HERE, "pareto.png")
ABLATION_PATH = os.path.join(_HERE, "ablation.json")

ABLATIONS = (
    ("full", MASS_PROBLEM_ROWS),
    ("no deflection", tuple(r for r in MASS_PROBLEM_ROWS if r != "deflection")),
    ("no stress", tuple(r for r in MASS_PROBLEM_ROWS if r != "stress")),
    ("no stress, no deflection", tuple(r for r in MASS_PROBLEM_ROWS
                                       if r not in ("stress", "deflection"))),
    ("no manufacturing (monotone + min chord)", tuple(r for r in MASS_PROBLEM_ROWS if r != "manufacturing")),
    ("no moment cap", tuple(r for r in MASS_PROBLEM_ROWS if r != "moment")),
)
FEASIBILITY_TOL = 1e-6
#: Two optima of the same problem agree in the objective to SLSQP's `ftol`
#: (1e-10 on the mass fraction, 1e-8 in %); the monotonicity statements use
#: this looser figure in % points so a solver-tolerance tie is a tie.
SAME_SAVING_TOL_PCT = 1e-6


def solve(delta, include, u_start, label, names):
    """One solve; returns the record `run_mass_slsqp.py` would (no files)."""

    problem, u0 = C.prepared_problem()
    start = np.asarray(problem.u0, dtype=float) if u_start is None else np.asarray(u_start)
    result, recorder, wall = C.run_mass_slsqp(problem, start, delta, include, verbose=False)
    u_m = np.array(result.x, dtype=float)
    record = C.blade_record(problem, u_m, delta, label, u0=u0)
    kkt = C.kkt_report(problem, u_m, delta, include, names)
    s = record["slacks"]
    worst = min(s["aep_floor"]["slack"], s["moment"]["slack"], s["stress"]["slack"],
                s["deflection"]["slack"], s["manufacturing"]["slack_min"])
    # Rows outside `include` are reported (their slack at this optimum says
    # what removing them bought) but do not count towards feasibility.
    counted = [s["aep_floor"]["slack"]]
    if "moment" in include:
        counted.append(s["moment"]["slack"])
    if "stress" in include:
        counted.append(s["stress"]["slack"])
    if "deflection" in include:
        counted.append(s["deflection"]["slack"])
    if "manufacturing" in include:
        counted.append(s["manufacturing"]["slack_min"])
    out = {
        "label": label, "delta": float(delta), "rows": list(include),
        "warm_started": u_start is not None,
        "success": bool(result.success), "status": int(result.status),
        "message": str(result.message), "nit": int(result.nit),
        "u_m": record["u"], "x_m": record["x"], "control_points": record["control_points"],
        "mass": record["mass_objective"],
        "shell_pct_vs_x0": record["shell_pct_vs_x0"], "solid_pct_vs_x0": record["solid_pct_vs_x0"],
        "aep_mwh_per_yr": record["aep_mwh_per_yr"], "aep_pct_vs_x0": record["aep_pct_vs_x0"],
        "energy_given_up_pct": -record["aep_pct_vs_x0"],
        "KS_over_KS0": record["KS_over_KS0"],
        "stress_ratio": record["stress_ratio"], "deflection_ratio": record["deflection_ratio"],
        "root_chord_m": record["root_chord_m"],
        "slacks": s, "worst_slack_all_rows": float(worst),
        "feasible_for_its_rows": bool(min(counted) >= -FEASIBILITY_TOL),
        "active_rows": kkt["active_rows"], "multipliers": kkt["multipliers"],
        "exchange_rate_pct_material_per_pct_energy": kkt["exchange_rate_pct_material_per_pct_energy"],
        "kkt_residual": kkt["residual_norm_after_projection"],
        "reynolds_min": record["reynolds_min"],
        "post_check_uncapped_within": (record["post_check"]["uncapped"]["within"]
                                       if record["post_check"]["uncapped"] else None),
        "n_evaluation_failures": len(recorder.failures),
        "wall_time_s": float(wall),
    }
    print(f"  [{label}] delta={delta:g} exit {result.status} nit {result.nit:3d}  "
          f"shell {out['shell_pct_vs_x0']:+.3f} %  solid {out['solid_pct_vs_x0']:+.3f} %  "
          f"AEP {out['aep_pct_vs_x0']:+.4f} %  stress {out['stress_ratio']:.4f}  "
          f"defl {out['deflection_ratio']:.4f}  KS/KS0 {out['KS_over_KS0']:.4f}  "
          f"rate {out['exchange_rate_pct_material_per_pct_energy']}  "
          f"active {kkt['active_rows']}  {wall:.0f}s", flush=True)
    return out


def plot(steps, ablation, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    # (a) the front: material saved vs energy given up, with the exchange
    # rate at each optimum drawn as a tangent: in these axes (energy given
    # up, material saved) its slope is the AEP row's multiplier itself.
    xs = [s["chosen"]["energy_given_up_pct"] for s in steps]
    ys = [-s["chosen"]["shell_pct_vs_x0"] for s in steps]
    axes[0].plot(xs, ys, "o-", color="#1f5fbf", ms=6, lw=1.6, label="shell material saved")
    axes[0].plot(xs, [-s["chosen"]["solid_pct_vs_x0"] for s in steps], "s--", color="#5aa9e6",
                 ms=5, lw=1.2, label="solid proxy saved (reported)")
    half = 0.35 * (max(xs) - min(xs)) / max(1, len(xs) - 1)
    for s, x, y in zip(steps, xs, ys):
        rate = s["chosen"]["exchange_rate_pct_material_per_pct_energy"]
        if rate is not None:
            axes[0].plot([x - half, x + half], [y - rate * half, y + rate * half], "-",
                         color="#d1495b", lw=1.0, alpha=0.8)
        axes[0].annotate(f"delta={s['delta']:g}", (x, y), textcoords="offset points",
                         xytext=(6, -12), fontsize=8)
        if not s["agrees_to_multistart_spread"]:
            for r in (s["cold"], s["warm"]):
                axes[0].plot([r["energy_given_up_pct"]], [-r["shell_pct_vs_x0"]], "x",
                             color="#e08a1e", ms=7)
    axes[0].plot([], [], "-", color="#d1495b", lw=1.0, label="KKT exchange rate (tangent)")
    axes[0].set_xlabel("energy given up vs x0 [%]")
    axes[0].set_ylabel("material saved vs x0 [%]")
    axes[0].set_title("material-energy front", fontsize=11)
    axes[0].grid(True, color="#dddddd", lw=0.6)
    axes[0].legend(fontsize=8, frameon=False)

    # (b) ablation at delta = 0.
    labels = [a["label"] for a in ablation]
    saved = [-a["shell_pct_vs_x0"] for a in ablation]
    axes[1].barh(range(len(ablation)), saved, color="#1f5fbf")
    axes[1].set_yticks(range(len(ablation)))
    axes[1].set_yticklabels(labels, fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("shell material saved at delta = 0 [%]")
    axes[1].set_title("ablation: rows removed one at a time", fontsize=11)
    axes[1].grid(True, axis="x", color="#dddddd", lw=0.6)

    # (c) chord along the front.
    probe, _ = C.prepared_problem()
    p = probe.parameterisation
    x0 = C.load_x0()
    axes[2].plot(p.radii, p.chord(x0), "-", color="#888888", lw=1.6, label="x0")
    for s in steps:
        axes[2].plot(p.radii, p.chord(np.array(s["chosen"]["x_m"])), "-", lw=1.4,
                     label=f"delta={s['delta']:g}")
    axes[2].set_xlabel("radius r [m]")
    axes[2].set_ylabel("chord [m]")
    axes[2].set_title("chord along the front", fontsize=11)
    axes[2].grid(True, color="#dddddd", lw=0.6)
    axes[2].legend(fontsize=8, frameon=False)

    fig.suptitle("Mass problem: energy-floor sweep and ablation -- " + C.BOUNDS_LABEL, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    spread = float(C.load_json(C.MULTISTART_PATH)["spread_of_optima_u_inf"])
    probe, _ = C.prepared_problem()
    names = C.variable_names(probe)
    started = time.perf_counter()

    steps = []
    previous = None
    for delta in C.DELTAS:
        print(f"\n== delta = {delta:g} ==", flush=True)
        cold = solve(delta, MASS_PROBLEM_ROWS, None, "cold", names)
        if previous is None:
            warm, du_inf, agrees, chosen, reason = None, 0.0, True, cold, None
        else:
            warm = solve(delta, MASS_PROBLEM_ROWS, previous, "warm", names)
            du_inf = float(np.max(np.abs(np.array(cold["u_m"]) - np.array(warm["u_m"]))))
            agrees = bool(du_inf <= spread)
            candidates = [r for r in (cold, warm) if r["feasible_for_its_rows"]] or [cold, warm]
            chosen = min(candidates, key=lambda r: r["mass"])
            reason = None if agrees else (
                f"cold and warm optima differ by {du_inf:.3e} in u, above the multi-start "
                f"spread {spread:.4f}; the lighter feasible optimum ({chosen['label']}) is kept")
        steps.append({"delta": float(delta), "cold": cold, "warm": warm, "du_inf": du_inf,
                      "agrees_to_multistart_spread": agrees, "chosen": chosen,
                      "chosen_label": chosen["label"], "reason": reason})
        previous = np.array(chosen["u_m"], dtype=float)

    savings = [-s["chosen"]["shell_pct_vs_x0"] for s in steps]
    saving_monotone = all(savings[k] <= savings[k + 1] + SAME_SAVING_TOL_PCT
                          for k in range(len(savings) - 1))
    # Secant slopes between neighbouring floors, to set against the KKT
    # exchange rates at the two ends of each segment.
    secants = []
    for a, b in zip(steps[:-1], steps[1:]):
        dx = b["chosen"]["energy_given_up_pct"] - a["chosen"]["energy_given_up_pct"]
        dy = -b["chosen"]["shell_pct_vs_x0"] + a["chosen"]["shell_pct_vs_x0"]
        secants.append({"from_delta": a["delta"], "to_delta": b["delta"],
                        "secant_pct_material_per_pct_energy": float(dy / dx) if dx else None,
                        "rate_at_from": a["chosen"]["exchange_rate_pct_material_per_pct_energy"],
                        "rate_at_to": b["chosen"]["exchange_rate_pct_material_per_pct_energy"]})

    print("\n== ablation at delta = 0 ==", flush=True)
    ablation = [solve(0.0, include, None, label, names) for label, include in ABLATIONS]
    full = ablation[0]["shell_pct_vs_x0"]
    for a in ablation:
        a["saving_gain_vs_full_pct_points"] = float(full - a["shell_pct_vs_x0"])
    ablation_monotone = all(a["shell_pct_vs_x0"] <= full + SAME_SAVING_TOL_PCT
                            for a in ablation[1:])
    both = next(a for a in ablation if a["label"] == "no stress, no deflection")
    nested_monotone = all(a["shell_pct_vs_x0"] >= both["shell_pct_vs_x0"] - SAME_SAVING_TOL_PCT
                          for a in ablation if a["label"] in ("no deflection", "no stress"))

    wall = time.perf_counter() - started
    plot(steps, ablation, FIGURE_PATH)

    pareto = {
        "problem": C.PROBLEM_LABEL,
        "bounds_label": C.BOUNDS_LABEL,
        "law_label": C.LAW_LABEL,
        "command": "python verification/mass_optimisation/run_sweep.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "deltas": list(C.DELTAS),
        "slsqp_options": {"ftol": C.SLSQP_FTOL, "maxiter": C.SLSQP_MAXITER},
        "multistart_u_spread_inf": spread,
        "steps": steps,
        "front": [{"delta": s["delta"],
                   "energy_given_up_pct": s["chosen"]["energy_given_up_pct"],
                   "shell_saved_pct": -s["chosen"]["shell_pct_vs_x0"],
                   "solid_saved_pct": -s["chosen"]["solid_pct_vs_x0"],
                   "exchange_rate_pct_material_per_pct_energy":
                       s["chosen"]["exchange_rate_pct_material_per_pct_energy"],
                   "active_rows": s["chosen"]["active_rows"]} for s in steps],
        "secants": secants,
        "shell_saving_monotone_in_delta": bool(saving_monotone),
        "all_agree_to_multistart_spread": bool(all(s["agrees_to_multistart_spread"] for s in steps)),
        "wall_time_s": wall,
    }
    with open(PARETO_PATH, "w", encoding="utf-8") as handle:
        json.dump(pareto, handle, indent=1)
    with open(ABLATION_PATH, "w", encoding="utf-8") as handle:
        json.dump({
            "problem": C.PROBLEM_LABEL,
            "command": "python verification/mass_optimisation/run_sweep.py",
            "generated": pareto["generated"], "src_commit": pareto["src_commit"],
            "delta": 0.0,
            "cases": ablation,
            "same_saving_tol_pct_points": SAME_SAVING_TOL_PCT,
            "saving_never_shrinks_when_a_row_is_removed": bool(ablation_monotone),
            "removing_both_proxies_saves_at_least_either": bool(nested_monotone),
        }, handle, indent=1)

    print(f"\nshell saving monotone in delta: {saving_monotone}; all agree to spread {spread:.4f}: "
          f"{pareto['all_agree_to_multistart_spread']}; ablation monotone: {ablation_monotone} "
          f"(nested {nested_monotone})  ({wall:.0f} s)")
    print(f"wrote {PARETO_PATH}\n      {ABLATION_PATH}\n      {FIGURE_PATH}")
    return pareto


if __name__ == "__main__":
    main()
