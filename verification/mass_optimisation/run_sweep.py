"""
The energy-floor sweep and the ablation of the mass
problem.

**Sweep.** For `delta in {0, 0.0025, 0.005, 0.01, 0.02}` the problem of
`run_mass_slsqp.py` is solved twice: **cold** from `x0` and **warm** from the
previous floor's optimum. The two must agree to the multi-start spread in
`u` (`verification/fd_optimisation_multistart/results.json`, 0.0166); where
they do not, the lighter feasible optimum is chosen and the README says why.
The KKT multiplier of the AEP row at each optimum is the exchange rate at
that floor (% material per % energy); it is drawn as the tangent of the
front in `pareto_front.png`, so the front's slope and the adjoint's
multiplier check each other.

**Ablation at delta = 0.** The full row set, then without the deflection
row, without the stress row, without both, without the manufacturing block
(the monotone rows and the 60 mm min-chord rows together -- the 45 mm box
result reappears there), and without the moment cap.
Removing rows can only enlarge the feasible set, so the saving must grow
(or stay) as rows are removed -- the check that says which row is doing the
work, and the table the report's "what makes the number defensible"
paragraph is built from.

Outputs, next to this script: `pareto.json`, `pareto_front.png`,
`ablation.json`, `pareto_ablation.png`.

Run from the repo root (about five minutes):

    python verification/mass_optimisation/run_sweep.py
    python verification/mass_optimisation/run_sweep.py --replot

`--replot` redraws `pareto_front.png` and `pareto_ablation.png` from the
committed `pareto.json` and `ablation.json`; no solve is run and no JSON is
written.

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
from gradients.problem import MASS_PROBLEM_ROWS  # noqa: E402

PARETO_PATH = os.path.join(_HERE, "pareto.json")
FRONT_FIGURE_PATH = os.path.join(_HERE, "pareto_front.png")
ABLATION_PATH = os.path.join(_HERE, "ablation.json")
ABLATION_FIGURE_PATH = os.path.join(_HERE, "pareto_ablation.png")

#: The ablation case labels as the figure reads them (the JSON keeps the
#: longer working labels).
ABLATION_CASE_LABELS = {
    "full": "full set",
    "no deflection": "no deflection row",
    "no stress": "no stress row",
    "no stress, no deflection": "no stress, no deflection",
    "no manufacturing (monotone + min chord)": "no manufacturing rows",
    "no moment cap": "no moment cap",
}

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
    worst = C.worst_slack(s)
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


def plot_pareto_front(parameterisation, x0, front, steps, path):
    """
    The energy-material front, two panels side by side, written through
    `figstyle.save`:

    * "Material-energy front" -- the shell (and the reported solid) material
      saved against the energy given up, with the KKT exchange rate at each
      optimum drawn as its tangent;
    * "Chord along the front" -- the chord distribution of the reference
      blade and of every optimum on the front.

    `front` and `steps` are `pareto.json`'s lists, `x0` the committed
    reference design vector ('x0.json'), `parameterisation` places both.
    Returns the path.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from plotting import figstyle

    figstyle.apply()

    blue, red = figstyle.PALETTE[0], figstyle.PALETTE[1]
    light = figstyle.LIGHT[blue]
    energy = np.array([f["energy_given_up_pct"] for f in front], dtype=float)
    shell = np.array([f["shell_saved_pct"] for f in front], dtype=float)
    solid = np.array([f["solid_saved_pct"] for f in front], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE)

    axes[0].plot(energy, shell, "o-", color=blue, lw=1.4, label="shell material saved")
    axes[0].plot(energy, solid, "s--", color=light, lw=1.2,
                 label="solid proxy saved")
    half = 0.35 * (energy.max() - energy.min()) / max(1, len(energy) - 1)
    for f, x, y in zip(front, energy, shell):
        rate = f["exchange_rate_pct_material_per_pct_energy"]
        if rate is not None:
            axes[0].plot([x - half, x + half], [y - rate * half, y + rate * half],
                         "-", color=red, lw=1.0, alpha=0.8)
    # Each optimum sits on its energy floor (energy given up = delta), so the
    # x-axis already reads delta; a per-point label would repeat it.
    axes[0].plot([], [], "-", color=red, lw=1.0, label="KKT exchange rate")
    axes[0].set_xlabel(figstyle.LABELS["energy_given_up"])
    # Shell and solid savings are both drawn: the axis names material.
    axes[0].set_ylabel(figstyle.label("Material saved", None, "%"))
    figstyle.title(axes[0], "Material-energy front")
    axes[0].legend(loc="lower right", fontsize=7)

    radii = parameterisation.radii
    x0_style = figstyle.BLADES["x0"]
    axes[1].plot(radii, parameterisation.chord(np.asarray(x0, dtype=float)) * 1e3,
                 x0_style["ls"], color=x0_style["color"], lw=1.4,
                 label=x0_style["label"])
    colours = plt.cm.Blues(np.linspace(0.35, 0.95, max(1, len(steps))))
    for step, colour in zip(steps, colours):
        axes[1].plot(radii,
                     parameterisation.chord(np.asarray(step["chosen"]["x_m"],
                                                       dtype=float)) * 1e3,
                     "-", color=colour, lw=1.2,
                     label=f"$\\delta$ = {step['delta'] * 1e2:g} %")
    axes[1].set_xlabel(figstyle.LABELS["radius"])
    axes[1].set_ylabel(figstyle.LABELS["chord"])
    figstyle.title(axes[1], "Chord along the front")
    axes[1].legend(loc="best", fontsize=7)

    fig.tight_layout()
    figstyle.save(fig, path)
    plt.close(fig)
    return path


def plot_pareto_ablation(cases, path):
    """
    The ablation bars at `delta = 0`, single panel, written through
    `figstyle.save`: the shell material saved when each row group of the
    full set is removed. `cases` is `ablation.json`'s `cases`. Returns the
    path.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from plotting import figstyle

    figstyle.apply()

    labels = [ABLATION_CASE_LABELS.get(c["label"], c["label"]) for c in cases]
    saved = [-c["shell_pct_vs_x0"] for c in cases]

    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    ax.bar(range(len(cases)), saved, color="#1f3b8b")
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=7)
    ax.set_ylabel(figstyle.LABELS["material_saved"])
    fig.tight_layout()
    figstyle.save(fig, path)
    plt.close(fig)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--replot", action="store_true",
                        help="redraw pareto_front.png and pareto_ablation.png "
                             "from the committed pareto.json and ablation.json; "
                             "no solve is run and no JSON is written")
    args = parser.parse_args(argv)

    if args.replot:
        pareto = C.load_json(PARETO_PATH)
        ablation = C.load_json(ABLATION_PATH)
        probe = C.build_problem()
        plot_pareto_front(probe.parameterisation, C.load_x0(), pareto["front"],
                          pareto["steps"], FRONT_FIGURE_PATH)
        plot_pareto_ablation(ablation["cases"], ABLATION_FIGURE_PATH)
        print(f"redrew {FRONT_FIGURE_PATH} and {ABLATION_FIGURE_PATH} from "
              f"{PARETO_PATH} and {ABLATION_PATH}; no solve", flush=True)
        return pareto

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
    ablation_record = {
        "problem": C.PROBLEM_LABEL,
        "command": "python verification/mass_optimisation/run_sweep.py",
        "generated": pareto["generated"], "src_commit": pareto["src_commit"],
        "delta": 0.0,
        "cases": ablation,
        "same_saving_tol_pct_points": SAME_SAVING_TOL_PCT,
        "saving_never_shrinks_when_a_row_is_removed": bool(ablation_monotone),
        "removing_both_proxies_saves_at_least_either": bool(nested_monotone),
    }
    with open(PARETO_PATH, "w", encoding="utf-8") as handle:
        json.dump(pareto, handle, indent=1)
    with open(ABLATION_PATH, "w", encoding="utf-8") as handle:
        json.dump(ablation_record, handle, indent=1)

    plot_pareto_front(probe.parameterisation, C.load_x0(), pareto["front"],
                      steps, FRONT_FIGURE_PATH)
    plot_pareto_ablation(ablation, ABLATION_FIGURE_PATH)

    print(f"\nshell saving monotone in delta: {saving_monotone}; all agree to spread {spread:.4f}: "
          f"{pareto['all_agree_to_multistart_spread']}; ablation monotone: {ablation_monotone} "
          f"(nested {nested_monotone})  ({wall:.0f} s)")
    print(f"wrote {PARETO_PATH}\n      {ABLATION_PATH}\n      {FRONT_FIGURE_PATH}\n      "
          f"{ABLATION_FIGURE_PATH}")
    return pareto


if __name__ == "__main__":
    main()
