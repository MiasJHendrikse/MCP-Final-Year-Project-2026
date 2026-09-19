"""
Phase 2, Stage A3: the central-FD step-size study at `x0`.

What it does
------------
At `u0 = bounds.to_scaled(x0)`, for all 10 scaled design variables, evaluates
the central-difference gradient of `fun(u) = J(u) / |J(u0)|` at 15 steps
`h = logspace(-2, -9, 15)` (half-decades) and keeps every value. From those:

  * `h*_j`   the step minimising the local flatness
             `|g_j(h) - g_j(h / sqrt10)| + |g_j(h) - g_j(h sqrt10)|`
  * `h*`     one global step: the grid step nearest the median of `log10 h*_j`
  * `eps_j`  `max(|g_j(h*_j) - g_j(h*_j / sqrt10)|, |g_j(h*_j) - g_j(h*_j sqrt10)|)`
             -- the FD **noise floor** the adjoint (Tier 3) is measured against.

Outputs, next to this script: `sweep.json` (every gradient at every step),
`v_curve.png` (`|g_j(h) - g_j(h*_j)|` vs `h`, log-log, with `h^2` and `1/h`
guide lines), `README.md` (written by hand; the numbers come from the JSON).

Why a plateau, not a V
-----------------------
`J` is C1 but not C2: the Buhl blend at `a = 0.4` is C1, the polar
interpolant has C2 breaks at every 0.5 deg alpha knot and every Reynolds row.
Central differences of a C1 function have an O(h) . delta f'' error whenever
the stencil straddles a break, so the truncation side of the curve flattens
into a plateau instead of falling as `h^2`.

Re-plotting against a reference
--------------------------------
`--reference gradient.json` re-plots the saved sweep (no re-run) with
`|g_j(h) - ref_j|` on the y-axis, where `ref_j` comes from the file. The file
needs one of `gradient_scaled` (units of `fun` per unit `u`, i.e. the same
units as the sweep) or `gradient_mwh_per_u` (MWh/yr per unit `u`; converted
with the sweep's own `J0`). This is how the Tier 3 adjoint is compared to the
whole sweep later without spending another 300 evaluations.

Provisional bounds: `chord_max_m = 0.45 m` is a placeholder pending the
hub-radius / root-attachment decision; the scaling `u = (d - lo)/span` and
therefore every gradient here is stated under provisional bounds.

Run from the repo root (about 65 s):

    python verification/fd_step_size/run_sweep.py
    python verification/fd_step_size/run_sweep.py --reference gradient.json

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import datetime
import json
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from design import BladeParameterisation, DesignBounds  # noqa: E402
from gradients import ScaledProblem, step_size_sweep  # noqa: E402
from gradients.finite_difference import OUT_OF_CACHE  # noqa: E402
from objective import WeibullResource  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
SWEEP_PATH = os.path.join(_HERE, "sweep.json")
FIGURE_PATH = os.path.join(_HERE, "v_curve.png")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")

#: 15 half-decade steps, 1e-2 down to 1e-9. Not larger than 1e-2: at h = 0.1
#: the tip chord control point (u0 ~ 0.054) leaves the polar-cache floor.
STEPS = np.logspace(-2, -9, 15)
SQRT10 = math.sqrt(10.0)


def load_x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"]
                    + artefact["twist_control_points_rad"])


def variable_names(parameterisation):
    return ([f"chord_{j}" for j in range(parameterisation.n_chord)]
            + [f"twist_{j}" for j in range(parameterisation.n_twist)])


def build_problem():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config())


# ---------------------------------------------------------------------------
# the sweep
# ---------------------------------------------------------------------------

def run_sweep(problem, u0):
    """`step_size_sweep` over STEPS, timed, as a list of records."""

    started = time.perf_counter()
    raw = step_size_sweep(problem.fun, u0, STEPS)
    wall = time.perf_counter() - started

    records = []
    for h in STEPS:
        grad = raw[float(h)]
        records.append({
            "h": float(h),
            "gradient": [OUT_OF_CACHE if math.isnan(g) else float(g) for g in grad],
        })
    return records, raw[OUT_OF_CACHE], wall


def gradient_matrix(records):
    """(n_steps, n) array of the sweep, `nan` where out of cache."""

    return np.array([[np.nan if g == OUT_OF_CACHE else g for g in r["gradient"]]
                     for r in records])


def choose_steps(steps, grads):
    """
    `h*_j`, the global `h*`, `eps_j` and the flatness table.

    `steps` descends, so index `k + 1` is `h / sqrt10` and `k - 1` is
    `h sqrt10`. Only interior steps have both neighbours; an `h` with a
    `nan` anywhere in its triple is not eligible.
    """

    n_steps, n = grads.shape
    flatness = np.full((n_steps, n), np.nan)
    for k in range(1, n_steps - 1):
        flatness[k] = (np.abs(grads[k] - grads[k + 1])
                       + np.abs(grads[k] - grads[k - 1]))

    best = np.empty(n, dtype=int)
    for j in range(n):
        column = flatness[:, j]
        if not np.any(np.isfinite(column)):
            raise RuntimeError(f"no eligible step for variable {j}")
        best[j] = int(np.nanargmin(column))

    h_star = steps[best]
    median_log = float(np.median(np.log10(h_star)))
    global_index = int(np.argmin(np.abs(np.log10(steps) - median_log)))

    eps = np.array([max(abs(grads[best[j], j] - grads[best[j] + 1, j]),
                        abs(grads[best[j], j] - grads[best[j] - 1, j]))
                    for j in range(n)])
    eps_global = np.array([max(abs(grads[global_index, j] - grads[global_index + 1, j]),
                               abs(grads[global_index, j] - grads[global_index - 1, j]))
                           for j in range(n)])

    return {
        "best_index": best,
        "h_star": h_star,
        "global_index": global_index,
        "h_star_global": float(steps[global_index]),
        "eps": eps,
        "eps_global": eps_global,
        "flatness": flatness,
    }


# ---------------------------------------------------------------------------
# the figure
# ---------------------------------------------------------------------------

def plot(steps, grads, reference, names, n_chord, h_star, h_star_global,
         reference_label, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Fixed hue order per panel, never cycled; one series per control point.
    colours = ["#1f5fbf", "#d1495b", "#2a9d8f", "#e9a03b", "#6f4bb3"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    blocks = [("chord control points", range(0, n_chord)),
              ("twist control points", range(n_chord, grads.shape[1]))]

    for ax, (title, block) in zip(axes, blocks):
        for colour, j in zip(colours, block):
            diff = np.abs(grads[:, j] - reference[j])
            mask = np.isfinite(diff) & (diff > 0.0)
            ax.plot(steps[mask], diff[mask], "-o", color=colour, lw=1.4, ms=4,
                    label=f"{names[j]}  (h*={h_star[j]:.0e})")

        # Guide lines through the geometric centre of the plotted data.
        finite = [np.abs(grads[:, j] - reference[j]) for j in block]
        finite = np.concatenate([d[np.isfinite(d) & (d > 0)] for d in finite])
        anchor = float(np.exp(np.mean(np.log(finite)))) if finite.size else 1e-6
        mid = steps[len(steps) // 2]
        ax.plot(steps, anchor * (steps / mid) ** 2, "--", color="#888888", lw=1,
                label="$h^2$ (truncation)")
        ax.plot(steps, anchor * (mid / steps), ":", color="#888888", lw=1,
                label="$1/h$ (round-off)")
        ax.axvline(h_star_global, color="#444444", lw=0.8, alpha=0.6)
        ax.text(h_star_global, ax.get_ylim()[0] if False else anchor * 1e-3,
                f" h*={h_star_global:.0e}", fontsize=8, color="#444444")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("step h (scaled variable u)")
        ax.set_title(title, fontsize=11)
        ax.grid(True, which="major", color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    axes[0].set_ylabel(f"|g_j(h) - {reference_label}|   [fun units per unit u]")
    fig.suptitle("Central-FD step-size study at x0 -- " + BOUNDS_LABEL, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# entry points
# ---------------------------------------------------------------------------

def load_reference(path, J0):
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    if "gradient_scaled" in data:
        return np.asarray(data["gradient_scaled"], dtype=float), data.get("label", "reference")
    if "gradient_mwh_per_u" in data:
        return (np.asarray(data["gradient_mwh_per_u"], dtype=float) / abs(J0),
                data.get("label", "reference"))
    raise ValueError(f"{path}: expected 'gradient_scaled' or 'gradient_mwh_per_u'")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--reference", metavar="gradient.json",
                        help="re-plot the saved sweep against this gradient (no re-run)")
    parser.add_argument("--replot", action="store_true",
                        help="re-plot the saved sweep.json without re-running")
    args = parser.parse_args(argv)

    if args.reference or args.replot:
        with open(SWEEP_PATH, encoding="utf-8") as handle:
            saved = json.load(handle)
        steps = np.array([r["h"] for r in saved["sweep"]])
        grads = gradient_matrix(saved["sweep"])
        names = saved["variables"]
        n_chord = saved["n_chord"]
        h_star = np.array(saved["h_star_per_variable"])
        if args.reference:
            reference, label = load_reference(args.reference, saved["J0_mwh_per_yr"])
        else:
            reference = np.array(saved["gradient_at_h_star_scaled"])
            label = "g_j(h*_j)"
        plot(steps, grads, reference, names, n_chord, h_star,
             saved["h_star_global"], label, FIGURE_PATH)
        print(f"re-plotted against {label}: {FIGURE_PATH}")
        return

    problem = build_problem()
    x0 = load_x0()
    u0 = problem.scaled(x0)
    names = variable_names(problem.parameterisation)

    started = time.perf_counter()
    J0 = problem.set_reference(u0)
    j_time = time.perf_counter() - started
    print(f"J(x0) = {J0:.6f} MWh/yr  ({j_time:.3f} s)")

    records, out_of_cache, wall = run_sweep(problem, u0)
    print(f"sweep: {len(STEPS)} steps x {2 * problem.n} evaluations in {wall:.1f} s; "
          f"{len(out_of_cache)} out-of-cache stencils")

    grads = gradient_matrix(records)
    chosen = choose_steps(STEPS, grads)
    g_star = np.array([grads[chosen["best_index"][j], j] for j in range(problem.n)])
    g_global = grads[chosen["global_index"]]

    summary = {
        "description": "Central-FD step-size study of fun(u) = J(u)/|J(u0)| at x0, "
                       "all 10 scaled design variables, 15 half-decade steps.",
        "bounds_label": BOUNDS_LABEL,
        "bounds": problem.bounds.as_record(),
        "command": "python verification/fd_step_size/run_sweep.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "variables": names,
        "n_chord": problem.parameterisation.n_chord,
        "x0": [float(v) for v in x0],
        "u0": [float(v) for v in u0],
        "J0_mwh_per_yr": J0,
        "objective_seconds": j_time,
        "sweep_seconds": wall,
        "n_evaluations": int(len(STEPS) * 2 * problem.n),
        "steps": [float(h) for h in STEPS],
        "sweep": records,
        "out_of_cache": [{"h": h, "variable": names[j]} for h, j in out_of_cache],
        "selection_rule": "h*_j minimises |g(h) - g(h/sqrt10)| + |g(h) - g(h sqrt10)|; "
                          "h* = grid step nearest median(log10 h*_j); "
                          "eps_j = max of those two neighbour differences at h*_j",
        "h_star_per_variable": [float(h) for h in chosen["h_star"]],
        "h_star_global": chosen["h_star_global"],
        "gradient_at_h_star_scaled": [float(g) for g in g_star],
        "gradient_at_h_star_mwh_per_u": [float(g) for g in problem.unscale(g_star)],
        "gradient_at_h_star_global_scaled": [float(g) for g in g_global],
        "gradient_at_h_star_global_mwh_per_u": [float(g) for g in problem.unscale(g_global)],
        "epsilon_scaled": [float(e) for e in chosen["eps"]],
        "epsilon_mwh_per_u": [float(e) for e in problem.unscale(chosen["eps"])],
        "epsilon_at_h_star_global_scaled": [float(e) for e in chosen["eps_global"]],
        "flatness": [[None if math.isnan(v) else float(v) for v in row]
                     for row in chosen["flatness"]],
    }
    with open(SWEEP_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)

    plot(STEPS, grads, g_star, names, problem.parameterisation.n_chord,
         chosen["h_star"], chosen["h_star_global"], "g_j(h*_j)", FIGURE_PATH)

    print(f"\n{'variable':<9} {'h*_j':>8} {'g_j(h*_j) [MWh/yr per u]':>26} {'eps_j':>12} {'eps_j/|g_j|':>12}")
    for j, name in enumerate(names):
        g = summary["gradient_at_h_star_mwh_per_u"][j]
        e = summary["epsilon_mwh_per_u"][j]
        print(f"{name:<9} {chosen['h_star'][j]:>8.0e} {g:>26.6f} {e:>12.3e} "
              f"{e / abs(g) if g else float('inf'):>12.2e}")
    print(f"\nglobal h* = {chosen['h_star_global']:.0e}")
    print(f"wrote {SWEEP_PATH}\n      {FIGURE_PATH}")


if __name__ == "__main__":
    main()
