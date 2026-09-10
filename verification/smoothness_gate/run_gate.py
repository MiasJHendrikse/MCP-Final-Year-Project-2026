"""
The objective smoothness gate (plan step 1.8).

The readiness test for all gradient work, and described in the Phase 1 brief as
the single most valuable step in Phase 1. For each design variable in turn,
hold the others at `x0` and sweep that one finely; plot `J` against it, and
plot the **first difference**, which is where the defects actually show.

    "Any of these defeats finite differences entirely and makes adjoint
     verification impossible to interpret. They are fixed here, where the
     cause is findable -- not after gradients exist, when the symptom is
     unexplained disagreement between two methods with no way to tell which
     one is wrong."

What is being swept
--------------------
`objective.energy_surrogate` -- the per-bin power summed with **unit weights**
rather than Weibull weights, because the wind resource is still `TODO` (see
`docs/OUTSTANDING-INPUTS.md` section 1). The Weibull weights are a fixed
convex combination over the bins and do not depend on `d`, so the surrogate
exercises the identical d-dependent chain:

    d -> control points -> chord/twist -> RotorGeometry -> BEM -> Cp -> P(V; d)

Every defect this gate looks for lives in that chain. Re-running against the
real objective once `k` and `c` arrive costs about 13 minutes.

The sweep range is provisional
-------------------------------
Plan step 1.8 says "across its feasible range", and the feasible range needs
the manufacturability bounds, which are also `TODO`. The ranges below are
stated here, locally, and deliberately not written into `config/` -- a
provisional bound in the configuration would silently constrain every later
optimisation result, which is what ground rule 3 forbids.

Run from the repo root (about 13 minutes at 300 points):

    python verification/smoothness_gate/run_gate.py [--points N]

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import json
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from design import BladeParameterisation  # noqa: E402
from objective import bin_powers, energy_surrogate  # noqa: E402

#: Provisional sweep half-widths. Chord relative to its own x0 value (control
#: points differ by 4x root to tip, so a single absolute range would be
#: meaningless); twist absolute, because twist control points pass through zero
#: near the tip and a relative range would collapse there.
CHORD_RELATIVE_HALF_WIDTH = 0.40
TWIST_ABSOLUTE_HALF_WIDTH_DEG = 8.0

DEFAULT_POINTS = 300


def load_x0():
    with open(os.path.join(REPO_ROOT, "verification", "baseline", "x0.json"),
              encoding="utf-8") as handle:
        artefact = json.load(handle)
    vector = np.array(artefact["chord_control_points_m"]
                      + artefact["twist_control_points_rad"])
    return vector, artefact


def sweep_range(index, x0, n_chord):
    """(lo, hi) for design variable `index`, and a label."""

    if index < n_chord:
        half = CHORD_RELATIVE_HALF_WIDTH * abs(x0[index])
        return x0[index] - half, x0[index] + half, f"chord CP {index}"
    half = math.radians(TWIST_ABSOLUTE_HALF_WIDTH_DEG)
    return x0[index] - half, x0[index] + half, f"twist CP {index - n_chord}"


def diagnose(values, j):
    """
    The four defects plan step 1.8 names, as numbers.

    Returns a dict. Every metric is scale-free so the same thresholds apply to
    a chord variable and a twist variable.
    """

    j = np.asarray(j, dtype=float)
    first = np.diff(j)
    second = np.diff(first)

    span = float(np.max(j) - np.min(j))
    distinct = len(np.unique(j))

    # Staircasing: a piecewise-constant interpolant somewhere in the chain
    # shows up as repeated identical values, or as a first difference that is
    # zero over stretches.
    flat_steps = int(np.sum(first == 0.0))

    # Discontinuous jumps and kinks: both appear in the difference sequences as
    # a single step far larger than typical. Measured against the MEDIAN, not
    # the mean, so one big spike does not inflate its own reference.
    median_first = float(np.median(np.abs(first))) or 1.0
    median_second = float(np.median(np.abs(second))) or 1.0

    return {
        "range": span,
        "distinct_values": distinct,
        "distinct_fraction": distinct / len(j),
        "flat_first_difference_steps": flat_steps,
        "max_first_difference_ratio": float(np.max(np.abs(first)) / median_first),
        "max_second_difference_ratio": float(np.max(np.abs(second)) / median_second),
        "first_difference_sign_changes": int(np.sum(np.diff(np.sign(first)) != 0)),
        "monotone": bool(np.all(first >= 0.0) or np.all(first <= 0.0)),
    }


def run(points):
    parameterisation = BladeParameterisation()
    x0, artefact = load_x0()
    n_chord = artefact["n_chord"]

    sweeps = []
    started = time.perf_counter()

    for index in range(len(x0)):
        lo, hi, label = sweep_range(index, x0, n_chord)
        grid = np.linspace(lo, hi, points)

        values, converged = [], []
        for value in grid:
            vector = x0.copy()
            vector[index] = value
            powers = bin_powers(vector, parameterisation=parameterisation)
            values.append(float(np.sum(powers["power_w"])))
            converged.append(bool(powers["all_converged"]))

        sweeps.append({
            "index": index,
            "label": label,
            "x0": float(x0[index]),
            "lo": float(lo),
            "hi": float(hi),
            "grid": [float(v) for v in grid],
            "J": values,
            "all_converged": all(converged),
            "n_not_converged": int(len(converged) - sum(converged)),
            "diagnostics": diagnose(grid, values),
        })

        elapsed = time.perf_counter() - started
        print(f"  [{index + 1}/{len(x0)}] {label:14} "
              f"range {min(values):,.0f}-{max(values):,.0f} W  "
              f"converged {all(converged)}  ({elapsed / 60:.1f} min elapsed)")

    return {
        "points_per_sweep": points,
        "objective": "energy_surrogate (unit-weighted per-bin power, watts)",
        "surrogate_reason": (
            "The Weibull parameters are TODO; the weights are a fixed convex "
            "combination independent of d, so the surrogate exercises the "
            "identical d-dependent chain. See docs/OUTSTANDING-INPUTS.md."
        ),
        "sweep_ranges": {
            "chord_relative_half_width": CHORD_RELATIVE_HALF_WIDTH,
            "twist_absolute_half_width_deg": TWIST_ABSOLUTE_HALF_WIDTH_DEG,
            "note": "Provisional: the real feasible range needs the bounds, "
                    "which are TODO. Deliberately not written into config/.",
        },
        "x0_source": "verification/baseline/x0.json",
        "sweeps": sweeps,
    }


def plot(results, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping the figure")
        return

    sweeps = results["sweeps"]
    blue, orange = "#0072B2", "#E69F00"
    n = len(sweeps)

    fig, axes = plt.subplots(2, n, figsize=(2.1 * n, 5.4), sharex="col")
    for column, sweep in enumerate(sweeps):
        grid = np.array(sweep["grid"])
        j = np.array(sweep["J"])
        is_twist = sweep["label"].startswith("twist")
        axis_values = np.degrees(grid) if is_twist else grid * 1000.0

        ax = axes[0, column]
        ax.plot(axis_values, j, color=blue, linewidth=1.0)
        ax.axvline(math.degrees(sweep["x0"]) if is_twist else sweep["x0"] * 1000.0,
                   color="0.6", linestyle=":", linewidth=0.9)
        ax.set_title(sweep["label"], fontsize=8)
        ax.tick_params(labelsize=6)
        if column == 0:
            ax.set_ylabel("$J$ [W]", fontsize=8)
        ax.grid(True, linewidth=0.25, alpha=0.4)

        ax = axes[1, column]
        centres = 0.5 * (axis_values[:-1] + axis_values[1:])
        ax.plot(centres, np.diff(j) / np.diff(axis_values), color=orange,
                linewidth=0.9)
        ax.tick_params(labelsize=6)
        ax.set_xlabel("deg" if is_twist else "mm", fontsize=7)
        if column == 0:
            ax.set_ylabel(r"$\Delta J/\Delta d_i$", fontsize=8)
        ax.grid(True, linewidth=0.25, alpha=0.4)

    fig.suptitle("Objective smoothness gate: $J$ and its first difference, "
                 f"{results['points_per_sweep']} points per design variable",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=DEFAULT_POINTS)
    args = parser.parse_args()

    print(f"Sweeping 10 design variables at {args.points} points each "
          f"({10 * args.points} objective evaluations)...")
    results = run(args.points)

    out = os.path.join(_HERE, "smoothness_gate.json")
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=1)
    print(f"wrote {out}")

    plot(results, os.path.join(_HERE, "smoothness_gate.png"))

    print()
    header = (f"{'variable':>14}{'converged':>11}{'distinct':>10}"
              f"{'flat steps':>12}{'max dJ/med':>12}{'max d2J/med':>13}{'monotone':>10}")
    print(header)
    for sweep in results["sweeps"]:
        d = sweep["diagnostics"]
        print(f"{sweep['label']:>14}{str(sweep['all_converged']):>11}"
              f"{d['distinct_fraction']:>10.3f}{d['flat_first_difference_steps']:>12}"
              f"{d['max_first_difference_ratio']:>12.2f}"
              f"{d['max_second_difference_ratio']:>13.2f}"
              f"{str(d['monotone']):>10}")


if __name__ == "__main__":
    main()
