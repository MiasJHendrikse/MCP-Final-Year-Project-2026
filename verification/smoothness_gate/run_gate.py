"""
The objective smoothness gate.

The readiness test for all gradient work. For each design variable in turn,
hold the others at `x0` and sweep that one finely; plot `J` against it, and
plot the **first difference**, which is where the defects actually show.

    "Any of these defeats finite differences entirely and makes adjoint
     verification impossible to interpret. They are fixed here, where the
     cause is findable -- not after gradients exist, when the symptom is
     unexplained disagreement between two methods with no way to tell which
     one is wrong."

What is being swept
-------------------
`objective.annual_energy_mwh` -- the committed objective: the Weibull-weighted
annual energy over the 17 wind-speed bins under the committed operating law
(fixed tip-speed-ratio tracking with the 300 rpm rotor-speed ceiling and the
configured rating). This is the functional the gradients differentiate.

The probe of the surrogate form
-------------------------------
An earlier pass of this gate swept `objective.energy_surrogate`, the same
per-bin powers summed with **unit weights**, while the wind resource was
unresolved. The argument for the substitution was that the Weibull weights
are a fixed convex combination over the bins which does not depend on `d`, so
both forms exercise the identical d-dependent chain --

    d -> control points -> chord/twist -> RotorGeometry -> BEM -> Cp -> P(V; d)

-- and a positive fixed weighting can neither create nor remove a
discontinuity in it. Every defect this gate looks for lives in that chain.
That argument was tested rather than merely asserted, sweep by sweep: every
sweep agreed -- 3000/3000 converged, `distinct = 1.000` throughout, zero flat
first-difference steps, and the difference ratios tracking within a few
percent with the same variables elevated. Had the verdict changed, the
weighting would have been coupled to `d` somewhere it should not be, and that
is a defect worth finding before any gradient work.

The sweep range
---------------
The stated range is the design box of `config/rotor_design.yaml`
(`parameterisation.bounds`, read with `DesignBounds.from_config`): chord
0.045-0.30 m, twist -2 to 35 degrees. Sweeping the whole box is not possible:
at the low-chord end a station Reynolds number falls below the polar cache's
40 000 floor and the objective is **undefined** there (`PolarDomainError`), and
at the high-twist end of the twist sweeps the BEM solver does not converge at
every grid point. Each variable is therefore swept across the largest interval
containing its `x0` value on which the objective is **defined and converged**;
the two extents are found by bisection and recorded in the JSON beside the box
together with the fraction of the box they cover. This is a property of the
model, not a choice: the polar cache and the solver's convergence bound the
reachable design space.

Run from the repo root (about 15 minutes at 300 points):

    python verification/smoothness_gate/run_gate.py [--points N]

The figure
----------
Two files, each two graphs side by side. `smoothness_gate_objective.png`
shows `J` along each chord sweep and along each twist sweep, each curve
normalised to its own range so that ten curves of different magnitude can be
read on one axis. `smoothness_gate_difference.png` shows the first difference,
the same split, which is where the C2 defects appear as single elevated
steps at the tip control points.

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

from design import BladeParameterisation, DesignBounds  # noqa: E402
from objective import bin_powers  # noqa: E402
from objective.objective import HOURS_PER_YEAR  # noqa: E402
from objective.power import wind_speed_bins  # noqa: E402
from objective.weibull import WeibullResource  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402

DEFAULT_POINTS = 300

#: Bisection tolerance for the reachable extent of one variable: 0.1 mm of
#: chord, 0.05 degrees of twist.
TOLERANCE_CHORD_M = 1e-4
TOLERANCE_TWIST_RAD = math.radians(0.05)


def load_x0():
    with open(os.path.join(REPO_ROOT, "verification", "baseline", "x0.json"),
              encoding="utf-8") as handle:
        artefact = json.load(handle)
    vector = np.array(artefact["chord_control_points_m"]
                      + artefact["twist_control_points_rad"])
    return vector, artefact


def box_interval(index, x0, n_chord, bounds):
    """(lo, hi, label) of design variable `index` in the design box."""

    if index < n_chord:
        return bounds.chord_min_m, bounds.chord_max_m, f"chord CP {index}"
    return bounds.twist_min_rad, bounds.twist_max_rad, f"twist CP {index - n_chord}"


def objective_at(index, value, x0, parameterisation, mass):
    """(J, converged) at `x0` with variable `index` moved to `value`."""

    vector = x0.copy()
    vector[index] = value
    powers = bin_powers(vector, parameterisation=parameterisation)
    return (float(np.sum(mass * powers["power_w"]) * HOURS_PER_YEAR / 1e6),
            bool(powers["all_converged"]))


def is_good(index, value, x0, parameterisation, mass):
    """True if the objective is defined **and** the rotor converged there."""

    try:
        _j, converged = objective_at(index, value, x0, parameterisation, mass)
    except PolarDomainError:
        return False
    return converged


def reachable_interval(index, x0, lo, hi, parameterisation, mass, tolerance):
    """
    `(lower, upper)` -- the widest interval containing `x0[index]`, inside
    `[lo, hi]`, on which the objective is defined and every station
    converges. Bisection on each side; `x0` is good, so both searches start
    from a known-good point.
    """

    x0_value = float(x0[index])

    def extent(limit):
        if is_good(index, limit, x0, parameterisation, mass):
            return float(limit)
        good, bad = x0_value, float(limit)
        while abs(bad - good) > tolerance:
            mid = 0.5 * (good + bad)
            if is_good(index, mid, x0, parameterisation, mass):
                good = mid
            else:
                bad = mid
        return good

    return extent(lo), extent(hi)


def diagnose(grid, j):
    """
    The four defects the gate names, as numbers, over the defined points.

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
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)

    # The objective's weights, computed ONCE. They depend only on the bin
    # scheme and the site resource -- not on `d`.
    resource = WeibullResource.from_config()
    edges, _midpoints, _width = wind_speed_bins()
    mass = resource.probability_between(edges[:-1], edges[1:])
    x0, artefact = load_x0()
    n_chord = artefact["n_chord"]

    sweeps = []
    started = time.perf_counter()

    for index in range(len(x0)):
        lo, hi, label = box_interval(index, x0, n_chord, bounds)
        tolerance = (TOLERANCE_CHORD_M if index < n_chord
                     else TOLERANCE_TWIST_RAD)
        lower, upper = reachable_interval(index, x0, lo, hi, parameterisation,
                                          mass, tolerance)
        grid = np.linspace(lower, upper, points)

        values, converged, undefined, unconverged = [], [], [], []
        for value in grid:
            try:
                j, ok = objective_at(index, value, x0, parameterisation, mass)
            except PolarDomainError:
                j, ok = None, False
                undefined.append(float(value))
            if not ok:
                j = None
                unconverged.append(float(value))
            values.append(j)
            converged.append(ok)

        defined = [j for j in values if j is not None]
        diagnostics = diagnose(np.array([v for v, j in zip(grid, values)
                                         if j is not None]), defined) \
            if len(defined) >= 3 else None

        sweeps.append({
            "index": index,
            "label": label,
            "x0": float(x0[index]),
            "box_lo": float(lo),
            "box_hi": float(hi),
            "lo": float(lower),
            "hi": float(upper),
            "box_fraction_reachable": float((upper - lower) / (hi - lo)),
            "grid": [float(v) for v in grid],
            "J": values,
            "all_converged": all(converged),
            "n_not_converged": int(len(converged) - sum(converged)),
            "n_undefined": len(undefined),
            "undefined": undefined,
            "unconverged": unconverged,
            "diagnostics": diagnostics,
        })

        elapsed = time.perf_counter() - started
        shown = f"{min(defined):.3f}-{max(defined):.3f}" if defined else "n/a"
        print(f"  [{index + 1}/{len(x0)}] {label:14} "
              f"box {lo:+.4f}..{hi:+.4f} reachable {lower:+.4f}..{upper:+.4f} "
              f"({100.0 * (upper - lower) / (hi - lo):.0f}% of the box)  "
              f"J {shown} MWh/yr  converged {all(converged)}  "
              f"({elapsed / 60:.1f} min elapsed)", flush=True)

    return {
        "points_per_sweep": points,
        "objective": "annual_energy_mwh (Weibull-weighted AEP, MWh/yr)",
        "resource": {
            "k": resource.k,
            "c_ms": resource.c,
            "mean_speed_ms": resource.mean_speed,
            "source": "config/site.yaml; see verification/wind_resource/",
        },
        "sweep_ranges": {
            "basis": "the design box (config/rotor_design.yaml "
                     "parameterisation.bounds, DesignBounds.from_config), "
                     "clipped to the part of each variable's interval on "
                     "which the objective is defined and converged: a station "
                     "Reynolds number below the polar cache floor of 40000 "
                     "makes the objective undefined (PolarDomainError), and "
                     "at the high-twist end of the twist sweeps the BEM "
                     "solver does not converge at every grid point",
            "bounds": bounds.as_record(),
            "chord_min_m": float(bounds.chord_min_m),
            "chord_max_m": float(bounds.chord_max_m),
            "twist_min_deg": math.degrees(bounds.twist_min_rad),
            "twist_max_deg": math.degrees(bounds.twist_max_rad),
            "extent_method": "bisection from the x0 value towards each box "
                             "bound, tolerance 0.1 mm of chord / 0.05 deg of "
                             "twist",
        },
        "n_chord": n_chord,
        "x0_source": "verification/baseline/x0.json",
        "sweeps": sweeps,
    }


def plot(results, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from plotting import figstyle
    except ImportError:
        print("matplotlib not available; skipping the figures")
        return

    figstyle.apply()
    sweeps = results["sweeps"]
    n_chord = int(results["n_chord"])

    groups = [
        dict(title="Chord control points",
             panel=sweeps[:n_chord], scale=1e3,
             xlabel=figstyle.LABELS["chord"],
             gradient_label=figstyle.label(
                 r"$\mathrm{d}J/\mathrm{d}c$", None,
                 r"MWh yr$^{-1}$ mm$^{-1}$"),
             symbols=[rf"$c_{{{k}}}$" for k in range(n_chord)]),
        dict(title="Twist control points",
             panel=sweeps[n_chord:], scale=float(np.degrees(1.0)),
             xlabel=figstyle.LABELS["twist"],
             gradient_label=figstyle.label(
                 r"$\mathrm{d}J/\mathrm{d}\theta$", None,
                 r"MWh yr$^{-1}$ deg$^{-1}$"),
             symbols=[rf"$\theta_{{{k}}}$"
                      for k in range(len(sweeps) - n_chord)]),
    ]

    objective_fig, objective_axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE)
    difference_fig, difference_axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE)

    for group, ax_objective, ax_difference in zip(groups, objective_axes,
                                                  difference_axes):
        for symbol, sweep in zip(group["symbols"], group["panel"]):
            x = group["scale"] * np.array(sweep["grid"])
            j = np.array([np.nan if value is None else value
                          for value in sweep["J"]], dtype=float)
            finite = np.isfinite(j)
            x, j = x[finite], j[finite]
            span = float(j.max() - j.min())
            ax_objective.plot(x, (j - j.min()) / span, label=symbol)
            centres = 0.5 * (x[:-1] + x[1:])
            ax_difference.plot(centres, np.diff(j) / np.diff(x), label=symbol)
        figstyle.title(ax_objective, group["title"])
        figstyle.title(ax_difference, group["title"])
        ax_objective.set_xlabel(group["xlabel"])
        ax_objective.set_ylabel(figstyle.label("Objective", "J", "--"))
        ax_objective.legend(ncol=2)
        ax_difference.set_xlabel(group["xlabel"])
        ax_difference.set_ylabel(group["gradient_label"])
        ax_difference.legend(ncol=2)

    objective_fig.tight_layout()
    difference_fig.tight_layout()

    objective_path = os.path.join(out_dir, "smoothness_gate_objective.png")
    difference_path = os.path.join(out_dir, "smoothness_gate_difference.png")
    figstyle.save(objective_fig, objective_path)
    figstyle.save(difference_fig, difference_path)
    plt.close(objective_fig)
    plt.close(difference_fig)
    print(f"wrote {objective_path}")
    print(f"wrote {difference_path}")


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

    plot(results, _HERE)

    print()
    header = (f"{'variable':>14}{'converged':>11}{'undefined':>11}{'distinct':>10}"
              f"{'flat steps':>12}{'max dJ/med':>12}{'max d2J/med':>13}{'monotone':>10}")
    print(header)
    for sweep in results["sweeps"]:
        d = sweep["diagnostics"]
        if d is None:
            print(f"{sweep['label']:>14}  no defined points")
            continue
        print(f"{sweep['label']:>14}{str(sweep['all_converged']):>11}"
              f"{sweep['n_undefined']:>11}"
              f"{d['distinct_fraction']:>10.3f}{d['flat_first_difference_steps']:>12}"
              f"{d['max_first_difference_ratio']:>12.2f}"
              f"{d['max_second_difference_ratio']:>13.2f}"
              f"{str(d['monotone']):>10}")


if __name__ == "__main__":
    main()
