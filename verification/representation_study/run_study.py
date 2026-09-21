"""
Representation study: how many control points does the blade need?

Plan section 4.3 and step 1.6's last bullet. The final control-point count is
to be "justified against ... a representation study comparing 6, 8, 10, 12 and
16 control points against a reference distribution", documenting the
comparison rather than asserting a number.

Method
-------
The reference is the analytic Schmitz blade for this rotor (R = 2.0 m, B = 3,
lambda = 6.5) at the SG6043 maximum-L/D point. That is the honest target: it is
the shape the optimiser will be started from, so what matters is whether the
parameterisation can hold it, not whether it can hold some generic taper
chosen to be easy.

Because the parameterisation is linear in its control points, fitting is a
single linear least-squares solve -- exact, no starting guess, no local minima.
The fitting error reported is therefore a property of the *basis*, not of an
optimiser that fitted it.

Two things are reported that a pure fitting-error table would miss:

  * **Conditioning** of the basis matrix. Plan 4.3 lists "numerical
    conditioning of the reduced Jacobian" first among the criteria, and it is
    the one that degrades as control points are added -- the basis functions
    start to overlap and the reduced Jacobian becomes harder to invert, which
    an optimiser feels long before a human notices.
  * **Sawtooth reproduction**, the property section 4.2 says the
    parameterisation exists to suppress. Reported as the fraction of a pure
    strip-to-strip oscillation the basis can reproduce -- 0 means the null
    space is invisible to the optimiser, 1 is what per-station design
    variables would give. More control points fit the reference better AND
    admit more oscillation, and choosing a count is exactly that trade, so
    both are measured rather than only the flattering one.

Run from the repo root:

    python verification/representation_study/run_study.py

Redraw the committed figure from the committed JSON, without recomputing the
study:

    python verification/representation_study/run_study.py --replot

The JSON holds the reference profile and each count's error record, not the
fitted curves; the replot rebuilds each fit with `BladeParameterisation` (the
same linear least-squares solve) and refuses to draw if the rebuilt RMS error
disagrees with the committed one.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from config import load_design_rotor  # noqa: E402
from design import BladeParameterisation  # noqa: E402
from design.schmitz import (  # noqa: E402
    DEFAULT_DESIGN_REYNOLDS,
    max_lift_to_drag_point,
    schmitz_distribution,
)
from polars.polar import interpolant_for  # noqa: E402

#: Total control-point counts to compare, split evenly between chord and
#: twist. Plan 4.3 names 6, 8, 10, 12 and 16.
TOTAL_COUNTS = [6, 8, 10, 12, 16]

#: Shared with the baseline (design/schmitz.py) so the study and x0 cannot
#: drift apart on the design point.
DESIGN_REYNOLDS = DEFAULT_DESIGN_REYNOLDS

#: Sawtooth amplitude used to probe oscillation rejection, metres.
SAWTOOTH_AMPLITUDE = 0.02


def reference_distribution(n_stations):
    """The analytic Schmitz blade, sampled at the BEM strip midpoints."""

    reference = BladeParameterisation(n_stations=n_stations)
    interpolant = interpolant_for("sg6043")
    alpha_design, design_cl, lift_to_drag = max_lift_to_drag_point(
        interpolant, DESIGN_REYNOLDS)

    r_over_R = reference.radii / reference.radius_m
    chord, twist = schmitz_distribution(r_over_R, design_cl, alpha_design)

    return reference, chord, twist, {
        "reynolds": DESIGN_REYNOLDS,
        "alpha_design_deg": math.degrees(alpha_design),
        "design_cl": design_cl,
        "max_lift_to_drag": lift_to_drag,
    }


def study():
    design = load_design_rotor()
    n_stations = design.parameterisation.n_bem_strips
    reference, chord_target, twist_target, airfoil_point = reference_distribution(n_stations)

    sawtooth = chord_target + SAWTOOTH_AMPLITUDE * np.array(
        [(-1.0) ** i for i in range(n_stations)])

    rows = []
    for total in TOTAL_COUNTS:
        n_chord = total // 2
        n_twist = total - n_chord
        parameterisation = BladeParameterisation(
            n_chord=n_chord, n_twist=n_twist, n_stations=n_stations)

        _, chord_rms, twist_rms = parameterisation.fit(chord_target, twist_target)

        # Fraction of a strip-to-strip oscillation the basis can reproduce.
        #
        # The fit is linear, so fitting (smooth + oscillation) and fitting
        # smooth alone differ by exactly the fit of the oscillation on its own
        # -- which means this is the amplitude the basis would hand an
        # optimiser looking for the sawtooth null space of section 4.2.
        # 0.0 = the oscillation is invisible to the parameterisation, which is
        # the property that removes the null space; 1.0 = reproduced exactly,
        # which is what per-station design variables would give.
        oscillation = sawtooth - chord_target
        fitted_oscillation = parameterisation.chord_basis @ np.linalg.lstsq(
            parameterisation.chord_basis, oscillation, rcond=None)[0]
        reproduced = float(np.sqrt(np.mean(fitted_oscillation ** 2))
                           / np.sqrt(np.mean(oscillation ** 2)))

        rows.append({
            "total_control_points": total,
            "n_chord": n_chord,
            "n_twist": n_twist,
            "degree": parameterisation.degree,
            "chord_rms_error_m": chord_rms,
            "chord_rms_error_pct_of_mean": 100.0 * chord_rms / float(np.mean(chord_target)),
            "twist_rms_error_deg": math.degrees(twist_rms),
            "chord_basis_condition": float(np.linalg.cond(parameterisation.chord_basis)),
            "twist_basis_condition": float(np.linalg.cond(parameterisation.twist_basis)),
            "sawtooth_reproduced_fraction": reproduced,
        })

    return {
        "rotor": {
            "radius_m": design.radius_m,
            "n_blades": design.n_blades,
            "design_tsr": design.design_tsr,
            "root_fraction": design.root_fraction,
            "n_bem_strips": n_stations,
        },
        "airfoil_design_point": airfoil_point,
        "reference": {
            "source": "analytic Schmitz at the SG6043 max-L/D point",
            "chord_m": [float(value) for value in chord_target],
            "twist_deg": [math.degrees(value) for value in twist_target],
            "r_over_R": [float(value) for value in reference.radii / reference.radius_m],
        },
        "counts": rows,
    }


def fit_curve(results, total, chord_target, twist_target):
    """
    Rebuild one control-point count's fitted curves from the committed study.

    The JSON stores the fit errors, not the fitted curves; the fit is the
    deterministic linear least-squares solve of `BladeParameterisation`, so
    the curves are rebuilt here and checked against the committed RMS values
    before they are drawn. No solver is involved.
    """

    row = next(record for record in results["counts"]
               if record["total_control_points"] == total)
    parameterisation = BladeParameterisation(
        n_chord=row["n_chord"], n_twist=row["n_twist"],
        n_stations=len(chord_target), degree=row["degree"],
        radius_m=results["rotor"]["radius_m"],
        root_fraction=results["rotor"]["root_fraction"])
    vector, chord_rms, twist_rms = parameterisation.fit(chord_target,
                                                        twist_target)

    if not math.isclose(chord_rms, row["chord_rms_error_m"], rel_tol=1e-9):
        raise SystemExit(
            f"refit at {total} control points gives chord RMS "
            f"{chord_rms:.6e} m, committed {row['chord_rms_error_m']:.6e} m; "
            "the figure would misrepresent the study")
    twist_rms_deg = math.degrees(twist_rms)
    if not math.isclose(twist_rms_deg, row["twist_rms_error_deg"], rel_tol=1e-9):
        raise SystemExit(
            f"refit at {total} control points gives twist RMS "
            f"{twist_rms_deg:.6e} deg, committed "
            f"{row['twist_rms_error_deg']:.6e} deg; the figure would "
            "misrepresent the study")
    return parameterisation, vector


def plot_study(results, out_dir):
    """Draw `representation_study.png`: the fits | the fitting error."""

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from plotting import figstyle
    except ImportError:
        print("matplotlib not available; skipping the figure")
        return

    figstyle.apply()
    fig, axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE,
                             layout="constrained")

    reference = results["reference"]
    radius = [results["rotor"]["radius_m"] * value
              for value in reference["r_over_R"]]
    chord_target = reference["chord_m"]
    twist_target = [math.radians(value) for value in reference["twist_deg"]]

    # Chord on the left axis (solid), twist on the right (dashed); the
    # control-point colours are the same for both blocks.
    ax = axes[0]
    ax.plot(radius, [1000.0 * value for value in chord_target], color="0.2",
            label="analytic")
    ax.set_xlabel(figstyle.LABELS["radius"])
    ax.set_ylabel(figstyle.LABELS["chord"])
    figstyle.title(ax, "Fit to the analytic blade")

    ax_twist = ax.twinx()
    ax_twist.plot(radius, reference["twist_deg"], color="0.2", linestyle="--")
    ax_twist.set_ylabel(figstyle.LABELS["twist"])
    ax_twist.grid(False)

    for colour, total in zip(("#0072B2", "#E69F00", "#009E73"), (6, 10, 16)):
        parameterisation, vector = fit_curve(results, total, chord_target,
                                             twist_target)
        ax.plot(radius,
                [1000.0 * value for value in parameterisation.chord(vector)],
                color=colour, label=f"{total} control points")
        ax_twist.plot(radius,
                      [math.degrees(value)
                       for value in parameterisation.twist(vector)],
                      color=colour, linestyle="--")

    ax.legend(loc="upper right")

    ax = axes[1]
    counts = [row["total_control_points"] for row in results["counts"]]
    ax.semilogy(counts,
                [1000.0 * row["chord_rms_error_m"]
                 for row in results["counts"]], "o-", color="#0072B2")
    ax.set_xlabel(figstyle.LABELS["design_variables"])
    ax.set_ylabel(figstyle.label("Chord RMS error", "e_c", "mm"),
                  color="#0072B2")
    ax.tick_params(axis="y", labelcolor="#0072B2")
    ax.set_xticks(counts)
    figstyle.title(ax, "Fitting error")

    ax_twist = ax.twinx()
    ax_twist.semilogy(counts,
                      [row["twist_rms_error_deg"] for row in results["counts"]],
                      "s--", color="#E69F00")
    ax_twist.set_ylabel(figstyle.label("Twist RMS error", r"e_\theta", "°"),
                        color="#E69F00")
    ax_twist.tick_params(axis="y", labelcolor="#E69F00")
    ax_twist.grid(False)

    path = os.path.join(out_dir, "representation_study.png")
    figstyle.save(fig, path)
    plt.close(fig)
    print(f"wrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replot", action="store_true",
                        help="redraw representation_study.png from the "
                             "committed JSON, without recomputing the study")
    args = parser.parse_args(argv)

    if args.replot:
        with open(os.path.join(_HERE, "representation_study.json"),
                  encoding="utf-8") as handle:
            results = json.load(handle)
        plot_study(results, _HERE)
        return

    results = study()

    out = os.path.join(_HERE, "representation_study.json")
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=1)
    print(f"wrote {out}")

    plot_study(results, _HERE)

    point = results["airfoil_design_point"]
    print(f"\nSG6043 design point at Re = {point['reynolds']:,.0f}: "
          f"max L/D = {point['max_lift_to_drag']:.1f} at "
          f"alpha = {point['alpha_design_deg']:.2f} deg, Cl = {point['design_cl']:.4f}\n")

    header = (f"{'total':>6}{'chord':>4}{'twist':>6}{'deg':>5}"
              f"{'chord RMS [m]':>15}{'[% mean]':>10}{'twist RMS [deg]':>17}"
              f"{'cond(N)':>10}{'sawtooth repr.':>16}")
    print(header)
    for row in results["counts"]:
        print(f"{row['total_control_points']:>6}{row['n_chord']:>4}"
              f"{row['n_twist']:>6}{row['degree']:>5}"
              f"{row['chord_rms_error_m']:>15.3e}"
              f"{row['chord_rms_error_pct_of_mean']:>10.3f}"
              f"{row['twist_rms_error_deg']:>17.4f}"
              f"{row['chord_basis_condition']:>10.2f}"
              f"{row['sawtooth_reproduced_fraction']:>15.2%}")


if __name__ == "__main__":
    main()
