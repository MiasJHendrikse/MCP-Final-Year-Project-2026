"""
Build the Schmitz baseline, evaluate it, and commit `x0`.

The evaluated performance of this blade gives the reference numbers for the
entire results chapter, so this script's output is a fixed artefact:
`x0.json` is what Phases 2-5 start from, and `baseline_reference.json` is what
every later result is compared against.

Regeneration is deliberate, not routine -- the same rule as the golden files.
If these numbers move, something in the solver or the polar layer moved, and
the commit message has to say what and why.

Run from the repo root:

    python verification/baseline/generate_baseline.py

Redraw the committed figures from the committed JSON, without recomputing
anything:

    python verification/baseline/generate_baseline.py --replot

Writes `x0.json`, `baseline_reference.json`, `baseline_geometry.png` and
`baseline_operating_line.png` next to this file. The earlier three-panel
`baseline.png` is superseded by the two figure files and is no longer written.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import argparse
import json
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from config import load_design_rotor, load_site  # noqa: E402
from design import build_schmitz_baseline, evaluate_baseline  # noqa: E402


def build():
    baseline = build_schmitz_baseline()
    performance = evaluate_baseline(baseline)
    design = load_design_rotor()
    site = load_site()

    x0 = {
        "description": (
            "Fitted-Schmitz control-point vector. The starting point for all "
            "gradient and optimisation work. Chord block first "
            "(metres), then twist block (radians)."
        ),
        "n_chord": baseline.parameterisation.n_chord,
        "n_twist": baseline.parameterisation.n_twist,
        "n_bem_strips": baseline.parameterisation.n_stations,
        "spline_degree": baseline.parameterisation.degree,
        "chord_control_points_m": [
            float(value) for value in baseline.design_vector[:baseline.parameterisation.n_chord]],
        "twist_control_points_rad": [
            float(value) for value in baseline.design_vector[baseline.parameterisation.n_chord:]],
        "twist_control_points_deg": [
            math.degrees(value)
            for value in baseline.design_vector[baseline.parameterisation.n_chord:]],
        "construction": {
            "source": "analytic Schmitz, least-squares fitted onto the spline",
            "design_reynolds": baseline.design_reynolds,
            "alpha_design_deg": math.degrees(baseline.alpha_design_rad),
            "design_cl": baseline.design_cl,
            "max_lift_to_drag": baseline.max_lift_to_drag,
            "chord_rms_fit_error_m": baseline.chord_rms_error_m,
            "twist_rms_fit_error_deg": math.degrees(baseline.twist_rms_error_rad),
        },
        "feasibility": baseline.feasibility,
        "rotor": {
            "radius_m": design.radius_m,
            "n_blades": design.n_blades,
            "root_fraction": design.root_fraction,
            "design_tsr": design.design_tsr,
            "rated_wind_speed_ms": design.rated_wind_speed_ms,
            "cut_in_wind_speed_ms": design.cut_in_wind_speed_ms,
        },
        "site": {
            "air_density": site.air_density,
            "kinematic_viscosity": site.kinematic_viscosity,
        },
    }

    return baseline, x0, performance


def plot_geometry(performance, x0, out_dir):
    """`baseline_geometry.png`: chord and twist | the C_P-lambda curve."""

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from plotting import figstyle
    except ImportError:
        print("matplotlib not available; skipping the figures")
        return

    figstyle.apply()
    fig, axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE,
                             layout="constrained")
    navy, rust, teal = figstyle.PALETTE[:3]

    spanwise = performance["spanwise"]
    radius = [row["r"] for row in spanwise]

    # Both curves are the one reference blade; the panel title names it and
    # the legend names the two quantities, one per axis.
    ax = axes[0]
    ax.plot(radius, [1000.0 * row["chord_m"] for row in spanwise], color=navy,
            label=r"chord $c$")
    ax.set_xlabel(figstyle.LABELS["radius"])
    ax.set_ylabel(figstyle.LABELS["chord"], color=navy)
    ax.tick_params(axis="y", labelcolor=navy)
    figstyle.title(ax, r"Schmitz reference $\mathbf{x}_0$")

    ax_twist = figstyle.twin(ax)
    ax_twist.plot(radius, [row["twist_deg"] for row in spanwise], color=rust,
                  linestyle="--", label=r"twist $\theta$")
    ax_twist.set_ylabel(figstyle.LABELS["twist"], color=rust)
    ax_twist.tick_params(axis="y", labelcolor=rust)
    figstyle.merged_legend(ax, ax_twist, loc="upper right")

    ax = axes[1]
    ax.plot([row["tsr"] for row in performance["cp_lambda"]],
            [row["Cp"] for row in performance["cp_lambda"]], color=navy,
            label=r"$C_P(\lambda)$")
    ax.axhline(16.0 / 27.0, color="0.45", linestyle="--", linewidth=1.0,
               label="Betz limit")

    design_tsr = float(x0["rotor"]["design_tsr"])
    design_row = next(row for row in performance["cp_lambda"]
                      if abs(row["tsr"] - design_tsr) < 1e-9)
    ax.plot([design_tsr], [design_row["Cp"]], marker="o", linestyle="none",
            color=teal, label=fr"design $\lambda$ = {design_tsr:g}")

    rated = performance["design_point"]
    ax.plot([rated["tsr"]], [rated["Cp"]], marker="s", linestyle="none",
            color=rust, label=fr"rated point $\lambda$ = {rated['tsr']:.2f}")

    ax.set_xlabel(figstyle.LABELS["tip_speed_ratio"])
    ax.set_ylabel(figstyle.LABELS["power_coefficient"])
    ax.set_ylim(0.0, 0.7)
    ax.legend(loc="lower right")
    figstyle.title(ax, "Power coefficient")

    path = os.path.join(out_dir, "baseline_geometry.png")
    figstyle.save(fig, path)
    plt.close(fig)
    print(f"wrote {path}")


def plot_operating_line(performance, out_dir):
    """`baseline_operating_line.png`: power and thrust against wind speed."""

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from plotting import figstyle
    except ImportError:
        print("matplotlib not available; skipping the figures")
        return

    figstyle.apply()
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    navy, rust = figstyle.PALETTE[:2]

    speeds = [row["v_inf"] for row in performance["operating_line"]]

    ax.plot(speeds, [row["power_w"] for row in performance["operating_line"]],
            color=navy, label=r"power $P$")
    ax.set_xlabel(figstyle.LABELS["wind_speed"])
    ax.set_ylabel(figstyle.label("Power", "P", "W"), color=navy)
    ax.tick_params(axis="y", labelcolor=navy)

    ax_thrust = figstyle.twin(ax)
    ax_thrust.plot(speeds,
                   [row["thrust_n"] for row in performance["operating_line"]],
                   color=rust, linestyle="--", label=r"thrust $T$")
    ax_thrust.set_ylabel(figstyle.label("Thrust", "T", "N"), color=rust)
    ax_thrust.tick_params(axis="y", labelcolor=rust)
    figstyle.merged_legend(ax, ax_thrust, loc="lower right")

    path = os.path.join(out_dir, "baseline_operating_line.png")
    figstyle.save(fig, path)
    plt.close(fig)
    print(f"wrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replot", action="store_true",
                        help="redraw the figures from the committed JSON, "
                             "without recomputing anything")
    args = parser.parse_args(argv)

    if args.replot:
        with open(os.path.join(_HERE, "baseline_reference.json"),
                  encoding="utf-8") as handle:
            performance = json.load(handle)
        with open(os.path.join(_HERE, "x0.json"), encoding="utf-8") as handle:
            x0 = json.load(handle)
        plot_geometry(performance, x0, _HERE)
        plot_operating_line(performance, _HERE)
        return

    _, x0, performance = build()

    for name, payload in (("x0.json", x0),
                          ("baseline_reference.json", performance)):
        path = os.path.join(_HERE, name)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1)
        print(f"wrote {path}")

    plot_geometry(performance, x0, _HERE)
    plot_operating_line(performance, _HERE)

    point = performance["design_point"]
    print()
    print(f"  design point      Cp = {point['Cp']:.4f}, Ct = {point['Ct']:.4f} "
          f"at lambda = {point['tsr']}, V = {point['v_inf']} m/s")
    best = max(performance["cp_lambda"], key=lambda row: row["Cp"])
    print(f"  peak Cp           {best['Cp']:.4f} at lambda = {best['tsr']}")
    print(f"  root bending      {performance['root_bending_moment_nm']:.1f} N.m per blade")
    print(f"  peak thrust       {performance['peak_thrust_n']:.1f} N at "
          f"{performance['peak_thrust_wind_speed_ms']} m/s")
    print(f"  fit error         chord {x0['construction']['chord_rms_fit_error_m'] * 1000:.3f} mm, "
          f"twist {x0['construction']['twist_rms_fit_error_deg']:.4f} deg")
    print(f"  feasibility       checked = {x0['feasibility']['checked']}")

    band = performance["aep_sanity_band_mwh_per_year"]
    verdict = ("in band" if performance["aep_in_sanity_band"]
               else "OUTSIDE BAND -- see config/rotor_design.yaml `sanity:`")
    print(f"  AEP               {performance['aep_mwh_per_year']:.3f} MWh/yr "
          f"(band {band[0]:.1f}-{band[1]:.1f}, {verdict})")
    for key, value in performance["outstanding"].items():
        print(f"  OUTSTANDING       {key}: {value}")


if __name__ == "__main__":
    main()
