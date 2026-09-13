"""
Build the Schmitz baseline, evaluate it, and commit `x0` (plan step 1.7).

The plan calls the evaluated performance of this blade "the reference numbers
for the entire results chapter", so this script's output is a fixed artefact:
`x0.json` is what Phases 2-5 start from, and `baseline_reference.json` is what
every later result is compared against.

Regeneration is deliberate, not routine -- the same rule as the golden files.
If these numbers move, something in the solver or the polar layer moved, and
the commit message has to say what and why.

Run from the repo root:

    python verification/baseline/generate_baseline.py

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os
import sys

import numpy as np

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
            "gradient and optimisation work (plan step 1.7). Chord block first "
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


def plot(baseline, performance, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping the figure")
        return

    blue, orange, green = "#0072B2", "#E69F00", "#009E73"
    parameterisation = baseline.parameterisation
    r_over_R = parameterisation.radii / parameterisation.radius_m

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

    ax = axes[0]
    ax.plot(r_over_R, baseline.chord_target_m * 1000, "k--", linewidth=1.6,
            label="analytic Schmitz")
    ax.plot(r_over_R, baseline.chord_m * 1000, color=blue, linewidth=1.6,
            label="fitted (x0)")
    ax.set_xlabel("$r/R$")
    ax.set_ylabel("chord [mm]")
    ax.set_title("Baseline chord")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, linewidth=0.3, alpha=0.4)

    ax2 = axes[0].twinx()
    ax2.plot(r_over_R, np.degrees(baseline.twist_rad), color=orange,
             linewidth=1.2, alpha=0.8)
    ax2.set_ylabel("twist [deg]", color=orange)
    ax2.tick_params(axis="y", labelcolor=orange)

    ax = axes[1]
    tsr = [row["tsr"] for row in performance["cp_lambda"]]
    cp = [row["Cp"] for row in performance["cp_lambda"]]
    ax.plot(tsr, cp, color=blue, linewidth=1.6)
    ax.axvline(performance["design_point"]["tsr"], color=green, linestyle=":",
               linewidth=1.2, label=f"design $\\lambda$ = {performance['design_point']['tsr']}")
    ax.axhline(16.0 / 27.0, color="0.5", linestyle="--", linewidth=1.0,
               label="Betz limit")
    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("$C_p$")
    ax.set_title("Baseline $C_p$-$\\lambda$")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, linewidth=0.3, alpha=0.4)

    ax = axes[2]
    speeds = [row["v_inf"] for row in performance["operating_line"]]
    power = [row["power_w"] for row in performance["operating_line"]]
    thrust = [row["thrust_n"] for row in performance["operating_line"]]
    ax.plot(speeds, power, color=blue, linewidth=1.6, label="power [W]")
    ax.set_xlabel("wind speed [m/s]")
    ax.set_ylabel("power [W]", color=blue)
    ax.tick_params(axis="y", labelcolor=blue)
    ax.set_title("Below-rated operating line")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax3 = ax.twinx()
    ax3.plot(speeds, thrust, color=orange, linewidth=1.4, label="thrust [N]")
    ax3.set_ylabel("thrust [N]", color=orange)
    ax3.tick_params(axis="y", labelcolor=orange)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


def main():
    baseline, x0, performance = build()

    for name, payload in (("x0.json", x0),
                          ("baseline_reference.json", performance)):
        path = os.path.join(_HERE, name)
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1)
        print(f"wrote {path}")

    plot(baseline, performance, os.path.join(_HERE, "baseline.png"))

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
               else "OUTSIDE BAND -- see the 2026-09-13 journal entry")
    print(f"  AEP               {performance['aep_mwh_per_year']:.3f} MWh/yr "
          f"(band {band[0]:.1f}-{band[1]:.1f}, {verdict})")
    for key, value in performance["outstanding"].items():
        print(f"  OUTSTANDING       {key}: {value}")


if __name__ == "__main__":
    main()
