"""
Spline-fit error of the Schmitz baseline (PROJECT_DIRECTION_v2 §7.3 point 4).

The baseline `x0` is the analytic Schmitz chord and twist least-squares
projected onto the optimiser's 5 + 5 cubic B-spline. §7.3 point 4 requires
the fitting error to be reported, because an optimisation gain measured
against a spline-represented baseline is confounded with the representation
difference unless that difference is known.

What this script measures
--------------------------
1. `x0` from `verification/baseline/x0.json` is rebuilt with
   `design.build_schmitz_baseline()` and checked to be the same vector.
2. Analytic Schmitz chord and twist (same design Cl and alpha as the
   baseline construction, `design_reynolds = 200 k`) versus the spline through
   `x0`, at the 25 BEM stations and on a fine grid of 401 span points: max
   and RMS error in chord (m) and twist (deg), signed error vs span.
3. AEP of the *analytic* blade (station values straight from the Schmitz
   formulas, same solver, same corrections, same resource) versus AEP of
   `x0` (the spline), versus AEP of A4's optimum. The first difference is the
   representation difference; the second is the optimisation gain.

Outputs next to this script: `fit_error.json`, `fit_error.png`, and the
README written from them.

Bounds: the configured set (`DesignBounds.from_config()`), grounded
2026-09-19 -- `chord_max_m = 0.30 m`; A4's gain is stated under those bounds.

Run from the repo root (~5 s):

    python verification/spline_fit_error/run_fit_error.py

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import datetime
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from bem.rotor import RotorGeometry  # noqa: E402
from design import BladeParameterisation, basis_matrix, build_schmitz_baseline  # noqa: E402
from design.schmitz import schmitz_distribution  # noqa: E402
from objective import HOURS_PER_YEAR, WeibullResource, annual_energy_mwh, power_per_bin, wind_speed_bins  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
A4_RESULT_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation", "result.json")
JSON_PATH = os.path.join(_HERE, "fit_error.json")
FIGURE_PATH = os.path.join(_HERE, "fit_error.png")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")
N_FINE = 401


def spline_on(parameterisation, s, design_vector):
    """Chord and twist of `design_vector` at arbitrary span fractions `s`."""

    p = parameterisation
    chord_points, twist_points = p.split(design_vector)
    n_c = basis_matrix(s, p.n_chord, p.degree)
    n_t = basis_matrix(s, p.n_twist, p.degree)
    return n_c @ chord_points, n_t @ twist_points


def aep_of_geometry(geometry, resource):
    """`annual_energy_mwh` for a geometry that is not a design vector."""

    edges, _mid, _w = wind_speed_bins()
    powers = power_per_bin(geometry)
    mass = resource.probability_between(edges[:-1], edges[1:])
    return float(np.sum(powers["power_w"] * mass) * HOURS_PER_YEAR) / 1e6, powers


def errors(chord_fit, chord_ref, twist_fit, twist_ref):
    dc = chord_fit - chord_ref
    dt = np.degrees(twist_fit - twist_ref)
    return {
        "chord_max_abs_m": float(np.max(np.abs(dc))),
        "chord_rms_m": float(np.sqrt(np.mean(dc ** 2))),
        "chord_max_rel": float(np.max(np.abs(dc) / chord_ref)),
        "twist_max_abs_deg": float(np.max(np.abs(dt))),
        "twist_rms_deg": float(np.sqrt(np.mean(dt ** 2))),
        "chord_error_m": dc.tolist(),
        "twist_error_deg": dt.tolist(),
    }


def main():
    parameterisation = BladeParameterisation()
    baseline = build_schmitz_baseline(parameterisation)
    resource = WeibullResource.from_config()

    with open(X0_PATH, encoding="utf-8") as handle:
        x0_artefact = json.load(handle)
    x0 = np.array(x0_artefact["chord_control_points_m"] + x0_artefact["twist_control_points_rad"])
    x0_rebuild_diff = float(np.max(np.abs(baseline.design_vector - x0)))

    with open(A4_RESULT_PATH, encoding="utf-8") as handle:
        a4 = json.load(handle)
    x_star = np.array(a4["x_star"])

    # -- 1. at the 25 stations -------------------------------------------
    r = parameterisation.radii
    chord_fit, twist_fit = parameterisation.evaluate(x0)
    station_errors = errors(chord_fit, baseline.chord_target_m, twist_fit, baseline.twist_target_rad)

    # -- 2. on a fine grid --------------------------------------------------
    s_fine = np.linspace(0.0, 1.0, N_FINE)
    r_fine = parameterisation.radius_m * (
        parameterisation.root_fraction + (1.0 - parameterisation.root_fraction) * s_fine)
    chord_ref_f, twist_ref_f = schmitz_distribution(
        r_fine / parameterisation.radius_m, baseline.design_cl, baseline.alpha_design_rad)
    chord_fit_f, twist_fit_f = spline_on(parameterisation, s_fine, x0)
    fine_errors = errors(chord_fit_f, chord_ref_f, twist_fit_f, twist_ref_f)

    # -- 3. AEP: analytic vs x0 vs optimum ----------------------------------
    analytic_geometry = RotorGeometry(
        r=[float(v) for v in r],
        chord=[float(v) for v in baseline.chord_target_m],
        twist=[float(v) for v in baseline.twist_target_rad],
        R=parameterisation.radius_m, polar_cache="sg6043",
        n_blades=baseline.to_geometry().n_blades,
        r_hub=parameterisation.root_fraction * parameterisation.radius_m)
    aep_analytic, powers_analytic = aep_of_geometry(analytic_geometry, resource)
    aep_x0 = annual_energy_mwh(x0, resource=resource, parameterisation=parameterisation)
    aep_star = annual_energy_mwh(x_star, resource=resource, parameterisation=parameterisation)

    representation = aep_x0 - aep_analytic
    gain = aep_star - aep_x0

    # how far the optimiser moved, in the same units as the fit error
    chord_star, twist_star = parameterisation.evaluate(x_star)
    move = {
        "chord_max_abs_m": float(np.max(np.abs(chord_star - chord_fit))),
        "chord_rms_m": float(np.sqrt(np.mean((chord_star - chord_fit) ** 2))),
        "twist_max_abs_deg": float(np.max(np.abs(np.degrees(twist_star - twist_fit)))),
        "twist_rms_deg": float(np.sqrt(np.mean(np.degrees(twist_star - twist_fit) ** 2))),
    }

    summary = {
        "description": "Fit error between the analytic Schmitz chord/twist and x0's "
                       "spline projection (PROJECT_DIRECTION_v2 §7.3 point 4), and "
                       "the AEP of the analytic blade vs x0 vs A4's optimum.",
        "bounds_label": BOUNDS_LABEL,
        "command": "python verification/spline_fit_error/run_fit_error.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "parameterisation": {"n_chord": parameterisation.n_chord, "n_twist": parameterisation.n_twist,
                             "n_stations": parameterisation.n_stations, "degree": parameterisation.degree},
        "schmitz_design_point": {"design_reynolds": baseline.design_reynolds,
                                 "alpha_design_deg": math.degrees(baseline.alpha_design_rad),
                                 "design_cl": baseline.design_cl,
                                 "max_lift_to_drag": baseline.max_lift_to_drag},
        "x0_rebuild_max_abs_diff": x0_rebuild_diff,
        "radii_m": r.tolist(),
        "analytic_chord_m": baseline.chord_target_m.tolist(),
        "analytic_twist_deg": np.degrees(baseline.twist_target_rad).tolist(),
        "spline_chord_m": chord_fit.tolist(),
        "spline_twist_deg": np.degrees(twist_fit).tolist(),
        "optimum_chord_m": chord_star.tolist(),
        "optimum_twist_deg": np.degrees(twist_star).tolist(),
        "errors_at_stations": station_errors,
        "errors_on_fine_grid": {"n_points": N_FINE, "r_m": r_fine.tolist(), **fine_errors},
        "optimiser_move_x0_to_optimum": move,
        "aep_mwh_per_yr": {
            "analytic_schmitz": aep_analytic,
            "x0_spline": aep_x0,
            "a4_optimum": aep_star,
            "representation_difference_x0_minus_analytic": representation,
            "representation_difference_pct_of_analytic": 100.0 * representation / aep_analytic,
            "optimisation_gain_optimum_minus_x0": gain,
            "optimisation_gain_pct_of_x0": 100.0 * gain / aep_x0,
            "gain_over_representation_ratio": (gain / abs(representation)
                                               if representation else float("inf")),
            "analytic_all_converged": bool(powers_analytic["all_converged"]),
            "analytic_rated_power_w": float(powers_analytic["rated_power_w"]),
        },
    }
    with open(JSON_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)

    plot(summary, FIGURE_PATH)

    e, f = station_errors, fine_errors
    print(f"x0 rebuild max |diff| = {x0_rebuild_diff:.2e}")
    print(f"stations : chord max {e['chord_max_abs_m']*1e3:.3f} mm, rms {e['chord_rms_m']*1e3:.3f} mm "
          f"(max rel {100*e['chord_max_rel']:.2f} %); twist max {e['twist_max_abs_deg']:.4f} deg, "
          f"rms {e['twist_rms_deg']:.4f} deg")
    print(f"fine grid: chord max {f['chord_max_abs_m']*1e3:.3f} mm, rms {f['chord_rms_m']*1e3:.3f} mm "
          f"(max rel {100*f['chord_max_rel']:.2f} %); twist max {f['twist_max_abs_deg']:.4f} deg, "
          f"rms {f['twist_rms_deg']:.4f} deg")
    print(f"optimiser move x0 -> x*: chord max {move['chord_max_abs_m']*1e3:.1f} mm, rms "
          f"{move['chord_rms_m']*1e3:.1f} mm; twist max {move['twist_max_abs_deg']:.2f} deg, rms "
          f"{move['twist_rms_deg']:.2f} deg")
    a = summary["aep_mwh_per_yr"]
    print(f"AEP analytic {a['analytic_schmitz']:.6f}  x0 {a['x0_spline']:.6f}  optimum {a['a4_optimum']:.6f}")
    print(f"representation difference {representation:+.6f} MWh/yr "
          f"({a['representation_difference_pct_of_analytic']:+.4f} %); optimisation gain {gain:+.6f} "
          f"({a['optimisation_gain_pct_of_x0']:+.4f} %); ratio {a['gain_over_representation_ratio']:.1f}")
    print(f"wrote {JSON_PATH}\n      {FIGURE_PATH}")


def plot(summary, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = np.array(summary["radii_m"])
    fine = summary["errors_on_fine_grid"]
    r_f = np.array(fine["r_m"])
    st = summary["errors_at_stations"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(r, summary["analytic_chord_m"], "-", color="#888888", lw=1.6, label="analytic Schmitz")
    ax.plot(r, summary["spline_chord_m"], "--", color="#1f5fbf", lw=1.6, label="x0 spline (5 control points)")
    ax.plot(r, summary["optimum_chord_m"], ":", color="#d1495b", lw=1.6, label="A4 optimum")
    ax.set_ylabel("chord [m]")
    ax.set_title("chord", fontsize=11)

    ax = axes[0, 1]
    ax.plot(r, summary["analytic_twist_deg"], "-", color="#888888", lw=1.6, label="analytic Schmitz")
    ax.plot(r, summary["spline_twist_deg"], "--", color="#1f5fbf", lw=1.6, label="x0 spline (5 control points)")
    ax.plot(r, summary["optimum_twist_deg"], ":", color="#d1495b", lw=1.6, label="A4 optimum")
    ax.set_ylabel("twist [deg]")
    ax.set_title("twist", fontsize=11)

    ax = axes[1, 0]
    ax.plot(r_f, 1e3 * np.array(fine["chord_error_m"]), "-", color="#1f5fbf", lw=1.2, label="spline − analytic (fine grid)")
    ax.plot(r, 1e3 * np.array(st["chord_error_m"]), "o", color="#1f5fbf", ms=4, mfc="white", label="at the 25 stations")
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.set_ylabel("chord fit error [mm]")
    ax.set_title(f"chord error: max {1e3*st['chord_max_abs_m']:.2f} mm, rms {1e3*st['chord_rms_m']:.2f} mm at stations", fontsize=10)

    ax = axes[1, 1]
    ax.plot(r_f, fine["twist_error_deg"], "-", color="#1f5fbf", lw=1.2, label="spline − analytic (fine grid)")
    ax.plot(r, st["twist_error_deg"], "o", color="#1f5fbf", ms=4, mfc="white", label="at the 25 stations")
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.set_ylabel("twist fit error [deg]")
    ax.set_title(f"twist error: max {st['twist_max_abs_deg']:.3f}°, rms {st['twist_rms_deg']:.3f}° at stations", fontsize=10)

    for ax in axes.flat:
        ax.set_xlabel("radius r [m]")
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
