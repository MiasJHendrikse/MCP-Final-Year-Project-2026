"""
The site wind resource at 20 m hub height -- the formal record.

Run:

    python verification/wind_resource/run_extrapolation.py

Writes `wind_resource_20m.json` next to this file. That artefact, not this
script's stdout, is what `config/site.yaml` was filled from, and the two are
checked against each other by `tests/test_wind_resource.py`.

What this is
------------
The extraction (`gasp-point-data-50m.png`, beside this file; originally a
screenshot, `Screenshot 2026-09-13 114138.png`)
reports
the Weibull pair at **50 m**. The design rotor's hub height is **20 m**. This
script performs and records the height extrapolation, per the project log
(2026-09-10):

    "not at GWA's default 10 m or 50 m, and if only those are available the
     height extrapolation must be done and recorded, not fudged."

Three things are deliberately recorded that a bare calculation would not:

1. **The input, verbatim, with its own self-check.** The tool displayed
   A = 8.8, k = 1.87 and, separately, U = 7.8 m/s. Those three are not
   independent: U must equal A * Gamma(1 + 1/k). It does, to 0.017 m/s. That
   is the evidence the screenshot was read correctly -- a misread digit in
   either parameter would break the identity.

2. **An independent cross-check with a different physical basis.** The log
   law bounds the extrapolated scale across the plausible roughness range for
   this terrain. Justus & Mikhail must land inside it; if it did not, one of
   the two would be wrong and that is a finding, not a rounding difference.

3. **The comparison against my prior expectation** -- including where it
   disagrees. I expected k ~ 1.8-2.4, c ~ 6-7 m/s, for sanity-checking an
   extraction only. The extraction disagrees on both. The
   disagreement is reported here; it is NOT tuned away, and the prior band is
   NOT widened to admit it.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))

from objective.height_extrapolation import (  # noqa: E402
    extrapolate_weibull, log_law_scale_ms)

HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.join(HERE, "..", "baseline", "x0.json")

# ---------------------------------------------------------------------------
# The extraction, exactly as displayed. Nothing here is derived.
# ---------------------------------------------------------------------------
# Source: GASP point data, 2026-09-13 (`gasp-point-data-50m.png`).
#
# The tool is GASP ("GASP: Point Data"), NOT globalwindatlas.info, which is
# what the original plan named. That difference is real and is
# recorded rather than smoothed over -- see the README beside this file.
SOURCE = {
    "tool": "GASP -- Point Data",
    "screenshot": "verification/wind_resource/gasp-point-data-50m.png",
    "screenshot_as_supplied": "misc/Screenshot 2026-09-13 114138.png",
    "supplied_by": "MJ Hendrikse",
    "supplied_on": "2026-09-13",
    "extraction_kind": "single point, not an area average",
    "sector": "Total (omni-directional)",
    "displayed_coordinates": "S22.428, W16.558",
}

REFERENCE_HEIGHT_M = 50.0    # the GASP panel's selected height
REFERENCE_K = 1.87           # "Omni A / k: 8.8 / 1.87"
REFERENCE_A_MS = 8.8         # "Omni A / k: 8.8 / 1.87"
DISPLAYED_MEAN_MS = 7.8      # "Avg. wind speed (U): 7.8 m/s"

# The rotor's hub height, from config/site.yaml:hub_height_m. Hard-coded here
# rather than loaded, because site.yaml is this script's OUTPUT: loading it
# would make the artefact depend on the file it is used to fill.
TARGET_HEIGHT_M = 20.0

# Plausible aerodynamic roughness lengths for Khomas Hochland highland bush
# savanna -- open thornbush over broken terrain, as visible in the supplied
# aerial screenshots. A RANGE, deliberately. No roughness survey exists for
# this site, so no single value may be stated.
ROUGHNESS_RANGE_M = (0.05, 0.10, 0.20, 0.30, 0.50)

# Plan 1.3's prior expectation, for comparison only. Never a substitute.
PRIOR_K_RANGE = (1.8, 2.4)
PRIOR_C_RANGE = (6.0, 7.0)


def held_k_penalty(extrapolated):
    """
    The energy cost of pinning `k` at its 50 m value.

    Both annual energies are `objective.annual_energy_mwh` for the reference
    blade `verification/baseline/x0.json`, under the committed operating law
    and rating (`config/rotor_design.yaml`): one with the extrapolated
    (configured) pair, one with `k` held at the 50 m value and the
    extrapolated scale. The only thing that changes between the two is the
    shape parameter, so the difference is what the height extrapolation's `k`
    reduction buys -- the penalty a fixed-`k` resource assumption would pay.

    The configured pair is built from this run's own extrapolation rather
    than read from `config/site.yaml`: that file is this script's output, and
    an artefact that reads what it fills cannot be regenerated from scratch.
    `tests/test_wind_resource.py` is what holds the two equal.
    """

    from objective.objective import annual_energy_mwh
    from objective.weibull import WeibullResource

    with open(X0_PATH, encoding="utf-8") as handle:
        record = json.load(handle)
    x0 = np.array(record["chord_control_points_m"]
                  + record["twist_control_points_rad"], dtype=float)

    configured = WeibullResource(k=extrapolated.k, c=extrapolated.scale_ms)
    held = WeibullResource(k=REFERENCE_K, c=extrapolated.scale_ms)
    aep_configured = annual_energy_mwh(x0, resource=configured)
    aep_held = annual_energy_mwh(x0, resource=held)

    return {
        "what": (
            "AEP of the reference blade under the committed operating law "
            "and rating, with the resource at the extrapolated (configured) "
            "shape parameter and with k pinned at its 50 m value, c held at "
            "the extrapolated scale. The difference is the energy the k "
            "reduction buys."
        ),
        "blade": "verification/baseline/x0.json",
        "law": "the committed operating law and rating (config/rotor_design.yaml)",
        "k_configured": float(configured.k),
        "c_configured_ms": float(configured.c),
        "k_held": float(held.k),
        "c_held_ms": float(held.c),
        "aep_configured": float(aep_configured),
        "aep_k_held": float(aep_held),
        "penalty_pct": float(100.0 * (aep_held / aep_configured - 1.0)),
    }


def main():
    # --- 1. Is the input self-consistent? -------------------------------
    implied_mean = REFERENCE_A_MS * math.gamma(1.0 + 1.0 / REFERENCE_K)
    mean_residual = implied_mean - DISPLAYED_MEAN_MS

    # --- 2. The extrapolation -------------------------------------------
    result = extrapolate_weibull(
        k=REFERENCE_K,
        scale_ms=REFERENCE_A_MS,
        reference_height_m=REFERENCE_HEIGHT_M,
        target_height_m=TARGET_HEIGHT_M,
    )

    # --- 3. The independent cross-check ---------------------------------
    log_law = {
        f"z0={z0:g}m": log_law_scale_ms(
            REFERENCE_A_MS, REFERENCE_HEIGHT_M, TARGET_HEIGHT_M, z0)
        for z0 in ROUGHNESS_RANGE_M
    }
    log_law_min, log_law_max = min(log_law.values()), max(log_law.values())
    inside = log_law_min <= result.scale_ms <= log_law_max

    # --- 4. Against the prior expectation --------------------------------
    prior = {
        "k_range": list(PRIOR_K_RANGE),
        "c_range_ms": list(PRIOR_C_RANGE),
        "k_inside": PRIOR_K_RANGE[0] <= result.k <= PRIOR_K_RANGE[1],
        "c_inside": PRIOR_C_RANGE[0] <= result.scale_ms <= PRIOR_C_RANGE[1],
    }

    artefact = {
        "what": (
            "Weibull k and A for the Khomas Hochland site at the 20 m hub "
            "height, extrapolated from a 50 m extraction. This artefact is "
            "the source of config/site.yaml:wind_resource."
        ),
        "source": SOURCE,
        "reference_level": {
            "height_m": REFERENCE_HEIGHT_M,
            "k": REFERENCE_K,
            "scale_ms": REFERENCE_A_MS,
            "displayed_mean_speed_ms": DISPLAYED_MEAN_MS,
            "mean_speed_implied_by_k_and_A_ms": implied_mean,
            "mean_speed_residual_ms": mean_residual,
            "self_consistent": abs(mean_residual) < 0.05,
        },
        "extrapolation": result.as_dict(),
        "cross_check_log_law": {
            "what": (
                "A(20 m) under the logarithmic profile across plausible "
                "roughness lengths. Bounds the central case; supplies no "
                "number of its own."
            ),
            "roughness_lengths_m": list(ROUGHNESS_RANGE_M),
            "scale_ms": log_law,
            "band_ms": [log_law_min, log_law_max],
            "central_case_inside_band": inside,
        },
        "against_plan_prior_expectation": prior,
        "held_k_penalty": held_k_penalty(result),
    }

    path = os.path.join(HERE, "wind_resource_20m.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(artefact, f, indent=2)
        f.write("\n")

    print(f"Reference, {REFERENCE_HEIGHT_M:g} m:  "
          f"k = {REFERENCE_K}, A = {REFERENCE_A_MS} m/s")
    print(f"  displayed U = {DISPLAYED_MEAN_MS} m/s, "
          f"implied = {implied_mean:.4f} m/s, "
          f"residual = {mean_residual:+.4f} m/s")
    print()
    print(f"Justus & Mikhail exponent n = {result.shear_exponent:.6f}")
    print(f"Target, {TARGET_HEIGHT_M:g} m:      "
          f"k = {result.k:.6f}, A = {result.scale_ms:.6f} m/s")
    print(f"  V_bar = {result.mean_speed_ms:.6f} m/s")
    print()
    print(f"Log-law band at {TARGET_HEIGHT_M:g} m: "
          f"[{log_law_min:.4f}, {log_law_max:.4f}] m/s  "
          f"-> central case {'INSIDE' if inside else 'OUTSIDE'}")
    print()
    print(f"Plan prior k {PRIOR_K_RANGE}: "
          f"{'inside' if prior['k_inside'] else 'OUTSIDE -- a finding'}")
    print(f"Plan prior c {PRIOR_C_RANGE}: "
          f"{'inside' if prior['c_inside'] else 'OUTSIDE -- a finding'}")
    print()
    penalty = artefact["held_k_penalty"]
    print(f"Held-k penalty (reference blade, committed law and rating):")
    print(f"  AEP(configured k = {penalty['k_configured']:.6f}, "
          f"c = {penalty['c_configured_ms']:.6f}) = "
          f"{penalty['aep_configured']:.6f} MWh/yr")
    print(f"  AEP(k held at {penalty['k_held']:.6f})            = "
          f"{penalty['aep_k_held']:.6f} MWh/yr  "
          f"({penalty['penalty_pct']:+.4f} %)")
    print()
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
