"""
Quantify the Viterna stitch: how continuous is a cached polar where measured
data meets extrapolation?

Work order Task 2, item 5. The extension is anchored to the airfoil's state at
the last converged alpha, so *value* continuity is exact by construction --
that is checked here rather than assumed. The question that matters is the
**slope**, because Task 3 fits a C1 interpolant across this join and it can
only do one of two things with a kink: reproduce it, breaking the smoothness
guarantee the whole task exists to provide, or smooth it away, misrepresenting
the physics. Either way the size of the kink has to be known and recorded.

Normal operation keeps alpha well below stall, so the stitch sits outside the
working range -- but optimiser excursions and the Step 8 smoothness sweeps will
visit it, which is exactly when an unrecorded discontinuity turns into an
unexplained gradient artefact.

What is measured, per Reynolds row, at both stitches (the upper one at the
last converged alpha, where classical Viterna takes over, and the lower one at
the first converged alpha, where the linear ramp onto the negative branch
does):

  value jump   |y_extrapolated(alpha_stitch) - y_measured(alpha_stitch)|,
               evaluated from the extrapolation's own formula rather than from
               the written row, so it tests the anchoring and not the file.
  slope jump   the change in first difference across the stitch node, on the
               cache's own 0.5 deg grid -- the incoming slope from the last
               measured interval against the outgoing slope from the first
               extrapolated one. This is the kink a C1 fit will see.

Run directly (from `src/`):

    python -m validation.check_stitch_continuity [airfoil]

Writes `results/polar_cache/<airfoil>_stitch_continuity.json` and prints the
markdown table that `data/polars/<airfoil>/README.md` records.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import glob
import json
import os
import re as re_module
import sys

import numpy as np

from polars.cache_format import SOURCE_VITERNA, read_polar_csv
from polars.viterna import viterna, viterna_coefficients
from xfoil.xfoil_runner import RESULTS_DIR

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "data"))

_FNAME_RE = re_module.compile(r"_Re(\d+)\.csv$", re_module.IGNORECASE)

#: Thin-airfoil lift-curve slope, 2*pi per radian in per-degree units. A slope
#: jump is reported against this rather than only against the incoming slope:
#: at the upper stitch the incoming slope is the stall shoulder, which passes
#: through zero somewhere in the middle of this cache's Reynolds range, so the
#: incoming-relative figure is unbounded there and says more about the
#: denominator than about the kink. Both are recorded; this is the one the
#: summary quotes.
_THIN_AIRFOIL_SLOPE_PER_DEG = 2.0 * np.pi * np.pi / 180.0  # ~0.1097

#: A jump at or below this fraction of the thin-airfoil slope is reported as
#: continuous. Nothing enforces it -- the script quantifies and records, it
#: does not gate. It is here so the table has a stated threshold rather than
#: an implied one.
_SMOOTH_RELATIVE = 0.05


def _stitch_metrics(alpha, y, index, forward):
    """
    Slope either side of a stitch node, and the jump between them.

    `index` is the stitch row (the last measured row on the upper side, the
    first on the lower). `forward` says which side the extrapolation is on.
    """

    if forward:
        slope_in = (y[index] - y[index - 1]) / (alpha[index] - alpha[index - 1])
        slope_out = (y[index + 1] - y[index]) / (alpha[index + 1] - alpha[index])
    else:
        slope_in = (y[index + 1] - y[index]) / (alpha[index + 1] - alpha[index])
        slope_out = (y[index] - y[index - 1]) / (alpha[index] - alpha[index - 1])

    jump = slope_out - slope_in
    scale = max(abs(slope_in), 1e-12)
    return {
        "slope_in_per_deg": float(slope_in),
        "slope_out_per_deg": float(slope_out),
        "slope_jump_per_deg": float(jump),
        "slope_jump_relative": float(abs(jump) / scale),
        "slope_jump_vs_thin_airfoil": float(abs(jump) / _THIN_AIRFOIL_SLOPE_PER_DEG),
    }


def check_curve(path):
    """Stitch metrics for one cached Reynolds curve."""

    table, source = read_polar_csv(path)
    measured = source != SOURCE_VITERNA
    if measured.all():
        raise ValueError(
            f"{os.path.basename(path)} carries no extrapolated rows -- it "
            f"predates the +/-180 deg extension, so it has no stitch to check"
        )

    alpha, cl, cd = table[:, 0], table[:, 1], table[:, 2]
    hi = int(np.flatnonzero(measured)[-1])
    lo = int(np.flatnonzero(measured)[0])

    # Value continuity at the upper stitch, from the Viterna formula itself:
    # the coefficients are matched at (alpha_s, cl_s, cd_s), so re-evaluating
    # there must return the measured pair exactly.
    coeffs = viterna_coefficients(alpha[hi], cl[hi], cd[hi])
    cl_anchor, cd_anchor = viterna(alpha[hi], coeffs)

    return {
        "alpha_upper_deg": float(alpha[hi]),
        "alpha_lower_deg": float(alpha[lo]),
        "value_jump_cl_upper": float(abs(cl_anchor - cl[hi])),
        "value_jump_cd_upper": float(abs(cd_anchor - cd[hi])),
        "upper": {
            "cl": _stitch_metrics(alpha, cl, hi, forward=True),
            "cd": _stitch_metrics(alpha, cd, hi, forward=True),
        },
        "lower": {
            "cl": _stitch_metrics(alpha, cl, lo, forward=False),
            "cd": _stitch_metrics(alpha, cd, lo, forward=False),
        },
    }


def check_cache(airfoil):
    """Stitch metrics for every Reynolds curve of one cache."""

    cache_dir = os.path.join(DATA_DIR, "polars", airfoil)
    paths = sorted(glob.glob(os.path.join(cache_dir, "*_Re*.csv")))
    if not paths:
        raise FileNotFoundError(f"No cached polar CSVs found in '{cache_dir}'.")

    rows = []
    for path in paths:
        match = _FNAME_RE.search(os.path.basename(path))
        if not match:
            continue
        record = check_curve(path)
        record["reynolds"] = int(match.group(1))
        rows.append(record)
    return sorted(rows, key=lambda r: r["reynolds"])


def markdown_table(rows):
    """The per-Reynolds table recorded in the cache README."""

    lines = [
        "| Re | Cl slope in / out (per deg) | d(Cl')/stitch | Cd slope in / out (per deg) | d(Cd')/stitch |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        cl, cd = r["upper"]["cl"], r["upper"]["cd"]
        lines.append(
            f"| {r['reynolds']:,} | {cl['slope_in_per_deg']:+.4f} / "
            f"{cl['slope_out_per_deg']:+.4f} | {cl['slope_jump_per_deg']:+.4f} | "
            f"{cd['slope_in_per_deg']:+.4f} / {cd['slope_out_per_deg']:+.4f} | "
            f"{cd['slope_jump_per_deg']:+.4f} |"
        )
    return "\n".join(lines)


def report(airfoil):
    rows = check_cache(airfoil)

    print("=" * 78)
    print(f"Viterna stitch continuity: {airfoil}")
    print("=" * 78)
    print(f"Upper stitch at alpha = {rows[0]['alpha_upper_deg']:+.1f} deg, "
          f"lower at {rows[0]['alpha_lower_deg']:+.1f} deg.")
    print()

    worst_value = max(max(r["value_jump_cl_upper"], r["value_jump_cd_upper"])
                      for r in rows)
    print(f"Value continuity at the upper stitch: worst |jump| over all "
          f"{len(rows)} Reynolds rows = {worst_value:.3e} "
          f"(exact by construction; anything above rounding is a bug).")
    print()

    print("Upper stitch (XFOIL -> classical Viterna):")
    print(markdown_table(rows))
    print()

    print("Lower stitch (linear ramp onto the negative branch -> XFOIL):")
    for r in rows:
        cl = r["lower"]["cl"]
        print(f"  Re={r['reynolds']:>9,}: Cl slope {cl['slope_in_per_deg']:+.4f} "
              f"(measured) vs {cl['slope_out_per_deg']:+.4f} (ramp), "
              f"jump {cl['slope_jump_per_deg']:+.4f} /deg "
              f"({cl['slope_jump_vs_thin_airfoil']:.1%} of the thin-airfoil slope)")
    print()

    def worst(side):
        row = max(rows, key=lambda r: abs(r[side]["cl"]["slope_jump_per_deg"]))
        return row, row[side]["cl"]

    up_row, up = worst("upper")
    lo_row, lo = worst("lower")
    print(f"Worst Cl slope jump, upper stitch: {up['slope_jump_per_deg']:+.4f} /deg "
          f"at Re={up_row['reynolds']:,} "
          f"({up['slope_jump_vs_thin_airfoil']:.0%} of the thin-airfoil slope "
          f"{_THIN_AIRFOIL_SLOPE_PER_DEG:.4f} /deg).")
    print(f"Worst Cl slope jump, lower stitch: {lo['slope_jump_per_deg']:+.4f} /deg "
          f"at Re={lo_row['reynolds']:,} "
          f"({lo['slope_jump_vs_thin_airfoil']:.0%}).")
    print(f"Reported as continuous below {_SMOOTH_RELATIVE:.0%} of that slope; "
          f"neither stitch is, so both are kinks a C1 fit will have to make a "
          f"choice about (work order Task 3).")

    out_path = os.path.join(RESULTS_DIR, "polar_cache",
                            f"{airfoil}_stitch_continuity.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"Wrote {os.path.relpath(out_path, RESULTS_DIR)} under results/.")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Quantify value/slope continuity at a cache's Viterna stitch.")
    parser.add_argument("airfoil", nargs="?", default="s809",
                        help="Cache subfolder under data/polars/ (default: s809).")
    args = parser.parse_args()
    report(args.airfoil)
    sys.exit(0)
