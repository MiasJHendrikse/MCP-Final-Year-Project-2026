"""
Validate the built SG6043 polar cache against the UIUC LSAT Vol. 3 clean-model
measurements (data/sg6043_uiuc_lsat/), at the six Reynolds numbers the two
share (100k-500k).

This is a different check from `validate_polars.py`'s five physical-
plausibility checks, and a stronger one for this cache specifically: S809 has
no equivalent independent experimental reference at its own build settings, so
this is a genuine measured-vs-computed comparison, not just an internal
consistency check. It is also different from the n_crit sensitivity study
(`ncrit_sensitivity.py`): that study compared *fresh, un-extended* XFOIL runs
at candidate n_crit values to pick one; this script compares the *finished,
gap-closed cache* -- built at the calibrated n_crit=9, gap-filled where XFOIL
would not converge -- against the same measurements, as the final acceptance
check before calling the cache done (work order Task 2 / plan 1.2).

Only the XFOIL-converged band (`cache_format.load_xfoil_band`) is compared --
the +/-180 deg Viterna extension is checked for stitch continuity, not
against measurement, by `check_stitch_continuity.py`. Points inside the
compared range that came from a gap retry or a local fit are counted and
reported separately, so a low RMSE achieved partly across filled points (not
independent XFOIL solutions) is visible rather than hidden inside one number.

Run from the repo root (no XFOIL: both paths read the committed SG6043
cache CSVs under `data/polars/sg6043/` and the committed UIUC tables under
`data/sg6043_uiuc_lsat/`):

    python src/validation/compare_sg6043_uiuc.py           # compute the RMSEs, write the JSON, plot
    python src/validation/compare_sg6043_uiuc.py --replot  # redraw the six PNGs only, no JSON write

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, ".."))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from polars.cache_format import load_xfoil_band, SOURCE_LABELS  # noqa: E402
from validation.uiuc_sg6043 import load_lift_clean, load_drag_clean, NOMINAL_RE  # noqa: E402
from xfoil.xfoil_runner import RESULTS_DIR  # noqa: E402

DATA_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "data"))
CACHE_DIR = os.path.join(DATA_DIR, "polars", "sg6043")
OUT_DIR = os.path.join(RESULTS_DIR, "sg6043_uiuc_validation")


def _rmse_against(cache_alpha, cache_value, measured, value_col):
    """RMSE between the cache's own curve and measured (alpha, value) points.

    Only measured points inside the cache's converged alpha range are used --
    comparing outside it would measure extrapolation reach, not agreement.
    Returns (rmse, n_points_compared); n_points_compared travels with every
    number so a thin overlap cannot read as a strong result.
    """
    lo, hi = cache_alpha.min(), cache_alpha.max()
    mask = (measured[:, 0] >= lo) & (measured[:, 0] <= hi)
    if mask.sum() == 0:
        return None, 0
    interp = np.interp(measured[mask, 0], cache_alpha, cache_value)
    err = interp - measured[mask, value_col]
    return float(np.sqrt(np.mean(err ** 2))), int(mask.sum())


def _provenance_in_range(alpha, source, measured_alpha):
    """Count cache rows by provenance, restricted to the measured alpha span.

    Tells the difference between "RMSE is low because XFOIL agrees with the
    tunnel" and "RMSE is low partly because a filled row was nudged toward
    the neighbours it was fit from" -- a local-fit row is not independent
    evidence of agreement.
    """
    if len(measured_alpha) == 0:
        return {}
    lo, hi = measured_alpha.min(), measured_alpha.max()
    in_range = (alpha >= lo) & (alpha <= hi)
    counts = {}
    for code in sorted(set(source[in_range].tolist())):
        counts[SOURCE_LABELS[code]] = int(np.count_nonzero(source[in_range] == code))
    return counts


def compare():
    lift = load_lift_clean()
    drag = load_drag_clean()

    rows = []
    for re in NOMINAL_RE:
        path = os.path.join(CACHE_DIR, f"SG6043_Re{re}.csv")
        if not os.path.exists(path):
            print(f"Re={re:,}: [MISSING] no cache file at {path}")
            rows.append(dict(re=re, missing=True))
            continue

        table, source = load_xfoil_band(path)
        alpha = table[:, 0]

        cl_rmse, cl_n = _rmse_against(alpha, table[:, 1], lift[re][:, [0, 1]], 1)
        cd_rmse, cd_n = _rmse_against(alpha, table[:, 2], drag[re][:, [0, 2]], 1)

        cl_prov = _provenance_in_range(alpha, source, lift[re][:, 0])
        cd_prov = _provenance_in_range(alpha, source, drag[re][:, 0])

        print(f"Re={re:>7,}: cl_rmse={cl_rmse:.4f} (n={cl_n}, provenance={cl_prov})  "
              f"cd_rmse={cd_rmse:.5f} (n={cd_n}, provenance={cd_prov})"
              if cl_rmse is not None and cd_rmse is not None else
              f"Re={re:>7,}: [WARN] no overlap with measured data "
              f"(cl_rmse={cl_rmse}, cd_rmse={cd_rmse})")

        rows.append(dict(
            re=re, missing=False,
            cl_rmse=cl_rmse, cl_n=cl_n, cl_provenance_in_range=cl_prov,
            cd_rmse=cd_rmse, cd_n=cd_n, cd_provenance_in_range=cd_prov,
            cache_alpha_range=[float(alpha.min()), float(alpha.max())],
            cache_n_rows=int(len(alpha)),
        ))

    return rows, lift, drag


def make_plots(lift, drag):
    """The committed two-panel figures, one per shared Reynolds number.

    The Reynolds number is not drawn: it goes in the report caption. The
    figure is written through `plotting.figstyle.save`, which is also the
    check that no project history leaks into the image.
    """

    from plotting import figstyle

    figstyle.apply()

    for re in NOMINAL_RE:
        path = os.path.join(CACHE_DIR, f"SG6043_Re{re}.csv")
        if not os.path.exists(path):
            continue
        table, _source = load_xfoil_band(path)
        alpha = table[:, 0]

        fig, (ax_cl, ax_cd) = plt.subplots(1, 2, figsize=figstyle.DOUBLE)

        for ax, cache_values, measured, value_col, title, ylabel in (
            (ax_cl, table[:, 1], lift[re], 1, "Lift",
             figstyle.LABELS["lift_coefficient"]),
            (ax_cd, table[:, 2], drag[re], 2, "Drag",
             figstyle.LABELS["drag_coefficient"]),
        ):
            ax.plot(alpha, cache_values, "o-", ms=3.0, color="#1f3b8b",
                    label=r"XFOIL, $N_{crit}$ = 9")
            ax.plot(measured[:, 0], measured[:, value_col], "s", ms=4, mfc="none",
                    color="#b3452a", label="UIUC measurement")
            ax.set_xlabel(figstyle.LABELS["angle_of_attack"])
            ax.set_ylabel(ylabel)
            figstyle.title(ax, title)
            ax.legend()

        out_path = os.path.join(OUT_DIR, f"sg6043_vs_uiuc_Re{re}.png")
        figstyle.save(fig, out_path)
        plt.close(fig)
        print(f"  wrote {os.path.relpath(out_path, DATA_DIR)}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--replot", action="store_true",
                        help="redraw the six comparison PNGs from the committed "
                             "polar cache and UIUC files; no JSON write")
    args = parser.parse_args(argv)

    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 78)
    print("SG6043 built cache vs UIUC LSAT Vol. 3 clean-model measurements")
    print("=" * 78)

    if args.replot:
        print("--replot: reading the committed SG6043 cache CSVs and "
              "data/sg6043_uiuc_lsat/; no XFOIL, no JSON write.")
        lift = load_lift_clean()
        drag = load_drag_clean()
        print()
        print("Writing comparison plots ...")
        make_plots(lift, drag)
        return 0

    rows, lift, drag = compare()

    valid = [r for r in rows if not r["missing"] and r["cl_rmse"] is not None]
    if valid:
        mean_cl = float(np.mean([r["cl_rmse"] for r in valid]))
        mean_cd = float(np.mean([r["cd_rmse"] for r in valid]))
        print()
        print(f"Mean over {len(valid)} Reynolds numbers: "
              f"Cl RMSE={mean_cl:.4f}, Cd RMSE={mean_cd:.5f}")

    with open(os.path.join(OUT_DIR, "sg6043_uiuc_comparison.json"), "w") as f:
        json.dump(rows, f, indent=2)
    print(f"Wrote {os.path.relpath(os.path.join(OUT_DIR, 'sg6043_uiuc_comparison.json'), DATA_DIR)}")

    print()
    print("Writing comparison plots ...")
    make_plots(lift, drag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
