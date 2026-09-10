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

Run directly (from `src/`): python -m validation.compare_sg6043_uiuc

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from polars.cache_format import load_xfoil_band, SOURCE_LABELS
from validation.uiuc_sg6043 import load_lift_clean, load_drag_clean, NOMINAL_RE
from xfoil.xfoil_runner import RESULTS_DIR

_HERE = os.path.dirname(os.path.abspath(__file__))
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
    for re in NOMINAL_RE:
        path = os.path.join(CACHE_DIR, f"SG6043_Re{re}.csv")
        if not os.path.exists(path):
            continue
        table, source = load_xfoil_band(path)
        alpha = table[:, 0]

        fig, (ax_cl, ax_cd) = plt.subplots(1, 2, figsize=(11, 4.5))

        from polars.cache_format import SOURCE_XFOIL, SOURCE_GAP_RETRY, SOURCE_LOCAL_FIT
        colors = {SOURCE_XFOIL: "C0", SOURCE_GAP_RETRY: "C1", SOURCE_LOCAL_FIT: "C3"}
        for code, label in SOURCE_LABELS.items():
            if code not in colors:
                continue
            m = source == code
            if not m.any():
                continue
            ax_cl.plot(alpha[m], table[m, 1], "o", ms=4, color=colors[code], label=f"cache ({label})")
            ax_cd.plot(alpha[m], table[m, 2], "o", ms=4, color=colors[code], label=f"cache ({label})")
        ax_cl.plot(alpha, table[:, 1], "-", lw=0.8, color="C0", alpha=0.5)
        ax_cd.plot(alpha, table[:, 2], "-", lw=0.8, color="C0", alpha=0.5)

        ax_cl.plot(lift[re][:, 0], lift[re][:, 1], "ks", ms=4, mfc="none", label="UIUC LSAT")
        ax_cd.plot(drag[re][:, 0], drag[re][:, 2], "ks", ms=4, mfc="none", label="UIUC LSAT")

        ax_cl.set_xlabel("alpha (deg)"); ax_cl.set_ylabel("Cl")
        ax_cl.set_title(f"Cl vs alpha, Re={re:,}")
        ax_cl.legend(fontsize=7); ax_cl.grid(alpha=0.3)

        ax_cd.set_xlabel("alpha (deg)"); ax_cd.set_ylabel("Cd")
        ax_cd.set_title(f"Cd vs alpha, Re={re:,}")
        ax_cd.legend(fontsize=7); ax_cd.grid(alpha=0.3)

        fig.tight_layout()
        out_path = os.path.join(OUT_DIR, f"sg6043_vs_uiuc_Re{re}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  wrote {os.path.relpath(out_path, DATA_DIR)}")


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 78)
    print("SG6043 built cache vs UIUC LSAT Vol. 3 clean-model measurements")
    print("=" * 78)

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
