"""
SG6043 n_crit sensitivity study.

`config/polars_sg6043.yaml` specifies the design-rotor polar cache but leaves
`build.ncrit` a raising sentinel: S809's Ncrit=5 was calibrated against that
section's own tunnel data at 21% thickness and low Re, and does not transfer
to a 10%-thick laminar-flow section (no cache inherits another cache's
calibration). This script runs that calibration for SG6043.

Method: XFOIL is swept at Re = 100k, 150k, 200k, 300k, 400k, 500k -- the
Reynolds numbers the UIUC LSAT Vol. 3 clean-model data actually covers
(data/sg6043_uiuc_lsat/) -- for each candidate n_crit in {5, 7, 9, 11, 13}.
Each XFOIL curve is compared against the measured Cl(alpha) (fine-resolution
lift file) and Cd(alpha) (drag file) at the same Re, and the n_crit
minimising the combined RMS error is selected.

This is a calibration study, not a cache build: nothing here writes to
data/polars/sg6043/. Its outputs are the evidence for the choice recorded in
config/polars_sg6043.yaml, written to results/ncrit_sensitivity/.

Run: python -m validation.ncrit_sensitivity   (from src/)

Author: MJ Hendrikse
Project: MCP820S — Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validation.uiuc_sg6043 import load_lift_clean, load_drag_clean, NOMINAL_RE
from xfoil.xfoil_runner import run_xfoil_polar, RESULTS_DIR

DATA_DIR = os.path.abspath(os.path.join(RESULTS_DIR, "..", "data"))
COORD_PATH = os.path.join(DATA_DIR, "airfoils", "sg6043.dat")
STUDY_DIR = os.path.join(RESULTS_DIR, "ncrit_sensitivity")
RAW_DIR = os.path.join(RESULTS_DIR, "_archive", "ncrit_sensitivity")

NCRIT_VALUES = [5.0, 7.0, 9.0, 11.0, 13.0]

# Same alpha range and panel count config/polars_sg6043.yaml specifies for the
# real cache build, so this study calibrates the settings that build will
# actually use rather than a different, easier sweep.
ALPHA_MIN, ALPHA_MAX, ALPHA_STEP = -8.0, 18.0, 0.5
N_PANEL = 240
N_ITER = 200
TIMEOUT = 400


def _run_one(re, ncrit):
    raw_path = os.path.join(
        "_archive", "ncrit_sensitivity", f"sg6043_re{re}_ncrit{int(ncrit)}.txt"
    )
    return run_xfoil_polar(
        airfoil_cmd=f"LOAD {COORD_PATH}",
        reynolds=re,
        alpha_min=ALPHA_MIN,
        alpha_max=ALPHA_MAX,
        alpha_step=ALPHA_STEP,
        polar_path=raw_path,
        n_iter=N_ITER,
        timeout=TIMEOUT,
        ncrit=ncrit,
        n_panel=N_PANEL,
        bidirectional=True,
    )


def _rmse_against(polar, measured, value_col):
    """RMS error between an XFOIL polar and measured (alpha, value) points.

    Only measured points inside the XFOIL polar's own converged alpha range
    are used -- comparing against an extrapolation would not be measuring
    n_crit, it would be measuring how far outside the sweep the point fell.
    Returns (rmse, n_points_compared) -- n_points_compared is reported
    alongside every number so a near-empty overlap cannot be mistaken for a
    good fit.
    """
    if polar is None or len(polar) < 2:
        return None, 0
    xfoil_alpha = polar[:, 0]
    xfoil_value = polar[:, value_col]
    lo, hi = xfoil_alpha.min(), xfoil_alpha.max()
    mask = (measured[:, 0] >= lo) & (measured[:, 0] <= hi)
    if mask.sum() == 0:
        return None, 0
    interp = np.interp(measured[mask, 0], xfoil_alpha, xfoil_value)
    err = interp - measured[mask, 1]
    return float(np.sqrt(np.mean(err ** 2))), int(mask.sum())


def run_study():
    lift = load_lift_clean()
    drag = load_drag_clean()

    os.makedirs(STUDY_DIR, exist_ok=True)
    os.makedirs(RAW_DIR, exist_ok=True)

    # polars[re][ncrit] = xfoil polar array or None
    polars = {re: {} for re in NOMINAL_RE}
    rows = []  # flat metric table, one row per (re, ncrit)

    for re in NOMINAL_RE:
        for ncrit in NCRIT_VALUES:
            print(f"Running SG6043 Re={re:,} ncrit={ncrit:.0f} ...")
            polar = _run_one(re, ncrit)
            polars[re][ncrit] = polar

            if polar is None:
                print(f"  [WARN] no converged polar")
                rows.append(dict(re=re, ncrit=ncrit, n_converged=0,
                                  cl_rmse=None, cl_n=0, cd_rmse=None, cd_n=0,
                                  cl_max_xfoil=None, cd_min_xfoil=None))
                continue

            cl_rmse, cl_n = _rmse_against(polar, lift[re][:, [0, 1]], value_col=1)
            cd_rmse, cd_n = _rmse_against(polar, drag[re][:, [0, 2]], value_col=2)

            print(f"  converged={len(polar):3d}  "
                  f"cl_rmse={cl_rmse if cl_rmse is None else f'{cl_rmse:.4f}'} (n={cl_n})  "
                  f"cd_rmse={cd_rmse if cd_rmse is None else f'{cd_rmse:.5f}'} (n={cd_n})")

            rows.append(dict(
                re=re, ncrit=ncrit, n_converged=len(polar),
                cl_rmse=cl_rmse, cl_n=cl_n, cd_rmse=cd_rmse, cd_n=cd_n,
                cl_max_xfoil=float(polar[:, 1].max()),
                cd_min_xfoil=float(polar[:, 2].min()),
            ))

    return rows, polars, lift, drag


def summarise(rows):
    """Aggregate per-ncrit scores across all six Reynolds numbers.

    Score is the mean Cl RMSE across Re where a comparison was possible, plus
    the mean Cd RMSE scaled by 10 (Cd is an order of magnitude smaller than
    Cl, so an unscaled sum would let Cl dominate the pick even where Cd
    disagreement is what actually distinguishes two n_crit values). Any
    (re, ncrit) with no converged polar, or no measured overlap, is a
    penalised outcome, not a missing one -- a candidate that fails to
    converge at some Re is disqualified by that failure, not scored on the
    Re where it happened to work.
    """
    per_ncrit = {}
    for ncrit in NCRIT_VALUES:
        sub = [r for r in rows if r["ncrit"] == ncrit]
        cl_vals = [r["cl_rmse"] for r in sub if r["cl_rmse"] is not None]
        cd_vals = [r["cd_rmse"] for r in sub if r["cd_rmse"] is not None]
        n_failed = sum(1 for r in sub if r["cl_rmse"] is None)
        per_ncrit[ncrit] = dict(
            mean_cl_rmse=float(np.mean(cl_vals)) if cl_vals else None,
            mean_cd_rmse=float(np.mean(cd_vals)) if cd_vals else None,
            n_re_failed=n_failed,
            n_re_ok=len(cl_vals),
        )
        p = per_ncrit[ncrit]
        if p["mean_cl_rmse"] is not None and n_failed == 0:
            p["score"] = p["mean_cl_rmse"] + 10.0 * (p["mean_cd_rmse"] or 0.0)
        else:
            p["score"] = float("inf")  # disqualified: incomplete coverage

    best_ncrit = min(per_ncrit, key=lambda k: per_ncrit[k]["score"])

    # Re-dependence check: does the same ncrit win at every individual Re?
    per_re_best = {}
    for re in NOMINAL_RE:
        sub = [r for r in rows if r["re"] == re and r["cl_rmse"] is not None]
        if not sub:
            per_re_best[re] = None
            continue
        best = min(sub, key=lambda r: r["cl_rmse"] + 10.0 * (r["cd_rmse"] or 0.0))
        per_re_best[re] = best["ncrit"]

    return dict(per_ncrit=per_ncrit, best_ncrit=best_ncrit, per_re_best=per_re_best)


def make_plots(polars, lift, drag):
    for re in NOMINAL_RE:
        fig, (ax_cl, ax_cd) = plt.subplots(1, 2, figsize=(11, 4.5))
        for ncrit in NCRIT_VALUES:
            polar = polars[re][ncrit]
            if polar is None:
                continue
            ax_cl.plot(polar[:, 0], polar[:, 1], label=f"n_crit={ncrit:.0f}", lw=1.2)
            ax_cd.plot(polar[:, 0], polar[:, 2], label=f"n_crit={ncrit:.0f}", lw=1.2)

        ax_cl.plot(lift[re][:, 0], lift[re][:, 1], "ko", ms=3, label="UIUC LSAT")
        ax_cd.plot(drag[re][:, 0], drag[re][:, 2], "ko", ms=3, label="UIUC LSAT")

        ax_cl.set_xlabel("alpha (deg)"); ax_cl.set_ylabel("Cl")
        ax_cl.set_title(f"Cl vs alpha, Re={re:,}")
        ax_cl.set_xlim(ALPHA_MIN, ALPHA_MAX)
        ax_cl.legend(fontsize=7); ax_cl.grid(alpha=0.3)

        ax_cd.set_xlabel("alpha (deg)"); ax_cd.set_ylabel("Cd")
        ax_cd.set_title(f"Cd vs alpha, Re={re:,}")
        ax_cd.set_xlim(ALPHA_MIN, ALPHA_MAX)
        ax_cd.legend(fontsize=7); ax_cd.grid(alpha=0.3)

        fig.tight_layout()
        out_path = os.path.join(STUDY_DIR, f"sg6043_ncrit_sweep_Re{re}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  wrote {os.path.relpath(out_path, DATA_DIR)}")

    # Aggregate: RMSE vs n_crit across all six Re, one line per metric.
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for re in NOMINAL_RE:
        cl_rmse_by_ncrit = []
        for ncrit in NCRIT_VALUES:
            polar = polars[re][ncrit]
            if polar is None:
                cl_rmse_by_ncrit.append(np.nan)
                continue
            rmse, _ = _rmse_against(polar, lift[re][:, [0, 1]], value_col=1)
            cl_rmse_by_ncrit.append(rmse if rmse is not None else np.nan)
        ax.plot(NCRIT_VALUES, cl_rmse_by_ncrit, "o-", label=f"Re={re:,}")
    ax.set_xlabel("n_crit"); ax.set_ylabel("Cl RMSE vs UIUC")
    ax.set_title("SG6043: Cl RMSE vs n_crit, by Reynolds number")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout()
    out_path = os.path.join(STUDY_DIR, "sg6043_ncrit_rmse_summary.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {os.path.relpath(out_path, DATA_DIR)}")


if __name__ == "__main__":
    print("=" * 70)
    print("SG6043 n_crit sensitivity study")
    print("=" * 70)

    rows, polars, lift, drag = run_study()
    summary = summarise(rows)

    print()
    print("Per-n_crit aggregate (mean over the 6 Reynolds numbers):")
    for ncrit, s in sorted(summary["per_ncrit"].items()):
        flag = "  <-- disqualified (incomplete coverage)" if s["score"] == float("inf") else ""
        print(f"  n_crit={ncrit:5.1f}  "
              f"mean_cl_rmse={s['mean_cl_rmse']}  mean_cd_rmse={s['mean_cd_rmse']}  "
              f"n_re_failed={s['n_re_failed']}{flag}")

    print()
    print(f"Selected n_crit = {summary['best_ncrit']:.0f}")
    print("Per-Reynolds best n_crit (Re-dependence check):")
    for re, ncrit in summary["per_re_best"].items():
        marker = "" if ncrit == summary["best_ncrit"] else "  <-- differs from overall pick"
        print(f"  Re={re:,}: n_crit={ncrit}{marker}")

    os.makedirs(STUDY_DIR, exist_ok=True)

    with open(os.path.join(STUDY_DIR, "ncrit_sensitivity_metrics.csv"), "w") as f:
        f.write("re,ncrit,n_converged,cl_rmse,cl_n,cd_rmse,cd_n,cl_max_xfoil,cd_min_xfoil\n")
        for r in rows:
            f.write(",".join(str(r[k]) for k in
                     ["re", "ncrit", "n_converged", "cl_rmse", "cl_n",
                      "cd_rmse", "cd_n", "cl_max_xfoil", "cd_min_xfoil"]) + "\n")

    with open(os.path.join(STUDY_DIR, "ncrit_sensitivity_summary.json"), "w") as f:
        json.dump(dict(
            alpha_range=[ALPHA_MIN, ALPHA_MAX, ALPHA_STEP],
            n_panel=N_PANEL,
            ncrit_values=NCRIT_VALUES,
            reynolds_values=NOMINAL_RE,
            per_ncrit=summary["per_ncrit"],
            selected_ncrit=summary["best_ncrit"],
            per_re_best_ncrit=summary["per_re_best"],
        ), f, indent=2)

    print()
    print("Writing comparison plots ...")
    make_plots(polars, lift, drag)

    print()
    print(f"Done. Metrics table, summary and plots written to "
          f"{os.path.relpath(STUDY_DIR, DATA_DIR)}/")
