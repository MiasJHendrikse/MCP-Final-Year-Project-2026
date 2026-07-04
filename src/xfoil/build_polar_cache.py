"""
XFOIL polar cache builder.

Runs an alpha/Reynolds sweep for a single airfoil through the existing XFOIL wrapper
(xfoil_runner.run_xfoil_polar) and writes one clean CSV per Reynolds number into
data/polars/<airfoil>/. This is the cache the BEM solver and adjoint FD verification
read from, instead of shelling out to XFOIL on every call.

CSV schema (written by save_polar_csv): one header row "alpha,cl,cd,cm", one row per
converged alpha, columns:
    alpha : angle of attack, degrees
    cl    : lift coefficient
    cd    : total drag coefficient
    cm    : quarter-chord pitching-moment coefficient

Author: MJ Hendrikse
Project: DSP810S — Inverse Design of Small Wind Turbine Blades
"""

import argparse
import os

import numpy as np

from xfoil_runner import run_xfoil_polar, RESULTS_DIR

DATA_DIR = os.path.abspath(os.path.join(RESULTS_DIR, "..", "data"))

# Same Re/alpha sweep logic applied to every airfoil cached so far: Re = 100k-500k
# spans the expected operating range for a small turbine at a ~5-6 m/s mean-wind,
# low-Re site (Windhoek-area); alpha -8 to 18 deg covers attached flow plus enough
# post-stall range for the Viterna extrapolation (Phase 1 BEM) to have a real
# XFOIL-derived baseline near the stall boundary.
_DEFAULT_REYNOLDS_LIST = [100_000, 150_000, 200_000, 300_000, 400_000, 500_000]
_DEFAULT_ALPHA_MIN, _DEFAULT_ALPHA_MAX, _DEFAULT_ALPHA_STEP = -8, 18, 0.5

# Registry of airfoils this pipeline knows how to cache. NACA 4412 was the
# original primary target (see PROJECT_PLAN.md); the primary target has since
# moved to S809 (NREL Phase VI's actual airfoil, with real low-Re tunnel data from
# Delft/OSU/CSU) but the NACA 4412 cache is retained as a secondary reference
# dataset — the pipeline and interpolation layer are airfoil-agnostic, so there is
# no reason to discard validated, working cache data.
AIRFOILS = {
    "naca4412": dict(airfoil_cmd="NACA 4412", label="naca4412"),
    "s809": dict(
        airfoil_cmd="LOAD " + os.path.join(DATA_DIR, "airfoils", "s809.dat"),
        label="s809",
    ),
}


def save_polar_csv(polar, csv_path):
    """
    Write a converged XFOIL polar array to a CSV with columns alpha,cl,cd,cm.

    Parameters
    ----------
    polar : numpy.ndarray
        Polar array as returned by run_xfoil_polar, columns
        [alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr].
    csv_path : str
        Destination CSV path. Parent directory is created if missing.
    """

    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    subset = polar[:, [0, 1, 2, 4]]  # alpha, cl, cd, cm
    np.savetxt(csv_path, subset, delimiter=",", header="alpha,cl,cd,cm",
               comments="", fmt="%.6f")


def build_polar_cache(airfoil_cmd, airfoil_label, reynolds_list,
                       alpha_min, alpha_max, alpha_step,
                       n_iter=200, timeout=120):
    """
    Sweep Reynolds numbers for one airfoil and cache each converged polar as a CSV.

    Raw XFOIL output is written to results/_archive/ as scratch (mirroring the
    convention already used for ad hoc wrapper runs); the cleaned alpha/cl/cd/cm
    table is written to data/polars/<airfoil_label>/<AIRFOIL_LABEL>_Re<re>.csv.

    Parameters
    ----------
    airfoil_cmd : str
        XFOIL airfoil-load command, e.g. "NACA 4412".
    airfoil_label : str
        Filesystem-safe label for the airfoil, e.g. "naca4412". Used for both the
        cache subfolder and the CSV filename stem.
    reynolds_list : list of float
        Reynolds numbers to sweep.
    alpha_min, alpha_max, alpha_step : float
        Alpha sweep bounds and step, in degrees.
    n_iter : int, optional
        Max viscous-solver iterations per alpha (default 200).
    timeout : float, optional
        Max wall-clock seconds per XFOIL run before it is killed (default 120).

    Returns
    -------
    dict
        Maps Reynolds number -> output CSV path for successfully cached polars.
        Reynolds numbers where XFOIL produced no converged polar are omitted.
    """

    output_dir = os.path.join(DATA_DIR, "polars", airfoil_label)
    cached = {}

    for re in reynolds_list:
        print(f"  Running {airfoil_cmd} at Re = {re:,} ...")
        raw_path = os.path.join("_archive", f"_raw_{airfoil_label}_re{re}.txt")
        polar = run_xfoil_polar(
            airfoil_cmd=airfoil_cmd,
            reynolds=re,
            alpha_min=alpha_min,
            alpha_max=alpha_max,
            alpha_step=alpha_step,
            polar_path=raw_path,
            n_iter=n_iter,
            timeout=timeout,
        )

        if polar is None:
            print(f"  [WARN] No converged polar for {airfoil_cmd} at Re={re:,} — skipping cache entry.")
            continue

        csv_path = os.path.join(output_dir, f"{airfoil_label.upper()}_Re{re}.csv")
        save_polar_csv(polar, csv_path)
        cached[re] = csv_path
        print(f"    Converged points: {len(polar)} "
              f"(alpha {polar[:, 0].min():+.1f} to {polar[:, 0].max():+.1f} deg) "
              f"-> {os.path.relpath(csv_path, DATA_DIR)}")

    return cached


# ----------------------------------------------------------------------------
# Run when this file is executed directly: build the cache for one registered
# airfoil (default: s809, the current primary target).
# ----------------------------------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Build an XFOIL polar cache.")
    parser.add_argument("airfoil", nargs="?", default="s809", choices=sorted(AIRFOILS),
                         help="Which registered airfoil to cache (default: s809).")
    args = parser.parse_args()

    config = AIRFOILS[args.airfoil]

    print("=" * 60)
    print(f"Building XFOIL polar cache: {config['label']}")
    print("=" * 60)

    cached = build_polar_cache(
        airfoil_cmd=config["airfoil_cmd"],
        airfoil_label=config["label"],
        reynolds_list=_DEFAULT_REYNOLDS_LIST,
        alpha_min=_DEFAULT_ALPHA_MIN,
        alpha_max=_DEFAULT_ALPHA_MAX,
        alpha_step=_DEFAULT_ALPHA_STEP,
        n_iter=300,
    )

    print()
    print(f"Done. Cached {len(cached)}/{len(_DEFAULT_REYNOLDS_LIST)} Reynolds numbers "
          f"to data/polars/{config['label']}/.")
