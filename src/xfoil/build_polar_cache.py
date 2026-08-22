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

from xfoil.xfoil_runner import run_xfoil_polar, RESULTS_DIR

DATA_DIR = os.path.abspath(os.path.join(RESULTS_DIR, "..", "data"))

# Re = 100k-500k was the original range: it spans the expected operating range
# for a small turbine at a ~5-6 m/s mean-wind, low-Re site (Windhoek-area).
# That range does NOT cover the real NREL Phase VI rotor used from Stage 4
# onward: bem.rotor.phase_vi_geometry's chord-based station Reynolds numbers
# run from ~525k (5 m/s, inboard) to ~1.31M (25 m/s, outboard) across the
# Sequence S wind-speed sweep (see compare_pybemt.py) -- every station at
# every wind speed was clamping to this cache's old 500k ceiling (S809Polar's
# out-of-range clamp in bem/airfoil.py), so there was zero real Re variation
# in any Phase VI run. 600k-1.3M in 100k steps was added (2026-07-28) to
# bracket that actual range with the same resolution as the original sweep;
# see the 2026-07-28 journal entry for the before/after and the convergence
# check at the new, higher-Re end.
# alpha -8 to 18 deg covers attached flow plus enough post-stall range for the
# Viterna extrapolation (Phase 1 BEM) to have a real XFOIL-derived baseline
# near the stall boundary.
_DEFAULT_REYNOLDS_LIST = [
    100_000, 150_000, 200_000, 300_000, 400_000, 500_000,
    600_000, 700_000, 800_000, 900_000, 1_000_000, 1_100_000, 1_200_000, 1_300_000,
]
_DEFAULT_ALPHA_MIN, _DEFAULT_ALPHA_MAX, _DEFAULT_ALPHA_STEP = -8, 18, 0.5

# Transition/paneling settings. These are NOT XFOIL's defaults, and the reason is
# specific to S809 at the low end of the Re range above (see the 2026-07-26
# journal entry for the full comparison):
#
#   ncrit=5 — XFOIL's default Ncrit=9 models a clean wind tunnel. On a section as
#     thick as the S809 (21% t/c) at Re <= 150k that assumption keeps the boundary
#     layer laminar far enough back to form a long separation bubble, and when the
#     bubble bursts XFOIL settles onto a separated solution branch: the cached
#     polars built that way had cl collapsing from 0.44 to -0.12 between alpha 1.5
#     and 3 deg, and cl(alpha=0) coming out at 0.43 against a published ~0.15.
#     Ncrit=5 represents the freestream turbulence a turbine actually operates in
#     (atmospheric boundary layer, Ncrit ~ 4-7, not a tunnel), transitions the
#     layer earlier, and suppresses the bubble. It yields a polar family that is
#     monotonic in Re for cl(0), cl_max, cd_min and L/D_max, which Ncrit=9 was not.
#   n_panel=240 — resolves the leading-edge suction peak finely enough that the
#     remaining near-stall behaviour is smooth rather than panel-noise.
#   bidirectional=True — sweep outward from alpha=0 rather than straight through
#     from alpha_min, so the continuation starts from attached, unambiguous flow.
#     This is what recovered full -8..+18 deg coverage at every Re; the previous
#     single-pass sweep silently lost up to 6 alphas per curve.
_DEFAULT_NCRIT = 5.0
_DEFAULT_N_PANEL = 240

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


def _expected_alphas(alpha_min, alpha_max, alpha_step):
    """The alpha grid a sweep is supposed to produce, as a rounded array."""
    n = int(round((alpha_max - alpha_min) / alpha_step)) + 1
    return np.round(alpha_min + alpha_step * np.arange(n), 6)


def fill_alpha_gaps(polar, airfoil_cmd, reynolds, alpha_min, alpha_max, alpha_step,
                    raw_path, n_iter, timeout, ncrit, n_panel, refine=4):
    """
    Recover individual alphas the main sweep failed to converge.

    XFOIL drops a non-converged alpha silently, which leaves interior holes in an
    otherwise complete polar. Each hole is bracketed by converged neighbours, so
    re-running just that neighbourhood at a finer step usually walks the viscous
    solver through it: the continuation takes smaller steps across the awkward
    region and lands on the missing alpha from a nearby converged state. Only rows
    landing exactly on the original alpha grid are kept, so the cached curve keeps
    its uniform spacing.

    Parameters
    ----------
    polar : numpy.ndarray
        Polar from the main sweep, columns as returned by run_xfoil_polar.
    airfoil_cmd, reynolds, alpha_min, alpha_max, alpha_step : see build_polar_cache.
    raw_path : str
        Base path for the scratch polar files these retry runs write.
    n_iter, timeout, ncrit, n_panel : see build_polar_cache.
    refine : int, optional
        Factor by which to shrink the alpha step within a gap (default 4).

    Returns
    -------
    tuple of (numpy.ndarray, int)
        The polar with any recovered rows merged in and re-sorted, and the number
        of alphas recovered.
    """

    expected = _expected_alphas(alpha_min, alpha_max, alpha_step)
    have = np.round(polar[:, 0], 6)
    missing = np.setdiff1d(expected, have)
    if len(missing) == 0:
        return polar, 0

    # Group consecutive missing alphas so one retry run covers a whole hole.
    groups = []
    for a in missing:
        if groups and abs(a - groups[-1][-1] - alpha_step) < 1e-6:
            groups[-1].append(a)
        else:
            groups.append([a])

    recovered = []
    stem, ext = os.path.splitext(raw_path)
    for gi, group in enumerate(groups):
        lo = max(group[0] - alpha_step, alpha_min)
        hi = min(group[-1] + alpha_step, alpha_max)
        segment = run_xfoil_polar(
            airfoil_cmd=airfoil_cmd,
            reynolds=reynolds,
            alpha_min=lo,
            alpha_max=hi,
            alpha_step=alpha_step / refine,
            polar_path=f"{stem}_gap{gi}{ext}",
            n_iter=n_iter,
            timeout=timeout,
            ncrit=ncrit,
            n_panel=n_panel,
            bidirectional=False,
        )
        if segment is None:
            continue
        for row in segment:
            if any(abs(row[0] - a) < 1e-6 for a in group):
                recovered.append(row)

    if not recovered:
        return polar, 0

    merged = np.vstack([polar, np.array(recovered)])
    return merged[np.argsort(merged[:, 0])], len(recovered)


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
                       n_iter=200, timeout=120,
                       ncrit=_DEFAULT_NCRIT, n_panel=_DEFAULT_N_PANEL,
                       bidirectional=True):
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
    ncrit : float, optional
        e^N transition criterion. Defaults to _DEFAULT_NCRIT (5.0) — see the note
        on that constant for why this deliberately differs from XFOIL's 9.0.
    n_panel : int or None, optional
        Repanelling count, default _DEFAULT_N_PANEL (240). None uses XFOIL's
        default 160-panel PANE.
    bidirectional : bool, optional
        Sweep outward from alpha = 0 in both directions (default True) rather
        than straight through from alpha_min.

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
            ncrit=ncrit,
            n_panel=n_panel,
            bidirectional=bidirectional,
        )

        if polar is None:
            print(f"  [WARN] No converged polar for {airfoil_cmd} at Re={re:,} — skipping cache entry.")
            continue

        polar, n_recovered = fill_alpha_gaps(
            polar, airfoil_cmd, re, alpha_min, alpha_max, alpha_step,
            raw_path=raw_path, n_iter=n_iter, timeout=timeout,
            ncrit=ncrit, n_panel=n_panel,
        )
        if n_recovered:
            print(f"    Gap-fill recovered {n_recovered} alpha(s) the main sweep dropped.")

        still_missing = np.setdiff1d(
            _expected_alphas(alpha_min, alpha_max, alpha_step),
            np.round(polar[:, 0], 6),
        )
        if len(still_missing):
            print(f"    [WARN] {len(still_missing)} alpha(s) still unconverged: "
                  f"{[float(a) for a in still_missing]}")

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
        n_iter=400,
        # A 240-panel bidirectional sweep is a few hundred viscous solves; the
        # 120 s default is not enough and silently truncates the polar (XFOIL is
        # killed mid-sweep, so the run looks "converged up to alpha X" rather
        # than failed). 900 s per Reynolds number leaves ample headroom.
        timeout=900,
    )

    print()
    print(f"Done. Cached {len(cached)}/{len(_DEFAULT_REYNOLDS_LIST)} Reynolds numbers "
          f"to data/polars/{config['label']}/.")
