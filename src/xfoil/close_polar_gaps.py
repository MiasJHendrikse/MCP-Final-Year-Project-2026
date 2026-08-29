"""
Close the interior holes in a committed polar cache, and write it out over the
full -180..180 deg circle.

Work order Task 2, items 1 and 4. This is a *repair* pass over an existing
cache, not a rebuild (ground rule 1): every converged row already in the cache
is written back bit-for-bit, at the same `%.6f` it was written with, and the
`Ncrit = 5` / 240-panel settings the cache was built at are read from
`config/polars_<name>.yaml` rather than restated, so a retry cannot silently
run at different settings from the rows it is filling between.

Three things happen per Reynolds curve:

1. **Missing alphas are re-run through XFOIL.**
   `build_polar_cache.fill_alpha_gaps` re-runs the hole's neighbourhood at a
   quarter step, which usually walks the viscous solver through it. Rows
   recovered this way are marked `gap_retry` -- real viscous solutions, just
   not ones the main sweep reached.

2. **Anything still missing is filled from a documented local fit** -- a
   quadratic through the four nearest converged points, which is the same
   estimator `validate_polars.py` check 1 already uses to *judge* whether a
   hole is harmless. Rows filled this way are marked `local_fit`, so a filled
   point stays distinguishable from a converged one in the data itself and not
   only in prose.

3. **The curve is extended to +/-180 deg** by Viterna (`polars.viterna`),
   marked `viterna`.

Why holes matter enough to repair. `PolarLookup` reindexes every curve onto the
union alpha axis; an interior hole is bridged by linear interpolation between
its neighbours, so it does not raise -- it quietly returns a chord across a
region where XFOIL could not find a solution, at exactly the near-stall alphas
where the curve has the most structure. Task 3 fits a C1 interpolant over this
grid, and a bridged hole becomes a fitted feature.

Run directly (from `src/`):

    python -m xfoil.close_polar_gaps s809            # repair + extend
    python -m xfoil.close_polar_gaps s809 --dry-run  # report only, no writes
    python -m xfoil.close_polar_gaps s809 --no-xfoil # extend only, no retries

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

from config import load_polar_cache
from polars import cache_format
from xfoil.build_polar_cache import AIRFOILS, DATA_DIR, fill_alpha_gaps
from xfoil.xfoil_runner import RESULTS_DIR

_FNAME_RE = re_module.compile(r"_Re(\d+)\.csv$", re_module.IGNORECASE)

#: XFOIL retry budget. The retry sweeps one hole's neighbourhood at a quarter
#: step -- a handful of viscous solves, not a full polar -- so it needs the
#: main build's iteration count but nothing like its timeout.
_RETRY_N_ITER = 400
_RETRY_TIMEOUT = 300

#: Points either side used by the local quadratic fill, matching the estimator
#: validate_polars.py check 1 judges an isolated dropout with.
_LOCAL_FIT_POINTS = 4


def _expected_alphas(cache_config):
    n = int(round((cache_config.alpha_max_deg - cache_config.alpha_min_deg)
                  / cache_config.alpha_step_deg)) + 1
    return np.round(cache_config.alpha_min_deg
                    + cache_config.alpha_step_deg * np.arange(n), 6)


def _as_xfoil_columns(table):
    """
    Pad a cached (alpha, cl, cd, cm) table into run_xfoil_polar's column layout.

    `fill_alpha_gaps` merges rows in that layout and takes columns
    [alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr] from XFOIL. Only 0, 1, 2 and 4
    are ever read back out, so the unused columns are zero-filled rather than
    invented.
    """

    padded = np.zeros((len(table), 7))
    padded[:, 0] = table[:, 0]
    padded[:, 1] = table[:, 1]
    padded[:, 2] = table[:, 2]
    padded[:, 4] = table[:, 3]
    return padded


def _local_fit(table, alpha_gap):
    """
    Fill one alpha from a quadratic through its nearest converged neighbours.

    Returns (cl, cd, cm). Used only where XFOIL declined to converge twice --
    the value is a documented interpolation of the surrounding curve, not a
    solution, and is marked as such in the cache.
    """

    alpha = table[:, 0]
    near = np.argsort(np.abs(alpha - alpha_gap))[:_LOCAL_FIT_POINTS]
    fitted = []
    for column in (1, 2, 3):
        coeffs = np.polyfit(alpha[near], table[near, column], 2)
        fitted.append(float(np.polyval(coeffs, alpha_gap)))
    return tuple(fitted)


def close_gaps(airfoil, run_xfoil=True, dry_run=False, cd_max=None,
               refine=4, n_iter=_RETRY_N_ITER):
    """
    Repair and extend every Reynolds curve of one committed cache.

    Parameters
    ----------
    airfoil : str
        Cache name -- the subfolder under data/polars/ and the suffix of the
        `config/polars_<name>.yaml` that records how it was built.
    run_xfoil : bool, optional
        Attempt the XFOIL retry for missing alphas (default True). False skips
        straight to the local fit, for a machine without XFOIL.
    dry_run : bool, optional
        Report what would change and write nothing.
    refine : int, optional
        Factor by which the retry shrinks the alpha step inside a hole
        (default 4, `fill_alpha_gaps`'s own default). A larger factor is worth
        a second pass over a cache whose holes fell through to a local fit:
        an already-fitted row is discarded and re-attempted on every run, so
        a fill can only ever be upgraded to a converged solution, never the
        other way round.
    n_iter : int, optional
        Viscous iteration budget per alpha in the retry.
    cd_max : float or None, optional
        Post-stall drag maximum for the extension. None uses the Viterna
        module default (1.8, the value this cache and the committed QBlade
        files were produced with); a design-rotor cache should pass
        `polars.viterna.cd_max_finite_blade(AR)`.

    Returns
    -------
    list of dict
        One record per Reynolds curve: the alphas retried, recovered, fitted
        and still missing.
    """

    cache_config = load_polar_cache(airfoil)
    cache_dir = os.path.join(DATA_DIR, "polars", airfoil)
    expected = _expected_alphas(cache_config)
    airfoil_cmd = AIRFOILS[airfoil]["airfoil_cmd"] if airfoil in AIRFOILS else None

    viterna_kwargs = {} if cd_max is None else {"cd_max": cd_max}
    records = []

    for path in sorted(glob.glob(os.path.join(cache_dir, "*_Re*.csv"))):
        match = _FNAME_RE.search(os.path.basename(path))
        if not match:
            continue
        reynolds = int(match.group(1))

        table, source = cache_format.load_xfoil_band(path)

        # A previously fitted row is not data: drop it and let this pass try
        # XFOIL again. If the retry fails once more it is refilled below from
        # the same neighbours, so nothing is lost by re-attempting it.
        keep = source != cache_format.SOURCE_LOCAL_FIT
        previously_fitted = [float(a) for a in table[~keep, 0]]
        table, source = table[keep], source[keep]

        missing = np.setdiff1d(expected, np.round(table[:, 0], 6))
        record = {
            "reynolds": reynolds,
            "missing_before": [float(a) for a in missing],
            "previously_fitted": previously_fitted,
            "recovered_by_retry": [],
            "filled_by_local_fit": [],
        }

        if len(missing) and run_xfoil:
            if airfoil_cmd is None:
                raise KeyError(
                    f"'{airfoil}' is not in build_polar_cache.AIRFOILS, so the "
                    f"retry cannot know which airfoil to load into XFOIL"
                )
            print(f"  Re={reynolds:,}: retrying {[float(a) for a in missing]} "
                  f"at Ncrit={cache_config.ncrit:g}, {cache_config.n_panel} panels ...")
            padded, n_recovered = fill_alpha_gaps(
                _as_xfoil_columns(table), airfoil_cmd, reynolds,
                cache_config.alpha_min_deg, cache_config.alpha_max_deg,
                cache_config.alpha_step_deg,
                raw_path=os.path.join("_archive", f"_gapfix_{airfoil}_re{reynolds}.txt"),
                n_iter=n_iter, timeout=_RETRY_TIMEOUT,
                ncrit=cache_config.ncrit, n_panel=cache_config.n_panel,
                refine=refine,
            )
            if n_recovered:
                merged = padded[:, [0, 1, 2, 4]]
                recovered_mask = ~np.isin(np.round(merged[:, 0], 6),
                                          np.round(table[:, 0], 6))
                record["recovered_by_retry"] = [
                    float(a) for a in merged[recovered_mask, 0]]
                new_source = np.where(recovered_mask, cache_format.SOURCE_GAP_RETRY,
                                      cache_format.SOURCE_XFOIL)
                # Rows already in the cache keep the provenance they came with.
                new_source[~recovered_mask] = source
                table, source = merged, new_source
                missing = np.setdiff1d(expected, np.round(table[:, 0], 6))

        for alpha_gap in missing:
            cl, cd, cm = _local_fit(table, float(alpha_gap))
            insert = np.searchsorted(table[:, 0], alpha_gap)
            table = np.insert(table, insert, [alpha_gap, cl, cd, cm], axis=0)
            source = np.insert(source, insert, cache_format.SOURCE_LOCAL_FIT)
            record["filled_by_local_fit"].append(float(alpha_gap))
            how = "on retry" if run_xfoil else "and no retry was attempted"
            print(f"  Re={reynolds:,}: XFOIL did not converge alpha={alpha_gap:+.1f} "
                  f"{how}; filled from a local quadratic "
                  f"(Cl={cl:.4f}, Cd={cd:.5f}, Cm={cm:.4f})")

        still_missing = np.setdiff1d(expected, np.round(table[:, 0], 6))
        assert not len(still_missing), still_missing  # local fit is unconditional
        record["n_rows_band"] = int(len(table))
        # The band's provenance as it now stands, not just what this pass did:
        # a row recovered by an earlier pass is carried, and the record should
        # describe the cache rather than the run.
        record["provenance"] = {
            cache_format.SOURCE_LABELS[code]: int(np.count_nonzero(source == code))
            for code in sorted(set(source.tolist()))
        }

        if not dry_run:
            full, full_source = cache_format.extend_to_full_range(
                table, source, step=cache_config.alpha_step_deg, **viterna_kwargs)
            cache_format.write_polar_csv(full, full_source, path)
            record["n_rows_written"] = int(len(full))

        records.append(record)

    return records


def _summarise(records, airfoil, dry_run):
    print()
    print("=" * 68)
    retried = sum(len(r["recovered_by_retry"]) for r in records)
    fitted = sum(len(r["filled_by_local_fit"]) for r in records)
    holes = sum(len(r["missing_before"]) for r in records)
    print(f"{airfoil}: {holes} hole(s) found; {retried} closed by XFOIL retry, "
          f"{fitted} filled from a local fit.")
    for record in records:
        if record["missing_before"]:
            print(f"  Re={record['reynolds']:>9,}: "
                  f"missing {record['missing_before']} -> "
                  f"retry {record['recovered_by_retry']}, "
                  f"local fit {record['filled_by_local_fit']}")
    if dry_run:
        print("Dry run: no files written.")
        return

    out_path = os.path.join(RESULTS_DIR, "polar_cache",
                            f"{airfoil}_gap_closure.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(records, f, indent=2)
    print(f"Wrote {os.path.relpath(out_path, RESULTS_DIR)} under results/.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Close interior gaps in a polar cache and extend it to +/-180 deg.")
    parser.add_argument("airfoil", nargs="?", default="s809",
                        help="Cache name under data/polars/ (default: s809).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would change; write nothing.")
    parser.add_argument("--no-xfoil", action="store_true",
                        help="Skip the XFOIL retry and fill any hole locally.")
    parser.add_argument("--refine", type=int, default=4,
                        help="Alpha-step refinement factor for the retry (default 4).")
    parser.add_argument("--n-iter", type=int, default=_RETRY_N_ITER,
                        help=f"Viscous iterations per alpha in the retry "
                             f"(default {_RETRY_N_ITER}).")
    args = parser.parse_args()

    print("=" * 68)
    print(f"Closing polar-cache gaps: {args.airfoil}")
    print("=" * 68)
    result = close_gaps(args.airfoil, run_xfoil=not args.no_xfoil,
                        dry_run=args.dry_run, refine=args.refine,
                        n_iter=args.n_iter)
    _summarise(result, args.airfoil, args.dry_run)
    sys.exit(0)
