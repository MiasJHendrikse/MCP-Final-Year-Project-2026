"""
The polar-cache CSV schema, and the one place its provenance column is defined.

A cache CSV is `alpha,cl,cd,cm,source`, one file per Reynolds number, ascending
in alpha. The first four columns are what they always were. The fifth is new
with work order Task 2 and says where each row came from:

    0  xfoil       converged in the main XFOIL sweep
    1  gap_retry   converged in the finer-step retry that closed a hole
                   (`xfoil.build_polar_cache.fill_alpha_gaps`) -- a real
                   viscous solution, just not one the first pass reached
    2  local_fit   NOT converged anywhere: filled from a documented local fit
                   through its converged neighbours, because the retry failed
    3  viterna     post-stall extrapolation (`polars.viterna`), outside the
                   XFOIL-converged band entirely

Why the column exists at all
-----------------------------
Task 2 requires two things that pull in opposite directions: the cache must
cover -180..180 with no holes (a C1 interpolant cannot be fitted over a grid
with a NaN in it, and a ragged edge is the same problem), and a filled or
extrapolated point must stay *distinguishable* from a converged measurement.
Writing the extrapolation into the same file without a provenance column would
satisfy the first and destroy the second: the alpha at which XFOIL stopped
converging would no longer be recoverable from the data, only from prose.

It is also load-bearing for reproducing existing artefacts. `export_qblade`
and the pyBEMT comparison tables must extrapolate from the *measured* band,
not from an already-extended table -- re-anchoring Viterna at alpha_s = 180
deg divides by cos(180 deg). They call `load_xfoil_band` for exactly that.

Legacy 4-column files (the NACA 4412 reference cache) load fine: every row is
reported as `xfoil`, which is what they are.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import os

import numpy as np

from polars.viterna import build_full_range_polar

#: Provenance codes, written as the fifth CSV column.
SOURCE_XFOIL = 0
SOURCE_GAP_RETRY = 1
SOURCE_LOCAL_FIT = 2
SOURCE_VITERNA = 3

SOURCE_LABELS = {
    SOURCE_XFOIL: "xfoil",
    SOURCE_GAP_RETRY: "gap_retry",
    SOURCE_LOCAL_FIT: "local_fit",
    SOURCE_VITERNA: "viterna",
}

#: The CSV header written by `write_polar_csv`.
HEADER = "alpha,cl,cd,cm,source"

#: Column format. alpha/cl/cd/cm keep the %.6f the cache has always used, so a
#: row that is rewritten unchanged is byte-identical to the row it replaces.
_FMT = ["%.6f", "%.6f", "%.6f", "%.6f", "%d"]


def read_polar_csv(path):
    """
    Load one cache CSV.

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        The (n, 4) alpha/cl/cd/cm table and the (n,) integer source column.
        A 4-column legacy file reports every row as `SOURCE_XFOIL`.
    """

    data = np.loadtxt(path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] >= 5:
        return data[:, :4], data[:, 4].astype(int)
    return data[:, :4], np.full(len(data), SOURCE_XFOIL, dtype=int)


def load_xfoil_band(path):
    """
    Load one cache CSV, keeping only rows inside the XFOIL-converged band.

    That is every row the extrapolation did not invent: converged points,
    gap-retry points and any documented local fill. This is what a caller
    wants whenever it is going to extrapolate, clamp, or validate the
    measured data -- see the module docstring.

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        The (m, 4) table and its (m,) source column.
    """

    data, source = read_polar_csv(path)
    keep = source != SOURCE_VITERNA
    return data[keep], source[keep]


def extend_to_full_range(table, source, step=0.5, **viterna_kwargs):
    """
    Extend one converged curve to -180..180 deg, carrying provenance with it.

    The extrapolated rows are marked `SOURCE_VITERNA`; the measured rows are
    copied back over the extension's own resampling of them, so a row that was
    already in the cache is written out bit-for-bit as it went in rather than
    round-tripped through an interpolation that happens to be exact.

    Parameters
    ----------
    table : numpy.ndarray
        (n, 4) converged rows only -- alpha, cl, cd, cm, ascending in alpha.
    source : array_like of int
        Provenance of those rows.
    step : float, optional
        Output alpha resolution, degrees (default 0.5). Every measured alpha
        must land on this grid; one that does not is an error rather than a
        silently dropped measurement.
    **viterna_kwargs
        Passed to `polars.viterna.build_full_range_polar` (`cd_max`,
        `negative_cl_scale`, `reversed_cd_floor`).

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        The (m, 4) full-range table and its (m,) source column.
    """

    table = np.asarray(table, dtype=float)
    source = np.asarray(source, dtype=int)

    full = build_full_range_polar(table[:, 0], table[:, 1], table[:, 2],
                                  table[:, 3], step=step, **viterna_kwargs)
    full_source = np.full(len(full), SOURCE_VITERNA, dtype=int)

    grid = full[:, 0]
    index = np.searchsorted(grid, table[:, 0])
    off_grid = np.abs(grid[index] - table[:, 0]) > 1e-6
    if np.any(off_grid):
        raise ValueError(
            f"measured alpha(s) {table[off_grid, 0].tolist()} do not land on "
            f"the {step} deg output grid; the extension would drop them"
        )
    full[index] = table
    full_source[index] = source

    return full, full_source


def write_polar_csv(table, source, csv_path):
    """
    Write an (n, 4) alpha/cl/cd/cm table plus its provenance column to CSV.

    Parameters
    ----------
    table : numpy.ndarray
        (n, 4) array, columns alpha, cl, cd, cm, ascending in alpha.
    source : array_like of int
        One provenance code per row (see the module constants).
    csv_path : str
        Destination path. Parent directory is created if missing.
    """

    table = np.asarray(table, dtype=float)
    source = np.asarray(source, dtype=int)
    if table.shape[1] != 4:
        raise ValueError("table must have 4 columns (alpha, cl, cd, cm)")
    if len(source) != len(table):
        raise ValueError("source must have one entry per row")

    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    out = np.column_stack([table, source])
    np.savetxt(csv_path, out, delimiter=",", header=HEADER, comments="", fmt=_FMT)
    return csv_path
