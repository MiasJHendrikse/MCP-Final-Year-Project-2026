"""
Load a cached polar set into a rectangular (Re, alpha) grid.

This is the CSV-load and ragged-to-rectangular reindexing logic from
`xfoil.polar_lookup.PolarLookup.__init__`, lifted out so the polar-*consuming*
layer (`interpolant.py`, and everything built on it) has no XFOIL dependency
-- see the package docstring. The reindexing algorithm itself is unchanged:
each curve is resampled by 1D linear interpolation onto the union of every
alpha seen anywhere in the cache, with NaN left wherever a curve's own
convergence range does not reach.

Since work order Task 2, every cache is written gap-free over the full
-180..180 deg circle with identical alpha coverage in every row, so in
practice no NaN survives this step for a Task-2-built cache. The reindexing
is kept anyway: it costs nothing when the grid is already rectangular, and it
is what turns a hole in a future or hand-edited cache into a NaN that a
consumer can detect, rather than a silently misaligned row.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import glob
import os
import re as re_module

import numpy as np

from polars.cache_format import SOURCE_VITERNA, read_polar_csv

_FNAME_RE = re_module.compile(r"_Re(\d+)\.csv$", re_module.IGNORECASE)


class PolarGridError(Exception):
    """Raised when a cache directory cannot be assembled into a valid grid."""


class PolarGrid:
    """
    A gap-checked, rectangular (Re, alpha) grid of cl/cd/cm for one airfoil.

    This is purely a data-loading and rectangularisation layer -- it exposes
    the grid arrays and the XFOIL-converged alpha band, and does no
    interpolation itself. `interpolant.PolarInterpolant` is what turns this
    into a callable, differentiable surface.

    Attributes
    ----------
    re_values : numpy.ndarray
        (M,) Reynolds numbers, ascending.
    alpha_values : numpy.ndarray
        (N,) angles of attack, degrees, ascending. Shared by every row.
    cl_grid, cd_grid, cm_grid : numpy.ndarray
        (M, N) coefficient grids. NaN where no cached curve reaches.
    xfoil_alpha_min, xfoil_alpha_max : float
        The alpha band every cached Reynolds number converged over -- the
        intersection, not the union, of each curve's own XFOIL-converged
        range (see `xfoil.polar_lookup.PolarLookup` for why this is the
        intersection). For a Task-2-built cache with no gaps this is simply
        the extrapolation boundary, which is the same at every Re.
    """

    def __init__(self, cache_dir):
        """
        Load every cached Reynolds-number CSV in `cache_dir` into a grid.

        Parameters
        ----------
        cache_dir : str
            Directory containing "<LABEL>_Re<value>.csv" files, e.g.
            data/polars/s809/.

        Raises
        ------
        FileNotFoundError
            If `cache_dir` contains no matching CSV files.
        """

        paths = sorted(glob.glob(os.path.join(cache_dir, "*_Re*.csv")))
        if not paths:
            raise FileNotFoundError(
                f"No cached polar CSVs found in '{cache_dir}'. "
                "Run build_polar_cache.py for this airfoil first."
            )

        re_values = []
        curves = []  # list of (alpha, cl, cd, cm) per Re, possibly ragged
        bands = []   # per-curve (alpha_min, alpha_max) of the XFOIL-converged rows

        for path in paths:
            match = _FNAME_RE.search(os.path.basename(path))
            if not match:
                continue  # skip any non-conforming file sitting in the folder
            reynolds = int(match.group(1))

            data, source = read_polar_csv(path)
            alpha, cl, cd, cm = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

            measured = alpha[source != SOURCE_VITERNA]
            if len(measured) == 0:
                raise PolarGridError(
                    f"'{path}' has no non-extrapolated rows -- nothing to "
                    "anchor the XFOIL-converged band to."
                )
            bands.append((float(measured.min()), float(measured.max())))

            re_values.append(reynolds)
            curves.append((alpha, cl, cd, cm))

        if len(curves) < 2:
            raise PolarGridError(
                f"'{cache_dir}' has only {len(curves)} Reynolds-number file(s); "
                "a Re-direction interpolant needs at least two."
            )

        # Common alpha axis = union of every alpha seen in any cached curve.
        alpha_ref = np.unique(np.concatenate([c[0] for c in curves]))

        order = np.argsort(re_values)
        self.re_values = np.array(re_values, dtype=float)[order]
        self.alpha_values = alpha_ref

        cl_grid = np.full((len(curves), len(alpha_ref)), np.nan)
        cd_grid = np.full((len(curves), len(alpha_ref)), np.nan)
        cm_grid = np.full((len(curves), len(alpha_ref)), np.nan)
        for row, i in enumerate(order):
            alpha, cl, cd, cm = curves[i]
            in_range = (alpha_ref >= alpha.min()) & (alpha_ref <= alpha.max())
            cl_grid[row, in_range] = np.interp(alpha_ref[in_range], alpha, cl)
            cd_grid[row, in_range] = np.interp(alpha_ref[in_range], alpha, cd)
            cm_grid[row, in_range] = np.interp(alpha_ref[in_range], alpha, cm)

        self.cl_grid = cl_grid
        self.cd_grid = cd_grid
        self.cm_grid = cm_grid

        self.xfoil_alpha_min = max(b[0] for b in bands)
        self.xfoil_alpha_max = min(b[1] for b in bands)

    @property
    def has_gaps(self):
        """True if the rectangularised grid contains any NaN."""

        return bool(
            np.isnan(self.cl_grid).any()
            or np.isnan(self.cd_grid).any()
            or np.isnan(self.cm_grid).any()
        )
