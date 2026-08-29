"""
Bilinear (alpha, Re) lookup over a cached XFOIL polar set.

Reads the CSVs produced by build_polar_cache.py (one file per Reynolds number,
columns alpha,cl,cd,cm) and exposes cl/cd/cm as a function of (alpha, Re) via
bilinear interpolation. This is what the BEM solver and adjoint FD verification
call instead of shelling out to XFOIL on every strip/iteration.

Airfoil-agnostic interface: BEM/adjoint code should call the module-level
get_polar(alpha, reynolds) and never reference an airfoil name directly. Which
airfoil's cache that resolves to is configured in exactly one place — the
ACTIVE_AIRFOIL constant below (or set_active_airfoil, e.g. for tests) — so
swapping the project's airfoil (NACA 4412 -> S809 -> a future DU-series section)
never requires a change on the BEM side.

Author: MJ Hendrikse
Project: DSP810S — Inverse Design of Small Wind Turbine Blades
"""

import glob
import os
import re as re_module

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from polars.cache_format import SOURCE_VITERNA, read_polar_csv

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "data"))

_FNAME_RE = re_module.compile(r"_Re(\d+)\.csv$", re_module.IGNORECASE)


class PolarCacheError(Exception):
    """Raised when a requested (alpha, Re) point falls outside the cached data."""


class PolarLookup:
    """
    Bilinear interpolator over a cached (alpha, Re) polar grid for one airfoil.

    The cache must be "rectangular": every Reynolds-number CSV must share the same
    alpha grid (this holds automatically if all cache CSVs came from the same
    build_polar_cache.py sweep, since alpha_min/alpha_max/alpha_step are shared).
    """

    def __init__(self, cache_dir):
        """
        Load every cached Reynolds-number CSV in cache_dir and build interpolators.

        XFOIL sometimes fails to converge one or more alphas near stall (it simply
        omits that row rather than erroring — see xfoil_runner.py), so different
        Reynolds numbers in a sweep can end up with different alpha coverage — at
        low Re this can include the sweep's extremes, not just interior points (as
        seen with S809 at Re=100k). To keep the (alpha, Re) grid rectangular for
        RegularGridInterpolator, each curve is reindexed onto the union of all
        alphas seen across the cache: alphas within that curve's own convergence
        range are filled by 1D linear interpolation (exact at points the curve
        already has); alphas outside that curve's own range are left as NaN rather
        than flat-extrapolated, so a query landing on one of those gaps surfaces
        as a NaN result and __call__ raises PolarCacheError instead of silently
        returning an extrapolated value.

        Since work order Task 2 the caches are written over the full -180..180
        deg circle with a provenance column, so a curve no longer *has* a
        ragged edge and no NaN survives this step -- but the reindexing is
        kept, both because the NACA 4412 reference cache predates the change
        and because "a hole here is a NaN, not a plausible number" is the
        behaviour the layer should keep if a future cache is written ragged.
        `xfoil_alpha_min`/`xfoil_alpha_max` report where the measured data
        ends and the extrapolation begins.

        Parameters
        ----------
        cache_dir : str
            Directory containing "<LABEL>_Re<value>.csv" files, e.g.
            data/polars/naca4412/.

        Raises
        ------
        FileNotFoundError
            If cache_dir contains no matching CSV files.
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
            bands.append((float(measured.min()), float(measured.max())))

            re_values.append(reynolds)
            curves.append((alpha, cl, cd, cm))

        # Common alpha axis = union of every alpha seen in any cached curve.
        alpha_ref = np.unique(np.concatenate([c[0] for c in curves]))

        order = np.argsort(re_values)
        self.re_values = np.array(re_values)[order]
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

        # The alpha band XFOIL actually converged at *every* cached Reynolds
        # number -- the intersection, not the union, so a caller that limits
        # itself to this range is on measured data whichever Re it lands on.
        # For a cache written before the provenance column existed this is
        # simply the whole table. Callers that must not step onto the Viterna
        # extrapolation (the QBlade export, the pyBEMT comparison tables, the
        # Stage 4 alpha clamp) read it from here rather than assuming -8..18.
        self.xfoil_alpha_min = max(b[0] for b in bands)
        self.xfoil_alpha_max = min(b[1] for b in bands)

        points = (self.re_values, self.alpha_values)
        self._cl_interp = RegularGridInterpolator(points, cl_grid)
        self._cd_interp = RegularGridInterpolator(points, cd_grid)
        self._cm_interp = RegularGridInterpolator(points, cm_grid)

    def __call__(self, alpha, reynolds):
        """
        Interpolate (cl, cd, cm) at the given angle of attack and Reynolds number.

        Parameters
        ----------
        alpha : float
            Angle of attack, degrees.
        reynolds : float
            Reynolds number.

        Returns
        -------
        tuple of float
            (cl, cd, cm) at the requested point.

        Raises
        ------
        PolarCacheError
            If alpha or reynolds falls outside the overall cached range, or if the
            point falls in a gap where one of the bracketing Reynolds curves never
            converged at that alpha (see __init__ — this is intentional: silently
            extrapolating a stalled/near-stall polar would produce plausible-
            looking but unvalidated numbers).
        """

        if not (self.alpha_values.min() <= alpha <= self.alpha_values.max()):
            raise PolarCacheError(
                f"alpha={alpha} deg is outside the cached range "
                f"[{self.alpha_values.min()}, {self.alpha_values.max()}] deg."
            )
        if not (self.re_values.min() <= reynolds <= self.re_values.max()):
            raise PolarCacheError(
                f"Re={reynolds:,.0f} is outside the cached range "
                f"[{self.re_values.min():,.0f}, {self.re_values.max():,.0f}]."
            )

        point = [(reynolds, alpha)]
        cl = float(self._cl_interp(point)[0])
        cd = float(self._cd_interp(point)[0])
        cm = float(self._cm_interp(point)[0])
        if np.isnan(cl) or np.isnan(cd) or np.isnan(cm):
            raise PolarCacheError(
                f"alpha={alpha} deg at Re={reynolds:,.0f} falls within the overall "
                f"cached range but at least one bracketing Reynolds curve never "
                f"converged at this alpha, so no reliable value is cached here."
            )
        return cl, cd, cm


# ----------------------------------------------------------------------------
# Airfoil-agnostic module-level interface
#
# This is the only part of the polar cache that BEM/adjoint code should import.
# The active airfoil is configured here, once, rather than by every caller
# supplying a cache directory or airfoil label — so a future airfoil swap is a
# one-line change in this module, not a search-and-replace across the BEM solver.
# ----------------------------------------------------------------------------

# Primary airfoil for the current phase of the project (see PROJECT_PLAN.md /
# the 2026-07-04 journal entry for why this is S809 rather than NACA 4412).
ACTIVE_AIRFOIL = "s809"

_lookup_cache = {}  # airfoil label -> PolarLookup, built lazily and reused


def set_active_airfoil(airfoil_label):
    """
    Change which cached airfoil get_polar() resolves to.

    Parameters
    ----------
    airfoil_label : str
        Subfolder name under data/polars/, e.g. "s809" or "naca4412".
    """

    global ACTIVE_AIRFOIL
    ACTIVE_AIRFOIL = airfoil_label


def get_polar(alpha, reynolds, airfoil=None):
    """
    Interpolate (cl, cd, cm) for the active (or explicitly given) airfoil.

    This is the airfoil-agnostic entry point BEM/adjoint code should use — it
    never needs to name an airfoil; the active one is resolved from
    ACTIVE_AIRFOIL (see set_active_airfoil to change it).

    Parameters
    ----------
    alpha : float
        Angle of attack, degrees.
    reynolds : float
        Reynolds number.
    airfoil : str or None, optional
        Override the active airfoil for this call only (subfolder name under
        data/polars/). Defaults to ACTIVE_AIRFOIL.

    Returns
    -------
    tuple of float
        (cl, cd, cm) at the requested point.

    Raises
    ------
    PolarCacheError
        If alpha or reynolds falls outside the cached range.
    FileNotFoundError
        If no cache exists yet for the resolved airfoil.
    """

    airfoil = airfoil or ACTIVE_AIRFOIL
    if airfoil not in _lookup_cache:
        cache_dir = os.path.join(DATA_DIR, "polars", airfoil)
        _lookup_cache[airfoil] = PolarLookup(cache_dir)

    return _lookup_cache[airfoil](alpha, reynolds)
