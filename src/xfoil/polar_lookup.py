"""
Bilinear (alpha, Re) lookup over a cached XFOIL polar set — RETIRED FROM THE
SOLVE PATH, kept only as the pre-remediation baseline.

Reads the CSVs produced by build_polar_cache.py (one file per Reynolds number,
columns alpha,cl,cd,cm) and exposes cl/cd/cm as a function of (alpha, Re) via
bilinear interpolation.

Nothing in the solver calls this any more. The interpolation core was
replaced with `polars.interpolant.PolarInterpolant` (C1, analytically
differentiable) because bilinear's alpha-derivative is piecewise constant with
a jump at every 0.5 deg knot, and the adapter above it was replaced with
`polars.polar.CachedPolar`. `tests/test_invariants.py` asserts that nothing
under `bem/` imports this package at all.

What this module is still for, and why it was not deleted with the rest:

  * `verification/polar_interpolant/generate_plots.py` draws the "before"
    curve of the committed staircase figure from it. That figure is the direct
    evidence for the interpolant's central claim, and it has to stay regenerable.
  * `tests/test_polar_cache.py` checks the gap-free guarantee through
    the same rectangularisation the audit measured, rather than through a
    re-derivation of it.

What didn't survive was removed: the module-level `ACTIVE_AIRFOIL` /
`set_active_airfoil()` pair, the `get_polar()` façade over them, and the
`_lookup_cache` dict. A single global "current airfoil" was reasonable with one
airfoil in the project; once two caches were in use at the same time it
made the answer depend on call order, a direct threat to the requirement
that the same vector via a different code path gives the same answer. The airfoil is
now a property of the blade (`bem.rotor.RotorGeometry.polar_cache`), stated
once where the blade is defined. Construct `PolarLookup(cache_dir)` explicitly.

Author: MJ Hendrikse
Project: MCP820S — Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
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

        The caches are written over the full -180..180
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
        # old alpha clamp) read it from here rather than assuming -8..18.
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
