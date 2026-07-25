"""
Airfoil polars for the BEM solver: a synthetic linear polar (Stage 1) and a
real-data S809 adapter over the Phase 0 XFOIL cache (Stage 4).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
import os
from dataclasses import dataclass

from xfoil.polar_lookup import DATA_DIR, PolarCacheError, PolarLookup

_s809_lookup = None  # module-level PolarLookup, built once and reused by
                      # every S809Polar instance (one per station) rather
                      # than re-parsing the cached CSVs per station.


def _get_s809_lookup():
    global _s809_lookup
    if _s809_lookup is None:
        _s809_lookup = PolarLookup(os.path.join(DATA_DIR, "polars", "s809"))
    return _s809_lookup


@dataclass
class LinearPolar:
    """
    Cl(alpha) = cl_slope * alpha, clamped to +/- cl_max beyond +/- stall_alpha.
    Cd(alpha) = cd0, constant (no stall drag rise modelled at this stage).

    Parameters
    ----------
    cl_slope : float
        Lift-curve slope, per radian (thin-airfoil default: 2*pi).
    stall_alpha : float
        Angle of attack, radians, at which Cl clamps to cl_max.
    cd0 : float
        Constant drag coefficient.
    """

    cl_slope: float = 6.2831853071795865  # 2*pi
    stall_alpha: float = 0.20943951023931953  # 12 deg, radians
    cd0: float = 0.01

    @property
    def cl_max(self):
        return self.cl_slope * self.stall_alpha

    def cl(self, alpha):
        """Lift coefficient at angle of attack alpha (radians)."""
        if alpha > self.stall_alpha:
            return self.cl_max
        if alpha < -self.stall_alpha:
            return -self.cl_max
        return self.cl_slope * alpha

    def cd(self, alpha):
        """Drag coefficient at angle of attack alpha (radians)."""
        return self.cd0

    def __call__(self, alpha):
        return self.cl(alpha), self.cd(alpha)


class S809Polar:
    """
    Real S809 polar, adapting the Phase 0 XFOIL cache (xfoil.polar_lookup)
    to station.py's duck-typed .cl(alpha)/.cd(alpha) interface (radians in,
    fixed Reynolds number per instance -- one instance per station, since
    each station in the rotor loop has its own local Reynolds estimate).

    Reuses the Phase 0 cache/interpolation machinery as-is (PolarLookup,
    bilinear over the cached (alpha, Re) grid) rather than inventing a new
    lookup -- see xfoil/polar_lookup.py.

    Alpha clamping, and why it matters here specifically
    ------------------------------------------------------
    The cached S809 sweep only covers alpha in [-8, 18] deg. station.py's
    bracket scan (`_select_bracket`) evaluates the residual -- hence
    Cl/Cd -- at ~2000 trial phi values spanning nearly all of (0, pi/2)
    per station, purely to *locate* the physical root; almost all of those
    trial points imply alpha far outside any physically realistic range.
    With Stage 1's synthetic polar this was harmless (linear Cl, defined
    everywhere). With a real, range-limited polar it is not: letting
    PolarCacheError propagate would abort the scan before it ever reaches
    the one physically meaningful alpha. So out-of-range alpha is clamped
    to the nearest cached boundary (flat extrapolation) for the *scan*;
    this only affects samples far from the true root, which the residual
    at those samples is never trusted for anyway. Reynolds is clamped the
    same way for the same reason.

    PolarLookup can also raise `PolarCacheError` for an alpha *inside* the
    overall cached range, if the bracketing Reynolds curves didn't both
    converge there (a real gap in the XFOIL sweep near stall at low Re --
    see PolarLookup's own docstring; hit in practice while scanning trial
    phi for Stage 4's rotor geometry). Since the scan needs *some* finite
    value at every trial phi to keep going, not necessarily the physically
    correct one, this is handled the same way as the out-of-range case:
    step alpha inward toward 0 deg (bounded, `_GAP_FALLBACK_STEPS` steps of
    `_GAP_FALLBACK_STEP_DEG`) until a converged point is found. This never
    triggers at the eventual solution for any of the Stage 4 geometries
    actually solved (verified in validate_stage4.py); it only smooths over
    scan samples that were never going to be the physical root anyway.
    """

    _GAP_FALLBACK_STEP_DEG = 0.25
    _GAP_FALLBACK_STEPS = 80  # up to 20 deg inward -- generous vs. the ~8 deg
                               # of documented gap width at the worst-case Re

    def __init__(self, reynolds):
        self._lookup = _get_s809_lookup()
        self._alpha_min = float(self._lookup.alpha_values.min())
        self._alpha_max = float(self._lookup.alpha_values.max())
        self._re_min = float(self._lookup.re_values.min())
        self._re_max = float(self._lookup.re_values.max())
        self.reynolds = min(max(reynolds, self._re_min), self._re_max)

    def _query(self, alpha_rad):
        alpha_deg = math.degrees(alpha_rad)
        alpha_deg = min(max(alpha_deg, self._alpha_min), self._alpha_max)
        try:
            return self._lookup(alpha_deg, self.reynolds)
        except PolarCacheError:
            pass

        step = self._GAP_FALLBACK_STEP_DEG if alpha_deg > 0 else -self._GAP_FALLBACK_STEP_DEG
        for i in range(1, self._GAP_FALLBACK_STEPS + 1):
            probe = alpha_deg - step * i
            if not (self._alpha_min <= probe <= self._alpha_max):
                break
            try:
                return self._lookup(probe, self.reynolds)
            except PolarCacheError:
                continue

        raise PolarCacheError(
            f"no converged S809 data found within {self._GAP_FALLBACK_STEPS * self._GAP_FALLBACK_STEP_DEG} "
            f"deg of alpha={alpha_deg} at Re={self.reynolds:,.0f}"
        )

    def cl(self, alpha):
        """Lift coefficient at angle of attack alpha (radians)."""
        cl, _cd, _cm = self._query(alpha)
        return cl

    def cd(self, alpha):
        """Drag coefficient at angle of attack alpha (radians)."""
        _cl, cd, _cm = self._query(alpha)
        return cd

    def __call__(self, alpha):
        cl, cd, _cm = self._query(alpha)
        return cl, cd
