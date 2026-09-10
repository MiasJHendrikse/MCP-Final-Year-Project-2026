"""
`CachedPolar`: one airfoil polar at one fixed Reynolds number, adapting the
C1 interpolant to the BEM solver's duck-typed .cl(alpha)/.cd(alpha) interface.

Work order Task 4. This replaces `bem.airfoil.S809Polar`, which carried three
behaviours that are fine for a lookup and fatal for anything differentiable:

  * **alpha clamped** to the XFOIL-converged band. Past the edge the
    derivative was identically zero and the value was whatever the boundary
    happened to be -- at the Phase VI 15 m/s point that meant reading
    Cl(18 deg) ~ 1.0 for a station actually sitting at 32 deg, i.e. lift that
    was never measured or extrapolated, only invented by the clamp.
  * **Reynolds clamped** the same way. This one had already caused a real,
    weeks-long defect: before the cache was extended to 1.3M every Phase VI
    station clamped to the old 500k ceiling, so there was zero Reynolds
    variation in any run and nothing said so (see
    `xfoil.build_polar_cache`'s own comment).
  * **gap substitution** -- on a `PolarCacheError` the old adapter stepped
    alpha inward in 0.25 deg increments, up to 80 times, and returned the
    value it found there *unflagged*. A station solved on Cl from a
    completely different angle of attack was indistinguishable from a correct
    one, inside the residual.

None of it is needed any more. Task 2 made both caches gap-free over the full
-180..180 deg circle, so there is no hole to step around and no edge to clamp
at within the angles a blade element can physically reach; Task 3 turned that
grid into a C1 surface. What is left for this module to do is convert radians
to degrees, hold the station's Reynolds number, and get out of the way.

Out of range raises
--------------------
`CachedPolar` clamps nothing and substitutes nothing: an alpha or Reynolds
outside the built grid raises `PolarDomainError`. That is the point of the
task, not a rough edge of it. For the design rotor, alpha leaving the table is
information the Step 8 smoothness gate needs to see; a plausible number
returned in its place is the failure mode this whole layer exists to remove.

Reynolds is checked once, at construction, rather than on first evaluation.
A station whose Reynolds number the cache does not cover is a cache-coverage
problem -- the envelope study of Task 2 step 0 is what sets those bounds --
and it should say so with the station's number and the cache's range in the
message, not surface later as an opaque failure from inside a residual
evaluation a root-finder happens to be in the middle of.

Complex-step safety
--------------------
The radian-to-degree conversion is a plain multiply by a scalar, not
`math.degrees` (which rejects complex input), so a complex perturbation to
alpha propagates into the interpolant untouched -- see `polars.interpolant`'s
module docstring for why the evaluation path below it is complex-safe too.
The alpha-derivatives this module exposes are per *radian*, chain-ruled from
the interpolant's per-degree values, because radians are what the BEM
residual differentiates with respect to.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
import os
from functools import lru_cache

from config import REPO_ROOT, load_polar_cache
from polars.cache import PolarGrid
from polars.interpolant import PolarDomainError, PolarInterpolant

#: Radians -> degrees as a plain multiplier. `math.degrees` would be the
#: obvious call, but it coerces to float and so silently drops the imaginary
#: part of a complex-step perturbation -- exactly the dtype trap the work
#: order warns about for the Jacobian arrays, in scalar form.
_DEG_PER_RAD = 180.0 / math.pi


class CachedPolar:
    """
    An airfoil polar at one fixed Reynolds number, over a C1 cached surface.

    One instance per blade station, matching how `bem.rotor.solve_rotor`
    estimates a per-station Reynolds number. The underlying `PolarInterpolant`
    is shared between them -- it is immutable and its construction costs about
    15 ms, so it is built once per cache and handed to every station (see
    `CachedPolarFactory` and `interpolant_for`).

    Parameters
    ----------
    interpolant : polars.interpolant.PolarInterpolant
        The C1 (alpha, Reynolds) -> (Cl, Cd, Cm) surface to read.
    reynolds : float
        This station's Reynolds number. Must lie within the interpolant's
        built Reynolds range; see the module docstring on why that is checked
        here rather than on first use.

    Raises
    ------
    polars.interpolant.PolarDomainError
        If `reynolds` is outside the cache's Reynolds range.
    """

    def __init__(self, interpolant, reynolds):
        re_lo = float(interpolant.re_values[0])
        re_hi = float(interpolant.re_values[-1])
        reynolds = float(reynolds)
        if not (re_lo <= reynolds <= re_hi):
            raise PolarDomainError(
                f"Reynolds={reynolds:,.0f} is outside the cached range "
                f"[{re_lo:,.0f}, {re_hi:,.0f}]. The cache does not cover this "
                "station's operating point; extend the cache's Reynolds "
                "bounds rather than clamping to the edge (work order Task 2 "
                "step 0, Task 4)."
            )

        self._interpolant = interpolant
        self.reynolds = reynolds

    def cl(self, alpha):
        """Lift coefficient at angle of attack `alpha` (radians)."""

        return self._interpolant.cl(alpha * _DEG_PER_RAD, self.reynolds)

    def cd(self, alpha):
        """Drag coefficient at angle of attack `alpha` (radians)."""

        return self._interpolant.cd(alpha * _DEG_PER_RAD, self.reynolds)

    def cm(self, alpha):
        """Moment coefficient at angle of attack `alpha` (radians)."""

        return self._interpolant.cm(alpha * _DEG_PER_RAD, self.reynolds)

    def dcl_dalpha(self, alpha):
        """dCl/dalpha at `alpha` (radians), per radian."""

        return self._interpolant.dcl_dalpha(
            alpha * _DEG_PER_RAD, self.reynolds) * _DEG_PER_RAD

    def dcd_dalpha(self, alpha):
        """dCd/dalpha at `alpha` (radians), per radian."""

        return self._interpolant.dcd_dalpha(
            alpha * _DEG_PER_RAD, self.reynolds) * _DEG_PER_RAD

    def dcl_dre(self, alpha):
        """dCl/dRe at `alpha` (radians), per unit Reynolds number."""

        return self._interpolant.dcl_dre(alpha * _DEG_PER_RAD, self.reynolds)

    def dcd_dre(self, alpha):
        """dCd/dRe at `alpha` (radians), per unit Reynolds number."""

        return self._interpolant.dcd_dre(alpha * _DEG_PER_RAD, self.reynolds)

    def __call__(self, alpha):
        """(Cl, Cd) at `alpha` (radians) -- the pair the residual needs."""

        return self.cl(alpha), self.cd(alpha)

    def __repr__(self):
        return f"CachedPolar(reynolds={self.reynolds:,.0f})"


class CachedPolarFactory:
    """
    Turns one shared `PolarInterpolant` into per-station `CachedPolar`s.

    This is the object `bem.rotor.solve_rotor` injects: it estimates a
    Reynolds number per station and needs a polar for it, without knowing
    which cache is behind it or paying to rebuild the surface each time.

    Parameters
    ----------
    interpolant : polars.interpolant.PolarInterpolant
    """

    def __init__(self, interpolant):
        self.interpolant = interpolant

    def __call__(self, reynolds):
        return CachedPolar(self.interpolant, reynolds)


@lru_cache(maxsize=None)
def interpolant_for(cache_name):
    """
    The C1 interpolant over a named polar cache, built once and reused.

    Why this memo is not the `ACTIVE_AIRFOIL` global back again
    -----------------------------------------------------------
    Task 4 retires `xfoil.polar_lookup`'s module-level "active airfoil"
    because it made the answer depend on call order: two caches are live from
    Phase 1.4 onward, and whichever ran last decided what the next caller got.
    This is a different thing. It is a pure memo keyed on the argument -- a
    given `cache_name` always returns the same surface, nothing mutates it,
    and no caller can change what another caller sees. `PolarInterpolant` and
    `PolarGrid` are both read-only once built, so sharing one between every
    station of every rotor is safe as well as ~15 ms cheaper per rotor.

    The cache's location comes from `config/polars_<name>.yaml`, so no data
    path is written down anywhere under `src/` (the Task 1 rule).

    Parameters
    ----------
    cache_name : str
        A cache named by a `config/polars_<name>.yaml` file, e.g. "s809" or
        "sg6043". No default -- see `bem.rotor.RotorGeometry.polar_cache` for
        why the airfoil is a property of the blade, stated once where the
        blade is defined.

    Returns
    -------
    polars.interpolant.PolarInterpolant
    """

    cache = load_polar_cache(cache_name)
    directory = os.path.join(REPO_ROOT, cache.directory)
    return PolarInterpolant(PolarGrid(directory))


def polar_factory_for(cache_name):
    """`CachedPolarFactory` over `interpolant_for(cache_name)`."""

    return CachedPolarFactory(interpolant_for(cache_name))
