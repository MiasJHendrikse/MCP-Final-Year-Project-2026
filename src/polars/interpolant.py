"""
A C1, analytically differentiable (alpha, Reynolds) -> (Cl, Cd, Cm) surface.

Work order Task 3. This replaces querying `RegularGridInterpolator` with
`method='linear'` (the default `xfoil.polar_lookup.PolarLookup` uses), whose
alpha-derivative is piecewise constant with a jump at every 0.5 deg knot --
measured in the audit as exactly four distinct values of dCl/dalpha across
alpha in [4.0, 5.6] deg at Re=600k. That staircase is what a gradient-based
optimiser (Phase 2/3) or a smoothness sweep (plan Step 8) sees as noise.

Interpolation scheme, and why it is two different schemes glued together
--------------------------------------------------------------------------
The alpha axis is dense and uniform (0.5 deg steps, ~700+ points): a global
cubic spline is exactly the right tool there, and gives the C1 (in fact C2)
continuity the task asks for.

The Reynolds axis is the opposite shape: on the order of ten knots, unevenly
spaced, spanning a decade or more. A cubic spline fitted through points that
sparse can overshoot between them and invent a non-physical Re trend that
was never in the data -- Cd is not free to wiggle just because the spline
wants a smooth second derivative. A shape-preserving scheme (PCHIP) does not
overshoot, and is C1, which is the only continuity this task actually
requires in that direction.

So: cubic spline in alpha, PCHIP in log(Re) -- log rather than raw Re because
the knots span roughly 40k to 1.3M, and PCHIP's slope-limiting behaves the
same way at every decade only if the axis it operates on is itself linear in
log-space.

Bivariate construction
-----------------------
This is not two independent 1D interpolations bolted together (that does not
compose into a valid 2D C1 surface). Instead:

  1. For each Reynolds row, fit a not-a-knot cubic spline over the (shared)
     alpha grid, and take from it that row's alpha-DERIVATIVE at every alpha
     node. The values at the nodes are the cached data itself.
  2. Blend two channels per alpha node across log(Re) with PCHIP: the value,
     and that alpha-derivative. `PchipInterpolator` accepts multi-dimensional
     `y` with an `axis` argument, so all nodes are fit in one vectorised call.
  3. At query time, evaluate both channels at the query Reynolds number for
     the two alpha nodes bracketing the query alpha, and build the **cubic
     Hermite** on that interval.

  Adjacent alpha intervals share the value AND the slope at their common node,
  so C1 in alpha holds identically, at every Reynolds number, by construction.
  PCHIP keeps it C1 in log(Re). Differentiating the Hermite basis and the
  channel polynomials in closed form gives d/dalpha and d/d(log Re)
  analytically, with no finite differencing anywhere in the evaluation path.
  On a cache row at a cache node the scheme returns the cached number bitwise,
  because it interpolates the data rather than reconstructing it.

Why not blend the spline coefficients (the scheme this replaced)
-----------------------------------------------------------------
The original construction fitted a cubic spline in alpha per Reynolds row and
PCHIP-blended the resulting *coefficients* across log(Re). That is subtly but
definitely wrong, and plan step 1.7 caught it: continuity at an alpha knot is
a **linear** constraint relating the coefficient channels of the two adjacent
intervals, and PCHIP is **nonlinear** in its data -- its slope limiter is a
harmonic mean. Blending channel-by-channel therefore does not preserve the
constraint, and the surface acquired a jump at every alpha knot for any
Reynolds number *between* rows.

It hid because on a cache row the blended coefficients are the originals, so
continuity was exactly restored there -- and because complex-step verification
evaluates a single polynomial piece and is blind by construction to a mismatch
between adjacent pieces. Measured on SG6043 before the fix: Cl jumped 5.4e-3
across alpha = 6.0 deg at Re = 70,995, while dCl/dalpha stayed continuous. A
BEM station at that Reynolds number failed to converge, which is how it
surfaced -- the residual really was discontinuous.

`tests/test_polar_interpolant.py` now checks continuity across every alpha
knot at the geometric midpoint of every Reynolds interval, on both caches.

Complex-step safety
--------------------
Every step past the one-off `CubicSpline`/`PchipInterpolator` fits (which
only ever see real data) is plain scalar polynomial evaluation -- Horner's
method, implemented in this module with ordinary Python arithmetic rather
than NumPy array ops (see the perf note below). Python's numeric tower
promotes float + complex -> complex automatically and has no fixed-width
buffer to truncate into, which is what sidesteps the failure mode the work
order calls out explicitly: a NumPy array pre-allocated with a real dtype
(e.g. `np.zeros(n)`) silently discards the imaginary part of anything
assigned into it, capping observed complex-step accuracy near 1e-6 while
looking entirely plausible. There is no such buffer here -- every
intermediate is either a Python float, a Python complex, or a NumPy scalar
produced by ordinary arithmetic on those, so a complex perturbation to
`alpha` or `reynolds` propagates through untouched. See
`tests/test_polar_interpolant.py` for the resulting ~1e-14 agreement against
complex-step differentiation, and `verify_interpolant.py` under
`verification/polar_interpolant/` for the numbers measured against the real
S809 cache.

Why plain Python arithmetic, not vectorised NumPy, inside `evaluate`
----------------------------------------------------------------------
Measured (see the Task 3 journal entry): a version of this evaluator written
with NumPy array slicing and elementwise ops costs roughly 30 microseconds
per call, dominated by NumPy's per-op dispatch overhead on 4-element arrays,
not by the arithmetic itself. Rewritten with plain Python floats/lists and
`bisect` for the interval search, the same evaluation costs roughly 2
microseconds -- an order of magnitude under the task's 20 microsecond target.
NumPy is still exactly what the one-off construction (`CubicSpline`,
`PchipInterpolator`) is built from; only the hot per-query path avoids it.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import bisect
import cmath
import math

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator


class PolarDomainError(Exception):
    """Raised when a query point falls outside the interpolant's built domain."""


def _horner(coeffs, x):
    """
    Evaluate a polynomial given highest-power-first coefficients.

    `coeffs` is a plain Python sequence, e.g. [c0, c1, c2, c3] representing
    c0*x**3 + c1*x**2 + c2*x + c3 -- the convention scipy's PPoly-family
    classes (`CubicSpline`, `PchipInterpolator`) use for `.c`.
    """

    acc = coeffs[0]
    for c in coeffs[1:]:
        acc = acc * x + c
    return acc


def _deriv_coeffs(coeffs):
    """
    Highest-power-first coefficients of d/dx of the polynomial `coeffs`.

    [c0, c1, c2, c3] (cubic) -> [3*c0, 2*c1, c2] (quadratic), i.e. the usual
    power-rule differentiation of a polynomial written in this convention.
    """

    n = len(coeffs)
    return [coeffs[i] * (n - 1 - i) for i in range(n - 1)]


def _log(x):
    """`math.log`, or `cmath.log` if `x` is complex-typed (complex-step)."""

    if isinstance(x, (complex, np.complexfloating)):
        return cmath.log(x)
    return math.log(x)


def _bracket(sorted_values, x, label):
    """
    Locate the half-open interval `[sorted_values[i], sorted_values[i+1])`
    containing `x`'s real part, clipped so the upper edge of the domain
    resolves to the last interval rather than raising a boundary indexing
    error.

    Parameters
    ----------
    sorted_values : list of float
        Ascending breakpoints.
    x : float or complex
        Query value. Only `x.real` decides which interval is used -- an
        infinitesimal complex-step perturbation must not change the branch,
        or the "closed-form derivative" and "complex-step derivative" would
        silently be evaluating two different polynomial pieces.
    label : str
        Name used in the out-of-domain error message.

    Returns
    -------
    int
        Interval index i such that sorted_values[i] <= x.real <= sorted_values[i+1]
        (the upper bound of the whole range included).
    """

    xr = x.real
    lo, hi = sorted_values[0], sorted_values[-1]
    if xr < lo or xr > hi:
        raise PolarDomainError(
            f"{label}={xr!r} is outside the interpolant's built range "
            f"[{lo!r}, {hi!r}]."
        )
    i = bisect.bisect_right(sorted_values, xr) - 1
    return min(max(i, 0), len(sorted_values) - 2)


class _CubicPchipSurface:
    """
    One scalar field (Cl, Cd or Cm) as a C1 function of (alpha, Reynolds).

    See the module docstring for the construction and evaluation scheme.
    Not part of the public API -- `PolarInterpolant` composes three of these,
    one per coefficient.
    """

    def __init__(self, alpha_values, re_values, field_grid):
        if np.isnan(field_grid).any():
            raise PolarDomainError(
                "cannot build a C1 interpolant over a grid containing NaN -- "
                "the gap-free, +/-180 deg extension (Task 2) is a "
                "precondition of this module, not something it enforces "
                "itself."
            )

        alpha_values = np.asarray(alpha_values, dtype=float)
        re_values = np.asarray(re_values, dtype=float)
        field_grid = np.asarray(field_grid, dtype=float)

        # Per Reynolds row: a not-a-knot cubic spline in alpha, used only to
        # obtain that row's alpha-DERIVATIVE at each node. The values are the
        # cached data itself, so nothing is smoothed on the way in.
        rows = [CubicSpline(alpha_values, field_grid[m], bc_type="not-a-knot")
                for m in range(len(re_values))]
        node_values = field_grid
        node_slopes = np.stack([row(alpha_values, 1) for row in rows], axis=0)

        logre = np.log(re_values)
        value_blend = PchipInterpolator(logre, node_values, axis=0)
        slope_blend = PchipInterpolator(logre, node_slopes, axis=0)

        # `.c` has shape (4, M-1, N): [poly degree, Re interval, alpha node].
        # Converted to nested plain lists once, so the hot per-query path never
        # pays NumPy's per-element dispatch overhead (module docstring).
        self._alpha = alpha_values.tolist()
        self._re = re_values.tolist()
        self._re_breaks = value_blend.x.tolist()  # == logre, ascending
        self._value_coeffs = value_blend.c.tolist()
        self._slope_coeffs = slope_blend.c.tolist()

    def evaluate(self, alpha, reynolds):
        """
        Return (value, d(value)/d(alpha), d(value)/d(reynolds)) at one point.

        `alpha` and `reynolds` may be Python/NumPy floats, or complex (for
        complex-step verification) -- see the module docstring's dtype note.

        Raises
        ------
        PolarDomainError
            If `alpha` or `reynolds` falls outside the grid this surface was
            built from.
        """

        m0 = _bracket(self._re, reynolds, "reynolds")
        k0 = _bracket(self._alpha, alpha, "alpha")

        logre = _log(reynolds)
        u = logre - self._re_breaks[m0]

        # Blend the two bracketing alpha nodes' value and alpha-slope across
        # log(Re). Four scalar channels, each a cubic in u.
        values = self._value_coeffs
        slopes = self._slope_coeffs

        v0_col = [values[p][m0][k0] for p in range(4)]
        v1_col = [values[p][m0][k0 + 1] for p in range(4)]
        s0_col = [slopes[p][m0][k0] for p in range(4)]
        s1_col = [slopes[p][m0][k0 + 1] for p in range(4)]

        v0 = _horner(v0_col, u)
        v1 = _horner(v1_col, u)
        s0 = _horner(s0_col, u)
        s1 = _horner(s1_col, u)

        dv0 = _horner(_deriv_coeffs(v0_col), u)
        dv1 = _horner(_deriv_coeffs(v1_col), u)
        ds0 = _horner(_deriv_coeffs(s0_col), u)
        ds1 = _horner(_deriv_coeffs(s1_col), u)

        # Cubic Hermite across the alpha interval, in the normalised
        # coordinate t = (alpha - a_k) / h.
        a_lo = self._alpha[k0]
        h = self._alpha[k0 + 1] - a_lo
        t = (alpha - a_lo) / h
        t2 = t * t
        t3 = t2 * t

        # Basis and its t-derivative.
        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2

        g00 = 6.0 * t2 - 6.0 * t
        g10 = 3.0 * t2 - 4.0 * t + 1.0
        g01 = -6.0 * t2 + 6.0 * t
        g11 = 3.0 * t2 - 2.0 * t

        m0_scaled = s0 * h
        m1_scaled = s1 * h

        value = h00 * v0 + h10 * m0_scaled + h01 * v1 + h11 * m1_scaled
        dvalue_dalpha = (g00 * v0 + g10 * m0_scaled
                         + g01 * v1 + g11 * m1_scaled) / h
        dvalue_dlogre = (h00 * dv0 + h10 * ds0 * h
                         + h01 * dv1 + h11 * ds1 * h)
        dvalue_dre = dvalue_dlogre / reynolds

        return value, dvalue_dalpha, dvalue_dre


class PolarInterpolant:
    """
    A C1, analytically differentiable (alpha, Reynolds) -> (Cl, Cd, Cm)
    surface built from a `polars.cache.PolarGrid`.

    Exposes value accessors (`cl`, `cd`, `cm`) and the four analytic partials
    Task 3 asks for (`dcl_dalpha`, `dcl_dre`, `dcd_dalpha`, `dcd_dre`). Cm has
    no partial exposed -- nothing downstream differentiates it (the BEM
    residual and its Jacobian need only Cl and Cd; Cm is carried through for
    the moment/pitching calculations that consume it as a value).

    Parameters
    ----------
    grid : polars.cache.PolarGrid
        A rectangular, gap-free polar grid (see `polars.cache.PolarGrid`).

    Raises
    ------
    PolarDomainError
        If `grid` still contains a NaN gap.
    """

    def __init__(self, grid):
        if grid.has_gaps:
            raise PolarDomainError(
                "cache grid has NaN gaps; a C1 interpolant needs the "
                "gap-free, +/-180 deg extension from Task 2 first."
            )

        self.alpha_values = grid.alpha_values
        self.re_values = grid.re_values
        self.xfoil_alpha_min = grid.xfoil_alpha_min
        self.xfoil_alpha_max = grid.xfoil_alpha_max

        self._cl = _CubicPchipSurface(grid.alpha_values, grid.re_values, grid.cl_grid)
        self._cd = _CubicPchipSurface(grid.alpha_values, grid.re_values, grid.cd_grid)
        self._cm = _CubicPchipSurface(grid.alpha_values, grid.re_values, grid.cm_grid)

    def __call__(self, alpha, reynolds):
        """(Cl, Cd, Cm) at one (alpha, Reynolds) point -- values only."""

        cl, _, _ = self._cl.evaluate(alpha, reynolds)
        cd, _, _ = self._cd.evaluate(alpha, reynolds)
        cm, _, _ = self._cm.evaluate(alpha, reynolds)
        return cl, cd, cm

    def cl(self, alpha, reynolds):
        value, _, _ = self._cl.evaluate(alpha, reynolds)
        return value

    def dcl_dalpha(self, alpha, reynolds):
        _, derivative, _ = self._cl.evaluate(alpha, reynolds)
        return derivative

    def dcl_dre(self, alpha, reynolds):
        _, _, derivative = self._cl.evaluate(alpha, reynolds)
        return derivative

    def cd(self, alpha, reynolds):
        value, _, _ = self._cd.evaluate(alpha, reynolds)
        return value

    def dcd_dalpha(self, alpha, reynolds):
        _, derivative, _ = self._cd.evaluate(alpha, reynolds)
        return derivative

    def dcd_dre(self, alpha, reynolds):
        _, _, derivative = self._cd.evaluate(alpha, reynolds)
        return derivative

    def cm(self, alpha, reynolds):
        value, _, _ = self._cm.evaluate(alpha, reynolds)
        return value
