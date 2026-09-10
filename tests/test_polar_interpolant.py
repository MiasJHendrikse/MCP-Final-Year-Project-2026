"""
Task 3 acceptance: the C1 polar interpolant's analytic derivatives agree with
complex-step differentiation, and the audit's staircase measurement is gone.

Complex-step reference: for a holomorphic f, f(x + ih).imag / h equals f'(x)
to machine precision as h -> 0, with none of finite differencing's
subtraction cancellation error. It is not usable as gradient code in the
solver (out of scope, ground rule 4) -- here it is purely a verification
oracle for the analytic partials `PolarInterpolant` exposes.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
import os

import numpy as np
import pytest

from polars.cache import PolarGrid
from polars.interpolant import PolarDomainError, PolarInterpolant

_HERE = os.path.dirname(os.path.abspath(__file__))
S809_DIR = os.path.abspath(os.path.join(_HERE, "..", "data", "polars", "s809"))

#: Complex-step size. Small enough that truncation error is negligible,
#: large enough to stay well clear of double-precision underflow in the
#: imaginary part (see e.g. Martins, Sturdza & Alonso 2003).
_H = 1e-30

#: (alpha, Reynolds) points spanning attached flow, near-stall and the
#: Reynolds range's interior and edges.
_TEST_POINTS = [
    (4.0, 6e5),        # exact grid point in both axes
    (5.6, 6e5),
    (10.3, 3.5e5),
    (-3.7, 1.05e6),
    (0.25, 1.5e5),
    (17.9, 9.9e5),
    (-20.0, 1.0e5),    # Reynolds range's lower edge
    (170.0, 1.3e6),    # Reynolds range's upper edge, deep in the Viterna extension
]


@pytest.fixture(scope="module")
def interpolant():
    grid = PolarGrid(S809_DIR)
    assert not grid.has_gaps, "S809 cache has gaps; Task 2 precondition violated"
    return PolarInterpolant(grid)


@pytest.mark.parametrize("alpha,reynolds", _TEST_POINTS)
def test_dcl_dalpha_matches_complex_step(interpolant, alpha, reynolds):
    analytic = interpolant.dcl_dalpha(alpha, reynolds)
    perturbed = interpolant.cl(alpha + _H * 1j, reynolds)
    complex_step = perturbed.imag / _H
    assert complex_step == pytest.approx(analytic, abs=1e-14, rel=1e-14)


@pytest.mark.parametrize("alpha,reynolds", _TEST_POINTS)
def test_dcl_dre_matches_complex_step(interpolant, alpha, reynolds):
    analytic = interpolant.dcl_dre(alpha, reynolds)
    perturbed = interpolant.cl(alpha, reynolds + _H * 1j)
    complex_step = perturbed.imag / _H
    assert complex_step == pytest.approx(analytic, abs=1e-14, rel=1e-14)


@pytest.mark.parametrize("alpha,reynolds", _TEST_POINTS)
def test_dcd_dalpha_matches_complex_step(interpolant, alpha, reynolds):
    analytic = interpolant.dcd_dalpha(alpha, reynolds)
    perturbed = interpolant.cd(alpha + _H * 1j, reynolds)
    complex_step = perturbed.imag / _H
    assert complex_step == pytest.approx(analytic, abs=1e-14, rel=1e-14)


@pytest.mark.parametrize("alpha,reynolds", _TEST_POINTS)
def test_dcd_dre_matches_complex_step(interpolant, alpha, reynolds):
    analytic = interpolant.dcd_dre(alpha, reynolds)
    perturbed = interpolant.cd(alpha, reynolds + _H * 1j)
    complex_step = perturbed.imag / _H
    assert complex_step == pytest.approx(analytic, abs=1e-14, rel=1e-14)


def test_exact_grid_point_reproduces_cached_row(interpolant):
    """
    At a point that is exactly on both axes' grids, the interpolant must
    reproduce the cached value exactly (mirrors the CSV row for
    S809_Re600000.csv, alpha=4.0: cl=0.618600, cd=0.010370, cm=-0.052200).
    """

    cl, cd, cm = interpolant(4.0, 6e5)
    assert cl == pytest.approx(0.618600, abs=1e-6)
    assert cd == pytest.approx(0.010370, abs=1e-6)
    assert cm == pytest.approx(-0.052200, abs=1e-6)


def test_no_staircase_in_dcl_dalpha():
    """
    Reproduces the audit's measurement directly: dCl/dalpha for alpha in
    [4.0, 5.6] deg at Re=600k took exactly four distinct values under the old
    bilinear lookup. The new interpolant must not collapse onto a handful of
    values across the same sweep.
    """

    grid = PolarGrid(S809_DIR)
    interp = PolarInterpolant(grid)

    alphas = np.linspace(4.0, 5.6, 200)
    derivative = np.array([interp.dcl_dalpha(a, 6e5) for a in alphas])
    distinct = len(np.unique(np.round(derivative, 8)))

    assert distinct > 20, (
        f"only {distinct} distinct dCl/dalpha values across 200 samples -- "
        "looks like the old bilinear staircase is back"
    )


def test_out_of_range_alpha_raises(interpolant):
    with pytest.raises(PolarDomainError):
        interpolant.cl(200.0, 6e5)


def test_out_of_range_reynolds_raises(interpolant):
    with pytest.raises(PolarDomainError):
        interpolant.cl(4.0, 1.0)


def test_grid_rejects_nan():
    """A grid with a gap must not silently produce a plausible-looking surface."""

    grid = PolarGrid(S809_DIR)
    grid.cl_grid[3, 100] = np.nan
    with pytest.raises(PolarDomainError):
        PolarInterpolant(grid)


# ---------------------------------------------------------------------------
# Continuity ACROSS the alpha knots, at Reynolds numbers between cache rows
# ---------------------------------------------------------------------------
#
# Added 2026-09-10, after plan step 1.7 found a real discontinuity that every
# test above passed straight through.
#
# The defect: the surface used to be built by fitting a cubic spline in alpha
# per Reynolds row and then PCHIP-blending the spline COEFFICIENTS across
# log(Re). Continuity at an alpha knot is a *linear* constraint on those
# coefficient channels, and PCHIP is *nonlinear* in its data -- its slope
# limiter is a harmonic mean -- so blending channel-by-channel does not
# preserve it. On a cache row the blended coefficients are the originals and
# continuity was restored, which is why it hid: every Reynolds number these
# tests happened to use was either a row or close enough not to matter.
#
# Measured on SG6043 before the fix: Cl jumped 5.4e-3 across alpha = 6.0 deg
# at Re = 70,995, while dCl/dalpha stayed continuous (0.1766 -> 0.1776). A
# station at that Reynolds number failed to converge in the BEM solver, which
# is how it was found -- the residual was genuinely discontinuous.
#
# Why nothing above caught it: complex-step differentiation evaluates ONE
# polynomial piece. It verifies the derivative within a piece and is blind by
# construction to a mismatch between adjacent pieces. `test_no_staircase_in_
# dcl_dalpha` samples the derivative, which was continuous. The value jump was
# in a quantity no test compared across a knot.
#
# The fix blends node VALUES and node ALPHA-SLOPES across log(Re) and builds a
# cubic Hermite in alpha. Adjacent intervals share both endpoint quantities, so
# C1 in alpha holds identically at any Reynolds number, by construction rather
# than by luck.

#: Both caches. The defect was found on SG6043 (the design airfoil), but it
#: was a property of the construction, not of the data, so S809 carried it too.
_CACHE_DIRS = {
    "s809": S809_DIR,
    "sg6043": os.path.abspath(os.path.join(_HERE, "..", "data", "polars", "sg6043")),
}


def _off_row_reynolds(interpolant):
    """
    Geometric midpoints between adjacent cache rows.

    Derived from the cache rather than hard-coded: mid-interval in log(Re) is
    where the old blending error peaked, and deriving it keeps every point in
    range for whichever cache is being tested.
    """

    rows = [float(value) for value in interpolant.re_values]
    return [math.sqrt(lo * hi) for lo, hi in zip(rows[:-1], rows[1:])]


@pytest.fixture(scope="module", params=sorted(_CACHE_DIRS))
def any_cache(request):
    """(interpolant, grid) for each committed cache."""

    grid = PolarGrid(_CACHE_DIRS[request.param])
    assert not grid.has_gaps
    return PolarInterpolant(grid), grid


@pytest.fixture(scope="module")
def any_interpolant(any_cache):
    return any_cache[0]


def test_value_is_continuous_across_alpha_knots_off_row(any_interpolant):
    """
    No jump in Cl, Cd or Cm across a cache alpha node, at any Reynolds number.

    Every alpha node in the working band is checked at every off-row Reynolds
    number, not a sample: the old defect's magnitude varied node to node and a
    sparse sample could miss the worst one.
    """

    for reynolds in _off_row_reynolds(any_interpolant):
        for alpha in any_interpolant.alpha_values:
            alpha = float(alpha)
            if not (-30.0 <= alpha <= 30.0):
                continue
            for name, function in (("cl", any_interpolant.cl),
                                   ("cd", any_interpolant.cd),
                                   ("cm", any_interpolant.cm)):
                jump = abs(function(alpha, reynolds)
                           - function(alpha - 1e-9, reynolds))
                assert jump < 1e-8, (
                    f"{name} jumps {jump:.3e} across the alpha = {alpha} deg "
                    f"knot at Re = {reynolds:,.0f}; the surface is not C0 there")


def test_derivative_is_continuous_across_alpha_knots_off_row(any_interpolant):
    """
    The C1 half of the same statement.

    A Hermite basis guarantees the slope matches at shared nodes, so this holds
    to round-off rather than merely to a tolerance.
    """

    for reynolds in _off_row_reynolds(any_interpolant):
        for alpha in any_interpolant.alpha_values:
            alpha = float(alpha)
            if not (-30.0 <= alpha <= 30.0):
                continue
            jump = abs(any_interpolant.dcl_dalpha(alpha, reynolds)
                       - any_interpolant.dcl_dalpha(alpha - 1e-9, reynolds))
            assert jump < 1e-6, (
                f"dCl/dalpha jumps {jump:.3e} across alpha = {alpha} deg at "
                f"Re = {reynolds:,.0f}")


def test_cached_rows_are_reproduced_bitwise_at_every_node(any_cache):
    """
    On a cache row, at a cache alpha node, the surface returns the cached
    number itself -- bitwise.

    The Hermite construction interpolates the cached values directly, rather
    than reconstructing them from blended spline coefficients. That is what
    makes this exact rather than merely close, and it is the property that
    proves the continuity fix did not quietly change what the cache says: the
    data is reproduced, only the path between data points moved.
    """

    interpolant, grid = any_cache

    for row, reynolds in enumerate(float(value) for value in grid.re_values):
        for index in range(0, len(grid.alpha_values), 11):
            alpha = float(grid.alpha_values[index])
            assert interpolant.cl(alpha, reynolds) == grid.cl_grid[row][index]
            assert interpolant.cd(alpha, reynolds) == grid.cd_grid[row][index]
            assert interpolant.cm(alpha, reynolds) == grid.cm_grid[row][index]
