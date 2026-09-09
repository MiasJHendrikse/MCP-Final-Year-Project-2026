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
