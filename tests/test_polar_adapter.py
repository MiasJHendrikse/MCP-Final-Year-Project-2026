"""
Task 4 acceptance: the polar adapter clamps nothing, substitutes nothing, and
holds no module-level airfoil state.

Each test here corresponds to one of the work order's "done when" clauses for
Task 4. The clause "Phase VI cross-check figures unchanged" is not testable
here -- it is what `test_golden_regression.py` covers, and Task 4 moves those
figures deliberately (see `tests/golden/README.md`'s change log).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

import pytest

from polars.interpolant import PolarDomainError
from polars.polar import (
    CachedPolar,
    CachedPolarFactory,
    interpolant_for,
    polar_factory_for,
)

#: Both caches, so the design airfoil is covered as well as the validation one.
#: SG6043 is the harder case: its Re floor is 40k and it carries a real
#: first-derivative discontinuity inside the converged band at Re = 40k (the
#: laminar-separation-bubble jump documented in data/polars/sg6043/README.md).
CACHES = ["s809", "sg6043"]

#: An angle of attack safely inside every cache's XFOIL-converged band.
ALPHA_ATTACHED_DEG = 5.0


def _mid_reynolds(interpolant):
    """A Reynolds number strictly inside the built range, log-centred."""

    lo = float(interpolant.re_values[0])
    hi = float(interpolant.re_values[-1])
    return math.exp(0.5 * (math.log(lo) + math.log(hi)))


# ---------------------------------------------------------------------------
# "An out-of-table alpha raises rather than returning a plausible number"
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cache", CACHES)
@pytest.mark.parametrize("alpha_deg", [-180.5, 181.0, 1e4])
def test_out_of_table_alpha_raises(cache, alpha_deg):
    """
    Past +/-180 deg there is no data and no defensible extrapolation, so the
    adapter must refuse rather than return the boundary value.

    The old `S809Polar` clamped to the XFOIL-converged band instead, which at
    the Phase VI 15 m/s point meant reading Cl(18 deg) for a station at 32 deg
    -- a plausible number, indistinguishable from a correct one, inside the
    residual.
    """

    polar = polar_factory_for(cache)(_mid_reynolds(interpolant_for(cache)))
    with pytest.raises(PolarDomainError):
        polar.cl(math.radians(alpha_deg))
    with pytest.raises(PolarDomainError):
        polar.cd(math.radians(alpha_deg))


@pytest.mark.parametrize("cache", CACHES)
def test_alpha_past_the_xfoil_band_is_served_not_clamped(cache):
    """
    Inside +/-180 but past the XFOIL-converged band, the Viterna extension is
    real data and must be returned as such.

    This is the positive half of the clause above: "raises out of range" is
    only correct if the range is the whole circle. If this were still clamping
    at the converged band, Cl at 40 and 60 deg would be identical.
    """

    interpolant = interpolant_for(cache)
    polar = CachedPolar(interpolant, _mid_reynolds(interpolant))

    beyond = interpolant.xfoil_alpha_max + 10.0
    assert beyond < 180.0, "test assumes the converged band is well inside the circle"

    cl_40 = polar.cl(math.radians(40.0))
    cl_60 = polar.cl(math.radians(60.0))
    cl_edge = polar.cl(math.radians(interpolant.xfoil_alpha_max))

    assert cl_40 != cl_60, "post-stall Cl is constant -- alpha is being clamped"
    assert cl_40 != cl_edge, "post-stall Cl equals the band edge -- alpha is being clamped"


# ---------------------------------------------------------------------------
# "No clamping ... remains in the polar path" -- the Reynolds half
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cache", CACHES)
def test_out_of_range_reynolds_raises_at_construction(cache):
    """
    Both ends. The ceiling is the one with history: before the 2026-07-28
    extension every Phase VI station clamped to the cache's old 500k ceiling,
    so there was zero real Reynolds variation in any run and nothing said so.
    """

    interpolant = interpolant_for(cache)
    lo = float(interpolant.re_values[0])
    hi = float(interpolant.re_values[-1])

    with pytest.raises(PolarDomainError):
        CachedPolar(interpolant, lo * 0.5)
    with pytest.raises(PolarDomainError):
        CachedPolar(interpolant, hi * 2.0)

    # ...and the message names the range, so a coverage failure is diagnosable
    # without reading the cache directory.
    with pytest.raises(PolarDomainError, match=r"outside the cached range"):
        CachedPolar(interpolant, hi * 2.0)


@pytest.mark.parametrize("cache", CACHES)
def test_reynolds_is_not_clamped_inside_the_range(cache):
    """Two Reynolds numbers inside the range must give two different answers."""

    interpolant = interpolant_for(cache)
    lo = float(interpolant.re_values[0])
    hi = float(interpolant.re_values[-1])
    alpha = math.radians(ALPHA_ATTACHED_DEG)

    near_floor = CachedPolar(interpolant, lo * 1.01).cl(alpha)
    near_ceiling = CachedPolar(interpolant, hi * 0.99).cl(alpha)

    assert near_floor != near_ceiling, "Cl does not vary with Re -- Re is being clamped"


# ---------------------------------------------------------------------------
# "No silent substitution"
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cache", CACHES)
def test_adapter_is_a_pure_view_of_the_interpolant(cache):
    """
    `CachedPolar` reads the surface and converts units. Nothing else.

    Asserted as exact equality, not to a tolerance: the adapter does no
    arithmetic of its own beyond the radian-to-degree scaling, so any
    difference at all would be a substitution, a clamp or a fallback that this
    task exists to remove.

    The reference angle is scaled by the adapter's own factor rather than by
    `math.degrees(math.radians(x))`. That round trip is not exact in binary
    floating point -- 12 deg comes back as 12.000000000000002 and lands a
    unit in the last place away on Cd -- which would make this test fail on
    an arithmetic artefact of its own construction rather than on anything
    the adapter did.
    """

    interpolant = interpolant_for(cache)
    reynolds = _mid_reynolds(interpolant)
    polar = CachedPolar(interpolant, reynolds)
    deg_per_rad = 180.0 / math.pi

    for alpha_deg in (-20.0, -5.0, 0.0, ALPHA_ATTACHED_DEG, 12.0, 45.0, 120.0):
        alpha = math.radians(alpha_deg)
        reference_deg = alpha * deg_per_rad
        assert polar.cl(alpha) == interpolant.cl(reference_deg, reynolds)
        assert polar.cd(alpha) == interpolant.cd(reference_deg, reynolds)
        assert polar.cm(alpha) == interpolant.cm(reference_deg, reynolds)
        assert polar(alpha) == (polar.cl(alpha), polar.cd(alpha))


# ---------------------------------------------------------------------------
# "No module-level mutable airfoil state"
# ---------------------------------------------------------------------------

def test_two_caches_are_live_simultaneously_and_order_independent():
    """
    The direct guard on what the `ACTIVE_AIRFOIL` global cost.

    From Phase 1.4 the S809 validation cache and the SG6043 design cache are
    both in use in the same process. Under a module-level "active airfoil",
    whichever was set last decided what the next caller got, so the answer
    depended on call order -- a direct threat to the brief's "same vector via
    a different code path" determinism requirement. Constructed in both
    orders, interleaved, the answers must be identical.
    """

    alpha = math.radians(ALPHA_ATTACHED_DEG)
    reynolds = 200000.0  # inside both caches' Reynolds ranges

    s809_first = CachedPolar(interpolant_for("s809"), reynolds).cl(alpha)
    sg6043_second = CachedPolar(interpolant_for("sg6043"), reynolds).cl(alpha)

    sg6043_first = CachedPolar(interpolant_for("sg6043"), reynolds).cl(alpha)
    s809_second = CachedPolar(interpolant_for("s809"), reynolds).cl(alpha)

    assert s809_first == s809_second
    assert sg6043_first == sg6043_second
    assert s809_first != sg6043_first, "two different airfoils gave the same Cl"


def test_interpolant_is_shared_not_rebuilt():
    """
    `interpolant_for` memoises per cache name.

    Not a performance nicety: `solve_rotor` builds one polar per station per
    operating point, and a rebuild there would put the ~15 ms grid-and-spline
    construction inside the sweep. Identity, not equality, is the check --
    the object must be *the same* object.
    """

    assert interpolant_for("s809") is interpolant_for("s809")
    assert interpolant_for("s809") is not interpolant_for("sg6043")


def test_factory_produces_independent_polars():
    """`CachedPolarFactory` holds no per-call state of its own."""

    factory = CachedPolarFactory(interpolant_for("s809"))
    a = factory(200000.0)
    b = factory(800000.0)

    assert a.reynolds == 200000.0
    assert b.reynolds == 800000.0
    assert a.cl(math.radians(ALPHA_ATTACHED_DEG)) != b.cl(math.radians(ALPHA_ATTACHED_DEG))


# ---------------------------------------------------------------------------
# Derivatives: the radian/degree chain rule
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cache", CACHES)
def test_alpha_derivatives_are_per_radian(cache):
    """
    The adapter takes radians, so its alpha-derivatives are per radian, while
    the interpolant underneath works in degrees. Getting that conversion
    backwards is a factor of 3283 and would be caught by nothing else here --
    the values would still look like a smooth, plausible polar.

    Checked against complex step through the adapter itself, which also
    confirms the conversion did not coerce the perturbation to float on the
    way in (see `polars.polar`'s note on `math.degrees`).
    """

    interpolant = interpolant_for(cache)
    polar = CachedPolar(interpolant, _mid_reynolds(interpolant))
    alpha = math.radians(ALPHA_ATTACHED_DEG)
    h = 1e-30

    for value, derivative in ((polar.cl, polar.dcl_dalpha),
                              (polar.cd, polar.dcd_dalpha)):
        complex_step = value(complex(alpha, h)).imag / h
        assert derivative(alpha) == pytest.approx(complex_step, rel=1e-12)


@pytest.mark.parametrize("cache", CACHES)
def test_reynolds_derivatives_pass_through_unscaled(cache):
    """dC/dRe needs no conversion -- Reynolds number is dimensionless."""

    interpolant = interpolant_for(cache)
    reynolds = _mid_reynolds(interpolant)
    polar = CachedPolar(interpolant, reynolds)
    alpha_deg = ALPHA_ATTACHED_DEG
    alpha = math.radians(alpha_deg)

    assert polar.dcl_dre(alpha) == interpolant.dcl_dre(alpha_deg, reynolds)
    assert polar.dcd_dre(alpha) == interpolant.dcd_dre(alpha_deg, reynolds)
