"""
Plan step 1.5 acceptance, for the parts that do not need the wind resource.

The Weibull parameters are still `TODO`, so `annual_energy_mwh` raises and the
"baseline returns 4-6 MWh/yr" exit criterion cannot be met yet. Everything
below it -- the bin scheme, the operating strategy, the rated-power limit, the
distribution itself, and determinism -- is testable now and is tested here.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os

import numpy as np
import pytest

from config import load_design_rotor
from config.unresolved import UnresolvedConfigError
from design import BladeParameterisation
from objective import (
    BIN_WIDTH_MS,
    HOURS_PER_YEAR,
    WeibullResource,
    annual_energy_mwh,
    bin_powers,
    energy_surrogate,
    sanity_band,
    wind_speed_bins,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.abspath(os.path.join(
    _HERE, "..", "verification", "baseline", "x0.json"))


@pytest.fixture(scope="module")
def x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"]
                    + artefact["twist_control_points_rad"])


@pytest.fixture(scope="module")
def parameterisation():
    return BladeParameterisation()


@pytest.fixture(scope="module")
def powers(x0, parameterisation):
    return bin_powers(x0, parameterisation=parameterisation)


# ---------------------------------------------------------------------------
# The bin scheme, fixed and documented
# ---------------------------------------------------------------------------

def test_bins_span_cut_in_to_cut_out():
    design = load_design_rotor()
    edges, midpoints, width = wind_speed_bins()

    assert edges[0] == design.cut_in_wind_speed_ms
    assert edges[-1] == design.cut_out_wind_speed_ms
    assert width == BIN_WIDTH_MS
    assert len(midpoints) == len(edges) - 1
    assert np.allclose(np.diff(edges), BIN_WIDTH_MS)


def test_midpoints_are_strictly_inside_their_bins():
    """
    Why the midpoint rule rather than trapezoid.

    Trapezoid evaluates at the bin edges, so adjacent bins share a node and the
    rated-power cut at exactly 11 m/s would land on one. The midpoint rule puts
    every evaluation strictly inside its own bin, so no single evaluation has
    to be simultaneously below and above rated.
    """

    design = load_design_rotor()
    edges, midpoints, _width = wind_speed_bins()

    assert np.all(midpoints > edges[:-1])
    assert np.all(midpoints < edges[1:])
    assert design.rated_wind_speed_ms not in set(midpoints)


# ---------------------------------------------------------------------------
# The operating strategy
# ---------------------------------------------------------------------------

def test_power_is_limited_above_rated_and_free_below(powers):
    """
    The rated-power cut: bins above rated are held, bins below are not.

    Checked against the *unlimited* series, so this is a test of the limiter
    rather than of the aerodynamics.
    """

    design = load_design_rotor()
    midpoints = powers["midpoints"]
    rated = powers["rated_power_w"]

    for speed, limited, unlimited in zip(midpoints, powers["power_w"],
                                         powers["power_unlimited_w"]):
        assert limited <= rated + 1e-9
        if speed < design.rated_wind_speed_ms:
            assert limited == pytest.approx(unlimited, rel=1e-12)


def test_above_rated_bins_all_hold_exactly_the_rated_power(powers):
    """
    The property that keeps J smooth in d above rated.

    Every limited bin takes the value P_aero(V_rated; d), which is smooth in
    d. If instead the limit were read off the bin grid, the objective would
    depend on where the bins happened to fall.
    """

    limited = powers["power_w"][powers["limited"]]
    assert len(limited) > 0, "no bin is limited; the test case is not exercising the cut"
    assert np.allclose(limited, powers["rated_power_w"], rtol=1e-12)


def test_power_rises_monotonically_below_rated(powers):
    """Aerodynamic power must increase with wind speed on the fixed-TSR line."""

    design = load_design_rotor()
    below = [power for speed, power in zip(powers["midpoints"],
                                           powers["power_unlimited_w"])
             if speed < design.rated_wind_speed_ms]

    assert all(below[i] < below[i + 1] for i in range(len(below) - 1))


def test_every_bin_converges(powers):
    assert powers["all_converged"], "a bin failed to converge at x0"


# ---------------------------------------------------------------------------
# The Weibull distribution
# ---------------------------------------------------------------------------

def test_weibull_from_config_still_raises():
    """
    Passing is the correct state today. When the GWA extraction lands this
    fails deliberately, which is the reminder to replace it.
    """

    with pytest.raises(UnresolvedConfigError, match="still TODO"):
        WeibullResource.from_config()


def test_weibull_mean_matches_the_closed_form():
    """V-bar = c * Gamma(1 + 1/k)."""

    for k, c in ((1.5, 5.0), (2.0, 6.0), (2.5, 7.2)):
        resource = WeibullResource(k=k, c=c)
        assert resource.mean_speed == pytest.approx(
            c * math.gamma(1.0 + 1.0 / k), rel=1e-12)


def test_weibull_pdf_integrates_to_one():
    """The distribution is normalised -- the cheapest check that it is one."""

    resource = WeibullResource(k=2.0, c=6.0)
    speeds = np.linspace(1e-6, 60.0, 400_001)
    assert np.trapezoid(resource.pdf(speeds), speeds) == pytest.approx(1.0, abs=1e-6)


def test_bin_mass_from_cdf_matches_the_integrated_pdf():
    """
    Bin masses come from the CDF, which is exact for any bin width.

    That keeps the bin scheme's accuracy a question about resolving the POWER
    curve, not the distribution -- two separate concerns that would otherwise
    be confounded in any bin-width study.
    """

    resource = WeibullResource(k=2.1, c=6.4)
    for lower, upper in ((3.0, 4.0), (7.0, 8.0), (19.0, 20.0)):
        speeds = np.linspace(lower, upper, 20_001)
        assert resource.probability_between(lower, upper) == pytest.approx(
            np.trapezoid(resource.pdf(speeds), speeds), rel=1e-9)


def test_weibull_rejects_invalid_parameters():
    for k, c in ((0.0, 6.0), (-1.0, 6.0), (2.0, 0.0), (2.0, -3.0)):
        with pytest.raises(ValueError):
            WeibullResource(k=k, c=c)


# ---------------------------------------------------------------------------
# AEP, and the surrogate that stands in for it
# ---------------------------------------------------------------------------

def test_aep_raises_while_the_resource_is_unresolved(x0, parameterisation):
    with pytest.raises(UnresolvedConfigError):
        annual_energy_mwh(x0, parameterisation=parameterisation)


def test_aep_works_with_an_explicit_resource(x0, parameterisation):
    """
    The mechanism is complete -- only the two numbers are missing.

    Uses an explicitly stated provisional resource, which is the supported way
    to run a study without writing provisional values into `config/`. The
    result is checked for plausibility only; it is NOT a project figure and the
    parameters are not the site's.
    """

    resource = WeibullResource(k=2.0, c=6.0)
    aep = annual_energy_mwh(x0, resource=resource, parameterisation=parameterisation)

    assert aep > 0.0
    # Sanity: below the energy a rotor held at rated power all year could make.
    rated = bin_powers(x0, parameterisation=parameterisation)["rated_power_w"]
    assert aep < rated * HOURS_PER_YEAR / 1e6


def test_surrogate_is_the_unit_weighted_sum(x0, parameterisation, powers):
    """
    The surrogate is exactly the per-bin powers summed -- no scaling, no
    hidden weighting.

    Asserted because the whole argument for using it in the smoothness gate
    rests on it sharing the d-dependent chain with AEP and differing only in
    the weights.
    """

    assert energy_surrogate(x0, parameterisation=parameterisation) == pytest.approx(
        float(np.sum(powers["power_w"])), rel=1e-12)


def test_surrogate_and_aep_share_the_same_bin_powers(x0, parameterisation):
    """
    Both go through `bin_powers`, so they cannot drift apart on the physics.

    A change that moved one and not the other would invalidate the gate's
    conclusions without any test failing, so this pins the shared path.
    """

    resource = WeibullResource(k=2.0, c=6.0)
    powers = bin_powers(x0, parameterisation=parameterisation)
    edges, _midpoints, _width = wind_speed_bins()
    mass = resource.probability_between(edges[:-1], edges[1:])

    expected = float(np.sum(powers["power_w"] * mass)) * HOURS_PER_YEAR / 1e6
    assert annual_energy_mwh(x0, resource=resource,
                             parameterisation=parameterisation) == pytest.approx(
        expected, rel=1e-12)


def test_hours_per_year_is_the_julian_year():
    """8766, not 8760 -- the latter biases every AEP low by about 0.07 %."""

    assert HOURS_PER_YEAR == 365.25 * 24.0


def test_sanity_band_comes_from_config():
    low, high = sanity_band()
    design = load_design_rotor()
    assert (low, high) == (design.aep_mwh_per_year_min, design.aep_mwh_per_year_max)


# ---------------------------------------------------------------------------
# Determinism, at the level the brief actually states it
# ---------------------------------------------------------------------------

def test_objective_is_bitwise_reproducible(x0, parameterisation):
    """
    The brief's determinism requirement, now at the AEP level rather than only
    the solver level (`tests/test_determinism.py` covers `solve_rotor`).

    Bitwise, repeated, and with an interleaved evaluation at a different design
    vector in between -- the guard against any accumulated state in the chain.
    """

    resource = WeibullResource(k=2.0, c=6.0)
    first = annual_energy_mwh(x0, resource=resource, parameterisation=parameterisation)

    perturbed = x0.copy()
    perturbed[0] *= 1.1
    annual_energy_mwh(perturbed, resource=resource, parameterisation=parameterisation)

    assert annual_energy_mwh(x0, resource=resource,
                             parameterisation=parameterisation) == first
    assert energy_surrogate(x0, parameterisation=parameterisation) == \
        energy_surrogate(x0, parameterisation=parameterisation)
