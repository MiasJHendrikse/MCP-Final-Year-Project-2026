"""
Plan step 1.5 acceptance.

The Weibull parameters resolved on 2026-09-13, so `annual_energy_mwh` now
returns a number and the whole of step 1.5 is testable: the bin scheme, the
operating strategy, the rated-power limit, the distribution itself,
determinism, and the config-backed resource.

The "baseline returns a figure inside the sanity band" exit criterion is met at
10.27 MWh/yr -- but only after that band was found to be defective and revised
from 4-6 to 8-12 MWh/yr on 2026-09-13. See `test_baseline.py`, which asserts
both the band and that it is still narrow enough to catch a bug.

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
WIND_RESOURCE_PATH = os.path.abspath(os.path.join(
    _HERE, "..", "verification", "wind_resource", "wind_resource_20m.json"))
X0_PATH = os.path.abspath(os.path.join(
    _HERE, "..", "verification", "baseline", "x0.json"))


def _wind_resource_artefact():
    with open(WIND_RESOURCE_PATH, encoding="utf-8") as handle:
        return json.load(handle)


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

def test_weibull_from_config_returns_the_site_resource():
    """
    Replaces `test_weibull_from_config_still_raises`, which was the placeholder
    guarding the TODO and whose job ended when the resource landed.

    The values are pinned against `verification/wind_resource/`, not restated:
    a test that repeated the literals would agree with a typo in site.yaml.

    Tolerance is 1e-6 rather than exact because site.yaml records six decimal
    places while the artefact carries full double precision. That rounding is
    deliberate -- the config file is read by people -- and 1e-6 is tight enough
    that only the rounding fits inside it. A transposed digit would not.
    """

    resource = WeibullResource.from_config()
    artefact = _wind_resource_artefact()["extrapolation"]

    assert resource.k == pytest.approx(artefact["k"], rel=1e-6)
    assert resource.c == pytest.approx(artefact["scale_ms"], rel=1e-6)


def test_weibull_from_config_mean_matches_the_closed_form():
    """
    `site.mean_wind_speed_ms` is derived, and must equal c * Gamma(1 + 1/k).

    The loader enforces this at 1e-4; this asserts the much tighter agreement
    that the recorded value actually has, so a future edit that merely squeaks
    past the loader's tolerance is still caught here.
    """

    from config import load_site

    resource = WeibullResource.from_config()
    recorded = float(load_site().mean_wind_speed_ms)

    assert resource.mean_speed == pytest.approx(recorded, rel=1e-6)
    assert resource.mean_speed == pytest.approx(
        resource.c * math.gamma(1.0 + 1.0 / resource.k), rel=1e-12)


def test_the_site_resource_is_an_extrapolation_not_an_extraction():
    """
    The one property of this resource most likely to be forgotten later.

    `k` and `c` are NOT the numbers on the screenshot -- those are at 50 m.
    This asserts the 20 m pair differs from the 50 m pair in the directions
    downward extrapolation requires: lower scale (less wind nearer the ground)
    and lower shape (a broader distribution). A future edit that pasted the
    50 m pair straight into site.yaml would fail here.
    """

    artefact = _wind_resource_artefact()
    reference = artefact["reference_level"]
    target = artefact["extrapolation"]
    resource = WeibullResource.from_config()

    assert target["target_height_m"] == 20.0
    assert reference["height_m"] == 50.0
    assert resource.c < reference["scale_ms"]
    assert resource.k < reference["k"]


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

def test_aep_defaults_to_the_config_resource(x0, parameterisation):
    """
    Replaces `test_aep_raises_while_the_resource_is_unresolved`.

    With no `resource=`, `annual_energy_mwh` must use the site's resource --
    and must give the identical value to passing that same resource
    explicitly. Those two paths going out of step would mean some caller is
    silently integrating against a different distribution than it thinks.
    """

    implicit = annual_energy_mwh(x0, parameterisation=parameterisation)
    explicit = annual_energy_mwh(x0, resource=WeibullResource.from_config(),
                                 parameterisation=parameterisation)

    assert implicit == explicit
    assert implicit > 0.0


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
