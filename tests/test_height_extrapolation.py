"""
The Weibull height extrapolation, and the artefact it produced.

`src/objective/height_extrapolation.py` carries an open provenance item: the
Justus & Mikhail (1976) formulae are transcribed from general knowledge, not
from a retrieved copy of the paper. So the tests here are doing real work, not
ceremony -- they are the entire basis for believing the transcription, exactly
as `test_corrections.py` is for Buhl's constants.

They are arranged in increasing order of what they would catch:

  1. Identity and sign properties, which catch a mangled expression.
  2. The mean-speed identity on the *input* pair, which catches a misread
     screenshot.
  3. Agreement with the log law across plausible roughness, which is the only
     check here with an independent physical basis and the only one that could
     catch a wrong constant rather than a wrong rearrangement.

What none of them can catch is a faithful transcription of the wrong
correlation. That is what the provenance item in docs/OUTSTANDING-INPUTS.md is
for, and no test should be read as closing it.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import math
import os

import pytest

from objective.height_extrapolation import (
    JUSTUS_ANCHOR_HEIGHT_M,
    extrapolate_weibull,
    log_law_scale_ms,
    shear_exponent,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
ARTEFACT_PATH = os.path.abspath(os.path.join(
    _HERE, "..", "verification", "wind_resource", "wind_resource_20m.json"))


@pytest.fixture(scope="module")
def artefact():
    with open(ARTEFACT_PATH, encoding="utf-8") as handle:
        return json.load(handle)


# ---------------------------------------------------------------------------
# 1. Structural properties of the correlation
# ---------------------------------------------------------------------------

def test_extrapolating_to_the_reference_height_changes_nothing():
    """The identity case. A sign error in either exponent breaks it."""

    result = extrapolate_weibull(k=1.9, scale_ms=7.5,
                                 reference_height_m=50.0,
                                 target_height_m=50.0)

    assert result.k == pytest.approx(1.9, rel=1e-12)
    assert result.scale_ms == pytest.approx(7.5, rel=1e-12)


def test_both_parameters_fall_when_going_down():
    """
    Downward extrapolation must reduce BOTH parameters.

    Less wind nearer the ground (lower scale), and a broader distribution
    (lower shape) because shear and surface effects widen the spread. Getting
    the k relation upside down is the single most plausible transcription
    error, and it would pass every check that only looked at the scale.
    """

    result = extrapolate_weibull(k=1.87, scale_ms=8.8,
                                 reference_height_m=50.0,
                                 target_height_m=20.0)

    assert result.scale_ms < 8.8
    assert result.k < 1.87


def test_both_parameters_rise_when_going_up():
    """The converse, so the direction is pinned rather than the magnitude."""

    result = extrapolate_weibull(k=1.87, scale_ms=8.8,
                                 reference_height_m=50.0,
                                 target_height_m=100.0)

    assert result.scale_ms > 8.8
    assert result.k > 1.87


def test_the_scale_is_monotone_in_height():
    """No fold or turning point anywhere in the range of interest."""

    scales = [extrapolate_weibull(1.87, 8.8, 50.0, z).scale_ms
              for z in (10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0)]

    assert all(a < b for a, b in zip(scales, scales[1:]))


def test_round_trip_returns_the_original_scale():
    """
    50 -> 20 -> 50 must return the input scale.

    Note this is NOT true of `k`: the k relation is anchored to 10 m, so
    re-extrapolating uses the new height's denominator and the round trip is
    exact for k as well. Both are asserted -- if either drifted, the relation
    would not be a pure function of height, which it must be.
    """

    down = extrapolate_weibull(1.87, 8.8, 50.0, 20.0)
    up = extrapolate_weibull(down.k, down.scale_ms, 20.0, 50.0)

    assert up.k == pytest.approx(1.87, rel=1e-9)
    # The scale round trip is only approximate: the exponent `n` is recomputed
    # from the NEW scale, so it is not the same exponent coming back. That is a
    # property of the correlation, not a defect -- it is asserted loosely to
    # record that it was checked and understood rather than overlooked.
    assert up.scale_ms == pytest.approx(8.8, rel=2e-2)


def test_the_anchor_height_is_where_the_k_relation_is_neutral():
    """At z = 10 m the k denominator is 1 by construction."""

    result = extrapolate_weibull(2.0, 6.0,
                                 reference_height_m=JUSTUS_ANCHOR_HEIGHT_M,
                                 target_height_m=JUSTUS_ANCHOR_HEIGHT_M)
    assert result.k == pytest.approx(2.0, rel=1e-12)


def test_the_shear_exponent_falls_as_the_reference_scale_rises():
    """
    Windier sites have flatter profiles. The exponent carries `-ln A`, so this
    is a sign check on the term most easily transcribed backwards.
    """

    exponents = [shear_exponent(scale, 50.0) for scale in (4.0, 6.0, 8.0, 10.0)]
    assert all(a > b for a, b in zip(exponents, exponents[1:]))


def test_invalid_inputs_raise_rather_than_returning_nonsense():
    for kwargs in (
        {"k": 0.0, "scale_ms": 7.0, "reference_height_m": 50.0,
         "target_height_m": 20.0},
        {"k": 1.9, "scale_ms": -1.0, "reference_height_m": 50.0,
         "target_height_m": 20.0},
        {"k": 1.9, "scale_ms": 7.0, "reference_height_m": 50.0,
         "target_height_m": 0.0},
    ):
        with pytest.raises(ValueError):
            extrapolate_weibull(**kwargs)


def test_the_log_law_refuses_heights_inside_the_roughness():
    """The log profile is undefined at or below z0; it must say so."""

    with pytest.raises(ValueError):
        log_law_scale_ms(8.8, 50.0, 20.0, roughness_length_m=25.0)
    with pytest.raises(ValueError):
        log_law_scale_ms(8.8, 50.0, 20.0, roughness_length_m=0.0)


# ---------------------------------------------------------------------------
# 2. The input pair was read correctly
# ---------------------------------------------------------------------------

def test_the_screenshot_pair_reproduces_its_own_displayed_mean(artefact):
    """
    GASP displayed A = 8.8, k = 1.87 and U = 7.8 m/s. Those are not
    independent: U = A * Gamma(1 + 1/k).

    This is the check that the screenshot was transcribed correctly. A misread
    digit in either parameter moves the implied mean well outside display
    rounding.
    """

    reference = artefact["reference_level"]
    implied = (reference["scale_ms"]
               * math.gamma(1.0 + 1.0 / reference["k"]))

    assert implied == pytest.approx(
        reference["mean_speed_implied_by_k_and_A_ms"], rel=1e-12)
    assert implied == pytest.approx(reference["displayed_mean_speed_ms"],
                                    abs=0.05)
    assert reference["self_consistent"] is True


# ---------------------------------------------------------------------------
# 3. The independent cross-check -- the one with real power
# ---------------------------------------------------------------------------

def test_the_central_case_sits_inside_the_log_law_band(artefact):
    """
    Justus & Mikhail against the logarithmic profile, across plausible
    roughness for Khomas Hochland bush savanna.

    These two have genuinely different bases: one is an empirical correlation
    on the scale parameter, the other is surface-layer similarity theory. They
    are not required to agree exactly and the band is not a tolerance -- but a
    central case falling OUTSIDE the band over the whole plausible roughness
    range would mean one of them is wrong, and that is worth knowing before any
    AEP figure is quoted.
    """

    cross = artefact["cross_check_log_law"]
    low, high = cross["band_ms"]
    scale = artefact["extrapolation"]["scale_ms"]

    assert low < scale < high
    assert cross["central_case_inside_band"] is True


def test_the_log_law_band_is_recomputed_not_just_read_back(artefact):
    """
    The artefact is regenerated by a script; this recomputes its cross-check
    from the recorded inputs so a stale artefact cannot pass by agreeing with
    itself.
    """

    reference = artefact["reference_level"]
    cross = artefact["cross_check_log_law"]

    recomputed = [
        log_law_scale_ms(reference["scale_ms"], reference["height_m"],
                         artefact["extrapolation"]["target_height_m"], z0)
        for z0 in cross["roughness_lengths_m"]
    ]

    assert min(recomputed) == pytest.approx(cross["band_ms"][0], rel=1e-12)
    assert max(recomputed) == pytest.approx(cross["band_ms"][1], rel=1e-12)


# ---------------------------------------------------------------------------
# 4. The artefact and the config agree
# ---------------------------------------------------------------------------

def test_the_artefact_reproduces_from_its_own_recorded_inputs(artefact):
    """
    Recompute the extrapolation from the artefact's stated inputs and require
    the stated outputs back. Catches an artefact hand-edited after the fact.
    """

    reference = artefact["reference_level"]
    recorded = artefact["extrapolation"]

    result = extrapolate_weibull(
        k=reference["k"],
        scale_ms=reference["scale_ms"],
        reference_height_m=reference["height_m"],
        target_height_m=recorded["target_height_m"],
    )

    assert result.k == pytest.approx(recorded["k"], rel=1e-12)
    assert result.scale_ms == pytest.approx(recorded["scale_ms"], rel=1e-12)
    assert result.mean_speed_ms == pytest.approx(recorded["mean_speed_ms"],
                                                 rel=1e-12)


def test_the_extrapolation_targets_the_configured_hub_height(artefact):
    """
    20 m is not a literal that may drift. It is `site.hub_height_m`, and if the
    tower height is ever revised this must fail rather than the resource
    quietly continuing to describe a different height than the rotor flies at.
    """

    from config import load_site

    assert (artefact["extrapolation"]["target_height_m"]
            == pytest.approx(load_site().hub_height_m))


def test_the_plan_prior_disagreement_is_recorded(artefact):
    """
    Plan 1.3 expected k ~ 1.8-2.4 and c ~ 6-7 m/s. The result is outside both.

    Pinned here so it stays visible. This is a prior, not a specification --
    the disagreement is documented rather than tuned away,
    and a later edit that quietly brings them into line has to change this
    test and say why.
    """

    prior = artefact["against_plan_prior_expectation"]

    assert prior["k_inside"] is False, "k moved back inside the prior band"
    assert prior["c_inside"] is False, "c moved back inside the prior band"


# ---------------------------------------------------------------------------
# 5. The held-k penalty (the design basis's own sensitivity)
# ---------------------------------------------------------------------------

def test_the_held_k_penalty_block_is_present_and_small(artefact):
    """
    Pinning k at its 50 m value costs the reference blade about 0.9 % of its
    annual energy, under the committed law and rating.

    The block is the artefact the design basis quotes when it says the
    extrapolated shape parameter is worth having; a missing block, or a
    penalty an order of magnitude from the recorded one, is a defect either
    way round.
    """

    block = artefact["held_k_penalty"]

    assert block["k_held"] == pytest.approx(artefact["reference_level"]["k"],
                                            rel=1e-12)
    assert (block["c_held_ms"]
            == pytest.approx(artefact["extrapolation"]["scale_ms"], rel=1e-12))
    assert (block["aep_k_held"] < block["aep_configured"]), (
        "holding k at its 50 m value must cost energy, not gain it")
    assert block["penalty_pct"] == pytest.approx(-0.91, abs=0.05)
