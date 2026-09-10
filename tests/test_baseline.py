"""
Plan step 1.7 acceptance: the Schmitz baseline and the committed `x0`.

The plan calls the evaluated performance of this blade "the reference numbers
for the entire results chapter". These tests guard the two things that would
quietly invalidate that: `x0` drifting away from the artefact Phases 2-5 read,
and the feasibility status reading as "passed" when it was never checked.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os

import numpy as np
import pytest

from config import load_design_rotor
from design import build_schmitz_baseline, evaluate_baseline
from design.schmitz import (
    DEFAULT_DESIGN_REYNOLDS,
    schmitz_chord,
    schmitz_inflow_angle,
    schmitz_twist,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
ARTEFACT_DIR = os.path.abspath(os.path.join(_HERE, "..", "verification", "baseline"))


@pytest.fixture(scope="module")
def baseline():
    return build_schmitz_baseline()


@pytest.fixture(scope="module")
def performance(baseline):
    return evaluate_baseline(baseline)


@pytest.fixture(scope="module")
def committed_x0():
    with open(os.path.join(ARTEFACT_DIR, "x0.json"), encoding="utf-8") as handle:
        return json.load(handle)


# ---------------------------------------------------------------------------
# The analytic formulae
# ---------------------------------------------------------------------------

def test_schmitz_inflow_angle_is_two_thirds_of_the_zero_swirl_angle():
    """
    The 2/3 is the entire content of the wake-rotation correction.

    The Betz construction gives arctan(1/lambda_r); the Schmitz optimum sits
    two-thirds of the way from the rotor plane toward it. Checking the ratio
    rather than the value makes this a test of the physics rather than a
    restatement of the code.
    """

    for r_over_R in (0.15, 0.4, 0.7, 1.0):
        for tsr in (4.0, 6.5, 9.0):
            phi = schmitz_inflow_angle(r_over_R, tsr)
            zero_swirl = math.atan(1.0 / (tsr * r_over_R))
            assert phi == pytest.approx((2.0 / 3.0) * zero_swirl, rel=1e-14)


def test_schmitz_chord_falls_monotonically_outboard():
    """A taper, and a steep one -- the shape the parameterisation must hold."""

    r_over_R = np.linspace(0.15, 1.0, 40)
    chord = schmitz_chord(r_over_R, 6.5, 2.0, 3, 1.2825)

    assert np.all(np.diff(chord) < 0.0)
    assert chord[0] > 4.0 * chord[-1], "Schmitz roots are wide; this one is not"


def test_schmitz_twist_is_inflow_minus_design_alpha():
    """
    theta = phi - alpha_design, in the convention station.py uses
    (alpha = phi - twist).

    Getting this sign wrong produces a blade that looks plausible and is wrong
    everywhere, so it is asserted rather than assumed.
    """

    alpha_design = math.radians(5.36)
    for r_over_R in (0.2, 0.6, 1.0):
        twist = schmitz_twist(r_over_R, 6.5, alpha_design)
        phi = schmitz_inflow_angle(r_over_R, 6.5)
        assert twist == pytest.approx(phi - alpha_design, rel=1e-14)


def test_chord_requires_positive_design_cl():
    with pytest.raises(ValueError):
        schmitz_chord(0.5, 6.5, 2.0, 3, 0.0)


# ---------------------------------------------------------------------------
# The baseline at the design point
# ---------------------------------------------------------------------------

def test_the_blade_operates_near_its_design_angle_of_attack(baseline, performance):
    """
    A Schmitz blade at its design tip-speed ratio should sit near
    alpha_design across the span -- that is what the construction is for.

    Inboard stations deviate most, which is expected: the wide root chord
    raises local solidity and induction there.
    """

    alphas = [station["alpha_deg"] for station in performance["spanwise"]]
    target = math.degrees(baseline.alpha_design_rad)

    assert max(alphas) == pytest.approx(target, abs=1.0)
    assert min(alphas) > target - 6.0


def test_peak_cp_lands_at_the_design_tip_speed_ratio(performance):
    """
    The strongest available check on the whole construction.

    The blade is designed for lambda = 6.5. If the Cp-lambda peak landed
    somewhere else, the Schmitz formulae, the twist convention or the solver
    wiring would be wrong -- and nothing else in this file would catch it.
    """

    design = load_design_rotor()
    best = max(performance["cp_lambda"], key=lambda row: row["Cp"])

    assert best["tsr"] == pytest.approx(design.design_tsr, abs=0.3)


def test_baseline_is_physically_plausible(performance):
    point = performance["design_point"]

    assert 0.35 < point["Cp"] < 16.0 / 27.0, "Cp implausible or above Betz"
    assert 0.0 < point["Ct"] < 1.2
    assert performance["root_bending_moment_nm"] > 0.0
    assert performance["peak_thrust_n"] > 0.0


def test_every_operating_point_converges(performance):
    """
    The whole below-rated line, and the Cp-lambda sweep.

    This is the test that would have caught the interpolant knot
    discontinuity found while building this step: a station at Re = 70,995
    failed here, with a genuinely discontinuous residual.
    """

    for row in performance["operating_line"]:
        assert row["converged"], (row["v_inf"], row["failed_stations"])
    for row in performance["cp_lambda"]:
        assert row["converged"], row["tsr"]


def test_aep_is_reported_as_outstanding_not_estimated(performance):
    """
    AEP needs the Weibull parameters, which are TODO.

    It must be absent and *named as absent*, not approximated -- an estimated
    AEP in the reference numbers would be indistinguishable from a real one to
    anyone reading the artefact later.
    """

    assert "aep_mwh_per_year" in performance["outstanding"]
    assert "BLOCKED" in performance["outstanding"]["aep_mwh_per_year"]
    assert "aep_mwh_per_year" not in performance


# ---------------------------------------------------------------------------
# Feasibility must not read as "passed" when it was never checked
# ---------------------------------------------------------------------------

def test_feasibility_reports_that_it_was_not_checked(baseline):
    """
    The bounds are still TODO, so the status must say so explicitly.

    A baseline that was never checked must not be indistinguishable from one
    that passed. When the manufacturability study lands this test fails
    deliberately, which is the reminder to replace it with a real check.
    """

    assert baseline.feasibility["checked"] is False
    assert "TODO" in baseline.feasibility["reason"]
    assert baseline.feasibility["violations"] == []


# ---------------------------------------------------------------------------
# The committed artefact
# ---------------------------------------------------------------------------

def test_committed_x0_matches_what_the_code_builds(baseline, committed_x0):
    """
    `verification/baseline/x0.json` is what Phases 2-5 start from.

    If the code and the artefact drift apart, every later result is measured
    against a baseline nothing reproduces.
    """

    n_chord = baseline.parameterisation.n_chord

    assert committed_x0["n_chord"] == n_chord
    assert committed_x0["n_twist"] == baseline.parameterisation.n_twist
    assert np.allclose(committed_x0["chord_control_points_m"],
                       baseline.design_vector[:n_chord], rtol=1e-12)
    assert np.allclose(committed_x0["twist_control_points_rad"],
                       baseline.design_vector[n_chord:], rtol=1e-12)


def test_committed_x0_records_its_construction(committed_x0):
    """The artefact has to be self-describing -- it outlives this session."""

    construction = committed_x0["construction"]

    assert construction["design_reynolds"] == DEFAULT_DESIGN_REYNOLDS
    assert construction["chord_rms_fit_error_m"] < 1e-3
    assert committed_x0["feasibility"]["checked"] is False
    assert committed_x0["rotor"]["design_tsr"] == load_design_rotor().design_tsr
