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
from bem.rotor import solve_rotor
from config import load_site
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

def test_the_blade_operates_near_its_design_angle_of_attack(baseline):
    """
    A Schmitz blade at its design tip-speed ratio should sit near
    alpha_design across the span -- that is what the construction is for.

    Inboard stations deviate most, which is expected: the wide root chord
    raises local solidity and induction there.

    THIS TEST NOW SOLVES AT `lambda_design` ITSELF, which it did not have to
    do before 2026-09-19. It used to read the artefact's `spanwise` block,
    because that block *was* the design point (11 m/s, lambda = 6.5, 341 rpm).
    Since B1 (300 rpm) the machine's design condition is the tip-speed ceiling
    -- `lambda(11) = 5.71`, 300 rpm -- so the artefact evaluates there and its
    alphas run 2.1..8.7 deg, off the construction's design point by design.
    Reading `spanwise` here would test the schedule, not the Schmitz
    construction, so the solve is explicit and the claim stays about the
    blade. The artefact's own design point is pinned separately in
    `test_the_design_point_is_the_rated_speed_on_the_schedule`.
    """

    design, site = load_design_rotor(), load_site()
    geometry = baseline.to_geometry()
    result = solve_rotor(geometry, tsr=design.design_tsr,
                         v_inf=design.rated_wind_speed_ms,
                         air_density=site.air_density,
                         kinematic_viscosity=site.kinematic_viscosity)
    alphas = [math.degrees(station["alpha"]) for station in result["stations"]]
    target = math.degrees(baseline.alpha_design_rad)

    assert max(alphas) == pytest.approx(target, abs=1.0)
    assert min(alphas) > target - 6.0


def test_the_design_point_is_the_rated_speed_on_the_schedule(performance):
    """
    The reference design point is the machine's, not the blade's.

    `baseline_reference.json` used to solve its design point at
    `design.design_tsr` (341 rpm at 11 m/s), a rotor speed the 300 rpm
    machine cannot reach. It is now the schedule's lambda at the rated wind
    speed, which for this machine is the tip-speed ceiling -- so this pins
    both the law and the reason the reference numbers moved.
    """

    design = load_design_rotor()
    expected = design.max_tip_speed_ms / design.rated_wind_speed_ms
    point = performance["design_point"]

    assert expected < design.design_tsr, "300 rpm must actually bind at 11 m/s"
    assert point["v_inf"] == pytest.approx(design.rated_wind_speed_ms)
    assert point["tsr"] == pytest.approx(expected, rel=1e-12)
    assert point["rpm"] == pytest.approx(design.max_rotor_speed_rpm, rel=1e-9)
    assert point["converged"]


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


def test_aep_is_reported_as_a_real_number(performance):
    """
    Replaces `test_aep_is_reported_as_outstanding_not_estimated`.

    WHY THE OLD TEST WAS NOT A GOOD REMINDER, recorded because the lesson
    generalises. It asserted that `outstanding["aep_mwh_per_year"]` contained
    the string "BLOCKED". But that string was a literal in `baseline.py`, which
    did not import the objective at all -- so the test and the code were both
    *describing* a blockage rather than being blocked by one. When the Weibull
    parameters landed on 2026-09-13 it stayed green, and nothing in the suite
    said the data had arrived. A placeholder guard has to be attached to the
    thing that actually changes: here, `WeibullResource.from_config()` raising.

    What is asserted now is that AEP is present, positive, finite, and no
    longer listed as outstanding.
    """

    assert "aep_mwh_per_year" not in performance["outstanding"]

    aep = performance["aep_mwh_per_year"]
    assert math.isfinite(aep)
    assert aep > 0.0


def test_the_baseline_aep_meets_the_revised_sanity_band(performance):
    """
    Plan 1.4's exit criterion, against the band as revised on 2026-09-13.

    THE BAND WAS WIDENED, 4-6 -> 8-12 MWh/yr, on MJ's explicit instruction, and
    that is a deliberate revision of a stated expectation rather than a
    tolerance tweak. The full basis is in `config/rotor_design.yaml` beside the
    numbers and in the 2026-09-13 journal entry; the short version is that the
    old band was defective in two independent ways, both predating the data
    that exposed them:

      1. It compared different quantities -- a capacity factor on a ~3.1 kW
         ELECTRICAL rating (Cp = 0.42, eta = 0.90) against a model that
         produces AERODYNAMIC shaft energy at the solver's actual Cp of 0.472
         with no drivetrain efficiency anywhere in the chain.
      2. It was never consistent with plan 1.3's own resource prior. At the
         calmest corner of that prior this rotor already returns 5.9 MWh/yr.

    So the band could not have been met by the design it was written for. It is
    now derived from the resource uncertainty this project has actually
    recorded (+/-10 % on `c` and `k` spans 8.21-12.01) rather than from a
    capacity factor borrowed from generic small turbines.
    """

    low, high = performance["aep_sanity_band_mwh_per_year"]
    assert (low, high) == (8.0, 12.0), (
        "plan 1.4's sanity band changed again. It was set to 8-12 on "
        "2026-09-13 with its derivation written out in "
        "config/rotor_design.yaml. Changing it is a deliberate act that needs "
        "the same treatment, not a tolerance tweak to make a test pass."
    )

    assert performance["aep_in_sanity_band"] is True
    assert performance["aep_mwh_per_year"] == pytest.approx(10.27, abs=0.1)


def test_the_sanity_band_is_still_narrow_enough_to_catch_a_bug(performance):
    """
    A band widened until everything fits is a note, not a check.

    This is the test that keeps the widening honest. It reconstructs the
    failure modes the band exists to catch and asserts each still lands
    outside it -- so if the band is ever widened again to accommodate some
    future result, the cost of that is visible here rather than silent.

    The modes are the ones an AEP chain actually gets wrong: forgetting the
    rated-power limit, using sea-level density at an 1800 m site, taking bin
    masses from the PDF without normalising, and any factor-of-two slip.
    """

    low, high = performance["aep_sanity_band_mwh_per_year"]
    aep = performance["aep_mwh_per_year"]

    # Ratios measured on 2026-09-13 against the correct 10.27 MWh/yr; see the
    # journal entry. Applied to the live AEP so this tracks the real chain.
    failure_modes = {
        "rated-power limiting not applied": aep * 1.424,
        "sea-level density (1.225) at an 1800 m site": aep * 1.265,
        "bin masses not normalised": aep / 0.79889,
        "factor of two low": aep / 2.0,
        "factor of two high": aep * 2.0,
    }

    for description, wrong in failure_modes.items():
        assert not (low <= wrong <= high), (
            f"the sanity band no longer catches: {description} "
            f"({wrong:.2f} MWh/yr sits inside {low}-{high}). The band has been "
            f"widened past the point of being a check."
        )


# ---------------------------------------------------------------------------
# Feasibility is checked, and x0 passes it unclipped
# ---------------------------------------------------------------------------

def test_x0_is_feasible_against_the_configured_bounds(baseline):
    """
    Plan step 1.7: feasibility *recorded*. Until 2026-09-19 the bounds were
    TODO and this test asserted `checked is False`; with `chord_max_m =
    0.30 m` in config the check runs, and the Schmitz root control point
    (276 mm) sits inside it -- the fairness argument of the audit's section
    3.4 depends on the baseline being unclipped, so that is asserted, not
    merely reported.
    """

    feasibility = baseline.feasibility
    assert feasibility["checked"] is True
    assert feasibility["reason"] is None
    assert feasibility["clipped"] is False
    assert feasibility["violations"] == []
    assert np.allclose(feasibility["clipped_vector"], baseline.design_vector)


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


def test_configured_rating_is_the_baselines_aerodynamic_rated_power(baseline):
    """
    `operating.rated_power_w` in `config/rotor_design.yaml` is provisional:
    it stands in for the generator nameplate (outstanding input B2) and its
    stated basis is `P_aero(V_rated; x0)` -- the Schmitz baseline's own
    aerodynamic power at the rated wind speed, at full precision so freezing
    the rating left the baseline's AEP unchanged.

    This pins the config to that basis, so the number cannot be edited
    without the basis being restated. **When B2 lands and the nameplate
    replaces it, this test is deleted** (not loosened): the rating then has a
    basis outside the code and nothing here can check it.
    """

    from objective.power import aerodynamic_power
    from config import load_site

    design, site = load_design_rotor(), load_site()
    geometry = baseline.parameterisation.to_geometry(baseline.design_vector)
    p_aero, result = aerodynamic_power(geometry, design.rated_wind_speed_ms,
                                       design.design_tsr, site.air_density,
                                       site.kinematic_viscosity)
    assert result["converged"]
    assert design.rated_power_w == pytest.approx(p_aero, rel=1e-12)


def test_committed_x0_records_its_construction(committed_x0):
    """The artefact has to be self-describing -- it outlives this session."""

    construction = committed_x0["construction"]

    assert construction["design_reynolds"] == DEFAULT_DESIGN_REYNOLDS
    assert construction["chord_rms_fit_error_m"] < 1e-3
    assert committed_x0["feasibility"]["checked"] is True
    assert committed_x0["feasibility"]["reason"] is None
    assert committed_x0["feasibility"]["clipped"] is False
    assert committed_x0["feasibility"]["violations"] == []
    assert committed_x0["rotor"]["design_tsr"] == load_design_rotor().design_tsr
