"""
The staged validators of the BEM solver, ported to pytest.

`src/validation/validate_stage1..4.py` were nine executable scripts in
print-and-assert style, run by hand from `src/`. The physics was good but the
harness wasn't: around nine acceptance criteria are machine-checkable claims, and those cannot live in a script
whose known-good state is a non-zero exit.

This file is the physics, unchanged. Every assertion below is the assertion
the corresponding script made, with the same tolerances and the same reference
values; what has gone is the printing, the manual invocation and the `main()`
that had to be remembered. Where a check's *reasoning* has moved on since it
was written -- the bracket scan and the branch selection have both been replaced --
the docstring says so rather than leaving a stale explanation attached to a
still-correct number.

Stage map:
    1  single-station solve vs the classical Hansen fixed-point iteration
    2  Prandtl tip/hub loss: the F -> 1 limits, bounds, and direction
    3  Glauert/Buhl high-thrust: Ct(a) shape, the low-induction identity,
       regression against stages 1-2, and a real turbulent-wake station
    4  the spanwise loop: smoothness, bounds, rotor-level plausibility, and
       that the loop reproduces a direct single-station solve exactly

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math

import pytest

from bem.airfoil import LinearPolar
from bem.corrections import (
    BUHL_AC,
    corrected_axial_induction,
    high_thrust_correction,
    hub_loss_factor,
    tip_loss_factor,
)
from bem.rotor import RotorGeometry, demo_rotor_geometry, solve_rotor
from bem.station import StationParams, residual, solve_station
from config import load_phase_vi_rotor
from polars.polar import CachedPolar, interpolant_for

#: The synthetic geometry stages 1-3 are defined on. R = 1000 puts r/R ~ 0.005,
#: far enough from the tip that F is 1.0 to double precision -- which is what
#: makes the uncorrected numbers a valid reference even after the tip-loss
#: correction made R a required parameter.
CHORD = 0.3
TWIST = math.radians(5.0)
TSR = 5.0
R_TIP = 5.5

#: The uncorrected result, carried forward as a regression anchor for the
#: corrected solver as well. These are this solver's own values, not an external
#: reference -- their role is to catch drift, not to validate physics.
STAGE_1_REFERENCE = {
    "phi": math.radians(10.013235),
    "a": 0.11399664,
    "a_prime": 0.00359930,
}

_AIR = load_phase_vi_rotor()


def _station(r, chord=CHORD, twist=TWIST, tsr=TSR, n_blades=3, R=R_TIP, r_hub=None):
    return StationParams(r=r, chord=chord, twist=twist, airfoil=LinearPolar(),
                         tsr=tsr, R=R, n_blades=n_blades, r_hub=r_hub)


# ===========================================================================
# Single-station solve
# ===========================================================================

def _hansen_fixed_point(station, tol=1e-10, max_iter=200):
    """
    Classical uncorrected BEM fixed-point iteration (Hansen 2008, Ch. 6).

    Deliberately a *different solution method* for the same equations: direct
    substitution on (a, a') rather than a root-find on a single residual in
    phi. With no corrections applied the two are algebraically equivalent, so
    agreement is a real cross-check of the Ning reduction rather than a
    restatement of it.
    """

    a, a_prime = 0.0, 0.0
    sigma = station.solidity
    for _ in range(max_iter):
        phi = math.atan((1 - a) / ((1 + a_prime) * station.tsr))
        alpha = phi - station.twist
        cl = station.airfoil.cl(alpha)
        cd = station.airfoil.cd(alpha)
        cn = cl * math.cos(phi) + cd * math.sin(phi)
        ct = cl * math.sin(phi) - cd * math.cos(phi)

        a_new = 1.0 / ((4 * math.sin(phi) ** 2) / (sigma * cn) + 1.0)
        a_prime_new = 1.0 / ((4 * math.sin(phi) * math.cos(phi)) / (sigma * ct) - 1.0)

        converged = abs(a_new - a) < tol and abs(a_prime_new - a_prime) < tol
        a, a_prime = a_new, a_prime_new
        if converged:
            break

    return {"phi": math.atan((1 - a) / ((1 + a_prime) * station.tsr)),
            "a": a, "a_prime": a_prime}


def test_stage1_induction_within_physical_bounds():
    """0 < a < 0.5 and a' small for normal operation."""

    result = solve_station(_station(r=5.0, R=1000.0))

    assert result["converged"]
    assert 0.0 < result["a"] < 0.5
    assert -0.5 < result["a_prime"] < 0.5


def test_stage1_agrees_with_the_classical_hansen_fixed_point():
    """
    The Ning reduction and the classical fixed-point iteration must agree.

    They are the same equations solved two different ways, so this is the
    check that the single-residual formulation did not quietly change the
    physics while changing the solution method.
    """

    station = _station(r=5.0, R=1000.0)
    ning = solve_station(station)
    hansen = _hansen_fixed_point(station)

    for key in ("phi", "a", "a_prime"):
        assert abs(ning[key] - hansen[key]) < 1e-6, key


def test_stage1_residual_is_smooth_over_the_physical_bracket():
    """
    No discontinuity in R(phi) in a +/-5 deg window around the root.

    The original comment here pointed at the Ct = 0 pole and noted the window
    was chosen to stay clear of it. With Ning's residual that is structural rather
    than a matter of choosing the window: Ning's residual has no a' pole at
    all, and the 'a' pole requires Cn < 0, which the momentum-region bracket
    excludes. The window is kept as written so the measurement is comparable.
    """

    station = _station(r=5.0, R=1000.0)
    phi_star = solve_station(station)["phi"]

    lo = max(1e-3, phi_star - math.radians(5))
    hi = phi_star + math.radians(5)
    samples = [residual(lo + (hi - lo) * i / 24, station) for i in range(25)]
    worst = max(abs(samples[i] - samples[i - 1]) for i in range(1, len(samples)))

    assert worst < 1.0


# ===========================================================================
# Prandtl tip/hub loss
# ===========================================================================

def test_stage2_many_blades_recovers_no_tip_loss():
    """More blades means less tip-vortex loss per blade, so F -> 1."""

    result = solve_station(_station(r=5.0, n_blades=200))
    assert abs(result["F"] - 1.0) < 1e-3


def test_stage2_midspan_limit_reproduces_stage1():
    """
    At r/R ~ 0.005, F is 1 to double precision and the corrected solver must
    return the uncorrected numbers.
    """

    result = solve_station(_station(r=5.0, R=1000.0))

    assert abs(result["F"] - 1.0) < 1e-9
    for key, expected in STAGE_1_REFERENCE.items():
        assert abs(result[key] - expected) < 1e-5, key


def test_stage2_loss_factors_are_bounded_and_monotonic():
    """F_tip and F_hub in (0, 1], decreasing toward the tip and away from it."""

    phi = math.radians(10.0)
    r_hub = 0.5
    radii = [r for r in (0.5 + 0.05 * i for i in range(1, 100)) if r < R_TIP]

    f_tip = [tip_loss_factor(r, R_TIP, 3, phi) for r in radii]
    f_hub = [hub_loss_factor(r, r_hub, 3, phi) for r in radii]

    assert all(0.0 < value <= 1.0 for value in f_tip)
    assert all(0.0 < value <= 1.0 for value in f_hub)
    assert all(f_tip[i] <= f_tip[i - 1] + 1e-12 for i in range(1, len(f_tip)))
    assert all(f_hub[i] >= f_hub[i - 1] - 1e-12 for i in range(1, len(f_hub)))


def test_stage2_tip_loss_shifts_induction_at_the_tip_not_midspan():
    """
    The correction has to act where it is supposed to act.

    A tip-loss model that changed mid-span answers by as much as tip answers
    would be wrong even if every individual number looked plausible, so this
    compares each station against its own F = 1 reference rather than against
    a fixed value.
    """

    near_tip = solve_station(_station(r=5.4))
    midspan = solve_station(_station(r=2.75))
    ref_tip = solve_station(_station(r=5.4, R=5000.0))
    ref_mid = solve_station(_station(r=2.75, R=5000.0))

    delta_tip = abs(near_tip["a"] - ref_tip["a"])
    delta_mid = abs(midspan["a"] - ref_mid["a"])

    assert near_tip["F"] < 0.99
    assert midspan["F"] > 0.999
    assert near_tip["a"] > ref_tip["a"], "tip loss must increase a near the tip"
    assert delta_tip > delta_mid
    assert delta_mid < 1e-3


# ===========================================================================
# Glauert/Buhl high-thrust correction
# ===========================================================================

def test_stage3_ct_curve_is_monotonic_across_the_blend():
    """Ct(a) rises throughout, including through the switch at a = 0.4."""

    F = 0.85
    values = []
    for i in range(400):
        a = 1e-4 + (0.999 - 1e-4) * i / 399
        values.append(high_thrust_correction(a, 4.0 * a * F * (1.0 - a), F))

    worst_backstep = max(values[i - 1] - values[i] for i in range(1, len(values)))
    assert worst_backstep < 1e-9


def test_stage3_low_induction_is_bit_identical_to_the_plain_relation():
    """
    Below the blend the corrected solve must be the tip-loss closed form
    exactly -- not almost.

    Bit-identity is the point: it is what makes the Buhl correction a pure
    addition rather than a change to low-induction behaviour, and the gamma form keeps
    `Y / (4F + Y)` written exactly that way (rather than the algebraically
    equal kappa/(1+kappa)) to preserve it through the gamma reparameterisation.
    """

    for F in (0.3, 0.6, 1.0):
        for Y in (0.01, 0.5, 1.0, 1.5):
            plain = Y / (4.0 * F + Y)
            if plain <= BUHL_AC:
                assert corrected_axial_induction(Y, F) == plain


def test_stage3_reproduces_the_stage1_and_stage2_references():
    """Adding the correction must not move stations that stay below the blend."""

    stage1 = solve_station(_station(r=5.0, R=1000.0))
    for key, expected in STAGE_1_REFERENCE.items():
        assert abs(stage1[key] - expected) < 1e-5, key

    near_tip = solve_station(_station(r=5.4))
    assert abs(near_tip["F"] - 0.3755) < 1e-3
    assert abs(near_tip["a"] - 0.236561) < 1e-4


def test_stage3_turbulent_wake_station_converges_and_stays_bounded():
    """
    A station actually in the turbulent-wake state (a > 0.4).

    High solidity drives the loading up until the naive momentum value exceeds
    the blend, which is the regime plain momentum theory gets wrong.
    """

    result = solve_station(
        _station(r=2.0, chord=1.4, twist=math.radians(2), tsr=6.0, R=8.0))

    assert result["converged"]
    assert result["a"] > BUHL_AC, "case did not reach the turbulent wake state"
    assert 0.0 < result["a"] < 1.0
    assert math.isfinite(result["phi"])


def test_stage3_residual_stays_smooth_around_a_turbulent_wake_root():
    """
    Ct(a) being C1 is not enough on its own; what the root-finder sees is the
    residual, so the smoothness has to survive composition into it.
    """

    station = _station(r=2.0, chord=1.4, twist=math.radians(2), tsr=6.0, R=8.0)
    phi_star = solve_station(station)["phi"]

    lo, hi = phi_star - math.radians(0.5), phi_star + math.radians(0.5)
    samples = [residual(lo + (hi - lo) * i / 39, station) for i in range(40)]
    worst = max(abs(samples[i] - samples[i - 1]) for i in range(1, len(samples)))

    assert worst < 0.05


# ===========================================================================
# The spanwise loop
# ===========================================================================

def _demo_solve():
    geometry = demo_rotor_geometry()
    return geometry, solve_rotor(
        geometry, tsr=5.0, air_density=_AIR.air_density,
        kinematic_viscosity=_AIR.kinematic_viscosity)


def test_stage4_spanwise_quantities_vary_smoothly():
    """
    No spike between adjacent stations in any tracked quantity.

    A discontinuity/spike check, not a tight tolerance: chord, twist and
    Reynolds all vary smoothly by construction, so a real solver bug shows up
    as an obvious jump. The Cd threshold is 0.25 rather than the original 0.1
    because the alpha clamp that used to flatten this blade's stalled root
    stations onto Cd(18 deg) has been removed. This blade's three innermost
    stations sit at alpha = 31.6, 25.6 and 20.3 deg at tsr = 5 and always did;
    the clamp flattened all three onto Cd(18 deg) = 0.079, so a 0.1 limit was
    unreachable by construction and the check was, in effect, testing the
    clamp. Real post-stall drag rises 0.437 -> 0.273 -> 0.149 across ~6 deg
    alpha steps per station -- steep, monotone, max step 0.165 -- against a
    ~0.41 full-span Cd range.
    """

    limits = {"alpha": math.radians(15), "Cl": 0.5, "Cd": 0.25,
              "a": 0.2, "a_prime": 0.1}
    _geometry, result = _demo_solve()
    stations = result["stations"]

    for key, limit in limits.items():
        worst = max(abs(stations[i][key] - stations[i - 1][key])
                    for i in range(1, len(stations)))
        assert worst < limit, f"spanwise discontinuity in {key}: {worst}"


def test_stage4_every_station_is_physically_bounded():
    _geometry, result = _demo_solve()

    for station in result["stations"]:
        assert station["converged"]
        assert 0.0 < station["a"] < 1.0
        assert 0.0 < station["F"] <= 1.0
        assert math.isfinite(station["phi"])


def test_stage4_integrated_rotor_is_plausible():
    """Cp below Betz, Ct in a sane range -- the coarsest possible sanity net."""

    _geometry, result = _demo_solve()

    assert 0.0 < result["Cp"] < 16.0 / 27.0
    assert 0.0 < result["Ct"] < 1.5
    assert result["converged"]


def test_stage4_loop_reproduces_a_direct_single_station_solve():
    """
    The spanwise loop must add nothing of its own.

    Exact equality, not a tolerance: `solve_rotor` is a loop over
    `solve_station` plus an integration, so a single-station geometry must
    come back bit-identical to calling the station solver directly. Any
    difference is wiring, not physics.
    """

    r, chord, twist_deg, R, tsr_tip, v_inf = 3.0, 0.3, 6.0, 8.0, 5.0, 7.0
    geometry = RotorGeometry(r=[r], chord=[chord], twist=[math.radians(twist_deg)],
                             R=R, polar_cache="s809", n_blades=3)

    via_loop = solve_rotor(
        geometry, tsr=tsr_tip, v_inf=v_inf, air_density=_AIR.air_density,
        kinematic_viscosity=_AIR.kinematic_viscosity)["stations"][0]

    omega = tsr_tip * v_inf / R
    reynolds = math.hypot(v_inf, omega * r) * chord / _AIR.kinematic_viscosity
    direct = solve_station(StationParams(
        r=r, chord=chord, twist=math.radians(twist_deg),
        airfoil=CachedPolar(interpolant_for("s809"), reynolds),
        tsr=omega * r / v_inf, R=R, n_blades=3))

    for key in ("phi", "a", "a_prime", "Cl", "Cd"):
        assert via_loop[key] == direct[key], key
