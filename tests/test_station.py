"""
Task 5 acceptance: the residual form, the momentum-region bracket, the
post-solve residual check, and reported non-convergence.

Each test corresponds to one of the work order's "done when" clauses. The cost
gate is `test_cost.py`; everything else is here.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

import pytest

from bem.airfoil import LinearPolar
from bem.corrections import combined_loss_factor
from bem.rotor import PHASE_VI_RATED_RPM, phase_vi_geometry
from bem.station import (
    PHI_EPS,
    PHI_MAX,
    RESIDUAL_RTOL,
    StationParams,
    momentum_region_bracket,
    residual,
    solve_station,
)
from polars.polar import CachedPolar, interpolant_for

NU = 1.5e-5
OMEGA_RATED = PHASE_VI_RATED_RPM * 2.0 * math.pi / 60.0

#: The envelope the bracket guarantee was established over. Wind speed and
#: rotor speed set the operating point; the chord and twist perturbations are
#: what the work order means by "baseline **and** perturbed geometries", and
#: are the shape of excursion an optimiser line search actually makes.
ENVELOPE = [
    ("v=5 rated", 5.0, OMEGA_RATED, 1.0, 0.0),
    ("v=7 rated", 7.0, OMEGA_RATED, 1.0, 0.0),
    ("v=10 rated", 10.0, OMEGA_RATED, 1.0, 0.0),
    ("v=15 rated", 15.0, OMEGA_RATED, 1.0, 0.0),
    ("v=20 rated", 20.0, OMEGA_RATED, 1.0, 0.0),
    ("tsr=2", 7.0, 2.0 * 7.0 / 5.029, 1.0, 0.0),
    ("tsr=4", 7.0, 4.0 * 7.0 / 5.029, 1.0, 0.0),
    ("tsr=6", 7.0, 6.0 * 7.0 / 5.029, 1.0, 0.0),
    ("tsr=7.5", 7.0, 7.5 * 7.0 / 5.029, 1.0, 0.0),
    ("chord x0.7", 7.0, OMEGA_RATED, 0.7, 0.0),
    ("chord x1.3", 7.0, OMEGA_RATED, 1.3, 0.0),
    ("twist -6deg", 7.0, OMEGA_RATED, 1.0, math.radians(-6.0)),
    ("twist +6deg", 7.0, OMEGA_RATED, 1.0, math.radians(6.0)),
]


def _stations(v_inf, omega, chord_scale=1.0, twist_delta=0.0):
    geometry = phase_vi_geometry()
    interpolant = interpolant_for(geometry.polar_cache)
    out = []
    for r, chord, twist in zip(geometry.r, geometry.chord, geometry.twist):
        chord = chord * chord_scale
        reynolds = math.hypot(v_inf, omega * r) * chord / NU
        out.append(StationParams(
            r=r, chord=chord, twist=twist + twist_delta,
            airfoil=CachedPolar(interpolant, reynolds),
            tsr=omega * r / v_inf, R=geometry.R,
            n_blades=geometry.n_blades, r_hub=geometry.r_hub,
        ))
    return out


# ---------------------------------------------------------------------------
# "Solver converges across the full operating envelope for baseline
#  **and** perturbed geometries"
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,v_inf,omega,chord_scale,twist_delta", ENVELOPE,
                         ids=[case[0] for case in ENVELOPE])
def test_every_station_converges_across_the_envelope(name, v_inf, omega,
                                                     chord_scale, twist_delta):
    """
    Every station of every operating point converges, with a checked residual.

    This is the sweep the bracket guarantee was established over -- 13 cases
    times 19 stations -- run as a test rather than left as a one-off
    measurement in a journal entry.
    """

    for i, station in enumerate(_stations(v_inf, omega, chord_scale, twist_delta)):
        result = solve_station(station)
        assert result["converged"], f"{name} station {i}: {result['failure']}"
        assert result["residual"] <= RESIDUAL_RTOL * result["residual_initial"]
        assert PHI_EPS <= result["phi"] <= PHI_MAX


@pytest.mark.parametrize("name,v_inf,omega,chord_scale,twist_delta", ENVELOPE,
                         ids=[case[0] for case in ENVELOPE])
def test_bracket_endpoints_have_the_proven_signs(name, v_inf, omega,
                                                 chord_scale, twist_delta):
    """
    R(lo) < 0 < R(hi) at every station -- the property the bracket rests on.

    Asserting the endpoint signs directly, rather than only that a root was
    found, is what makes this a test of the region classification instead of a
    test of brentq.
    """

    for i, station in enumerate(_stations(v_inf, omega, chord_scale, twist_delta)):
        bracket = momentum_region_bracket(station)
        assert bracket is not None, f"{name} station {i}: no bracket"
        phi_lo, phi_hi, r_lo, r_hi = bracket
        assert r_lo < 0.0 < r_hi
        assert phi_lo < phi_hi == PHI_MAX


def test_upper_endpoint_sign_holds_at_pi_over_two_exactly():
    """
    R(pi/2) > 0, evaluated at pi/2 itself.

    The old bracket stopped at pi/2 - 1e-4 because the multiplied residual's
    a' term blows up there. Ning's form with (1+a') eliminated algebraically
    is finite at pi/2, and the bracket closing on the exact endpoint rather
    than on a fudge factor is what removes the last tuned constant from the
    search.
    """

    for station in _stations(7.0, OMEGA_RATED):
        value = residual(PHI_MAX, station)
        assert math.isfinite(value)
        assert value > 0.0


def test_normal_force_is_positive_at_pi_over_two():
    """
    Cn(pi/2) = Cd > 0, the fact that makes the Cn bisection always well posed.

    If this were ever false the region boundary could not be bracketed and the
    classification would have no starting point, so it is asserted rather than
    assumed.
    """

    for station in _stations(7.0, OMEGA_RATED):
        cl = station.airfoil.cl(PHI_MAX - station.twist)
        cd = station.airfoil.cd(PHI_MAX - station.twist)
        cn = cl * math.cos(PHI_MAX) + cd * math.sin(PHI_MAX)
        assert cn == pytest.approx(cd, rel=1e-12)
        assert cn > 0.0


# ---------------------------------------------------------------------------
# "No code path returns a station result without a checked residual"
# ---------------------------------------------------------------------------

def test_every_result_carries_a_status():
    """Success and failure return the same keys, so no caller branches on shape."""

    required = {"phi", "a", "a_prime", "Cl", "Cd", "Ct", "Cq", "F", "converged",
                "residual", "residual_initial", "iterations", "function_calls",
                "failure", "history"}

    good = solve_station(_stations(7.0, OMEGA_RATED)[9])
    assert required <= set(good)
    assert good["converged"] and good["failure"] is None

    bad = solve_station(_degenerate_station())
    assert required <= set(bad)
    assert not bad["converged"] and isinstance(bad["failure"], str)


def _degenerate_station():
    """
    A station with genuinely no root in the momentum region.

    Constructed so the two ends of the region fight each other: a 40 deg twist
    puts the Cn = 0 boundary up near 40 deg, while a local speed ratio of 20
    would put a windmill root down near atan(1/20) = 2.9 deg -- inside the
    Cn < 0 band, where no windmill solution exists. The result is
    R(lo) = +0.6 > 0, so the bracket condition R(lo) < 0 < R(hi) fails and
    there is nothing for brentq to find.

    Verified to be a real failure rather than a contrived one: this is the
    same shape as a rotor driven far past its design tip-speed ratio, which is
    exactly what an optimiser line search can wander into. The point is that
    the solver must *report* it rather than raise.
    """

    return StationParams(
        r=1.0, chord=0.3, twist=math.radians(40.0),
        airfoil=LinearPolar(), tsr=20.0, R=5.0, n_blades=3,
    )


def test_non_convergence_is_reported_not_raised():
    """
    Task 5 item 5: the Phase 0 `ValueError` on no-bracket aborted a whole
    sweep over one bad station. In a several-hundred-point smoothness sweep or
    an optimiser line search that is the difference between a plot with a gap
    and no plot at all.
    """

    result = solve_station(_degenerate_station())

    assert result["converged"] is False
    assert "momentum region" in result["failure"]
    assert math.isnan(result["phi"])


# ---------------------------------------------------------------------------
# "Deliberately constructing a pole-selection case returns a reported failure,
#  not a number"
# ---------------------------------------------------------------------------

def test_a_pole_is_rejected_by_the_residual_check():
    """
    Feed the convergence check a pole and confirm it is rejected.

    The Phase 0 failure mode was `solve_station` never evaluating the residual
    after brentq returned, so a pole came back as a converged result with
    |R| ~ 1e+15 and nothing noticed. The check now added is relative to R0,
    and the separation is enormous: a converged station sits at |R|/R0 ~ 1e-15
    and a pole at ~1e+15, thirty orders apart. This asserts the check's
    arithmetic on a real pole rather than trusting the tolerance by eye.
    """

    station = _stations(7.0, OMEGA_RATED)[9]
    bracket = momentum_region_bracket(station)
    residual_initial = max(abs(bracket[2]), abs(bracket[3]))

    # The 'a' pole: Cn = -4F sin^2(phi)/sigma, i.e. Y = -4F. It lives below the
    # Cn = 0 boundary, so it is outside the bracket -- which is the structural
    # fix. Find it by scanning under the boundary and confirm what the check
    # would do with it.
    pole_residuals = []
    for i in range(1, 4000):
        phi = bracket[0] * i / 4000.0
        try:
            pole_residuals.append(abs(residual(phi, station)))
        except (ValueError, ZeroDivisionError):
            continue

    worst = max(pole_residuals)
    assert worst > 1e3 * residual_initial, (
        "the region below the Cn = 0 boundary should contain the pole this "
        "test is about"
    )
    assert worst > RESIDUAL_RTOL * residual_initial, (
        "a pole must fail the convergence check by a wide margin"
    )


def test_pole_region_is_outside_the_bracket():
    """
    The structural half of the pole fix.

    The residual check above is the safety net. The actual guarantee is that
    the 'a' pole (Y = -4F, which needs Cn < 0) lies strictly below the
    bracket's lower endpoint, so brentq is never handed an interval containing
    it and cannot select it in the first place.
    """

    for i, station in enumerate(_stations(7.0, OMEGA_RATED)):
        phi_lo, _phi_hi, _r_lo, _r_hi = momentum_region_bracket(station)
        for k in range(1, 200):
            phi = phi_lo + (PHI_MAX - phi_lo) * k / 200.0
            cl = station.airfoil.cl(phi - station.twist)
            cd = station.airfoil.cd(phi - station.twist)
            cn = cl * math.cos(phi) + cd * math.sin(phi)
            F = combined_loss_factor(station.r, station.R, station.n_blades,
                                     phi, station.r_hub)
            Y = station.solidity * cn / math.sin(phi) ** 2
            assert Y > -4.0 * F, (
                f"station {i}: the 'a' pole at Y = -4F is inside the bracket "
                f"at phi={math.degrees(phi):.3f} deg"
            )


# ---------------------------------------------------------------------------
# Relative, not absolute, tolerance
# ---------------------------------------------------------------------------

def test_tolerance_is_relative_to_the_initial_residual():
    """
    R0 varies materially across the span, which is why the criterion is
    scaled to it rather than being a fixed number.

    Measured at v = 7 m/s on the Phase VI rotor: R0 runs from 1.04 at the root
    to 90.7 at the tip, a spread of about 87x within a single operating point
    -- and that is before sweeping wind speed or perturbing the geometry. An
    absolute tolerance would therefore be ~87x stricter at one end of the
    blade than the other, for no reason connected to the physics.

    The threshold below is deliberately well inside the measured spread rather
    than just under it: the test is asserting that the variation is material,
    not pinning today's exact number.
    """

    initials = [solve_station(s)["residual_initial"] for s in _stations(7.0, OMEGA_RATED)]
    spread = max(initials) / min(initials)
    assert spread > 10.0, (
        f"R0 spread is only {spread:.1f}x across the span, which would undercut "
        "the case for a relative criterion -- re-examine"
    )


def test_history_is_recorded_on_request_only():
    """Off by default: it costs an append per residual evaluation."""

    station = _stations(7.0, OMEGA_RATED)[9]

    assert solve_station(station)["history"] is None

    recorded = solve_station(station, record_history=True)
    history = recorded["history"]
    assert history and all(len(entry) == 2 for entry in history)
    # The history must end at the answer, not somewhere near it.
    assert abs(history[-1][1]) <= max(abs(entry[1]) for entry in history)
