"""
Task 5's hard gate: one operating point under 0.1 s.

This is the criterion the work order calls hard, and it is hard for a reason
that is about feasibility rather than speed: the Step 8 objective-smoothness
sweep is roughly 3000 evaluations x 18 wind-speed bins. At the audit's
measured 9.0 s per operating point that is ~5.5 days of compute per pass, and
the brief names Step 8 the single most valuable step in Phase 1. At that cost
it gets skipped or truncated, which is how a gradient-based project ends up
unable to say whether its objective is smooth.

    9.0 s    audit baseline (2000-point bracket scan x 120 us polar lookup)
    0.56 s   after Task 4   (CachedPolar over the C1 interpolant, ~6 us)
    ~0.015 s after Task 5   (Ning bracket: ~48 residual evaluations/station)

So the gate is met with roughly 6.7x of headroom. The work order is explicit
that this is not to be optimised past -- "the gate exists to make Step 8
feasible, not to be fast for its own sake" -- so this test asserts the gate,
not the current number.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
import time

import pytest

from bem.rotor import PHASE_VI_RATED_RPM, phase_vi_geometry, solve_rotor
from config import load_phase_vi_rotor

#: The gate, seconds per operating point.
COST_TARGET_S = 0.1

#: Repeats. Best-of, not mean: a wall-clock measurement on a shared machine
#: has a hard floor and an unbounded tail, so the minimum is the robust
#: estimator of the cost and the mean mostly measures what else is running.
REPEATS = 5

_PHASE_VI = load_phase_vi_rotor()


def _operating_point():
    geometry = phase_vi_geometry()
    omega = PHASE_VI_RATED_RPM * 2.0 * math.pi / 60.0
    return solve_rotor(
        geometry, tsr=omega * geometry.R / 7.0, v_inf=7.0,
        air_density=_PHASE_VI.air_density,
        kinematic_viscosity=_PHASE_VI.kinematic_viscosity,
    )


def test_operating_point_is_under_the_cost_target():
    """
    One Phase VI operating point, 19 stations, in under 0.1 s.

    The first call is discarded: it pays for the polar cache load and the
    spline construction (~15 ms, memoised per cache by
    `polars.polar.interpolant_for`), which a sweep pays once and not per
    point.
    """

    _operating_point()  # warm the memoised interpolant

    best = math.inf
    for _ in range(REPEATS):
        start = time.perf_counter()
        result = _operating_point()
        best = min(best, time.perf_counter() - start)

    assert result["converged"], result["failed_stations"]
    assert best < COST_TARGET_S, (
        f"operating point took {best * 1000:.1f} ms, over the {COST_TARGET_S * 1000:.0f} ms "
        f"gate. The Step 8 smoothness sweep is ~3000 evaluations x 18 bins, so "
        f"this multiplies straight into whether that step is feasible."
    )


def test_residual_evaluations_per_station_are_bounded():
    """
    The structural reason the gate is met, asserted directly.

    Wall-clock is machine-dependent; the evaluation count is not. The Phase 0
    scan cost 2000 residual evaluations per station before brentq even
    started. The bracket now costs a bounded Cn bisection plus ~10 brentq
    iterations, and this pins that so a future change that reintroduces a scan
    fails here with a clear reason rather than as a slow test somewhere else.
    """

    result = _operating_point()
    worst = max(s["iterations"] for s in result["stations"])
    assert worst < 60, f"{worst} root-finder iterations at some station"
