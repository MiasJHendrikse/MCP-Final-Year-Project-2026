"""
`validate_powercurve.py`'s six checks, ported to pytest (work order Task 7).

`bem/powercurve.py` adds no aerodynamics of its own -- it sweeps
`rotor.solve_rotor` and converts to dimensional units -- so everything here is
about the sweep layer being transparent: that the dimensional conversions
invert exactly, that the two machine conventions stay distinguishable, and
that nothing drifts between calling the sweep and calling the solver directly.

Two of the swept ranges are narrower than the original script used, because
Task 4 removed the Reynolds clamp and a sweep now stops where the S809 cache's
coverage actually stops instead of running on a clamped Reynolds number. The
`PHASE_VI_MAX_CACHED_WIND_SPEED_MS` constant below and the comments on the
check 3 and check 4 tests carry the arithmetic.

`PHASE_VI_SEQUENCE_S_REFERENCE` and `PHASE_VI_MAX_CACHED_WIND_SPEED_MS` live
here rather than in `src/` deliberately: they are regression anchors for this
suite, not inputs to the solver, and `validate_powercurve.py` -- their previous
home -- was retired by this task.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

import pytest

from bem.powercurve import (
    BETZ_LIMIT,
    cp_lambda_curve,
    omega_to_rpm,
    operating_point,
    peak_cp,
    power_curve,
    rpm_to_omega,
    swept_area,
)
from bem.rotor import (
    PHASE_VI_RATED_RPM,
    demo_rotor_geometry,
    phase_vi_geometry,
    solve_rotor,
)
from config import load_phase_vi_rotor

#: NREL Phase VI (Sequence S) at 71.63 RPM. v_inf -> (Cp, Ct). These are this
#: project's own solver's values, not experimental data -- a regression
#: anchor, not a validation against ground truth.
#:
#: Re-anchored 2026-09-10 by work order Task 4, which deleted the alpha clamp
#: in the polar adapter. The values below 10 m/s barely moved (the whole span
#: is attached, so only the bilinear -> C1 interpolant change shows); above
#: it they moved a lot, because that is where stations run past the S809
#: cache's +18 deg converged band and the old adapter returned Cl(18 deg) for
#: them instead of the Viterna extrapolation. Previous anchor, for the record:
#:
#:      v      old Cp    new Cp      old Ct    new Ct
#:      5.0    0.4062    0.4062      0.6777    0.6777
#:      7.0    0.3670    0.3671      0.5415    0.5416
#:     10.0    0.2587    0.2566      0.3653    0.3648
#:     13.0    0.1830    0.1582      0.2518    0.2505
#:     15.0    0.1501    0.1123      0.2016    0.2025
#:     20.0    0.1017    0.0533      0.1298    0.1414
#:
#: The old post-stall figures were high because clamped lift does not fall
#: off and clamped drag does not rise. See tests/golden/README.md's change
#: log for the same movement measured station by station.
PHASE_VI_SEQUENCE_S_REFERENCE = {
    5.0: (0.4062, 0.6777),
    7.0: (0.3671, 0.5416),
    10.0: (0.2566, 0.3648),
    13.0: (0.1582, 0.2505),
    15.0: (0.1123, 0.2025),
    20.0: (0.0533, 0.1414),
    25.0: (0.0760, 0.0948),  # stale, and not swept -- see
                             # PHASE_VI_MAX_CACHED_WIND_SPEED_MS. Kept as the
                             # record of a Sequence S point this cache cannot
                             # serve, not as a value to compare against.
}

#: Highest Sequence S wind speed the S809 cache can actually serve at
#: 71.63 RPM, and the reason the 25 m/s reference point above is recorded but
#: not swept.
#:
#: The root station (r = 1.2575 m, chord 0.737 m) has the largest chord, so at
#: a fixed rotor speed it is the *root* that runs out of cache first as wind
#: speed rises: at 25 m/s its Reynolds estimate is 1,312,857 against the
#: cache's 1,300,000 ceiling -- over by 1 %. At 20 m/s the whole span sits at
#: 1,014,447-1,125,817 and is comfortably inside.
#:
#: Before work order Task 4 that station silently clamped to the 1.3M curve
#: and the sweep ran to 25 m/s on a Reynolds number that was not the
#: station's. It now raises, so the sweep stops where the cache's coverage
#: stops. Extending the cache above 1.3M would be a fresh XFOIL build for a
#: band the Phase VI validation does not otherwise need; the deliberate
#: choice is to record the limit rather than paper over it.
PHASE_VI_MAX_CACHED_WIND_SPEED_MS = 20.0


_AIR = load_phase_vi_rotor()
AIR = {"air_density": _AIR.air_density,
       "kinematic_viscosity": _AIR.kinematic_viscosity}


# ---------------------------------------------------------------------------
# Check 1 -- dimensional round-trip
# ---------------------------------------------------------------------------

def test_dimensional_conversions_invert_exactly():
    """
    Power/thrust/torque invert back to the Cp/Ct they were derived from.

    Exact equality: Cp and Ct are *defined* against these references, so
    inverting them is an identity rather than an independent estimate. A
    tolerance here would hide a genuine algebra error.
    """

    geometry = phase_vi_geometry()
    area = swept_area(geometry)
    curve = power_curve(geometry, [5.0, 10.0, 20.0], rpm=PHASE_VI_RATED_RPM, **AIR)

    for point in curve["points"]:
        v = point["v_inf"]
        cp = point["power"] / (0.5 * AIR["air_density"] * v ** 3 * area)
        ct = point["thrust"] / (0.5 * AIR["air_density"] * v ** 2 * area)
        assert cp == pytest.approx(point["Cp"], rel=1e-12)
        assert ct == pytest.approx(point["Ct"], rel=1e-12)
        assert point["torque"] * point["omega"] == pytest.approx(point["power"], rel=1e-12)


def test_rpm_and_omega_round_trip():
    for rpm in (1.0, 71.63, 600.0):
        assert abs(omega_to_rpm(rpm_to_omega(rpm)) - rpm) < 1e-12


# ---------------------------------------------------------------------------
# Check 2 -- fixed-speed regression against the Sequence S anchor
# ---------------------------------------------------------------------------

def test_fixed_speed_reproduces_the_phase_vi_anchor():
    """
    The pinned Phase VI values, re-anchored by Task 4 when the alpha clamp was
    deleted. These are this solver's own numbers -- a drift detector, not a
    validation against experiment.
    """

    geometry = phase_vi_geometry()
    speeds = sorted(v for v in PHASE_VI_SEQUENCE_S_REFERENCE
                    if v <= PHASE_VI_MAX_CACHED_WIND_SPEED_MS)
    curve = power_curve(geometry, speeds, rpm=PHASE_VI_RATED_RPM, **AIR)

    assert curve["mode"] == "fixed-speed"
    for point in curve["points"]:
        cp_ref, ct_ref = PHASE_VI_SEQUENCE_S_REFERENCE[point["v_inf"]]
        assert abs(point["Cp"] - cp_ref) < 5e-5, point["v_inf"]
        assert abs(point["Ct"] - ct_ref) < 5e-5, point["v_inf"]
        assert abs(point["rpm"] - PHASE_VI_RATED_RPM) < 1e-9


def test_fixed_speed_holds_rotor_speed_so_tsr_falls():
    """
    The distinguishing property of the fixed-speed convention.

    Mixing the two conventions up is the single easiest way to produce a
    meaningless comparison -- the 2026-07-28 entry records exactly that
    mistake faking a 40-degree AoA disagreement against QBlade -- so each
    convention gets a test that fails if the modes were swapped.
    """

    geometry = phase_vi_geometry()
    speeds = sorted(v for v in PHASE_VI_SEQUENCE_S_REFERENCE
                    if v <= PHASE_VI_MAX_CACHED_WIND_SPEED_MS)
    tsrs = power_curve(geometry, speeds, rpm=PHASE_VI_RATED_RPM, **AIR)["tsr"]

    assert all(tsrs[i] > tsrs[i + 1] for i in range(len(tsrs) - 1))


# ---------------------------------------------------------------------------
# Check 3 -- variable-speed mode
# ---------------------------------------------------------------------------

def test_variable_speed_holds_tsr_and_cp_is_near_constant():
    """
    At fixed TSR, Cp varies only through Reynolds drift, and power scales v^3.

    The 3-7 m/s span is set by the cache, not the physics: at fixed TSR the
    Reynolds number scales linearly with wind speed, so the original 5-20 m/s
    sweep left the S809 cache's 1.3M ceiling at 10 m/s and was 2.3x past it at
    20. Three of its five points were therefore clamped -- in the check whose
    entire purpose is measuring Reynolds drift.
    """

    geometry = phase_vi_geometry()
    curve = power_curve(geometry, [3.0, 4.0, 5.0, 6.0, 7.0], tsr=6.0, **AIR)

    assert curve["mode"] == "variable-speed"
    assert all(abs(p["tsr"] - 6.0) < 1e-12 for p in curve["points"])

    spread = (max(curve["Cp"]) - min(curve["Cp"])) / max(curve["Cp"])
    assert spread < 0.05, f"Cp varies {spread:.1%} at fixed TSR -- too much for Re drift"

    first, last = curve["points"][0], curve["points"][-1]
    ratio = last["power"] / first["power"]
    cubic = (last["v_inf"] / first["v_inf"]) ** 3
    assert ratio == pytest.approx(cubic, rel=0.05)


# ---------------------------------------------------------------------------
# Check 4 -- Cp-lambda curve shape
# ---------------------------------------------------------------------------

def test_cp_lambda_curve_is_unimodal_and_sub_betz():
    """
    A single interior maximum, rising then falling, everywhere below Betz.

    Capped at lambda = 7.5 by the same cache limit as check 3: lambda = 8.0
    puts the tip station at Re = 1,384,108 against a 1,300,000 ceiling. The
    peak (~lambda = 7) is still bracketed, which is what the unimodality
    assertions need.
    """

    geometry = phase_vi_geometry()
    tsr_values = [round(1.0 + 0.5 * i, 2) for i in range(14)]
    curve = cp_lambda_curve(geometry, tsr_values, v_inf=7.0, **AIR)

    assert curve["mode"] == "cp-lambda"
    cps = curve["Cp"]
    assert max(cps) < BETZ_LIMIT
    assert all(value > 0.0 for value in cps)

    peak = peak_cp(curve)
    assert 0 < peak["index"] < len(cps) - 1, "the swept range does not bracket the peak"
    assert all(cps[i] < cps[i + 1] for i in range(peak["index"]))
    assert all(cps[i] > cps[i + 1] for i in range(peak["index"], len(cps) - 1))


# ---------------------------------------------------------------------------
# Check 5 -- the sweep layer adds nothing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("v_inf,tsr", [(6.0, 4.0), (9.0, 7.5)])
def test_sweep_layer_introduces_no_drift(v_inf, tsr):
    """
    `power_curve` and `operating_point` must be bit-identical to `solve_rotor`.

    Exact equality again, and for the same reason as the stage-4 loop test:
    this layer calls straight through, so anything other than equality is the
    wrapper doing arithmetic it should not be doing.
    """

    geometry = demo_rotor_geometry()

    direct = solve_rotor(geometry, tsr=tsr, v_inf=v_inf, **AIR)
    swept = power_curve(geometry, [v_inf], tsr=tsr, **AIR)["points"][0]
    single = operating_point(geometry, v_inf=v_inf, tsr=tsr, **AIR)

    assert swept["Cp"] == direct["Cp"]
    assert swept["Ct"] == direct["Ct"]
    assert single["Cp"] == direct["Cp"]


# ---------------------------------------------------------------------------
# Check 6 -- guard rails
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("label,call", [
    ("neither rpm nor tsr", lambda g: power_curve(g, [7.0], **AIR)),
    ("both rpm and tsr", lambda g: power_curve(g, [7.0], rpm=60.0, tsr=6.0, **AIR)),
    ("empty wind_speeds", lambda g: power_curve(g, [], tsr=6.0, **AIR)),
    ("zero wind speed", lambda g: power_curve(g, [0.0], tsr=6.0, **AIR)),
    ("negative wind speed", lambda g: power_curve(g, [-3.0], tsr=6.0, **AIR)),
    ("non-positive tsr", lambda g: operating_point(g, v_inf=7.0, tsr=0.0, **AIR)),
    ("empty tsr_values", lambda g: cp_lambda_curve(g, [], **AIR)),
])
def test_guard_rails_raise_rather_than_guessing(label, call):
    """
    Ambiguous or invalid input is refused, not interpreted.

    `rpm=` xor `tsr=` is the important one: it is what keeps the two machine
    conventions from being silently conflated.
    """

    with pytest.raises(ValueError):
        call(demo_rotor_geometry())


def test_air_properties_are_still_required():
    """
    Task 1's rule, checked at this layer too: omitting density is a TypeError,
    not a number 25 % high that still lands in the plausibility band.
    """

    with pytest.raises(TypeError):
        power_curve(demo_rotor_geometry(), [7.0], tsr=6.0)
