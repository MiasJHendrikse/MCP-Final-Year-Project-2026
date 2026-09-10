"""
Validation for the power-curve / Cp-lambda sweep layer (bem/powercurve.py).

powercurve.py adds no aerodynamics of its own -- it sweeps rotor.solve_rotor
and converts the dimensionless result to watts/newtons/newton-metres. These
checks are therefore aimed at exactly the two things that layer can get
wrong: the dimensional conversion, and the operating-point convention
(fixed-speed vs variable-speed, wind speed vs TSR).

Six checks:

1. Dimensional round-trip: power/thrust/torque invert back to solve_rotor's
   own Cp/Ct, and torque * omega reproduces power, to machine precision.
2. Regression against the real NREL Phase VI (Sequence S) numbers recorded in
   the 2026-07-28 journal entry -- 71.63 RPM fixed-speed at the seven
   published Sequence S wind-speed test points. Also confirms TSR falls
   monotonically with wind speed at fixed rotor speed.
3. Variable-speed mode: at constant TSR, Cp is constant with wind speed
   apart from Reynolds drift. A large spread here means the fixed-TSR path
   is not actually holding TSR, or the polar cache has a Reynolds defect.
4. Cp-lambda curve shape: every Cp below the Betz limit, a single interior
   maximum (unimodal), rising strictly before the peak and falling strictly
   after -- the shape a Cp-lambda curve must have for any real rotor.
5. Consistency: a power_curve point equals solve_rotor called directly at the
   same (v_inf, tsr), exactly -- confirms the sweep adds no drift.
6. Guard rails: the mutually-exclusive rpm=/tsr= contract and the positivity
   requirements raise rather than silently guessing.

Run directly: `python -m validation.validate_powercurve` (from src/).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

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
    demo_rotor_geometry,
    phase_vi_geometry,
    solve_rotor,
)
from config import load_phase_vi_rotor

# Air properties and rotor speed come from config/rotor_phase_vi.yaml.
# AIR_DENSITY was a module constant here (and in three sibling scripts)
# restating the same sea-level value the solver also defaulted to; both are
# gone. AIR_KINEMATIC_VISCOSITY is new to this script because solve_rotor now
# requires it explicitly.
#
# demo_rotor_geometry() in checks 5-6 is a synthetic pipeline-exercising blade
# with no site of its own; it is run at the same sea-level condition, which is
# what its recorded behaviour was produced at.
PHASE_VI = load_phase_vi_rotor()
AIR_DENSITY = PHASE_VI.air_density
AIR_KINEMATIC_VISCOSITY = PHASE_VI.kinematic_viscosity
PHASE_VI_RATED_RPM = PHASE_VI.rated_rpm

#: NREL Phase VI (Sequence S) at 71.63 RPM, from the 2026-07-28 journal entry
#: (the run that confirmed station Reynolds numbers were no longer clamped to
#: the old 500k cache ceiling). v_inf -> (Cp, Ct). These are this project's
#: own solver's values, not experimental data -- a regression anchor, not a
#: validation against ground truth.
PHASE_VI_SEQUENCE_S_REFERENCE = {
    5.0: (0.4062, 0.6777),
    7.0: (0.3670, 0.5415),
    10.0: (0.2587, 0.3653),
    13.0: (0.1830, 0.2518),
    15.0: (0.1501, 0.2016),
    20.0: (0.1017, 0.1298),
    25.0: (0.0760, 0.0948),  # not swept -- see PHASE_VI_MAX_CACHED_WIND_SPEED_MS
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


def check_1_dimensional_round_trip():
    print("=== Check 1: dimensional round-trip (power/thrust/torque <-> Cp/Ct) ===")
    geometry = phase_vi_geometry()
    curve = power_curve(geometry, [5.0, 10.0, 20.0], rpm=PHASE_VI_RATED_RPM,
                        air_density=AIR_DENSITY,
                        kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)
    area = swept_area(geometry)

    for p in curve["points"]:
        v = p["v_inf"]
        cp_back = p["power"] / (0.5 * AIR_DENSITY * v ** 3 * area)
        ct_back = p["thrust"] / (0.5 * AIR_DENSITY * v ** 2 * area)
        power_back = p["torque"] * p["omega"]

        print(f"  V={v:5.1f}  Cp={p['Cp']:.6f} (round-trip {cp_back:.6f})  "
              f"Ct={p['Ct']:.6f} (round-trip {ct_back:.6f})  "
              f"P={p['power']:10.3f} W (Q*omega {power_back:10.3f} W)")

        assert abs(cp_back - p["Cp"]) < 1e-12, f"Cp round-trip failed at V={v}"
        assert abs(ct_back - p["Ct"]) < 1e-12, f"Ct round-trip failed at V={v}"
        assert abs(power_back - p["power"]) < 1e-9, f"Q*omega != P at V={v}"

    # RPM <-> rad/s helpers are each other's inverse.
    for rpm in (1.0, 71.63, 1500.0):
        assert abs(omega_to_rpm(rpm_to_omega(rpm)) - rpm) < 1e-12
    print("PASS: dimensional conversions are exact and self-consistent.\n")


def check_2_phase_vi_fixed_speed_regression():
    print("=== Check 2: NREL Phase VI fixed-speed regression (71.63 RPM) ===")
    geometry = phase_vi_geometry()
    speeds = sorted(v for v in PHASE_VI_SEQUENCE_S_REFERENCE
                    if v <= PHASE_VI_MAX_CACHED_WIND_SPEED_MS)
    curve = power_curve(geometry, speeds, rpm=PHASE_VI_RATED_RPM,
                        air_density=AIR_DENSITY,
                        kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)

    assert curve["mode"] == "fixed-speed", curve["mode"]

    for p in curve["points"]:
        cp_ref, ct_ref = PHASE_VI_SEQUENCE_S_REFERENCE[p["v_inf"]]
        print(f"  V={p['v_inf']:5.1f}  tsr={p['tsr']:6.3f}  "
              f"Cp={p['Cp']:.4f} (ref {cp_ref:.4f})  "
              f"Ct={p['Ct']:.4f} (ref {ct_ref:.4f})  P={p['power']:9.1f} W")
        assert abs(p["Cp"] - cp_ref) < 5e-5, f"Cp regression at V={p['v_inf']}"
        assert abs(p["Ct"] - ct_ref) < 5e-5, f"Ct regression at V={p['v_inf']}"
        # Rotor speed is what is being held fixed, not TSR.
        assert abs(p["rpm"] - PHASE_VI_RATED_RPM) < 1e-9

    tsrs = curve["tsr"]
    assert all(tsrs[i] > tsrs[i + 1] for i in range(len(tsrs) - 1)), \
        "TSR must fall monotonically with wind speed at fixed rotor speed"
    print(f"  TSR falls monotonically {tsrs[0]:.3f} -> {tsrs[-1]:.3f} across "
          f"{speeds[0]:.0f}-{speeds[-1]:.0f} m/s, as it must at fixed RPM")
    print("PASS: reproduces the recorded Phase VI Sequence S values exactly.\n")


def check_3_variable_speed_cp_near_constant():
    print("=== Check 3: variable-speed mode holds Cp constant but for Reynolds ===")
    geometry = phase_vi_geometry()
    curve = power_curve(geometry, [5.0, 7.0, 10.0, 15.0, 20.0], tsr=6.0,
                        air_density=AIR_DENSITY,
                        kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)

    assert curve["mode"] == "variable-speed", curve["mode"]

    for p in curve["points"]:
        print(f"  V={p['v_inf']:5.1f}  tsr={p['tsr']:.3f}  rpm={p['rpm']:7.2f}  "
              f"Cp={p['Cp']:.4f}  P={p['power']:9.1f} W")
        assert abs(p["tsr"] - 6.0) < 1e-12, "variable-speed mode must hold TSR fixed"

    spread = max(curve["Cp"]) - min(curve["Cp"])
    relative = spread / max(curve["Cp"])
    print(f"  Cp spread across the sweep: {spread:.4f} ({relative * 100:.2f}% "
          f"of peak) -- Reynolds drift only")
    assert relative < 0.05, (
        f"Cp varies {relative * 100:.1f}% at fixed TSR -- too much to be "
        f"Reynolds drift alone"
    )

    # At fixed TSR and near-constant Cp, power must scale ~v^3.
    p0, p1 = curve["points"][0], curve["points"][-1]
    ratio = p1["power"] / p0["power"]
    cubic = (p1["v_inf"] / p0["v_inf"]) ** 3
    print(f"  power ratio {p0['v_inf']:.0f}->{p1['v_inf']:.0f} m/s: "
          f"{ratio:.3f} vs v^3 prediction {cubic:.3f}")
    assert abs(ratio - cubic) / cubic < 0.05, "power does not scale ~v^3 at fixed TSR"
    print("PASS: TSR held exactly, Cp near-constant, power scales as v^3.\n")


def check_4_cp_lambda_curve_shape():
    print("=== Check 4: Cp-lambda curve shape ===")
    geometry = phase_vi_geometry()
    # 1.0 .. 7.5. The upper end is set by the S809 cache, not by the physics.
    # At v_inf = 7 m/s the tip station's Reynolds estimate is
    # hypot(7, omega*r)*chord/nu; tsr = 7.5 puts it at 1,299,666 against the
    # cache's 1,300,000 ceiling, and tsr = 8.0 at 1,384,108 -- outside it.
    # Before Task 4 those points silently clamped to the ceiling curve and the
    # sweep ran to tsr = 12 on a Reynolds number that was not the station's.
    # Now they raise, correctly, so the sweep stops where the cache's coverage
    # actually stops. The Cp peak (~tsr = 7) is still bracketed, which is what
    # the unimodality assertions below need.
    tsr_values = [round(1.0 + 0.5 * i, 2) for i in range(14)]  # 1.0 .. 7.5
    curve = cp_lambda_curve(geometry, tsr_values, v_inf=7.0,
                            air_density=AIR_DENSITY,
                            kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)

    assert curve["mode"] == "cp-lambda", curve["mode"]

    for p in curve["points"]:
        print(f"  tsr={p['tsr']:5.2f}  Cp={p['Cp']:.4f}  Ct={p['Ct']:.4f}")

    cps = curve["Cp"]
    worst = max(cps)
    assert worst < BETZ_LIMIT, f"Cp={worst:.4f} exceeds the Betz limit {BETZ_LIMIT:.4f}"
    assert all(c > 0.0 for c in cps), "Cp must be positive across the swept range"

    pk = peak_cp(curve)
    print(f"  peak Cp={pk['Cp']:.4f} at tsr={pk['tsr']:.2f} "
          f"(Betz limit {BETZ_LIMIT:.4f})")
    assert 0 < pk["index"] < len(cps) - 1, \
        "peak Cp sits at a swept endpoint -- the range does not bracket the peak"

    rising = all(cps[i] < cps[i + 1] for i in range(pk["index"]))
    falling = all(cps[i] > cps[i + 1] for i in range(pk["index"], len(cps) - 1))
    assert rising, "Cp is not strictly rising below the peak -- curve is not unimodal"
    assert falling, "Cp is not strictly falling above the peak -- curve is not unimodal"
    print("PASS: single interior maximum, strictly unimodal, everywhere sub-Betz.\n")


def check_5_matches_solve_rotor_directly():
    print("=== Check 5: sweep agrees with solve_rotor called directly ===")
    geometry = demo_rotor_geometry()

    for v_inf, tsr in ((6.0, 4.0), (9.0, 7.5)):
        direct = solve_rotor(geometry, tsr=tsr, v_inf=v_inf,
                             air_density=AIR_DENSITY,
                             kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)
        swept = power_curve(geometry, [v_inf], tsr=tsr,
                            air_density=AIR_DENSITY,
                            kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)["points"][0]
        single = operating_point(geometry, v_inf=v_inf, tsr=tsr,
                                 air_density=AIR_DENSITY,
                                 kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)

        print(f"  V={v_inf:4.1f} tsr={tsr:4.1f}  solve_rotor Cp={direct['Cp']:.10f}  "
              f"power_curve Cp={swept['Cp']:.10f}  operating_point Cp={single['Cp']:.10f}")
        assert swept["Cp"] == direct["Cp"], "power_curve drifted from solve_rotor"
        assert swept["Ct"] == direct["Ct"], "power_curve drifted from solve_rotor"
        assert single["Cp"] == direct["Cp"], "operating_point drifted from solve_rotor"
    print("PASS: no numerical drift introduced by the sweep layer.\n")


def check_6_guard_rails():
    print("=== Check 6: guard rails reject ambiguous or invalid input ===")
    geometry = demo_rotor_geometry()

    # Air properties are supplied on every case: they are required arguments
    # now, and a missing one raises TypeError, which would mask the ValueError
    # each guard rail is actually being tested for.
    air = {"air_density": AIR_DENSITY,
           "kinematic_viscosity": AIR_KINEMATIC_VISCOSITY}

    cases = [
        ("neither rpm nor tsr",
         lambda: power_curve(geometry, [7.0], **air)),
        ("both rpm and tsr",
         lambda: power_curve(geometry, [7.0], rpm=60.0, tsr=6.0, **air)),
        ("empty wind_speeds",
         lambda: power_curve(geometry, [], tsr=6.0, **air)),
        ("zero wind speed",
         lambda: power_curve(geometry, [0.0], tsr=6.0, **air)),
        ("negative wind speed",
         lambda: power_curve(geometry, [-3.0], tsr=6.0, **air)),
        ("non-positive tsr",
         lambda: operating_point(geometry, v_inf=7.0, tsr=0.0, **air)),
        ("empty tsr_values",
         lambda: cp_lambda_curve(geometry, [], **air)),
    ]

    for label, call in cases:
        try:
            call()
        except ValueError as exc:
            print(f"  {label:24s} -> correctly raised ValueError: {exc}")
        else:
            raise AssertionError(f"{label} should have raised ValueError")
    print("PASS: all guard rails raise rather than guessing.\n")


def main():
    print("=" * 68)
    print("Power curve / Cp-lambda validation (bem/powercurve.py)")
    print("=" * 68 + "\n")

    check_1_dimensional_round_trip()
    check_2_phase_vi_fixed_speed_regression()
    check_3_variable_speed_cp_near_constant()
    check_4_cp_lambda_curve_shape()
    check_5_matches_solve_rotor_directly()
    check_6_guard_rails()

    print("Power curve validation: ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
