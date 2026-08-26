"""
Stage 4 validation for the multi-station spanwise loop over real S809
polar data (rotor.py, airfoil.S809Polar).

Four checks:

1. alpha, Cl, Cd, a, a_prime distributions across span are smooth (no
   station-to-station discontinuities/spikes) at a reasonable tsr.
2. Induction factors stay within physical bounds at every station, at both
   a moderate tsr (turbulent-wake correction inactive) and a higher tsr
   (correction actively engaged near the tip, a > 0.4).
3. Integrated rotor Ct, Cp are physically plausible (Cp below the Betz
   limit ~0.593, Ct in a sane range) at the reasonable operating point.
4. Regression: a single-station RotorGeometry, run through solve_rotor,
   reproduces station.solve_station called directly with the same
   StationParams/airfoil -- confirms the loop's station construction adds
   no discrepancy of its own.

Run directly: `python -m validation.validate_stage4` (from src/).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

from bem.airfoil import S809Polar
from bem.rotor import RotorGeometry, demo_rotor_geometry, solve_rotor
from bem.station import StationParams, solve_station
from config import load_phase_vi_rotor

# `bem.rotor.AIR_KINEMATIC_VISCOSITY` is gone: air properties are now required
# arguments on solve_rotor, supplied from config/. The stage-4 checks below
# exercise demo_rotor_geometry(), a synthetic pipeline blade with no site of
# its own, against the S809 cache -- so they run at the same sea-level
# condition their recorded behaviour was produced at, taken from
# config/rotor_phase_vi.yaml rather than restated here.
_AIR = load_phase_vi_rotor()
AIR_DENSITY = _AIR.air_density
AIR_KINEMATIC_VISCOSITY = _AIR.kinematic_viscosity


def check_1_spanwise_smoothness():
    print("=== Check 1: spanwise smoothness at a reasonable tsr ===")
    geometry = demo_rotor_geometry()
    result = solve_rotor(geometry, tsr=5.0, air_density=AIR_DENSITY,
                         kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)
    stations = result["stations"]

    for s in stations:
        print(f"  r={s['r']:5.2f} alpha={math.degrees(s['alpha']):7.3f} deg  "
              f"Cl={s['Cl']:7.4f}  Cd={s['Cd']:7.4f}  a={s['a']:7.4f}  "
              f"a'={s['a_prime']:8.5f}  F={s['F']:.4f}")

    # No large jump between adjacent stations for any tracked quantity --
    # thresholds are generous (this is a discontinuity/spike check, not a
    # tight tolerance -- taper/twist/Reynolds all vary smoothly by
    # construction, so a real solver bug would show up as an obvious spike,
    # not a borderline violation).
    thresholds = {"alpha": math.radians(15), "Cl": 0.5, "Cd": 0.1, "a": 0.2, "a_prime": 0.1}
    for key, limit in thresholds.items():
        steps = [abs(stations[i][key] - stations[i - 1][key]) for i in range(1, len(stations))]
        max_step = max(steps)
        print(f"  max adjacent-station |delta {key}| = {max_step:.4f} (limit {limit})")
        assert max_step < limit, f"spanwise discontinuity in {key}: {max_step}"
    print("PASS: alpha, Cl, Cd, a, a' vary smoothly across span, no spikes.\n")


def check_2_physical_bounds():
    print("=== Check 2: induction factors physically bounded (moderate and high tsr) ===")
    geometry = demo_rotor_geometry()

    for tsr, expect_correction in ((5.0, False), (9.0, True)):
        result = solve_rotor(geometry, tsr=tsr, air_density=AIR_DENSITY,
                             kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)
        stations = result["stations"]
        a_values = [s["a"] for s in stations]
        ap_values = [s["a_prime"] for s in stations]

        for s in stations:
            assert 0.0 < s["a"] < 1.0, f"a={s['a']} out of (0,1) at r={s['r']} (tsr={tsr})"
            assert -0.5 < s["a_prime"] < 0.5, f"a'={s['a_prime']} out of range at r={s['r']}"
            assert math.isfinite(s["a"]) and math.isfinite(s["a_prime"])

        max_a = max(a_values)
        engaged = max_a > 0.4
        print(f"  tsr={tsr}: a in [{min(a_values):.4f}, {max_a:.4f}], "
              f"a' in [{min(ap_values):.5f}, {max(ap_values):.5f}] "
              f"(turbulent-wake correction {'engaged' if engaged else 'not engaged'})")
        assert engaged == expect_correction, (
            f"expected turbulent-wake correction {'engaged' if expect_correction else 'inactive'} "
            f"at tsr={tsr}, got max a={max_a}"
        )
    print("PASS: a, a' stay physically bounded at every station, including with the "
          "turbulent-wake correction actively engaged near the tip at high tsr.\n")


def check_3_integrated_rotor_plausible():
    print("=== Check 3: integrated rotor Ct, Cp physically plausible ===")
    geometry = demo_rotor_geometry()
    result = solve_rotor(geometry, tsr=5.0, air_density=AIR_DENSITY,
                         kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)
    ct, cp = result["Ct"], result["Cp"]
    print(f"  Ct={ct:.4f}, Cp={cp:.4f} (Betz limit Cp <= 0.593)")
    assert 0.0 < cp < 0.593, f"Cp={cp} not physically plausible"
    assert 0.0 < ct < 1.5, f"Ct={ct} not in a sane range"
    print("PASS: Cp is below the Betz limit and Ct is in a sane range.\n")


def check_4_regression_single_station():
    print("=== Check 4: single-station regression (loop wiring vs. direct solve_station) ===")
    r, chord, twist_deg, R, tsr_tip, v_inf = 3.0, 0.3, 6.0, 8.0, 5.0, 7.0
    geometry = RotorGeometry(r=[r], chord=[chord], twist=[math.radians(twist_deg)], R=R, n_blades=3)

    result = solve_rotor(geometry, tsr=tsr_tip, v_inf=v_inf,
                         air_density=AIR_DENSITY,
                         kinematic_viscosity=AIR_KINEMATIC_VISCOSITY)
    via_loop = result["stations"][0]

    omega = tsr_tip * v_inf / R
    w_approx = math.hypot(v_inf, omega * r)
    reynolds = w_approx * chord / AIR_KINEMATIC_VISCOSITY
    tsr_local = omega * r / v_inf
    station = StationParams(
        r=r, chord=chord, twist=math.radians(twist_deg), airfoil=S809Polar(reynolds),
        tsr=tsr_local, R=R, n_blades=3,
    )
    direct = solve_station(station)

    print(f"  via solve_rotor : phi={math.degrees(via_loop['phi']):.6f} deg, "
          f"a={via_loop['a']:.8f}, a'={via_loop['a_prime']:.8f}")
    print(f"  direct solve_station: phi={math.degrees(direct['phi']):.6f} deg, "
          f"a={direct['a']:.8f}, a'={direct['a_prime']:.8f}")

    for key in ("phi", "a", "a_prime", "Cl", "Cd"):
        assert via_loop[key] == direct[key], f"{key} mismatch between loop and direct solve"
    print("PASS: solve_rotor reproduces solve_station exactly for an equivalent single station.\n")


def main():
    check_1_spanwise_smoothness()
    check_2_physical_bounds()
    check_3_integrated_rotor_plausible()
    check_4_regression_single_station()
    print("Stage 4 validation: ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
