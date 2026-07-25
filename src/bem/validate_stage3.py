"""
Stage 3 validation for the Glauert/Buhl high-thrust correction
(corrections.py: high_thrust_correction, corrected_axial_induction;
station.py: residual/_blade_element_and_induction).

Four checks:

1. Ct(a) is monotonically increasing over a in [0, ~1], with no
   discontinuity or kink at the blend point ac=0.4 (value and numerical
   derivative checked either side).
2. For a <= ac, corrected_axial_induction reproduces the plain Stage 2
   closed form a = Y / (4F + Y) bit-for-bit (this correction is a pure
   addition at low induction, not a change to it).
3. Regression: re-run the exact Stage 1/2 test stations and confirm
   identical results (those stations stay below ac).
4. A high-thrust test case (heavily loaded station) pushes a station's
   naive induction above ac; confirm the corrected solve still converges
   to a bounded, sensible a rather than failing, diverging, or exceeding 1.

Run directly: `python -m bem.validate_stage3` (from src/).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

from bem.airfoil import LinearPolar
from bem.corrections import BUHL_AC, corrected_axial_induction, high_thrust_correction
from bem.station import StationParams, residual, solve_station, _blade_element_and_induction

R = 5.5
CHORD = 0.3
TWIST = math.radians(5)
TSR = 5.0


def _station(r, chord=CHORD, twist=TWIST, tsr=TSR, n_blades=3, R_=R, r_hub=None):
    return StationParams(
        r=r, chord=chord, twist=twist, airfoil=LinearPolar(), tsr=tsr,
        R=R_, n_blades=n_blades, r_hub=r_hub,
    )


def check_1_ct_curve_monotonic_and_continuous():
    print("=== Check 1: Ct(a) monotonic, continuous, no kink at ac ===")
    F = 0.85
    n = 400
    a_values = [1e-4 + (0.999 - 1e-4) * i / (n - 1) for i in range(n)]

    ct_values = []
    for a in a_values:
        ct_naive = 4.0 * a * F * (1.0 - a)
        ct_values.append(high_thrust_correction(a, ct_naive, F))

    max_backstep = 0.0
    for i in range(1, n):
        max_backstep = max(max_backstep, ct_values[i - 1] - ct_values[i])
    print(f"F={F}: Ct sweep a=[{a_values[0]:.4f}, {a_values[-1]:.4f}], "
          f"max backward step in Ct = {max_backstep:.3e}")
    assert max_backstep < 1e-9, "Ct(a) is not monotonically increasing"

    eps = 1e-6
    ct_naive_lo = 4.0 * (BUHL_AC - eps) * F * (1.0 - (BUHL_AC - eps))
    ct_naive_hi = 4.0 * (BUHL_AC + eps) * F * (1.0 - (BUHL_AC + eps))
    ct_lo = high_thrust_correction(BUHL_AC - eps, ct_naive_lo, F)
    ct_hi = high_thrust_correction(BUHL_AC + eps, ct_naive_hi, F)
    value_jump = abs(ct_hi - ct_lo)

    slope_lo = (high_thrust_correction(BUHL_AC, 4 * BUHL_AC * F * (1 - BUHL_AC), F)
                - high_thrust_correction(BUHL_AC - 2 * eps, 4 * (BUHL_AC - 2 * eps) * F * (1 - (BUHL_AC - 2 * eps)), F)) / (2 * eps)
    slope_hi = (high_thrust_correction(BUHL_AC + 2 * eps, 4 * (BUHL_AC + 2 * eps) * F * (1 - (BUHL_AC + 2 * eps)), F)
                - high_thrust_correction(BUHL_AC, 4 * BUHL_AC * F * (1 - BUHL_AC), F)) / (2 * eps)
    slope_jump = abs(slope_hi - slope_lo)

    print(f"at ac={BUHL_AC}: |Ct jump| = {value_jump:.3e}, |slope jump| = {slope_jump:.3e}")
    # Thresholds scaled for the eps=1e-6 finite-difference step above (float64
    # roundoff, not a real discontinuity -- confirmed algebraically in
    # corrections.py's docstring derivation: the value/slope match exactly).
    assert value_jump < 1e-4, "Ct(a) has a discontinuity at the blend point"
    assert slope_jump < 1e-2, "Ct(a) has a visible kink (slope discontinuity) at the blend point"
    print("PASS: Ct(a) is monotonic and C0/C1-continuous at the blend point.\n")


def check_2_low_induction_unchanged():
    print("=== Check 2: a <= ac reproduces the plain Stage 2 closed form exactly ===")
    for F in (0.3, 0.6, 1.0):
        for Y in (0.01, 0.5, 1.0, 1.5):
            a_plain = Y / (4.0 * F + Y)
            a_corrected = corrected_axial_induction(Y, F)
            if a_plain <= BUHL_AC:
                assert a_corrected == a_plain, (
                    f"F={F} Y={Y}: corrected a={a_corrected} != plain a={a_plain}"
                )
    print("PASS: for a <= ac, corrected_axial_induction is bit-identical to the "
          "Stage 2 closed form (pure addition, no change to low-induction behaviour).\n")


def check_3_regression_stage1_stage2():
    print("=== Check 3: regression against Stage 1/2 test stations ===")

    station_stage1 = _station(r=5.0, R_=1000.0)  # F~1, a well below ac
    result = solve_station(station_stage1)
    stage1_reference = {"phi": math.radians(10.013235), "a": 0.11399664, "a_prime": 0.00359930}
    print(f"Stage 1 station: phi={math.degrees(result['phi']):.6f} deg, "
          f"a={result['a']:.8f}, a'={result['a_prime']:.8f}")
    for key in ("phi", "a", "a_prime"):
        diff = abs(result[key] - stage1_reference[key])
        assert diff < 1e-5, f"{key} regressed vs. Stage 1 reference, diff={diff:.3e}"

    station_stage2_near_tip = _station(r=5.4)  # r/R=0.982, Stage 2's near-tip case
    result_tip = solve_station(station_stage2_near_tip)
    print(f"Stage 2 near-tip station: F={result_tip['F']:.4f}, a={result_tip['a']:.6f} "
          f"(Stage 2 reference: F=0.3755, a=0.236561)")
    assert abs(result_tip["F"] - 0.3755) < 1e-3
    assert abs(result_tip["a"] - 0.236561) < 1e-4
    print("PASS: Stage 1 and Stage 2 reference results reproduced (both stay below ac).\n")


def check_4_turbulent_wake_case():
    print("=== Check 4: high-thrust station pushed into the turbulent wake state ===")
    # Large chord (high solidity) drives up Cn/thrust loading enough to
    # push the naive momentum-theory a above ac at a normal tsr.
    station_high_thrust = _station(r=2.0, chord=1.4, twist=math.radians(2), tsr=6.0, R_=8.0, n_blades=3)
    result = solve_station(station_high_thrust)
    print(f"phi={math.degrees(result['phi']):.4f} deg, a={result['a']:.6f}, "
          f"a'={result['a_prime']:.6f}, F={result['F']:.4f}")

    assert result["a"] > BUHL_AC, "test case did not actually reach the turbulent wake regime"
    assert 0.0 < result["a"] < 1.0, f"a={result['a']} not bounded in (0, 1)"
    assert math.isfinite(result["a"]) and math.isfinite(result["phi"])
    print(f"PASS: station converged with a={result['a']:.4f} > ac={BUHL_AC}, "
          f"bounded and finite (not the unphysical/failed value plain momentum theory would give).")

    # Stronger check: confirm the *residual* (not just Ct(a) in isolation)
    # is smooth in the neighbourhood of the actual solved root, where a(phi)
    # is above ac -- this is the property that actually matters for the
    # root-finder and a future adjoint (a Buhl kink far from any root the
    # solver would ever land on would be harmless).
    phi_star = result["phi"]
    lo, hi, n = phi_star - math.radians(0.8), phi_star + math.radians(0.8), 200
    phis = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    a_of_phi = [_blade_element_and_induction(p, station_high_thrust)[0] for p in phis]
    res_values = [residual(p, station_high_thrust) for p in phis]
    max_step = max(abs(res_values[i] - res_values[i - 1]) for i in range(1, n))
    assert min(a_of_phi) > BUHL_AC, "window did not stay in the corrected (a > ac) regime"
    print(f"phi in [{math.degrees(lo):.2f}, {math.degrees(hi):.2f}] deg (a > ac throughout): "
          f"max adjacent residual step = {max_step:.3e}")
    assert max_step < 0.05, "residual shows a visible kink near the solved root"
    print("PASS: residual itself stays smooth near the solved (a > ac) root, not just Ct(a) "
          "in isolation.")

    # Separately: flag (do not silently accept) a *pre-existing*, Buhl-
    # unrelated pole found while building this test case. The plain
    # momentum-theory 'a' formula (a = Y / (4F + Y), Y = sigma*Cn/sin^2(phi))
    # has its own singularity whenever Y = -4F, i.e. Cn = -4F*sin^2(phi)/sigma
    # -- mirroring the Ct=0 pole in a' documented back in Stage 1/2
    # (station._select_bracket), just never triggered by those stations'
    # geometry. This test case's high solidity happens to sweep through a
    # small negative Cn near phi~1.9 deg, well below the solved root at
    # phi~3.0 deg, and hits it. solve_station's existing bracket-scan
    # strategy (pick the sign change nearest the zero-induction guess) is
    # unaffected -- it already avoids this region -- but it is a real,
    # separate discontinuity in the raw residual, not a Buhl blend defect.
    # Confirmed here rather than left for a future session to rediscover.
    phi_pole_lo, phi_pole_hi = math.radians(1.7), math.radians(2.1)
    n_pole = 50
    phis_pole = [phi_pole_lo + (phi_pole_hi - phi_pole_lo) * i / (n_pole - 1) for i in range(n_pole)]
    res_pole = [residual(p, station_high_thrust) for p in phis_pole]
    max_pole_step = max(abs(res_pole[i] - res_pole[i - 1]) for i in range(1, n_pole))
    print(f"[flagged, not a Buhl defect] pre-existing 'a' pole near phi~1.9 deg "
          f"(Cn crosses -4F*sin^2(phi)/sigma): max adjacent residual step = {max_pole_step:.3e} "
          f"(vs {max_step:.3e} near the actual solved root) -- solve_station's bracket scan "
          f"avoids this region, but it is a genuine discontinuity, distinct from Stage 3's "
          f"Buhl blend, worth tracking before this solver is differentiated.\n")


def main():
    check_1_ct_curve_monotonic_and_continuous()
    check_2_low_induction_unchanged()
    check_3_regression_stage1_stage2()
    check_4_turbulent_wake_case()
    print("Stage 3 validation: ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
