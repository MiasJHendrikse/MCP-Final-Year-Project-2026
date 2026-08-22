"""
Stage 2 validation for the Prandtl tip/hub-loss correction (corrections.py,
station.py).

Four checks:

1. F -> 1 as B grows large, and F -> 1 at mid-span stations far from tip/hub;
   in both limits solve_station's (phi, a, a') converge to Stage 1's
   uncorrected values.
2. F is physically bounded: decreasing toward the tip and hub, in (0, 1]
   across the expected r/R range, never negative or > 1.
3. Re-solve the Stage 1 test case with realistic tip-loss (near-tip station)
   and confirm a/a' differ from the uncorrected result there, while a
   mid-span station stays close to the uncorrected result.
4. Residual smoothness with F included, printed over the physical bracket
   for a near-tip station.

Run directly: `python -m validation.validate_stage2` (from src/).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

from bem.airfoil import LinearPolar
from bem.corrections import tip_loss_factor, hub_loss_factor
from bem.station import StationParams, residual, solve_station

R = 5.5  # rotor radius, m -- close enough to the Stage 1 r=5.0 to make r/R meaningful
CHORD = 0.3
TWIST = math.radians(5)
TSR = 5.0


def _station(r, n_blades=3, R_=R, r_hub=None):
    return StationParams(
        r=r, chord=CHORD, twist=TWIST, airfoil=LinearPolar(), tsr=TSR,
        R=R_, n_blades=n_blades, r_hub=r_hub,
    )


def check_1_uncorrected_limit():
    print("=== Check 1: F -> 1 (many blades / mid-span) recovers Stage 1 ===")

    # (a) many blades at a moderate r/R -- F should be ~1 even fairly close
    # to the tip, since more blades means less tip-vortex loss per blade.
    station_many_blades = _station(r=5.0, n_blades=200)
    result_many_blades = solve_station(station_many_blades)
    print(f"B=200, r/R={5.0/R:.3f}: F={result_many_blades['F']:.8f}")
    assert abs(result_many_blades["F"] - 1.0) < 1e-3, "F should be ~1 for very large B"

    # (b) mid-span station, normal blade count, large R so r/R is small ->
    # F should be 1 to double precision, matching Stage 1's r=5, R=1000 case.
    station_midspan = _station(r=5.0, n_blades=3, R_=1000.0)
    result_midspan = solve_station(station_midspan)
    stage1_reference = {"phi": math.radians(10.013235), "a": 0.11399664, "a_prime": 0.00359930}
    print(
        f"B=3, r/R={5.0/1000.0:.4f}: F={result_midspan['F']:.10f}, "
        f"phi={math.degrees(result_midspan['phi']):.6f} deg, "
        f"a={result_midspan['a']:.8f}, a'={result_midspan['a_prime']:.8f}"
    )
    assert abs(result_midspan["F"] - 1.0) < 1e-9
    for key in ("phi", "a", "a_prime"):
        diff = abs(result_midspan[key] - stage1_reference[key])
        assert diff < 1e-5, f"{key} should match Stage 1 uncorrected result, diff={diff:.3e}"
    print("PASS: F -> 1 in both limits, and results converge to Stage 1's uncorrected values.\n")


def check_2_physical_bounds():
    print("=== Check 2: F physically bounded, decreasing toward tip/hub ===")
    phi_typical = math.radians(10.0)
    r_hub = 0.5

    r_values = [0.5 + 0.05 * i for i in range(1, 100)]  # r_hub < r < R, excluding endpoints
    r_values = [r for r in r_values if r < R]

    f_tip_values = [tip_loss_factor(r, R, 3, phi_typical) for r in r_values]
    f_hub_values = [hub_loss_factor(r, r_hub, 3, phi_typical) for r in r_values]

    for r, f_tip, f_hub in zip(r_values, f_tip_values, f_hub_values):
        assert 0.0 < f_tip <= 1.0, f"F_tip={f_tip} out of (0,1] at r={r}"
        assert 0.0 < f_hub <= 1.0, f"F_hub={f_hub} out of (0,1] at r={r}"

    # Monotonicity: F_tip strictly decreases as r increases toward R;
    # F_hub strictly increases as r increases away from r_hub.
    for i in range(1, len(f_tip_values)):
        assert f_tip_values[i] <= f_tip_values[i - 1] + 1e-12, "F_tip not monotonically decreasing toward tip"
        assert f_hub_values[i] >= f_hub_values[i - 1] - 1e-12, "F_hub not monotonically increasing away from hub"

    print(f"r/R sweep ({r_values[0]/R:.3f} to {r_values[-1]/R:.3f}): "
          f"F_tip {f_tip_values[0]:.4f} -> {f_tip_values[-1]:.4f} (decreasing toward tip)")
    print(f"                                                    "
          f"F_hub {f_hub_values[0]:.4f} -> {f_hub_values[-1]:.4f} (increasing away from hub)")
    print("PASS: F_tip and F_hub stay in (0, 1], monotonic in the expected directions.\n")


def check_3_tip_loss_direction():
    print("=== Check 3: tip-loss shifts induction near the tip, not mid-span ===")

    station_near_tip = _station(r=5.4, n_blades=3)  # r/R = 0.982, close to tip
    station_midspan = _station(r=2.75, n_blades=3)  # r/R = 0.5

    result_near_tip = solve_station(station_near_tip)
    result_midspan = solve_station(station_midspan)

    # Uncorrected (F=1) reference at each same r, via a huge R.
    ref_near_tip = solve_station(_station(r=5.4, n_blades=3, R_=5000.0))
    ref_midspan = solve_station(_station(r=2.75, n_blades=3, R_=5000.0))

    da_tip = abs(result_near_tip["a"] - ref_near_tip["a"])
    da_midspan = abs(result_midspan["a"] - ref_midspan["a"])

    print(f"near-tip (r/R=0.982): F={result_near_tip['F']:.4f}, "
          f"a={result_near_tip['a']:.6f} vs uncorrected a={ref_near_tip['a']:.6f} "
          f"(|delta a|={da_tip:.6f})")
    print(f"mid-span (r/R=0.500): F={result_midspan['F']:.4f}, "
          f"a={result_midspan['a']:.6f} vs uncorrected a={ref_midspan['a']:.6f} "
          f"(|delta a|={da_midspan:.6f})")

    assert result_near_tip["F"] < 0.99, "expected a meaningfully reduced F near the tip"
    assert result_midspan["F"] > 0.999, "expected F ~ 1 at mid-span"
    # Tip-loss (F<1) increases the momentum-consistent induction relative to
    # the uncorrected case -- see corrections.py / station.py r->R discussion.
    assert result_near_tip["a"] > ref_near_tip["a"], "tip-loss should increase a near the tip"
    assert da_tip > da_midspan, "near-tip correction should exceed mid-span correction"
    assert da_midspan < 1e-3, "mid-span station should stay close to the uncorrected result"
    print("PASS: induction shifts near the tip as expected and converges to Stage 1 mid-span.\n")


def check_4_residual_smoothness():
    print("=== Check 4: residual smoothness with F included (near-tip station) ===")
    station_near_tip = _station(r=5.4, n_blades=3)
    result = solve_station(station_near_tip)
    phi_star = result["phi"]

    lo = max(1e-3, phi_star - math.radians(5))
    hi = phi_star + math.radians(5)
    n = 25
    prev = None
    max_jump = 0.0
    for i in range(n):
        phi = lo + (hi - lo) * i / (n - 1)
        val = residual(phi, station_near_tip)
        print(f"  phi={math.degrees(phi):7.3f} deg  R(phi)={val: .6f}")
        if prev is not None:
            max_jump = max(max_jump, abs(val - prev))
        prev = val
    print(f"max |delta R| between adjacent samples: {max_jump:.6f}")
    assert max_jump < 1.0, "residual is not smooth over the physical bracket with F included"
    print("PASS: residual stays smooth and single-valued with F included.\n")


def main():
    check_1_uncorrected_limit()
    check_2_physical_bounds()
    check_3_tip_loss_direction()
    check_4_residual_smoothness()
    print("Stage 2 validation: ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
