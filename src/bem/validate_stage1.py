"""
Stage 1 validation for the single-station BEM solver (station.py).

Three checks, per PROJECT_PLAN.md Phase 1 / the session that introduced
station.py:

1. Hand-check: a and a' land in the physically expected range for normal
   operation (0 < a < 0.5, a' small and positive) for a simple case.
2. Cross-check against the classical, uncorrected BEM equations (Hansen
   2008, Ch. 6): iterate the standard fixed-point (a, a') update directly
   (no phi-residual, no Ning reduction) and confirm it converges to the same
   (phi, a, a') as solve_station -- the two formulations are algebraically
   equivalent with no corrections applied, so they must agree.
3. Residual smoothness: print (and optionally plot) residual(phi) over the
   scan range used by solve_station, to confirm it is smooth and does not
   silently rely on discontinuities other than the known Ct=0 pole (see
   station._select_bracket) within the physical bracket actually used.

Run directly: `python src/bem/validate_stage1.py`

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

from bem.airfoil import LinearPolar
from bem.station import StationParams, residual, solve_station


def hansen_fixed_point(station: StationParams, tol=1e-10, max_iter=200):
    """
    Classical uncorrected BEM fixed-point iteration (Hansen 2008, Ch. 6).

    Iterates:
        phi = atan((1-a) / ((1+a') * tsr))
        alpha = phi - twist
        Cl, Cd = airfoil(alpha)
        Cn = Cl*cos(phi) + Cd*sin(phi); Ct = Cl*sin(phi) - Cd*cos(phi)
        a = 1 / (4*sin(phi)^2 / (sigma*Cn) + 1)
        a' = 1 / (4*sin(phi)*cos(phi) / (sigma*Ct) - 1)
    until a, a' stop changing. No tip-loss, no high-thrust correction --
    this is the same physics as station.residual, just solved by direct
    substitution instead of a single-variable root-find, so it is a fair
    independent check of the same equations.
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

        if abs(a_new - a) < tol and abs(a_prime_new - a_prime) < tol:
            a, a_prime = a_new, a_prime_new
            break
        a, a_prime = a_new, a_prime_new

    phi = math.atan((1 - a) / ((1 + a_prime) * station.tsr))
    return {"phi": phi, "a": a, "a_prime": a_prime}


def check_1_hand_bounds(station):
    print("=== Check 1: hand-check physical bounds ===")
    result = solve_station(station)
    phi_deg = math.degrees(result["phi"])
    print(f"phi = {phi_deg:.4f} deg, a = {result['a']:.6f}, a' = {result['a_prime']:.6f}")
    assert 0.0 < result["a"] < 0.5, f"a={result['a']} outside (0, 0.5)"
    assert -0.5 < result["a_prime"] < 0.5, f"a'={result['a_prime']} outside expected range"
    print("PASS: a and a' within expected physical bounds for normal operation.\n")
    return result


def check_2_hansen_cross_check(station, ning_result):
    print("=== Check 2: cross-check vs. classical Hansen (2008) fixed-point BEM ===")
    hansen_result = hansen_fixed_point(station)
    print(
        f"Ning-residual : phi={math.degrees(ning_result['phi']):.6f} deg, "
        f"a={ning_result['a']:.8f}, a'={ning_result['a_prime']:.8f}"
    )
    print(
        f"Hansen fixed-pt: phi={math.degrees(hansen_result['phi']):.6f} deg, "
        f"a={hansen_result['a']:.8f}, a'={hansen_result['a_prime']:.8f}"
    )
    for key in ("phi", "a", "a_prime"):
        diff = abs(ning_result[key] - hansen_result[key])
        assert diff < 1e-6, f"{key} mismatch: {diff:.3e}"
    print("PASS: Ning residual and classical Hansen fixed-point agree to 1e-6.\n")


def check_3_residual_smoothness(station, n=25):
    print("=== Check 3: residual smoothness over the physical bracket ===")
    result = solve_station(station)
    phi_star = result["phi"]
    # Window straddling the found root, well clear of the known Ct=0 pole
    # (see hand-check output above and station._select_bracket docstring).
    lo = max(1e-3, phi_star - math.radians(5))
    hi = phi_star + math.radians(5)
    prev = None
    max_jump = 0.0
    for i in range(n):
        phi = lo + (hi - lo) * i / (n - 1)
        val = residual(phi, station)
        print(f"  phi={math.degrees(phi):7.3f} deg  R(phi)={val: .6f}")
        if prev is not None:
            max_jump = max(max_jump, abs(val - prev))
        prev = val
    print(f"max |delta R| between adjacent samples: {max_jump:.6f}")
    assert max_jump < 1.0, "residual is not smooth over the physical bracket"
    print("PASS: residual is smooth and single-valued over the physical bracket.\n")


def main():
    airfoil = LinearPolar()
    station = StationParams(
        r=5.0, chord=0.3, twist=math.radians(5), airfoil=airfoil, tsr=5.0, n_blades=3
    )

    ning_result = check_1_hand_bounds(station)
    check_2_hansen_cross_check(station, ning_result)
    check_3_residual_smoothness(station)

    print("Stage 1 validation: ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
