"""
Stage 1: single blade-element station BEM solver, no corrections applied.

Implements Ning (2014)'s reduction of the coupled axial/tangential induction
equations to a single residual equation in the inflow angle phi, solved with
a bracketed root-finder (Brent's method) rather than the classical
fixed-point iteration on (a, a') -- this is what keeps convergence robust
near high induction and keeps the solve a single scalar equation, which
matters later when this needs a discrete adjoint (one clean residual, no
iterative state, no branching on convergence history).

No tip-loss, no Glauert/Buhl high-thrust correction, no multi-station loop.
Those are explicitly deferred to later sessions (see PROJECT_PLAN.md Phase 1).

Reference
---------
Ning, S.A. (2014). "A simple solution method for the blade element momentum
equations with guaranteed convergence." Wind Energy, 17(9), 1327-1345.
Hansen, M.O.L. (2008). "Aerodynamics of Wind Turbines" (2nd ed.) -- used here
as the classical, uncorrected BEM formulation to cross-check against, since
with no corrections applied the two formulations must agree analytically.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
from dataclasses import dataclass

from scipy.optimize import brentq


@dataclass
class StationParams:
    """
    Inputs for one fixed blade-element station.

    Parameters
    ----------
    r : float
        Local radius, m.
    chord : float
        Local chord length, m.
    twist : float
        Local twist (+ pitch), radians, measured such that
        alpha = phi - twist.
    airfoil : object
        Callable/duck-typed polar with .cl(alpha) and .cd(alpha) in radians
        (e.g. bem.airfoil.LinearPolar). Stage 1 uses a synthetic polar only.
    tsr : float
        Local speed ratio at this station, Omega * r / Vinf.
    n_blades : int
        Number of blades, B.
    """

    r: float
    chord: float
    twist: float
    airfoil: object
    tsr: float
    n_blades: int = 3

    @property
    def solidity(self):
        """Local solidity sigma' = B*c / (2*pi*r)."""
        return self.n_blades * self.chord / (2.0 * math.pi * self.r)


def residual(phi, station: StationParams):
    """
    Ning-style BEM residual in the inflow angle phi, no corrections (F=1).

    Derivation: the blade-element normal/tangential force coefficients give
    momentum-consistent induction factors a(phi), a'(phi) directly (no
    fixed-point iteration needed); the residual is the mismatch between the
    input phi and the phi implied by the kinematic relation
    tan(phi) = (1-a) / ((1+a') * lambda_r), rearranged to avoid a division
    (so it stays finite and smooth even as phi -> 0 or a, a' -> singular
    values -- the property this formulation is chosen for).

    Parameters
    ----------
    phi : float
        Inflow angle, radians. Must lie strictly in (0, pi/2) for this
        uncorrected, non-turbulent-wake-state formulation.
    station : StationParams

    Returns
    -------
    float
        Residual value; solve_station finds the phi where this is zero.
    """

    alpha = phi - station.twist
    cl = station.airfoil.cl(alpha)
    cd = station.airfoil.cd(alpha)

    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    cn = cl * cos_phi + cd * sin_phi
    ct = cl * sin_phi - cd * cos_phi

    sigma = station.solidity

    # Momentum-consistent induction factors implied by this phi (F=1, no
    # high-thrust correction) -- see module docstring / Hansen (2008) Ch. 6.
    a = 1.0 / ((4.0 * sin_phi * sin_phi) / (sigma * cn) + 1.0)
    a_prime = 1.0 / ((4.0 * sin_phi * cos_phi) / (sigma * ct) - 1.0)

    return sin_phi * (1.0 + a_prime) * station.tsr - cos_phi * (1.0 - a)


def _induction_factors(phi, station: StationParams):
    """Recompute (a, a', Cl, Cd, Cn, Ct) at the solved phi (no residual)."""

    alpha = phi - station.twist
    cl = station.airfoil.cl(alpha)
    cd = station.airfoil.cd(alpha)

    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    cn = cl * cos_phi + cd * sin_phi
    ct = cl * sin_phi - cd * cos_phi

    sigma = station.solidity
    a = 1.0 / ((4.0 * sin_phi * sin_phi) / (sigma * cn) + 1.0)
    a_prime = 1.0 / ((4.0 * sin_phi * cos_phi) / (sigma * ct) - 1.0)

    return a, a_prime, cl, cd, cn, ct


def _select_bracket(station: StationParams, phi_range, n_scan):
    """
    Scan phi_range for the sign change nearest the zero-induction guess.

    The momentum-consistent a'(phi) used in `residual` has a genuine pole
    wherever Ct(phi) = 0 (a known BEM pathology, not specific to this
    implementation -- see Ning 2014 Sec. 3 for the full case analysis this
    motivates). A pole crosses the residual's sign just like a true root
    does, so naively bracketing on the endpoints of a wide range can hand
    brentq a pole instead of the physical root. Scanning first and picking
    the sign change closest to the classical zero-induction inflow angle
    (phi0 = atan(1/tsr)) avoids that failure mode without needing the full
    region-classification machinery Ning uses for the corrected solver --
    that is out of scope for this uncorrected Stage 1 solver.
    """

    phi_lo, phi_hi = phi_range
    phis = [phi_lo + (phi_hi - phi_lo) * i / (n_scan - 1) for i in range(n_scan)]
    vals = [residual(p, station) for p in phis]

    phi_guess = math.atan(1.0 / station.tsr)
    best = None
    best_dist = math.inf
    for i in range(n_scan - 1):
        if vals[i] == 0.0:
            return phis[i], phis[i]
        if vals[i] * vals[i + 1] < 0:
            midpoint = 0.5 * (phis[i] + phis[i + 1])
            dist = abs(midpoint - phi_guess)
            if dist < best_dist:
                best_dist = dist
                best = (phis[i], phis[i + 1])

    if best is None:
        raise ValueError(
            f"no sign change found for the residual over phi in {phi_range}; "
            "check station inputs."
        )
    return best


def solve_station(station: StationParams, bracket=(1e-4, math.pi / 2 - 1e-4), n_scan=2000):
    """
    Solve for the inflow angle phi at one station via bracketed root-find.

    Parameters
    ----------
    station : StationParams
    bracket : tuple of float
        (phi_lo, phi_hi), radians, the outer range to scan for candidate
        roots. Default spans the full normal-operation range (0, pi/2),
        excluding the singular endpoints.
    n_scan : int
        Number of samples used to isolate the physical root from the
        residual's spurious poles before refining with Brent's method (see
        `_select_bracket`).

    Returns
    -------
    dict
        phi, a, a_prime, Cl, Cd, Ct, Cq (all floats). Ct/Cq here are the
        tangential-force and torque coefficients at this station (Cq = Ct,
        reported separately for clarity when this is later integrated over
        the blade -- torque = force x radius, not done at the single-station
        stage).

    Raises
    ------
    ValueError
        If no physical root is found over the given range (e.g. physically
        inconsistent inputs).
    """

    phi_lo, phi_hi = _select_bracket(station, bracket, n_scan)
    if phi_lo == phi_hi:
        phi = phi_lo
    else:
        phi = brentq(residual, phi_lo, phi_hi, args=(station,), xtol=1e-12, rtol=1e-12)

    a, a_prime, cl, cd, cn, ct = _induction_factors(phi, station)

    return {
        "phi": phi,
        "a": a,
        "a_prime": a_prime,
        "Cl": cl,
        "Cd": cd,
        "Ct": ct,
        "Cq": ct,
    }
