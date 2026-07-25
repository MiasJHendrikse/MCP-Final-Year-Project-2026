"""
Stage 3: single blade-element station BEM solver with Prandtl tip/hub loss
and the Glauert/Buhl high-thrust (turbulent wake state) correction.

Implements Ning (2014)'s reduction of the coupled axial/tangential induction
equations to a single residual equation in the inflow angle phi, solved with
a bracketed root-finder (Brent's method) rather than the classical
fixed-point iteration on (a, a') -- this is what keeps convergence robust
near high induction and keeps the solve a single scalar equation, which
matters later when this needs a discrete adjoint (one clean residual, no
iterative state, no branching on convergence history).

Stage 2 added the Prandtl tip-loss factor F (and optional hub-loss) on the
momentum-theory side of the induction equations. Stage 3 adds the Buhl
(2005) empirical correction that replaces the plain momentum Ct(a) relation
once a exceeds ~0.4 (turbulent wake state), where plain momentum theory
predicts an unphysical decrease in thrust -- see corrections.py for the
closed-form used and why it is C0/C1-continuous with the plain relation at
the switch, so it does not introduce a kink into this residual. Still no
multi-station loop or real airfoil data; those are deferred to later
sessions (see PROJECT_PLAN.md Phase 1).

r -> R (tip) limit
------------------
F -> 0 exactly at r = R makes the momentum-consistent a(phi) and a'(phi)
collapse to the phi-independent constants 1 and -1 respectively (see
corrections.tip_loss_factor docstring), which makes `residual` identically
zero for *every* phi -- a completely degenerate, non-unique root. Station
inputs must therefore keep r strictly < R; `StationParams.__post_init__`
raises if r >= R rather than letting the solver silently fail on a
degenerate residual.

Reference
---------
Ning, S.A. (2014). "A simple solution method for the blade element momentum
equations with guaranteed convergence." Wind Energy, 17(9), 1327-1345.
Hansen, M.O.L. (2008). "Aerodynamics of Wind Turbines" (2nd ed.) -- used here
as the classical, uncorrected BEM formulation to cross-check Stage 1 against,
since with no corrections applied the two formulations must agree
analytically (Stage 2's F=1 limit reduces to the same case).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
from dataclasses import dataclass

from scipy.optimize import brentq

from bem.corrections import combined_loss_factor, corrected_axial_induction


@dataclass
class StationParams:
    """
    Inputs for one fixed blade-element station.

    Parameters
    ----------
    r : float
        Local radius, m. Must be strictly less than R (see module docstring
        on the r -> R limit).
    chord : float
        Local chord length, m.
    twist : float
        Local twist (+ pitch), radians, measured such that
        alpha = phi - twist.
    airfoil : object
        Callable/duck-typed polar with .cl(alpha) and .cd(alpha) in radians
        (e.g. bem.airfoil.LinearPolar). Stage 1/2 use a synthetic polar only.
    tsr : float
        Local speed ratio at this station, Omega * r / Vinf.
    R : float
        Rotor (blade tip) radius, m. Used only for the Prandtl tip-loss
        factor F (see corrections.tip_loss_factor).
    n_blades : int
        Number of blades, B.
    r_hub : float or None
        Hub radius, m. Optional -- if None (default), hub-loss is disabled
        (F_hub = 1) rather than requiring every station to define hub
        geometry (see corrections.hub_loss_factor).
    """

    r: float
    chord: float
    twist: float
    airfoil: object
    tsr: float
    R: float
    n_blades: int = 3
    r_hub: float = None

    def __post_init__(self):
        if self.r >= self.R:
            raise ValueError(
                f"r={self.r} must be strictly less than R={self.R}; F_tip -> 0 "
                "exactly at/beyond the tip degenerates the residual (see "
                "module docstring)."
            )

    @property
    def solidity(self):
        """Local solidity sigma' = B*c / (2*pi*r)."""
        return self.n_blades * self.chord / (2.0 * math.pi * self.r)


def _blade_element_and_induction(phi, station: StationParams):
    """
    Shared core for residual() and _induction_factors(): blade-element force
    coefficients (untouched by F or the Buhl correction) and the
    momentum-consistent induction factors a(phi), a'(phi) -- see
    corrections.py for why F only enters here, on the momentum side, and
    not into cn/ct.

    Axial induction a uses corrections.corrected_axial_induction, which
    transparently switches to the Buhl (2005) high-thrust relation once the
    naive momentum-theory value would exceed a=0.4 (turbulent wake state) --
    see corrections.py for why this stays a single closed-form expression
    (no inner iteration) and is C0/C1-continuous at the switch, so it does
    not introduce a kink into this residual. Tangential induction a' is
    unaffected -- the Buhl correction only addresses the axial thrust
    relation's known breakdown at high a; standard practice leaves a'
    on the plain momentum relation.
    """

    alpha = phi - station.twist
    cl = station.airfoil.cl(alpha)
    cd = station.airfoil.cd(alpha)

    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    cn = cl * cos_phi + cd * sin_phi
    ct = cl * sin_phi - cd * cos_phi

    sigma = station.solidity
    F = combined_loss_factor(station.r, station.R, station.n_blades, phi, station.r_hub)

    Y = sigma * cn / (sin_phi * sin_phi)
    a = corrected_axial_induction(Y, F)
    a_prime = 1.0 / ((4.0 * F * sin_phi * cos_phi) / (sigma * ct) - 1.0)

    return a, a_prime, cl, cd, cn, ct, F


def residual(phi, station: StationParams):
    """
    Ning-style BEM residual in the inflow angle phi, with Prandtl tip/hub
    loss (F) applied on the momentum side (see corrections.py).

    Derivation: the blade-element normal/tangential force coefficients give
    F-corrected, momentum-consistent induction factors a(phi), a'(phi)
    directly (no fixed-point iteration needed); the residual is the
    mismatch between the input phi and the phi implied by the kinematic
    relation tan(phi) = (1-a) / ((1+a') * lambda_r), rearranged to avoid a
    division (so it stays finite and smooth even as phi -> 0 or a, a' ->
    singular values -- the property this formulation is chosen for). This
    kinematic relation itself does not involve F -- F only enters through
    a(phi) and a'(phi).

    Parameters
    ----------
    phi : float
        Inflow angle, radians. Must lie strictly in (0, pi/2) for this
        non-turbulent-wake-state formulation.
    station : StationParams

    Returns
    -------
    float
        Residual value; solve_station finds the phi where this is zero.
    """

    a, a_prime, _, _, _, _, _ = _blade_element_and_induction(phi, station)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    return sin_phi * (1.0 + a_prime) * station.tsr - cos_phi * (1.0 - a)


def _induction_factors(phi, station: StationParams):
    """Recompute (a, a', Cl, Cd, Cn, Ct, F) at the solved phi (no residual)."""

    return _blade_element_and_induction(phi, station)


def _select_bracket(station: StationParams, phi_range, n_scan):
    """
    Scan phi_range for the sign change nearest the zero-induction guess.

    The momentum-consistent a'(phi) used in `residual` has a genuine pole
    wherever Ct(phi) = 0 (a known BEM pathology, not specific to this
    implementation -- see Ning 2014 Sec. 3 for the full case analysis this
    motivates). A symmetric pole exists on the axial side too: the plain
    momentum relation a = Y / (4F + Y), Y = sigma*Cn/sin^2(phi), blows up
    wherever Cn crosses -4*F*sin^2(phi)/sigma (found while building Stage
    3's turbulent-wake test case, which happens to sweep through a small
    negative Cn near stall-free zero-lift; flagged here rather than left
    for a later session to rediscover -- see validate_stage3.py's explicit
    check of it). Both poles cross the residual's sign just like a true
    root does, so naively bracketing on the endpoints of a wide range can
    hand brentq a pole instead of the physical root. Scanning first and
    picking the sign change closest to the classical zero-induction inflow
    angle (phi0 = atan(1/tsr)) avoids that failure mode without needing the
    full region-classification machinery Ning uses for the fully corrected
    solver -- that full case analysis is out of scope for this solver.
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
        phi, a, a_prime, Cl, Cd, Ct, Cq, F (all floats). Ct/Cq here are the
        tangential-force and torque coefficients at this station (Cq = Ct,
        reported separately for clarity when this is later integrated over
        the blade -- torque = force x radius, not done at the single-station
        stage). F is the combined Prandtl tip/hub loss factor at the solved
        phi (F=1 recovers Stage 1's uncorrected result).

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

    a, a_prime, cl, cd, cn, ct, F = _induction_factors(phi, station)

    return {
        "phi": phi,
        "a": a,
        "a_prime": a_prime,
        "Cl": cl,
        "Cd": cd,
        "Ct": ct,
        "Cq": ct,
        "F": F,
    }
