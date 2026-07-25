"""
Stage 2: Prandtl tip-loss (and optional hub-loss) correction factor.

Convention used, and why it matters
------------------------------------
BEM codes are inconsistent about *where* the loss factor F enters the
induction equations, and the two conventions give different numbers:

  (a) F multiplies the momentum-theory side: the annulus momentum balance
      assumes uniform induced velocity across a full annulus, which breaks
      down near a finite number of blade tips (trailed-vortex effects), so
      F corrects the *momentum* equation, e.g.
          4 F sin^2(phi) a (1-a) = sigma' Cn
      i.e. F appears alongside the "4 sin^2(phi)" momentum term.

  (b) F scales the blade-element force coefficients (Cn, Ct) directly before
      they enter the momentum balance.

This module (and station.py) uses convention (a) -- F enters the momentum
side only, exactly as in Ning (2014) and the standard Prandtl/Glauert
derivation: blade-element theory (Cl, Cd, hence Cn, Ct) is a statement about
aerodynamic forces on the actual blade and is not touched by F; momentum
theory (the annulus-averaged induced velocity) is what F is correcting for.
Folding F into Cn/Ct instead (convention (b)) would double-count the
correction if later combined with a Glauert/Buhl high-thrust correction
(Stage 3), since that correction is also formulated in terms of momentum
theory's a(phi) relation.

Reference
---------
Ning, S.A. (2014). "A simple solution method for the blade element momentum
equations with guaranteed convergence." Wind Energy, 17(9), 1327-1345.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math


def tip_loss_factor(r, R, B, phi):
    """
    Prandtl tip-loss factor F_tip in (0, 1].

    F_tip = (2/pi) * arccos(exp(-f)),  f = (B/2) * (R - r) / (r * sin(phi))

    Parameters
    ----------
    r : float
        Local radius, m. Must be < R (see module/station.py notes on the
        r -> R limit -- F_tip -> 0 there and the residual degenerates).
    R : float
        Rotor (blade tip) radius, m.
    B : int
        Number of blades.
    phi : float
        Inflow angle, radians.

    Returns
    -------
    float
        F_tip, bounded in (0, 1]. Returns 0.0 exactly at/beyond the tip
        (r >= R) rather than raising, since `station._select_bracket` scans
        many phi values and must not crash mid-scan; callers that need a
        physically meaningful solve must keep r strictly < R (see
        station.py's r >= R guard).
    """

    if r >= R:
        return 0.0

    sin_phi = abs(math.sin(phi))
    if sin_phi < 1e-8:
        # phi -> 0 limit: f -> +inf, F_tip -> 1 (no correction) for any r < R.
        return 1.0

    f = (B / 2.0) * (R - r) / (r * sin_phi)
    return (2.0 / math.pi) * math.acos(min(1.0, math.exp(-f)))


def hub_loss_factor(r, r_hub, B, phi):
    """
    Prandtl hub-loss factor F_hub in (0, 1], mirroring tip_loss_factor.

    F_hub = (2/pi) * arccos(exp(-f)),  f = (B/2) * (r - r_hub) / (r_hub * sin(phi))

    Parameters
    ----------
    r : float
        Local radius, m.
    r_hub : float or None
        Hub radius, m. If None or <= 0, hub-loss is disabled and this
        returns 1.0 unconditionally (most Stage 1/2 test cases have no hub
        geometry defined -- this keeps hub-loss strictly opt-in).
    B : int
        Number of blades.
    phi : float
        Inflow angle, radians.

    Returns
    -------
    float
        F_hub, bounded in (0, 1].
    """

    if r_hub is None or r_hub <= 0:
        return 1.0
    if r <= r_hub:
        return 0.0

    sin_phi = abs(math.sin(phi))
    if sin_phi < 1e-8:
        return 1.0

    f = (B / 2.0) * (r - r_hub) / (r_hub * sin_phi)
    return (2.0 / math.pi) * math.acos(min(1.0, math.exp(-f)))


def combined_loss_factor(r, R, B, phi, r_hub=None):
    """F = F_tip * F_hub (F_hub = 1.0 if r_hub is None/disabled)."""

    return tip_loss_factor(r, R, B, phi) * hub_loss_factor(r, r_hub, B, phi)
