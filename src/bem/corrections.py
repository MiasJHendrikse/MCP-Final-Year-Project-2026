"""
Stage 2/3: Prandtl tip/hub-loss and Glauert/Buhl high-thrust corrections.

Convention used, and why it matters (tip/hub loss, Stage 2)
-------------------------------------------------------------
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
correction when combined with the Glauert/Buhl high-thrust correction below,
since that correction is also formulated in terms of momentum theory's
a(phi) relation.

High-thrust (turbulent wake state) correction, Stage 3
--------------------------------------------------------
Momentum theory's local thrust relation, Ct(a) = 4 a F (1-a), peaks at a=0.5
and *decreases* for a > 0.5 -- physically wrong; measured rotors keep
loading up past that point (turbulent wake state). Buhl (2005) replaced the
classical Glauert empirical curve with a closed-form quadratic in a, chosen
specifically so it (and its slope) matches the momentum relation exactly at
a = ac = 0.4 for *any* F -- confirmed algebraically and numerically in this
module's tests, see BUHL_AC. Buhl is used here rather than classic Glauert
(which is typically defined as an iterative/graphical correction) because
its closed form inverts to an explicit quadratic-in-a solve, which is what
`station.py` needs to keep the residual a single closed-form expression per
trial phi -- no inner iteration, and no discontinuity for the root-finder or
a future adjoint to trip on.

Reference
---------
Ning, S.A. (2014). "A simple solution method for the blade element momentum
equations with guaranteed convergence." Wind Energy, 17(9), 1327-1345.
Buhl, M.L. (2005). "A New Empirical Relationship between Thrust Coefficient
and Induction Factor for the Turbulent Windmill State." NREL/TP-500-36834.

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


# ----------------------------------------------------------------------------
# Stage 3: Glauert/Buhl high-thrust (turbulent wake state) correction.
# ----------------------------------------------------------------------------

BUHL_AC = 0.4  # empirical blend point; see module docstring for why this
               # specific value gives exact C0/C1 continuity for any F.


def high_thrust_correction(a, Ct, F, ac=BUHL_AC):
    """
    Local thrust coefficient Ct(a), switching to Buhl (2005)'s empirical
    quadratic once a exceeds the turbulent-wake threshold ac.

    This is the forward a -> Ct relation only (see `corrected_axial_induction`
    for the inverse, i.e. what `station.py` actually calls to solve for a).
    Exposed standalone so it can be plotted/tested in isolation, e.g. to
    confirm continuity across the blend at a=ac.

    Parameters
    ----------
    a : float
        Axial induction factor.
    Ct : float
        Naive (uncorrected) local momentum thrust coefficient,
        Ct = 4*a*F*(1-a) -- the caller supplies this rather than having it
        recomputed here, since station.py already has it on hand.
    F : float
        Combined tip/hub loss factor (see combined_loss_factor).
    ac : float
        Blend threshold on a. Default 0.4 (Buhl 2005).

    Returns
    -------
    float
        Ct unchanged if a <= ac; the Buhl-corrected value otherwise.
    """

    if a <= ac:
        return Ct
    return (8.0 / 9.0) + (4.0 * F - 40.0 / 9.0) * a + (50.0 / 9.0 - 4.0 * F) * a * a


def corrected_axial_induction(Y, F, ac=BUHL_AC):
    """
    Solve Y*(1-a)^2 = Ct(a) for the axial induction factor a, where Ct(a) is
    `high_thrust_correction`'s piecewise Ct(a) and
    Y = sigma' * Cn / sin(phi)^2 is the blade-element-implied loading at
    this phi (see station.py -- this is the algebraic inverse the residual
    needs, not just the forward Ct(a) curve above).

    Below ac this reduces to the Stage 2 closed form a = Y / (4F + Y)
    exactly (solve the linear case directly, rather than route it through
    the ac<a branch's quadratic, to guarantee bit-identical Stage 2 results
    for a <= ac). Above ac, solves the quadratic obtained by substituting
    Buhl's Ct(a) into Y*(1-a)^2 = Ct(a):

        [Y - 50/9 + 4F] a^2 + [-2Y - 4F + 40/9] a + [Y - 8/9] = 0

    and takes the root >= ac (the branch continuous with the a <= ac side;
    the other root is a spurious artifact of squaring/rearranging and lies
    below ac).

    Parameters
    ----------
    Y : float
        sigma' * Cn / sin(phi)^2, the blade-element-implied local loading.
    F : float
        Combined tip/hub loss factor.
    ac : float
        Blend threshold (default 0.4).

    Returns
    -------
    float
        Axial induction factor a.
    """

    a_naive = Y / (4.0 * F + Y)
    if a_naive <= ac:
        return a_naive

    A = Y - 50.0 / 9.0 + 4.0 * F
    B = -2.0 * Y - 4.0 * F + 40.0 / 9.0
    C = Y - 8.0 / 9.0

    if abs(A) < 1e-14:
        # Degenerate (linear) case -- not expected for physical F in (0,1]
        # and Y > 0, but handled rather than dividing by ~0.
        return -C / B

    discriminant = B * B - 4.0 * A * C
    sqrt_disc = math.sqrt(max(0.0, discriminant))
    root_plus = (-B + sqrt_disc) / (2.0 * A)
    root_minus = (-B - sqrt_disc) / (2.0 * A)

    candidates = [r for r in (root_plus, root_minus) if r >= ac]
    if candidates:
        # Both roots >= ac only in edge cases; the one closest to a_naive
        # is the branch continuous with the low-induction solution.
        return min(candidates, key=lambda r: abs(r - a_naive))
    # No root >= ac found (can happen numerically right at the blend point,
    # or for pathological inputs) -- fall back to whichever root is closer
    # to ac rather than returning a value from the wrong branch.
    return max(root_plus, root_minus, key=lambda r: -abs(r - ac))
