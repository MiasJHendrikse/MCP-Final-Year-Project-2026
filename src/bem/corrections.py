"""
Prandtl tip/hub-loss and Glauert/Buhl high-thrust corrections.

Convention used, and why it matters (tip/hub loss)
--------------------------------------------------
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

High-thrust (turbulent wake state) correction
---------------------------------------------
Momentum theory's local thrust relation, Ct(a) = 4 a F (1-a), peaks at a=0.5
and *decreases* for a > 0.5 -- physically wrong; measured rotors keep
loading up past that point (turbulent wake state). Buhl (2005) replaced the
classical Glauert empirical curve with a closed-form quadratic in a, chosen
specifically so it (and its slope) matches the momentum relation exactly at
a = ac = 0.4 for *any* F -- confirmed algebraically and numerically in
`tests/test_corrections.py`, see BUHL_AC. Buhl is used here rather than classic
Glauert (which is typically defined as an iterative/graphical correction)
because its closed form inverts to an explicit quadratic-in-a solve, which is
what `station.py` needs to keep the residual a single closed-form expression
per trial phi -- no inner iteration, and no discontinuity for the root-finder
or a future adjoint to trip on.

The gamma1/gamma2/gamma3 form, and why it is used
-------------------------------------------------
Those gammas are *Ning's* (2014) reparameterisation of Buhl's correction, not
Buhl's own notation. They are the same equation: substituting Buhl's Ct(a)
into Y(1-a)^2 = Ct(a) gives a quadratic whose coefficients are exactly
A = 2*gamma3, B = -4*gamma1, with B^2 - 4AC = 16*gamma2. The choice between
writing it one way or the other is therefore not a modelling decision at all;
it is a decision about which algebraically identical expression to
differentiate and to hand a root-finder.

Adopted, for three reasons, in order of weight:

  1. It selects the physical root **analytically** -- the minus sign. The
     previous quadratic solve chose between the two roots with a
     nearest-value heuristic, which is not differentiable and can switch
     branch discontinuously. The adjoint is this project's main contribution
     and it differentiates this function.
  2. The solver's momentum-region classification is expressed in these same
     quantities, so the solver and the adjoint derivation are
     written in the notation of the cited method rather than translated into
     it afterwards.
  3. CCBlade -- the 0.51 %-agreement cross-validation reference -- uses this
     form, so a discrepancy against it is a discrepancy in the physics rather
     than in the algebra.

Verified equivalent to the previous implementation to 1.0e-13 worst case over
F in [0.05, 1.0] across the whole corrected region.

Reference
---------
Ning, S.A. (2014). "A simple solution method for the blade element momentum
equations with guaranteed convergence." Wind Energy, 17(9), 1327-1345.
Buhl, M.L. (2005). "A New Empirical Relationship between Thrust Coefficient
and Induction Factor for the Turbulent Windmill State." NREL/TP-500-36834.

    PROVENANCE NOT YET CLOSED. These constants should be checked against
    Buhl's paper directly --
    not against secondary sources and not against prior notes in this repo.
    That check has NOT been done: as of 2026-09-10 the report could not be
    retrieved (nrel.gov and docs.nrel.gov do not resolve, web.archive.org is
    unreachable, UNT's copy is behind a scripted interstitial, and OSTI's
    full-text link redirects to a look-alike domain that must not be trusted
    as a source). What HAS been established here is internal consistency: the
    three coefficients are uniquely pinned by C0 and C1 continuity with
    momentum theory at a = 0.4 plus the F-independent anchor Ct(1) = 2, all
    three asserted in tests/test_corrections.py. That proves they are
    *self-consistent*, not that they are *Buhl's*. The tests alone don't close
    this (docs/OUTSTANDING-INPUTS.md section 3).

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import cmath
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

    if isinstance(phi, complex):
        # Complex-step path. Every branch decision is taken on
        # the real part so an infinitesimal imaginary perturbation cannot
        # flip it; `abs` would drop the imaginary part, so the sign is
        # applied by hand. The real path below is untouched.
        return _prandtl_factor_complex((B / 2.0) * (R - r) / r, phi)

    sin_phi = abs(math.sin(phi))
    if sin_phi < 1e-8:
        # phi -> 0 limit: f -> +inf, F_tip -> 1 (no correction) for any r < R.
        return 1.0

    f = (B / 2.0) * (R - r) / (r * sin_phi)
    return (2.0 / math.pi) * math.acos(min(1.0, math.exp(-f)))


def _prandtl_factor_complex(numerator, phi):
    """
    `(2/pi) acos(exp(-numerator / |sin phi|))` for complex-typed `phi`.

    The complex-step twin of the real path in `tip_loss_factor` and
    `hub_loss_factor`, with `numerator = (B/2)(R - r)/r` (tip) or
    `(B/2)(r - r_hub)/r_hub` (hub). Branch decisions on `.real`; the
    `min(1, exp(-f))` clamp is inert for any station strictly inside
    (r_hub, R) since `f > 0` there, and is kept on `.real` for the same
    reason the guards are.
    """

    sin_phi = cmath.sin(phi)
    if sin_phi.real < 0.0:
        sin_phi = -sin_phi
    if sin_phi.real < 1e-8:
        return 1.0

    f = numerator / sin_phi
    exp_minus_f = cmath.exp(-f)
    if exp_minus_f.real > 1.0:
        return 0.0
    return (2.0 / math.pi) * cmath.acos(exp_minus_f)


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
        returns 1.0 unconditionally (most single-station test cases have no hub
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

    if isinstance(phi, complex):
        # Complex-step path; see `_prandtl_factor_complex`.
        return _prandtl_factor_complex((B / 2.0) * (r - r_hub) / r_hub, phi)

    sin_phi = abs(math.sin(phi))
    if sin_phi < 1e-8:
        return 1.0

    f = (B / 2.0) * (r - r_hub) / (r_hub * sin_phi)
    return (2.0 / math.pi) * math.acos(min(1.0, math.exp(-f)))


def combined_loss_factor(r, R, B, phi, r_hub=None):
    """F = F_tip * F_hub (F_hub = 1.0 if r_hub is None/disabled)."""

    return tip_loss_factor(r, R, B, phi) * hub_loss_factor(r, r_hub, B, phi)


# ----------------------------------------------------------------------------
# Glauert/Buhl high-thrust (turbulent wake state) correction.
# ----------------------------------------------------------------------------

#: Buhl's empirical blend point, the axial induction at which the momentum
#: relation is handed over to the empirical quadratic.
#:
#: This is a *constant*, not a parameter, and deliberately so. Both functions
#: below used to accept `ac=BUHL_AC` as an
#: argument, which was not honourable: the quadratic's coefficients (8/9, 40/9,
#: 50/9) are *derived from* ac = 0.4 and do not move with it. Passing ac = 0.5
#: would have switched branches at 0.5 while the quadratic still matched
#: momentum theory at 0.4, opening both a value and a slope discontinuity --
#: silently, since nothing checked. No caller ever passed a non-default value,
#: so it never bit. A knob that produces a wrong answer if turned is worse than
#: no knob; the same reasoning removed the old air-density default.
BUHL_AC = 0.4


def _sqrt(x):
    """`math.sqrt`, or `cmath.sqrt` if `x` is complex-typed (complex-step)."""

    if isinstance(x, complex):
        return cmath.sqrt(x)
    return math.sqrt(x)


def high_thrust_correction(a, Ct, F):
    """
    Local thrust coefficient Ct(a), switching to Buhl (2005)'s empirical
    quadratic once a exceeds the turbulent-wake threshold BUHL_AC.

        Ct = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2

    This is the forward a -> Ct relation only (see `corrected_axial_induction`
    for the inverse, i.e. what `station.py` actually calls to solve for a).
    Exposed standalone so it can be plotted/tested in isolation, e.g. to
    confirm continuity across the blend at a = BUHL_AC.

    The three coefficients are pinned uniquely by three conditions, which is
    what makes them self-verifying: value and slope both match the momentum
    relation 4aF(1-a) at a = 0.4 for *any* F, and Ct(1) = 2. Note the last one
    is F-independent -- the F terms cancel at a = 1 -- which is a
    transcription check costing one line of arithmetic. See
    `tests/test_corrections.py`, where all three are asserted.

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

    Returns
    -------
    float
        Ct unchanged if a <= BUHL_AC; the Buhl-corrected value otherwise.
    """

    if a.real <= BUHL_AC:
        return Ct
    return (8.0 / 9.0) + (4.0 * F - 40.0 / 9.0) * a + (50.0 / 9.0 - 4.0 * F) * a * a


def buhl_gammas(kappa, F):
    """
    Ning (2014)'s gamma1, gamma2, gamma3 for the empirical (turbulent-wake)
    region, given kappa = sigma'*Cn / (4 F sin^2 phi) and the loss factor F.

        gamma1 = 2 F kappa - (10/9 - F)
        gamma2 = 2 F kappa - F (4/3 - F)
        gamma3 = 2 F kappa - (25/9 - 2 F)

    Exposed rather than inlined so the tests can assert the two identities
    the closed form below rests on, independently of the code that uses them:

        gamma1^2 - gamma2 = gamma3 * (Y - 8/9) / 2          (Y = 4 F kappa)
        at kappa = 2/3:  sqrt(gamma2) = F, exactly, for any F

    The second is what makes the blend C0 by construction rather than by
    numerical luck -- substituting it gives a = 2/5 identically, for every F.
    """

    return (2.0 * F * kappa - (10.0 / 9.0 - F),
            2.0 * F * kappa - F * (4.0 / 3.0 - F),
            2.0 * F * kappa - (25.0 / 9.0 - 2.0 * F))


def corrected_axial_induction(Y, F):
    """
    Solve Y*(1-a)^2 = Ct(a) for the axial induction factor a, where Ct(a) is
    `high_thrust_correction`'s piecewise Ct(a) and
    Y = sigma' * Cn / sin(phi)^2 is the blade-element-implied loading at
    this phi (see station.py -- this is the algebraic inverse the residual
    needs, not just the forward Ct(a) curve above).

    Below BUHL_AC this is the tip-loss-corrected closed form a = Y / (4F + Y), written
    exactly that way (rather than as kappa/(1+kappa), which is algebraically
    identical but not bit-identical) so results in the momentum region are
    unchanged to the last bit by this reparameterisation.

    Form, and why it changed
    ------------------------
    Above the blend this used to solve the quadratic obtained by substituting
    Buhl's Ct(a) into Y(1-a)^2 = Ct(a),

        [Y - 50/9 + 4F] a^2 + [-2Y - 4F + 40/9] a + [Y - 8/9] = 0

    with `math.sqrt`, and then *chose between the two roots by heuristic* --
    "the one closest to a_naive", with a fallback of "whichever is closest to
    ac". Those choices were dead code for physical inputs, but a nearest-value
    heuristic has no place in a function that the adjoint differentiates: it can
    switch branch discontinuously, and nothing about it is differentiable.

    It is now written in Ning (2014)'s gamma1/gamma2/gamma3 reparameterisation
    of the same equation, which selects the branch *analytically* -- the minus
    sign -- so there is nothing left to choose:

        a = (gamma1 - sqrt(gamma2)) / gamma3

    This is the identical function, not an approximation: A = 2*gamma3,
    B = -4*gamma1, and B^2 - 4AC = 16*gamma2 exactly. Verified against the
    previous implementation to 1.0e-13 worst case over F in [0.05, 1.0] across
    the whole corrected region, and against a 60-digit Decimal reference in
    `tests/test_corrections.py`. Adopting it also means the adjoint derivation
    and the solver's region classification are expressed in the same quantities as
    the cited method, and it matches CCBlade -- the 0.51 %-agreement reference.

    Conditioning: two forms, one function
    --------------------------------------
    gamma3 vanishes at kappa = (25/9 - 2F)/(2F), which lies *inside* the
    corrected region whenever F <~ 0.86 -- and F falls toward 0 at the tip,
    which is exactly where this branch is live (measured: 1-3 stations per
    Phase VI operating point, all at F = 0.10-0.29). It is a removable
    singularity, since gamma1 = sqrt(gamma2) there too, so the published form
    does not blow up so much as cancel: measured relative error 1.9e-11 at
    gamma3 = 1.4e-5, 1.1e-9 at 1.4e-7, and ZeroDivisionError at 0.

    Rationalising by (gamma1 + sqrt(gamma2)) gives an algebraically identical
    expression with gamma3 removed entirely:

        a = (Y/2 - 4/9) / (gamma1 + sqrt(gamma2))

    exact to machine precision through gamma3 = 0. It has its own degeneracy
    where gamma1 + sqrt(gamma2) = 0, which happens only on Y = 8/9 (and only
    for F <= 1/3) -- and gamma3 is never zero there, for any physical F <= 1,
    so the two are never ill-conditioned at once. Taking whichever denominator
    is larger in magnitude is therefore always safe.

    That choice is a *branch*, and this docstring has just finished explaining
    why branches were removed, so the distinction matters: the two expressions
    are the same function, so the value and every derivative are continuous
    across the switch. It selects between two ways of computing one number,
    not between two different numbers -- unlike the root-selection heuristic
    it replaces. It is taken on `.real` so an infinitesimal complex-step
    perturbation cannot flip it (same rule as `polars.interpolant._bracket`).

    Parameters
    ----------
    Y : float or complex
        sigma' * Cn / sin(phi)^2, the blade-element-implied local loading.
    F : float
        Combined tip/hub loss factor.

    Returns
    -------
    float or complex
        Axial induction factor a.
    """

    a_momentum = Y / (4.0 * F + Y)
    if a_momentum.real <= BUHL_AC:
        return a_momentum

    kappa = Y / (4.0 * F)
    gamma1, gamma2, gamma3 = buhl_gammas(kappa, F)

    if gamma2.real < 0.0:
        # Past the a-pole: no real root, and no physical meaning either.
        #
        # gamma2 < 0 is reachable from here for exactly one set of inputs, and
        # it is worth naming precisely because the previous implementation hid
        # it behind `math.sqrt(max(0.0, discriminant))`. Reaching this branch
        # at all needs a_momentum = Y/(4F + Y) > 0.4. For Y > 0 that forces
        # Y > (8/3)F, hence kappa > 2/3, hence gamma2 > F^2 > 0 -- so a
        # negative gamma2 is impossible for any positive loading. It requires
        # Y < -4F: negative loading, past the point where 4F + Y changes sign
        # and a(Y) has its pole (see station.py's module docstring, which
        # documents that pole and why the residual's sign changes there).
        #
        # station._select_bracket evaluates the residual at ~2000 trial phi
        # per station, and some of those land beyond that pole. They are not
        # solutions and are never returned as one -- they exist only so the
        # scan can keep walking. What is returned here is the same value the
        # previous implementation produced for them (gamma1/gamma3, which is
        # -B/2A, the twice-repeated root of the clamped quadratic), so the
        # scan sees an unchanged residual landscape.
        #
        # The solver now uses Ning's region classification instead of a scan,
        # at which point nothing evaluates this function past the pole and
        # this branch can go with it.
        return gamma1 / gamma3

    root = _sqrt(gamma2)

    if abs(gamma3.real) >= abs((gamma1 + root).real):
        return (gamma1 - root) / gamma3
    return (0.5 * Y - 4.0 / 9.0) / (gamma1 + root)
