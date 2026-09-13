"""
Station kernels: the BEM residual, the power integrand, and every first
partial of both, at one station (Phase 3, B1).

This is `docs/adjoint_derivation.md` §2-§6 transcribed, one scalar station at
a time. The *values* of `R` and `q` are computed with the same calls, in the
same order, as `bem.station.residual` and `bem.rotor.solve_rotor`, so they
match those to the last bit; the partials are the hand-derived chain, and
they are verified against a complex step of the code's own residual in
`tests/test_adjoint_partials.py` (Tier 1) -- never against themselves.

Notation (subscripts dropped; `s = sin phi`, `k = cos phi`):

    sigma = B c / (2 pi r)                      Re = W c / nu  (dRe/dc = Re/c)
    alpha = phi - theta                         Cl, Cd = polar(alpha) at Re
    Cn = Cl k + Cd s                            Ct = Cl s - Cd k
    F = F_tip(phi) F_hub(phi)                   Y = sigma Cn / s^2
    a = a(Y, F)   momentum a = Y/(4F + Y) if that is <= 0.4, else Buhl
    T = sigma Ct / (4 F s)
    R = s/(1 - a) - (k - T)/lambda_r
    w = V (1 - a)/s                             q = 1/2 rho w^2 B c Ct r

The state is `phi` alone; `a` is an explicit function of `(phi, c, theta)`
and is differentiated as one. The Buhl branch is differentiated implicitly
from `Y (1 - a)^2 = Ct_B(a; F)`, which is branch-free within the region and
independent of which of the two algebraically equal expressions
`corrected_axial_induction` selected on conditioning.

Complex-safe: every transcendental dispatches on type, and no array with a
real dtype is allocated. The kernel is scalar; `adjoint.system` loops it
over (operating point, station).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import cmath
import math
from dataclasses import dataclass

from bem.corrections import (
    BUHL_AC,
    buhl_gammas,
    combined_loss_factor,
    corrected_axial_induction,
)


def _sin(x):
    return cmath.sin(x) if isinstance(x, complex) else math.sin(x)


def _cos(x):
    return cmath.cos(x) if isinstance(x, complex) else math.cos(x)


def _exp(x):
    return cmath.exp(x) if isinstance(x, complex) else math.exp(x)


def _sqrt(x):
    return cmath.sqrt(x) if isinstance(x, complex) else math.sqrt(x)


# ---------------------------------------------------------------------------
# the Prandtl loss factor and dF/dphi
# ---------------------------------------------------------------------------

def _prandtl_derivative(numerator, phi):
    """
    `dF/dphi` for `F = (2/pi) acos(exp(-f))`, `f = numerator / sin phi`.

        df/dphi   = -f cos phi / sin phi
        dF/df     = (2/pi) e^{-f} / sqrt(1 - e^{-2f})
        dF/dphi   = -(2/pi) f e^{-f} cos phi / (sin phi sqrt(1 - e^{-2f}))

    Mirrors the guards of `corrections.tip_loss_factor`: where the code
    returns the constant 1 (`sin phi < 1e-8`) or the clamp is live, the
    derivative is 0. On the momentum region `sin phi > 0`, so the code's
    `abs` is inert; it is honoured on `.real` anyway.
    """

    s = _sin(phi)
    if s.real < 0.0:
        s = -s
    if s.real < 1e-8:
        return 0.0

    f = numerator / s
    e = _exp(-f)
    if e.real >= 1.0:
        return 0.0
    k = _cos(phi)
    return -(2.0 / math.pi) * f * e * k / (s * _sqrt(1.0 - e * e))


def _prandtl_value(numerator, phi):
    """The factor itself, on the same guards, for the product rule."""

    s = _sin(phi)
    if s.real < 0.0:
        s = -s
    if s.real < 1e-8:
        return 1.0
    e = _exp(-numerator / s)
    if e.real >= 1.0:
        return 0.0
    acos = cmath.acos if isinstance(e, complex) else math.acos
    return (2.0 / math.pi) * acos(e)


def loss_factor_and_derivative(r, R, B, phi, r_hub):
    """
    `(F, dF/dphi)` with `F = F_tip F_hub`.

    The value is `corrections.combined_loss_factor` itself -- the code's own
    call, so it is bit-identical to what the residual used. The derivative is
    the product rule over the two factors, each from `_prandtl_derivative`.
    `F` has no chord or twist dependence.
    """

    F = combined_loss_factor(r, R, B, phi, r_hub)

    if r >= R:
        return F, 0.0
    tip_num = (B / 2.0) * (R - r) / r
    F_tip, dF_tip = _prandtl_value(tip_num, phi), _prandtl_derivative(tip_num, phi)

    if r_hub is None or r_hub <= 0 or r <= r_hub:
        F_hub, dF_hub = (1.0, 0.0) if (r_hub is None or r_hub <= 0) else (0.0, 0.0)
    else:
        hub_num = (B / 2.0) * (r - r_hub) / r_hub
        F_hub, dF_hub = _prandtl_value(hub_num, phi), _prandtl_derivative(hub_num, phi)

    return F, dF_tip * F_hub + F_tip * dF_hub


# ---------------------------------------------------------------------------
# the axial induction and its partials in (Y, F)
# ---------------------------------------------------------------------------

def induction_and_partials(Y, F):
    """
    `(a, da/dY, da/dF, buhl)`.

    The value is `corrections.corrected_axial_induction(Y, F)` -- the code's
    own call. The branch is decided exactly as the code decides it, on
    `a_momentum.real <= BUHL_AC`.

    Momentum branch, explicitly:

        a = Y/(4F + Y)      a_Y = 4F/(4F + Y)^2      a_F = -4Y/(4F + Y)^2

    Buhl branch, by implicit differentiation of the defining equation
    `G(a; Y, F) = Y (1 - a)^2 - Ct_B(a; F) = 0` with
    `Ct_B = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2`:

        G_a = -2Y(1 - a) - (4F - 40/9) - 2(50/9 - 4F) a
        G_Y = (1 - a)^2
        G_F = -4a(1 - a)
        a_Y = -G_Y/G_a      a_F = -G_F/G_a

    The `gamma2 < 0` branch of the code is unreachable for positive loading
    and is not modelled: reaching it here raises rather than differentiating
    a non-physical value.
    """

    a = corrected_axial_induction(Y, F)
    a_momentum = Y / (4.0 * F + Y)

    if a_momentum.real <= BUHL_AC:
        denominator = (4.0 * F + Y) * (4.0 * F + Y)
        return a, 4.0 * F / denominator, -4.0 * Y / denominator, False

    _gamma1, gamma2, _gamma3 = buhl_gammas(Y / (4.0 * F), F)
    if gamma2.real < 0.0:
        raise ValueError(
            f"gamma2 = {gamma2} < 0 at Y = {Y}, F = {F}: past the a-pole, a branch "
            "the adjoint does not model (see corrections.corrected_axial_induction)")

    one_minus_a = 1.0 - a
    G_a = (-2.0 * Y * one_minus_a - (4.0 * F - 40.0 / 9.0)
           - 2.0 * (50.0 / 9.0 - 4.0 * F) * a)
    G_Y = one_minus_a * one_minus_a
    G_F = -4.0 * a * one_minus_a
    return a, -G_Y / G_a, -G_F / G_a, True


# ---------------------------------------------------------------------------
# the station kernel
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class StationPartials:
    """
    Everything the adjoint needs from one station at one trial `phi`.

    `residual` and `q` are the values; the `d*` fields are the partials with
    respect to `phi`, the local chord `c` and the local twist `theta` (per
    station -- the constant spline Jacobians take them to the design vector
    in `adjoint.system`). The remaining fields are diagnostics.
    """

    residual: object
    dR_dphi: object
    dR_dc: object
    dR_dtheta: object
    q: object
    dq_dphi: object
    dq_dc: object
    dq_dtheta: object
    a: object
    F: object
    dF_dphi: object
    alpha: object
    cl: object
    cd: object
    cn: object
    ct: object
    buhl: bool


def station_partials(phi, chord, twist, reynolds, r, lam_r, R, B, r_hub, polar,
                     v_inf, air_density, derivatives=True):
    """
    The residual `R`, the power integrand `q`, and their partials in
    `(phi, c, theta)` at one station.

    Parameters
    ----------
    phi : float or complex
        Inflow angle, rad.
    chord, twist : float or complex
        Local chord (m) and twist (rad) at this station.
    reynolds : float or complex
        The station's Reynolds number, `W c / nu` -- the same zero-induction
        estimate `bem.rotor.solve_rotor` makes, so `dRe/dc = reynolds/chord`.
    r, lam_r, R, B, r_hub : float
        Radius, local speed ratio `Omega r / V`, tip radius, blade count, hub
        radius (or None). Constants of the station.
    polar : polars.polar.CachedPolar
        The polar at `reynolds` -- the object the residual read, so the value
        of `R` is the residual's own.
    v_inf, air_density : float
        Wind speed and density of the operating point, for `q`.
    derivatives : bool
        With `False`, only the values `R`, `q` and the diagnostics are
        computed (the four polar-slope lookups are skipped) and every
        partial field is `None`. Same value code either way; this is what
        `BEMSystem.J` uses so a complex step of `J` costs a value, not a
        gradient.

    Returns
    -------
    StationPartials
    """

    s = _sin(phi)
    k = _cos(phi)
    alpha = phi - twist

    # -- polar and its four partials (per rad, per unit Re) -----------------
    cl = polar.cl(alpha)
    cd = polar.cd(alpha)

    # -- force coefficients, same expressions as station._blade_element ------
    cn = cl * k + cd * s
    ct = cl * s - cd * k

    # -- loss factor ---------------------------------------------------------
    F, F_phi = loss_factor_and_derivative(r, R, B, phi, r_hub)

    # -- loading and induction -----------------------------------------------
    sigma = B * chord / (2.0 * math.pi * r)          # StationParams.solidity
    sigma_c = B / (2.0 * math.pi * r)
    s2 = s * s
    Y = sigma * cn / s2                               # as coded in residual()
    a, a_Y, a_F, buhl = induction_and_partials(Y, F)

    # -- the residual, in station.residual's own operation order -------------
    axial = s / (1.0 - a)
    tangential = (k - sigma * ct / (4.0 * F * s)) / lam_r
    residual = axial - tangential

    # -- the power integrand, in solve_rotor's own operation order -----------
    w = v_inf * (1.0 - a) / s
    q = 0.5 * air_density * w ** 2 * B * chord * ct * r

    if not derivatives:
        return StationPartials(
            residual=residual, dR_dphi=None, dR_dc=None, dR_dtheta=None,
            q=q, dq_dphi=None, dq_dc=None, dq_dtheta=None,
            a=a, F=F, dF_dphi=F_phi, alpha=alpha, cl=cl, cd=cd, cn=cn, ct=ct, buhl=buhl,
        )

    cl_alpha = polar.dcl_dalpha(alpha)
    cd_alpha = polar.dcd_dalpha(alpha)
    cl_re = polar.dcl_dre(alpha)
    cd_re = polar.dcd_dre(alpha)
    dre_dc = reynolds / chord

    cn_phi = cl_alpha * k + cd_alpha * s - ct
    ct_phi = cl_alpha * s - cd_alpha * k + cn
    cn_theta = -(cl_alpha * k + cd_alpha * s)
    ct_theta = -(cl_alpha * s - cd_alpha * k)
    cn_c = (cl_re * k + cd_re * s) * dre_dc
    ct_c = (cl_re * s - cd_re * k) * dre_dc

    # -- loading and induction partials -----------------------------------------
    Y_phi = sigma * cn_phi / s2 - 2.0 * sigma * cn * k / (s2 * s)
    Y_c = sigma_c * cn / s2 + sigma * cn_c / s2
    Y_theta = sigma * cn_theta / s2

    a_phi = a_Y * Y_phi + a_F * F_phi
    a_c = a_Y * Y_c
    a_theta = a_Y * Y_theta

    one_minus_a = 1.0 - a
    one_minus_a2 = one_minus_a * one_minus_a
    four_F_s = 4.0 * F * s
    T_phi = (sigma / 4.0) * (ct_phi / (F * s) - ct * F_phi / (F * F * s) - ct * k / (F * s2))
    T_c = (sigma_c * ct + sigma * ct_c) / four_F_s
    T_theta = sigma * ct_theta / four_F_s

    dR_dphi = k / one_minus_a + s * a_phi / one_minus_a2 + (s + T_phi) / lam_r
    dR_dc = s * a_c / one_minus_a2 + T_c / lam_r
    dR_dtheta = s * a_theta / one_minus_a2 + T_theta / lam_r

    # -- power-integrand partials -----------------------------------------------
    w_phi = -v_inf * (a_phi / s + one_minus_a * k / s2)
    w_c = -v_inf * a_c / s
    w_theta = -v_inf * a_theta / s
    prefactor = 0.5 * air_density * B * r
    dq_dphi = prefactor * chord * (2.0 * w * w_phi * ct + w * w * ct_phi)
    dq_dc = prefactor * (w * w * ct + chord * (2.0 * w * w_c * ct + w * w * ct_c))
    dq_dtheta = prefactor * chord * (2.0 * w * w_theta * ct + w * w * ct_theta)

    return StationPartials(
        residual=residual, dR_dphi=dR_dphi, dR_dc=dR_dc, dR_dtheta=dR_dtheta,
        q=q, dq_dphi=dq_dphi, dq_dc=dq_dc, dq_dtheta=dq_dtheta,
        a=a, F=F, dF_dphi=F_phi, alpha=alpha, cl=cl, cd=cd, cn=cn, ct=ct, buhl=buhl,
    )
