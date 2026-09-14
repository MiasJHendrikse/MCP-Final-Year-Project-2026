"""
Single blade-element station BEM solver: Prandtl tip/hub loss (Stage 2) and
the Glauert/Buhl high-thrust (turbulent wake state) correction (Stage 3).
Stage 4 (rotor.py) loops this over stations spanwise with real polar data;
this module itself is unchanged in Stage 4 beyond the r_hub guard below.

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
the switch, so it does not introduce a kink into this residual.

r -> R (tip) and r -> r_hub limits
-----------------------------------
F -> 0 exactly at r = R makes the momentum-consistent a(phi) and a'(phi)
collapse to the phi-independent constants 1 and -1 respectively (see
corrections.tip_loss_factor docstring), which makes `residual` identically
zero for *every* phi -- a completely degenerate, non-unique root. Station
inputs must therefore keep r strictly < R; `StationParams.__post_init__`
raises if r >= R rather than letting the solver silently fail on a
degenerate residual. The same thing happens at r = r_hub (F_hub -> 0) --
found while building Stage 4's rotor geometry, where the innermost station
was initially placed exactly at r_hub; `__post_init__` raises for that too.

What Task 5 changed
--------------------
Phase 0 had the right shape -- one scalar residual in phi, a(phi) and a'(phi)
in closed form, a bracketed root-find -- but three things stood between it and
a solver that can be differentiated and swept:

  * The residual was the *multiplied* form, sin(phi)(1+a')lambda_r -
    cos(phi)(1-a), not Ning's. Same zeros, different function, and the
    guarantees are properties of his. Now Ning's -- see `residual`.
  * The bracket came from a 2000-point scan plus "take the sign change nearest
    atan(1/lambda_r)". It could return a pole as a converged answer (the
    module used to document that both poles cross the residual's sign exactly
    as a root does, with the nearest-to-phi0 heuristic as the only defence),
    and it raised ValueError when it found nothing. Now derived from the
    momentum region's boundary -- see `momentum_region_bracket`.
  * Nothing checked the residual after brentq returned, and the tolerances
    were on the root location phi, not on the residual, and not relative to
    anything. Now every station carries a checked |R|/R0 -- see
    `solve_station` and RESIDUAL_RTOL.

Cost fell from 9.0 s per operating point (audit baseline) to about 15 ms, of
which Task 4 contributed the polar layer and this task the bracket: roughly
48 residual evaluations per station instead of 2000+.

    PROVENANCE: the region classification in `momentum_region_bracket` was
    DERIVED here from the residual's structure and verified numerically over
    247 stations (see verification/phase_vi/). It is not transcribed from
    Ning (2014), which could not be obtained -- it is a paywalled journal
    article and no copy was reachable. The mathematics is checked; what is
    not established is that it matches Ning's own region definitions. Ning's
    propeller-brake region on phi in (-pi/4, 0) is deliberately not
    implemented; a station that would need it is reported as non-convergent
    rather than solved by an untested branch. Treat "the Ning formulation"
    in any writeup as needing this closed first, alongside the Buhl
    provenance note in corrections.py.

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

import cmath
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
        if self.r_hub is not None and self.r_hub > 0 and self.r <= self.r_hub:
            # Mirrors the r >= R guard above: F_hub -> 0 exactly at/inside the
            # hub degenerates the residual the same way F_tip -> 0 does at
            # the tip (see corrections.hub_loss_factor) -- found while
            # building Stage 4's rotor geometry, where the innermost station
            # was placed exactly at r_hub.
            raise ValueError(
                f"r={self.r} must be strictly greater than r_hub={self.r_hub}; "
                "F_hub -> 0 exactly at/inside the hub degenerates the residual."
            )

    @property
    def solidity(self):
        """Local solidity sigma' = B*c / (2*pi*r)."""
        return self.n_blades * self.chord / (2.0 * math.pi * self.r)
#: Lower edge of the momentum region. phi = 0 is singular (sin phi -> 0 puts
#: the loading Y = sigma*Cn/sin^2(phi) and kappa' through the roof), so the
#: region is open at 0 and this is how far in the search starts.
PHI_EPS = 1e-7

#: Upper edge, included. Unlike the old scan's pi/2 - 1e-4, this is pi/2
#: exactly: the residual form below is finite there (see `residual`), and the
#: sign of R(pi/2) is one of the two analytic facts the bracket rests on.
PHI_MAX = math.pi / 2

#: Convergence criterion: |R(phi*)| <= RESIDUAL_RTOL * R0, where R0 is the
#: residual norm at the bracket endpoints.
#:
#: Relative to R0, never absolute, per the brief's item 5: "an absolute
#: tolerance below machine round-off for the residual's natural magnitude will
#: silently fail to converge while appearing to iterate". R0 varies by three
#: orders of magnitude across the span here, so any single absolute number
#: would be both too loose somewhere and unreachable somewhere else.
#:
#: 1e-9 is nine orders looser than what is actually achieved (measured
#: |R|/R0 ~ 1e-15 at every station of every operating point tested) and six
#: orders tighter than what a mis-solve produces: returning a pole gives
#: |R| ~ 1e+15. The gap either side is what makes this a real check rather
#: than a threshold that has to be tuned.
RESIDUAL_RTOL = 1e-9

#: Floor for the degenerate case where R0 itself is ~0, so the relative test
#: has nothing to scale against. Not the convergence criterion -- a guard on it.
RESIDUAL_ATOL = 1e-14


def _sin(x):
    """`math.sin`, or `cmath.sin` if `x` is complex-typed (complex-step)."""

    if isinstance(x, complex):
        return cmath.sin(x)
    return math.sin(x)


def _cos(x):
    """`math.cos`, or `cmath.cos` if `x` is complex-typed (complex-step)."""

    if isinstance(x, complex):
        return cmath.cos(x)
    return math.cos(x)


def _blade_element(phi, station: StationParams):
    """
    Blade-element force coefficients at one trial phi, plus the loss factor.

    Returns (cl, cd, cn, ct, F, sin_phi, cos_phi). Untouched by F, which
    enters only on the momentum side -- see corrections.py on why folding it
    into cn/ct instead would double-count against the Buhl correction.
    """

    alpha = phi - station.twist
    cl = station.airfoil.cl(alpha)
    cd = station.airfoil.cd(alpha)

    # Complex-safe (Phase 3, B0): a complex-step perturbation to phi must
    # pass through here with its imaginary part intact. `math.sin` rejects
    # complex input; the helpers dispatch on type and leave the real path
    # bit-identical (`_sin(float)` is `math.sin(float)`).
    sin_phi = _sin(phi)
    cos_phi = _cos(phi)

    cn = cl * cos_phi + cd * sin_phi
    ct = cl * sin_phi - cd * cos_phi
    F = combined_loss_factor(station.r, station.R, station.n_blades, phi, station.r_hub)

    return cl, cd, cn, ct, F, sin_phi, cos_phi


def _blade_element_and_induction(phi, station: StationParams):
    """
    (a, a', Cl, Cd, Cn, Ct, F) at one trial phi.

    Kept with this name and return shape because it is the decomposition the
    adjoint wants and because `tests/test_bem_stages.py` exercises the
    Stage 3 turbulent-wake path through it. `a'` is reported here but is deliberately *not* what `residual`
    divides by -- see that function on why.
    """

    cl, cd, cn, ct, F, sin_phi, cos_phi = _blade_element(phi, station)

    sigma = station.solidity
    Y = sigma * cn / (sin_phi * sin_phi)
    a = corrected_axial_induction(Y, F)
    a_prime = 1.0 / ((4.0 * F * sin_phi * cos_phi) / (sigma * ct) - 1.0)

    return a, a_prime, cl, cd, cn, ct, F


def residual(phi, station: StationParams):
    """
    Ning (2014)'s single-residual BEM equation in the inflow angle phi.

        R(phi) = sin(phi) / (1 - a)  -  [cos(phi) - sigma*Ct / (4 F sin phi)] / lambda_r

    Both forms, and why this one
    -----------------------------
    The kinematic relation is tan(phi) = (1-a) / ((1+a') lambda_r). It can be
    cleared of denominators two ways, and they are not equally well behaved:

      multiplied   sin(phi) (1+a') lambda_r - cos(phi) (1-a)     [used up to Task 5]
      divided      sin(phi)/(1-a) - cos(phi)/(lambda_r (1+a'))   [Ning's, used here]

    They share zeros -- one is the other times (1-a)(1+a') lambda_r -- but they
    are different functions, and Ning's convergence guarantees are properties
    of his. The multiplied form was chosen in Phase 0 to avoid a division; the
    cost of that trade is measurable. Scanning phi in (0, pi/2) at 4000 points
    for every Phase VI station:

        multiplied form   3 sign changes   (the root, plus both poles)
        divided form      2 sign changes   (the root, plus one)

    The pole that disappears is a'. In the multiplied form a' -> infinity
    wherever Ct(phi) = 0, and it crosses the residual's sign exactly as a true
    root does. In the divided form (1+a') sits in a denominator, so the same
    point is a *zero* of that term, not a pole. Better still, the division can
    be removed algebraically rather than merely surviving numerically: with
    kappa' = sigma*Ct / (4 F sin phi cos phi) the momentum relation gives
    a' = kappa'/(1 - kappa'), hence 1 + a' = 1/(1 - kappa'), hence

        cos(phi) / (lambda_r (1 + a')) = cos(phi) (1 - kappa') / lambda_r
                                       = [cos(phi) - sigma*Ct/(4 F sin phi)] / lambda_r

    which is what is coded above. There is no division by (1 + a') anywhere in
    this function, so a' having a pole is not something this residual can
    notice. It also makes R finite at phi = pi/2 exactly, where cos(phi) = 0
    and kappa' is infinite -- the two conspire to leave -sigma*Ct/(4 F lambda_r),
    and that is what lets the bracket close on pi/2 itself instead of on
    pi/2 minus a fudge factor.

    What this form does NOT fix
    ----------------------------
    The second sign change survives, and it is not a pole. It sits where
    Cn < 0 -- negative axial loading, where momentum theory has no windmill
    solution and `corrections.corrected_axial_induction` is running on the
    gamma2 < 0 branch it documents as non-physical. Excluding it is the job of
    the region classification in `momentum_region_bracket`, not of the
    residual form. Adopting Ning's residual is necessary for a guaranteed
    bracket; it is not sufficient.

    Parameters
    ----------
    phi : float
        Inflow angle, radians, in (0, pi/2].
    station : StationParams

    Returns
    -------
    float
        Residual; the solved phi is where this is zero.
    """

    _cl, _cd, cn, ct, F, sin_phi, cos_phi = _blade_element(phi, station)

    sigma = station.solidity
    a = corrected_axial_induction(sigma * cn / (sin_phi * sin_phi), F)

    axial = sin_phi / (1.0 - a)
    tangential = (cos_phi - sigma * ct / (4.0 * F * sin_phi)) / station.tsr

    return axial - tangential


def _induction_factors(phi, station: StationParams):
    """(a, a', Cl, Cd, Cn, Ct, F) at the solved phi (no residual)."""

    return _blade_element_and_induction(phi, station)


def _normal_force_coefficient(phi, station: StationParams):
    """Cn at one trial phi -- the quantity the region boundary is defined by."""

    return _blade_element(phi, station)[2]


def momentum_region_bracket(station: StationParams, eps=PHI_EPS):
    """
    The analytically guaranteed bracket for the momentum region, or None.

    Replaces the Phase 0 `_select_bracket`, which sampled the residual at 2000
    trial phi and took "the sign change nearest atan(1/lambda_r)". That was a
    heuristic standing in for Ning's region analysis, and it had two failure
    modes its own docstring admitted: it could select a pole instead of the
    root (both cross the residual's sign), and it raised ValueError when it
    found nothing, aborting a whole sweep over one station.

    The classification
    -------------------
    The momentum region is the part of phi in (0, pi/2] where the momentum
    relation actually has a windmill solution, which is where the axial
    loading is non-negative: Cn >= 0. Below that boundary the residual has a
    second sign change with no physical meaning (see `residual`), and it is
    the only other one there is.

    Two facts make the endpoints analytic rather than empirical:

      * **Cn(pi/2) = Cd > 0, always.** At phi = pi/2, Cn = Cl*cos(phi) +
        Cd*sin(phi) = Cd, and drag is positive for every airfoil at every
        angle. So if Cn(eps) < 0 the boundary is bracketed by (eps, pi/2) and
        bisection on it is always well posed -- there is no search that can
        fail to find a starting point.

      * **R(pi/2) > 0 whenever a < 1.** At phi = pi/2 the residual reduces to
        1/(1 - a) + sigma*Ct/(4 F lambda_r), and a < 1 holds in both branches
        of `corrected_axial_induction`: the momentum branch is Y/(4F + Y) < 1
        for Y >= 0, and the Buhl branch approaches 1 only asymptotically,
        since a = 1 would require Ct(1) = 0 where Buhl's quadratic gives 2.

    The lower endpoint's sign, R(phi_Cn0+) < 0, is the one that is argued
    rather than proven. At Cn = 0 the residual is sin(phi) - cos(phi)(1 -
    kappa')/lambda_r with kappa' < 0 strictly (Cn = 0 forces Cl = -Cd tan phi,
    hence Ct = -Cd/cos phi < 0), so it is negative whenever
    tan(phi) lambda_r < 1 - kappa', and the margin grows as phi_Cn0 -> 0.
    Verified over 247 stations spanning v = 5-20 m/s, lambda = 2-7.5, chord
    scaled by 0.7 and 1.3, and twist offset by +/-6 deg: exactly one sign
    change in the region, always with R(lo) < 0 < R(hi), always the physical
    root. It is nonetheless *checked* at run time rather than assumed -- if
    either endpoint has the wrong sign this returns None and the caller
    reports non-convergence.

    Provenance
    -----------
    Derived here from the residual's own structure and verified numerically as
    described. It is NOT transcribed from Ning (2014), which could not be
    obtained -- see the module docstring's PROVENANCE note. Ning also defines
    a propeller-brake region on phi in (-pi/4, 0); it is deliberately not
    implemented, because Cn(pi/2) > 0 guarantees the momentum region always has
    a valid upper endpoint and no case tested has ever failed to bracket there.
    A station that does reach that state is reported as non-convergent naming
    the suspected cause, rather than silently solved by an untested branch.

    Returns
    -------
    tuple or None
        (phi_lo, phi_hi, R_lo, R_hi) with R_lo < 0 < R_hi, or None if the
        momentum region does not bracket a root.
    """

    if _normal_force_coefficient(eps, station) >= 0.0:
        phi_lo = eps
    else:
        # Bisect the Cn sign change. Always well posed: Cn(eps) < 0 here and
        # Cn(pi/2) = Cd > 0 above.
        lo, hi = eps, PHI_MAX
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _normal_force_coefficient(mid, station) < 0.0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-12:
                break
        phi_lo = hi

    if phi_lo >= PHI_MAX:
        return None

    r_lo = residual(phi_lo, station)
    r_hi = residual(PHI_MAX, station)

    if not (r_lo < 0.0 < r_hi):
        return None

    return phi_lo, PHI_MAX, r_lo, r_hi


def solve_station(station: StationParams, *, record_history=False):
    """
    Solve for the inflow angle phi at one station, with a checked residual.

    Never raises on a solver failure. A station that cannot be solved comes
    back with ``converged`` False and a ``failure`` string; an optimiser line
    search or a several-hundred-point smoothness sweep gets a reportable point
    instead of an aborted run, which is what the Phase 0 `ValueError` cost.

    Parameters
    ----------
    station : StationParams
    record_history : bool
        Record every (phi, residual) the root-finder evaluates, under
        ``history``. Off by default: it costs a list append per residual
        evaluation and nothing in a normal sweep reads it. Used to generate
        the committed histories under `verification/phase_vi/`.

    Returns
    -------
    dict
        On success: phi, a, a_prime, Cl, Cd, Ct, Cq, F as before, plus
        ``converged`` True, ``residual`` (the checked |R| at the returned
        phi), ``residual_initial`` (R0, what it is measured against),
        ``iterations``, ``function_calls``, ``failure`` None, and ``history``.

        On failure: the same keys, ``converged`` False, ``failure`` naming the
        cause, and phi/a/... carrying the best available values (or NaN where
        there is none) so a caller can still plot the point.
    """

    history = [] if record_history else None

    def f(phi):
        value = residual(phi, station)
        if history is not None:
            history.append((phi, value))
        return value

    bracket = momentum_region_bracket(station)
    if bracket is None:
        return _failed_result(
            station,
            "no sign change in the momentum region (0, pi/2]: R(lo) < 0 < R(hi) "
            "does not hold. The rotor may be in the propeller-brake state, "
            "which this solver detects but does not solve -- see "
            "momentum_region_bracket.",
            history,
        )

    phi_lo, phi_hi, r_lo, r_hi = bracket
    residual_initial = max(abs(r_lo), abs(r_hi))

    phi, info = brentq(f, phi_lo, phi_hi, xtol=1e-14, rtol=8.9e-16,
                       maxiter=200, full_output=True)

    final_residual = abs(residual(phi, station))
    tolerance = max(RESIDUAL_RTOL * residual_initial, RESIDUAL_ATOL)

    if not info.converged:
        return _failed_result(
            station,
            f"brentq did not converge in {info.iterations} iterations",
            history, phi=phi, residual=final_residual,
            residual_initial=residual_initial, iterations=info.iterations,
            function_calls=info.function_calls,
        )

    if final_residual > tolerance:
        return _failed_result(
            station,
            f"residual not reduced: |R| = {final_residual:.3e} at the returned "
            f"phi exceeds {RESIDUAL_RTOL:g} * R0 = {tolerance:.3e}. A returned "
            "pole shows up here, six orders clear of the tolerance.",
            history, phi=phi, residual=final_residual,
            residual_initial=residual_initial, iterations=info.iterations,
            function_calls=info.function_calls,
        )

    a, a_prime, cl, cd, _cn, ct, F = _induction_factors(phi, station)

    return {
        "phi": phi,
        "a": a,
        "a_prime": a_prime,
        "Cl": cl,
        "Cd": cd,
        "Ct": ct,
        "Cq": ct,
        "F": F,
        "converged": True,
        "residual": final_residual,
        "residual_initial": residual_initial,
        "iterations": info.iterations,
        "function_calls": info.function_calls,
        "failure": None,
        "history": history,
    }


def _failed_result(station, reason, history, phi=float("nan"),
                   residual=float("nan"), residual_initial=float("nan"),
                   iterations=0, function_calls=0):
    """
    A reported non-convergence, shaped exactly like a successful result.

    Same keys either way, so no caller has to branch on the shape of what it
    got back -- only on ``converged``. Per-station aerodynamic values are
    recomputed at whatever phi was reached when there is one, so a failed
    point is still plottable by the smoothness gate.
    """

    values = dict.fromkeys(("a", "a_prime", "Cl", "Cd", "Ct", "Cq", "F"), float("nan"))
    if phi == phi:  # not NaN
        try:
            a, a_prime, cl, cd, _cn, ct, F = _induction_factors(phi, station)
            values.update(a=a, a_prime=a_prime, Cl=cl, Cd=cd, Ct=ct, Cq=ct, F=F)
        except (ValueError, ZeroDivisionError, ArithmeticError):
            pass

    return {
        "phi": phi,
        **values,
        "converged": False,
        "residual": residual,
        "residual_initial": residual_initial,
        "iterations": iterations,
        "function_calls": function_calls,
        "failure": reason,
        "history": history,
    }
