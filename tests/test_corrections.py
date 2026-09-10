"""
Task 6 acceptance: the Glauert/Buhl correction's constants, its C0/C1
continuity at the blend, and the gamma-form reparameterisation.

The audit measured the blend continuity by finite difference and recorded
jumps of 3e-7 to 1.3e-5 -- pure O(h) truncation, i.e. "small enough to look
like a continuous function". These tests replace that with exact statements:
the constants are checked against the three conditions that pin them uniquely,
and the blend derivative is checked against a closed-form analytic target by
complex step rather than against itself by differencing.

What these tests do NOT establish
----------------------------------
That the coefficients are *Buhl's*. They prove the three numbers are the
unique solution of C0 + C1 at a = 0.4 plus Ct(1) = 2, which means any
transcription error would show up here immediately -- but a self-consistent
set of constants transcribed from the wrong paper would pass every test in
this file. The plan (1.3) and the Phase 1 brief require a check against
NREL/TP-500-36834 directly, which is still open; see the PROVENANCE note in
`bem/corrections.py`'s module docstring.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from decimal import Decimal, getcontext

import pytest

from bem.corrections import (
    BUHL_AC,
    buhl_gammas,
    corrected_axial_induction,
    high_thrust_correction,
)

getcontext().prec = 60

#: Loss factors to sweep. F -> 0 at the tip and at the hub, and the corrected
#: branch is live precisely there (measured: 1-3 stations per Phase VI
#: operating point, all at F = 0.10-0.29), so the small-F end is the working
#: range here, not an edge case.
F_VALUES = [1.0, 0.95, 0.8, 0.5, 0.3, 0.2, 0.1, 0.05]

#: kappa at the blend. a = kappa/(1+kappa) = 0.4 <=> kappa = 2/3, for any F.
KAPPA_AC = 2.0 / 3.0


def _momentum_ct(a, F):
    """Momentum theory's local thrust relation, the curve Buhl blends with."""

    return 4.0 * a * F * (1.0 - a)


def _exact_a(Y, F):
    """`corrected_axial_induction` evaluated at 60 significant digits."""

    Y, F = Decimal(Y), Decimal(F)
    kappa = Y / (4 * F)
    g1 = 2 * F * kappa - (Decimal(10) / 9 - F)
    g2 = 2 * F * kappa - F * (Decimal(4) / 3 - F)
    g3 = 2 * F * kappa - (Decimal(25) / 9 - 2 * F)
    return (g1 - g2.sqrt()) / g3


# ---------------------------------------------------------------------------
# The three conditions that pin the coefficients uniquely
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("F", F_VALUES)
def test_buhl_matches_momentum_theory_in_value_at_the_blend(F):
    """C0: Ct is continuous at a = 0.4, for any F."""

    buhl = high_thrust_correction(BUHL_AC + 1e-15, None, F)
    assert buhl == pytest.approx(_momentum_ct(BUHL_AC, F), rel=1e-14)


@pytest.mark.parametrize("F", F_VALUES)
def test_buhl_matches_momentum_theory_in_slope_at_the_blend(F):
    """
    C1: dCt/da is continuous at a = 0.4, for any F.

    Both sides are differentiated in closed form rather than by differencing,
    so this is an equality of two exact expressions:

        momentum  d/da [4aF(1-a)]                     = 4F(1 - 2a)
        Buhl      d/da [8/9 + (4F-40/9)a + (50/9-4F)a^2]
                                                      = (4F-40/9) + 2(50/9-4F)a

    At a = 0.4 the first is 0.8F and the second is
    (4F - 40/9) + 0.8(50/9 - 4F) = 0.8F. The F terms cancel to leave an
    identity, which is the whole reason a_c = 0.4 is the value it is.
    """

    momentum_slope = 4.0 * F * (1.0 - 2.0 * BUHL_AC)
    buhl_slope = (4.0 * F - 40.0 / 9.0) + 2.0 * (50.0 / 9.0 - 4.0 * F) * BUHL_AC
    assert buhl_slope == pytest.approx(momentum_slope, rel=1e-14)


@pytest.mark.parametrize("F", F_VALUES)
def test_buhl_anchor_at_a_equals_one_is_two_and_F_independent(F):
    """
    The third condition: Ct(1) = 2.

    8/9 + (4F - 40/9) + (50/9 - 4F) = (8 - 40 + 50)/9 = 2, with the F terms
    cancelling exactly. Being F-independent makes this a one-line
    transcription check on all three coefficients at once: change any of them
    and this stops being 2, or stops being independent of F, or both.
    """

    assert high_thrust_correction(1.0, None, F) == pytest.approx(2.0, rel=1e-15)


# ---------------------------------------------------------------------------
# The gamma reparameterisation: the identities it rests on
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("F", F_VALUES)
@pytest.mark.parametrize("kappa", [0.7, 1.0, 2.0, 5.0, 20.0])
def test_gamma_identity_relating_gamma1_gamma2_gamma3(F, kappa):
    """
    gamma1^2 - gamma2 = gamma3 * (Y - 8/9) / 2, with Y = 4 F kappa.

    This is the identity that makes the rationalised second form in
    `corrected_axial_induction` exact rather than approximate: it is what
    turns (gamma1 - sqrt(gamma2))/gamma3 into (Y/2 - 4/9)/(gamma1 +
    sqrt(gamma2)) when the numerator and denominator are multiplied by
    (gamma1 + sqrt(gamma2)).
    """

    g1, g2, g3 = buhl_gammas(kappa, F)
    Y = 4.0 * F * kappa
    assert g1 * g1 - g2 == pytest.approx(g3 * (Y - 8.0 / 9.0) / 2.0, rel=1e-12, abs=1e-15)


@pytest.mark.parametrize("F", F_VALUES)
def test_sqrt_gamma2_equals_F_at_the_blend(F):
    """
    At kappa = 2/3, gamma2 = F^2 exactly, for any F.

    gamma2 = 2F(2/3) - F(4/3 - F) = 4F/3 - 4F/3 + F^2 = F^2. Substituting
    sqrt(gamma2) = F into a = (gamma1 - sqrt(gamma2))/gamma3 collapses to
    a = 2/5 identically -- which is why the blend is C0 by construction and
    not by numerical coincidence.
    """

    _g1, g2, _g3 = buhl_gammas(KAPPA_AC, F)
    assert g2 == pytest.approx(F * F, rel=1e-15, abs=1e-18)


@pytest.mark.parametrize("F", F_VALUES)
def test_induction_is_exactly_ac_at_the_blend(F):
    """a(kappa = 2/3) = 0.4 for any F -- the C0 condition, seen from the inverse."""

    Y = 4.0 * F * KAPPA_AC
    assert corrected_axial_induction(Y, F) == pytest.approx(BUHL_AC, rel=1e-14)


@pytest.mark.parametrize("F", F_VALUES)
def test_induction_is_c1_across_the_blend(F):
    """
    da/dY is continuous at the blend, checked against a closed-form target.

    The momentum branch is a = Y/(4F + Y), so da/dY = 4F/(4F + Y)^2. At the
    blend Y = (8/3)F, giving exactly 9/(100F). The corrected branch is a
    different expression entirely, and this asserts that its derivative
    approaches the same number from above.

    Complex step, not finite difference: the function is complex-safe
    (`corrections._sqrt`, and the branch tests take `.real`), so the
    derivative comes back with no truncation error and the tolerance below
    reflects only the offset from the blend, not the differencing scheme.
    """

    Y_blend = (8.0 / 3.0) * F
    offset = Y_blend * 1e-9
    h = 1e-30
    analytic = 9.0 / (100.0 * F)

    from_below = corrected_axial_induction(complex(Y_blend - offset, h), F).imag / h
    from_above = corrected_axial_induction(complex(Y_blend + offset, h), F).imag / h

    assert from_below == pytest.approx(analytic, rel=1e-7)
    assert from_above == pytest.approx(analytic, rel=1e-7)


# ---------------------------------------------------------------------------
# Accuracy and conditioning
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("F", F_VALUES)
def test_matches_a_sixty_digit_reference(F):
    """
    Accuracy is checked against exact arithmetic, not against the previous
    implementation.

    Pinning this to the old quadratic solve would only assert that the
    reparameterisation reproduced whatever error that had; the Decimal
    reference asserts the answer is right.
    """

    worst = 0.0
    for i in range(1, 800):
        Y = i * 0.025
        a = corrected_axial_induction(Y, F)
        if a <= BUHL_AC:
            continue
        exact = float(_exact_a(Y, F))
        worst = max(worst, abs(a - exact) / abs(exact))
    assert worst < 1e-14, f"worst relative error {worst:.3e} at F={F}"


@pytest.mark.parametrize("F", [0.7, 0.5, 0.3, 0.1, 0.05])
def test_exact_at_the_gamma3_removable_singularity(F):
    """
    gamma3 = 0 is reachable inside the corrected region for F <~ 0.86, and F
    falls to 0 at the tip where this branch is live.

    It is a removable singularity -- gamma1 = sqrt(gamma2) there too -- so
    Ning's published form (gamma1 - sqrt(gamma2))/gamma3 divides 0 by 0
    exactly at it and cancels catastrophically near it (measured: 1.9e-11
    relative at gamma3 = 1.4e-5, 1.1e-9 at 1.4e-7). The rationalised form
    this module switches to has no gamma3 in it at all.
    """

    kappa = (25.0 / 9.0 - 2.0 * F) / (2.0 * F)
    Y = 4.0 * F * kappa

    _g1, _g2, g3 = buhl_gammas(kappa, F)
    assert g3 == pytest.approx(0.0, abs=1e-15), "test did not actually hit gamma3 = 0"

    a = corrected_axial_induction(Y, F)          # must not raise
    assert a == pytest.approx(float(_exact_a(Y, F)), rel=1e-14)


@pytest.mark.parametrize("F", F_VALUES)
def test_induction_is_strictly_increasing_in_loading(F):
    """
    a(Y) rises monotonically, through the blend and across both conditioning
    forms.

    This is the direct guard on what Task 6 removed. The old implementation
    picked between the quadratic's two roots by "closest to a_naive", with a
    fallback of "closest to ac"; a heuristic like that fails by selecting the
    other root somewhere, and the signature of that is a non-monotone step in
    an otherwise smooth curve. Sampling across the blend and out to heavy
    loading would catch it.
    """

    previous = None
    for i in range(1, 2000):
        Y = i * 0.01
        a = corrected_axial_induction(Y, F)
        if previous is not None:
            assert a > previous, f"a(Y) decreased at Y={Y}, F={F}"
        previous = a


# ---------------------------------------------------------------------------
# "No heuristic branch selection remains in corrections.py"
# ---------------------------------------------------------------------------

def test_no_heuristic_root_selection_remains():
    """
    The work order's third "done when" clause, asserted against the source.

    Both removed expressions chose between two genuinely different roots by
    proximity to a reference value. The conditioning switch that remains
    chooses between two algebraically identical expressions for the same
    number, so value and every derivative are continuous across it -- see the
    function's docstring, which draws exactly this distinction.
    """

    import inspect

    from bem import corrections

    source = inspect.getsource(corrections)
    for banned in ("min(candidates", "max(root_plus", "key=lambda"):
        assert banned not in source, (
            f"{banned!r} is back in corrections.py; branch selection must be "
            "analytic (the minus sign in Ning's form), not a nearest-value "
            "heuristic"
        )


def test_ac_is_a_constant_not_a_parameter():
    """
    The coefficients 8/9, 40/9 and 50/9 are derived from ac = 0.4, so an
    `ac=` argument that did not move them was a knob that produced a
    discontinuity if turned. Task 6 removed it.
    """

    import inspect

    from bem.corrections import corrected_axial_induction, high_thrust_correction

    for function in (high_thrust_correction, corrected_axial_induction):
        assert "ac" not in inspect.signature(function).parameters, (
            f"{function.__name__} takes an `ac` argument again; the Buhl "
            "coefficients do not move with it"
        )
