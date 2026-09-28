"""
Viterna post-stall extrapolation: extend a converged XFOIL polar to +/-180 deg.

Promoted (behaviour unchanged) out of `validation/export_qblade.py`, where it
was reachable only by running a QBlade export. Two things now need it that have
nothing to do with QBlade:

  * the polar caches themselves, which are written over the full circle so the
    (alpha, Re) grid has no ragged edge for a C1 interpolant to fit across
    (see `polars/interpolant.py`);
  * the pyBEMT/QBlade comparison tables, which already reached into
    `export_qblade` across the `validation` package boundary to get at it.

The classical Viterna block transfers unchanged -- A1 = B1/2, B1 = CD_MAX,
B2 = (Cd_s - CD_MAX sin^2 a_s)/cos a_s, A2 = (Cl_s - CD_MAX sin a_s cos a_s)
sin a_s / cos^2 a_s -- and `export_qblade` produces byte-identical `.plr`
output through this module (asserted by tests/test_polar_cache.py).

The three tuned constants
-------------------------
`CD_MAX`, the reversed-flow lift amplitude and the negative-side scale were
fitted to match QBlade rather than derived. Reviewed here before they can enter
a design-rotor cache:

* **CD_MAX.** Now a parameter, not a constant. `cd_max_finite_blade()` gives
  Viterna & Corrigan's finite-blade relation `1.11 + 0.018 AR`, a *derived*
  value for a blade of known aspect ratio; that is what a design-rotor cache
  should be built with. The module default stays 1.8 because the S809 cache,
  the committed QBlade `.plr` pair and the pyBEMT comparison tables were all
  produced with it -- see `CD_MAX_QBLADE_MATCHED`.

* **Reversed-flow amplitude (0.63).** Not independent after all: the flat-plate
  lift amplitude in the Viterna form is A1 = CD_MAX/2 = 0.9 at CD_MAX = 1.8,
  and the negative branch is scaled by 0.7, giving 0.7 x 0.9 = 0.63 exactly.
  It is therefore *derived* here as `negative_cl_scale * cd_max / 2` rather
  than carried as a fitted number: it reproduces 0.63 at the S809 defaults and
  stays self-consistent when CD_MAX changes for another blade.

* **Negative-side scale (0.7).** Retained as a documented convention, not
  derived. Past negative stall a cambered section's pressure side becomes the
  suction side and produces less lift than the mirrored positive branch; a
  0.7-0.8 reduction is the usual BEM-code treatment, and there is no measured
  S809 or SG6043 data at those angles to fit it against. It is a named
  uncertainty. It acts only at alpha < -alpha_s, outside any converged
  operating state of either rotor.

None of the three changes a value inside the XFOIL-converged band: the
extrapolation is anchored to the measured endpoint and reproduces it exactly
there (verified per Reynolds row by `validation/check_stitch_continuity.py`).

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math

import numpy as np

#: Post-stall drag maximum the S809 cache, the committed QBlade `.plr` files
#: and the pyBEMT comparison tables were all built with. Retained as the
#: default so those artefacts keep reproducing byte-for-byte; it is a fitted
#: number, not a derived one (see the module docstring and
#: `cd_max_finite_blade`).
CD_MAX_QBLADE_MATCHED = 1.8

#: Name the original `export_qblade` constant was known by.
CD_MAX = CD_MAX_QBLADE_MATCHED

#: Cd floor on the reversed-flow branch -- a reversed airfoil aligned with the
#: flow still has profile drag.
REVERSED_CD_FLOOR = 0.020

#: Scaling applied to the mirrored Viterna curve on the negative-alpha side.
#: A retained convention; see the module docstring.
NEGATIVE_CL_SCALE = 0.7

#: Aspect ratio beyond which the finite-blade relation is capped, and the value
#: it caps at: above AR ~ 50 the section behaves as 2D and CD_MAX saturates.
_AR_CAP = 50.0
_CD_MAX_2D = 2.01


def cd_max_finite_blade(aspect_ratio):
    """
    Viterna & Corrigan's finite-blade post-stall drag maximum, `1.11 + 0.018 AR`.

    The derived alternative to the fitted 1.8. Use
    it for a cache built for a blade whose aspect ratio is known: the design
    rotor's Schmitz baseline gives AR ~ 13 and hence CD_MAX ~ 1.35, well below
    the fitted value.

    Parameters
    ----------
    aspect_ratio : float
        Blade aspect ratio, R / mean chord. Capped at 50, above which the
        section behaves as 2D and CD_MAX saturates at 2.01.

    Returns
    -------
    float
    """

    if aspect_ratio <= 0.0:
        raise ValueError("aspect_ratio must be positive")
    if aspect_ratio >= _AR_CAP:
        return _CD_MAX_2D
    return 1.11 + 0.018 * aspect_ratio


def reversed_cl_amplitude(cd_max=CD_MAX_QBLADE_MATCHED,
                          negative_cl_scale=NEGATIVE_CL_SCALE):
    """
    Reversed-flow (|alpha| > 90 deg) lift amplitude, `negative_cl_scale * cd_max / 2`.

    The flat-plate amplitude in the Viterna form is A1 = CD_MAX/2; the reversed
    branch is the negative branch seen from behind, so it carries the same
    reduction. Returns 0.63 at the S809 defaults -- exactly the fitted number
    this replaces, so nothing built with it moves.
    """

    return negative_cl_scale * cd_max / 2.0


def viterna_coefficients(alpha_s_deg, cl_s, cd_s, cd_max=CD_MAX_QBLADE_MATCHED):
    """Classical Viterna A1/A2/B1/B2, matched to the airfoil state at alpha_s."""
    a = math.radians(alpha_s_deg)
    sin_a, cos_a = math.sin(a), math.cos(a)
    b1 = cd_max
    b2 = (cd_s - cd_max * sin_a ** 2) / cos_a
    a1 = b1 / 2.0
    a2 = (cl_s - cd_max * sin_a * cos_a) * sin_a / cos_a ** 2
    return a1, a2, b1, b2


def viterna(alpha_deg, coeffs):
    """Viterna Cl, Cd at alpha (deg), valid for 0 < alpha <= 90."""
    a1, a2, b1, b2 = coeffs
    r = math.radians(alpha_deg)
    sin_r, cos_r = math.sin(r), math.cos(r)
    cl = a1 * math.sin(2 * r) + a2 * cos_r ** 2 / sin_r
    cd = b1 * sin_r ** 2 + b2 * cos_r
    return cl, cd


def build_full_range_polar(alpha_xf, cl_xf, cd_xf, cm_xf, step=0.5,
                           cd_max=CD_MAX_QBLADE_MATCHED,
                           negative_cl_scale=NEGATIVE_CL_SCALE,
                           reversed_cd_floor=REVERSED_CD_FLOOR):
    """
    Extend a converged XFOIL polar to the full -180..180 deg range.

    The scheme, unchanged from the original QBlade export:

      alpha_n <= a <= alpha_s   the cached XFOIL data verbatim (alpha_n,
                                alpha_s are the most negative / most positive
                                converged alphas).
      alpha_s < a <= 90         classical Viterna, matched to the airfoil's
                                state at alpha_s.
      90 < a <= 180             flat-plate reversed-flow branch,
                                Cl = -A sin(2(180-a)) with A from
                                `reversed_cl_amplitude`, Cd mirrored from the
                                Viterna branch and floored.
      -alpha_s <= a < alpha_n   linear ramp in both Cl and Cd, from the XFOIL
                                endpoint at alpha_n onto the negative branch's
                                value at -alpha_s. (XFOIL rarely reaches
                                negative stall, so there is no measured anchor
                                to Viterna against on this side.)
      a < -alpha_s              negative branch: Cl = -negative_cl_scale *
                                Cl_viterna(|a|), Cd = Cd_viterna(|a|); beyond
                                -90 deg it is the negation of the positive
                                branch.

      Cm is held constant outside the converged band at whichever endpoint
      value is nearer -- a placeholder rather than a model.

    Parameters
    ----------
    alpha_xf, cl_xf, cd_xf, cm_xf : array_like
        The cached XFOIL polar, ascending in alpha. Pass the converged rows
        only: feeding an already-extended table back in would re-anchor the
        extrapolation at +/-180 deg, where cos(alpha_s) = -1.
    step : float, optional
        Output alpha resolution in degrees (default 0.5).
    cd_max : float, optional
        Post-stall drag maximum (default `CD_MAX_QBLADE_MATCHED` = 1.8; see
        `cd_max_finite_blade` for the derived alternative).
    negative_cl_scale : float, optional
        Reduction applied to the mirrored curve at negative alpha.
    reversed_cd_floor : float, optional
        Minimum Cd on the reversed-flow branch.

    Returns
    -------
    numpy.ndarray
        Columns [alpha, cl, cd, cm] over -180..180 deg.
    """

    alpha_xf = np.asarray(alpha_xf, dtype=float)
    cl_xf = np.asarray(cl_xf, dtype=float)
    cd_xf = np.asarray(cd_xf, dtype=float)
    cm_xf = np.asarray(cm_xf, dtype=float)

    alpha_s, cl_s, cd_s = alpha_xf[-1], cl_xf[-1], cd_xf[-1]
    alpha_n, cl_n, cd_n = alpha_xf[0], cl_xf[0], cd_xf[0]
    cm_hi, cm_lo = cm_xf[-1], cm_xf[0]

    coeffs = viterna_coefficients(alpha_s, cl_s, cd_s, cd_max=cd_max)
    amplitude = reversed_cl_amplitude(cd_max=cd_max,
                                      negative_cl_scale=negative_cl_scale)

    # Anchor for the negative-side ramp: the negative branch evaluated at
    # -alpha_s, which is where that branch takes over.
    cl_anchor_pos, cd_anchor = viterna(alpha_s, coeffs)
    cl_anchor = -negative_cl_scale * cl_anchor_pos

    n = int(round(360.0 / step)) + 1
    out = np.empty((n, 4))

    for i in range(n):
        a = -180.0 + i * step
        a = round(a, 6)

        if alpha_n - 1e-9 <= a <= alpha_s + 1e-9:
            # Inside the converged band: use the XFOIL data as-is. np.interp is
            # exact at the sample points and covers the odd non-converged alpha
            # the cache may be missing internally.
            cl = float(np.interp(a, alpha_xf, cl_xf))
            cd = float(np.interp(a, alpha_xf, cd_xf))
            cm = float(np.interp(a, alpha_xf, cm_xf))

        elif a > alpha_s:
            cm = cm_hi
            if a <= 90.0:
                cl, cd = viterna(a, coeffs)
            else:
                mirror = 180.0 - a
                cl = -amplitude * math.sin(2 * math.radians(mirror))
                if mirror <= 1e-9:
                    cd = reversed_cd_floor
                else:
                    cd = max(viterna(mirror, coeffs)[1], reversed_cd_floor)

        else:  # a < alpha_n
            cm = cm_lo
            if a >= -alpha_s:
                # Linear ramp from the XFOIL endpoint onto the negative branch.
                t = (a - alpha_n) / (-alpha_s - alpha_n)
                cl = cl_n + t * (cl_anchor - cl_n)
                cd = cd_n + t * (cd_anchor - cd_n)
            else:
                mag = -a
                if mag <= 90.0:
                    cl_p, cd = viterna(mag, coeffs)
                    cl = -negative_cl_scale * cl_p
                else:
                    mirror = 180.0 - mag
                    cl = amplitude * math.sin(2 * math.radians(mirror))
                    if mirror <= 1e-9:
                        cd = reversed_cd_floor
                    else:
                        cd = max(viterna(mirror, coeffs)[1], reversed_cd_floor)

        out[i] = (a, cl, cd, cm)

    return out
