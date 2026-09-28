"""
The analytic Schmitz optimum-rotor blade.

Schmitz's design gives the chord and twist that maximise power extraction at a
chosen tip-speed ratio, *including* wake rotation -- which is what separates it
from the simpler Betz construction that ignores swirl and consequently
under-chords the inboard span badly. For a rotor at lambda = 6.5 with a
substantial root cut-out, that difference is not academic.

    phi(r)   = (2/3) * arctan(1 / lambda_r),      lambda_r = lambda * r / R
    c(r)     = (16 * pi * r) / (B * C_L) * sin^2(phi / 2)
    theta(r) = phi(r) - alpha_design

`C_L` and `alpha_design` are not free parameters: they are fixed at
the airfoil's **maximum lift-to-drag point** at a representative Reynolds
number, read from the polar interpolant rather than quoted from a datasheet --
so the baseline is consistent with the same polar data the solver uses.

What this is for
-----------------
Two things, and it is worth keeping them apart:

  1. **The reference distribution for the representation study**.
     A parameterisation is judged by what it can represent, and the honest
     target is the shape the optimiser will actually be asked to start from
     and move away from -- not a generic taper chosen to be easy to fit.

  2. **The baseline `x0`**, whose evaluated performance is
     the reference numbers for the entire results chapter. That construction lives in `baseline.py`; this module is only
     the analytic shape.

Schmitz is a *design* formula, not a performance model. It says nothing about
what Cp the resulting blade achieves -- that comes from running it through the
BEM solver like any other geometry, which is exactly what step 1.7 does.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math

import numpy as np

from config import load_design_rotor


#: Representative Reynolds number for the design point ("at a representative
#: Re", which makes the choice a documented input rather than
#: a detail). The design rotor's computed envelope is roughly 49k-854k; 200k
#: sits in the band the blade spends most of its energy-producing hours in, and
#: is inside the range where the SG6043 cache has independent UIUC experimental
#: support (100k-500k). Defined here so the baseline and the representation
#: study cannot drift apart on it.
DEFAULT_DESIGN_REYNOLDS = 200_000.0


def schmitz_inflow_angle(r_over_R, tip_speed_ratio):
    """
    phi(r) = (2/3) * arctan(1 / lambda_r), radians.

    The 2/3 is Schmitz's, and it is the whole content of the wake-rotation
    correction: Betz's construction gives arctan(1 / lambda_r) times 1, i.e.
    the zero-swirl inflow angle, and the optimum with swirl sits two-thirds of
    the way from the rotor plane toward it.
    """

    r_over_R = np.asarray(r_over_R, dtype=float)
    if np.any(r_over_R <= 0.0):
        raise ValueError("r/R must be positive; the formula is singular at the axis")

    local_speed_ratio = tip_speed_ratio * r_over_R
    return (2.0 / 3.0) * np.arctan(1.0 / local_speed_ratio)


def schmitz_chord(r_over_R, tip_speed_ratio, radius_m, n_blades, design_cl):
    """
    c(r) = (16 pi r) / (B C_L) * sin^2(phi / 2), metres.

    Note the shape this produces inboard: as r -> 0 the inflow angle tends to
    90 deg and the chord to a large value. Schmitz blades are famously wide at
    the root, which is why the result has to be checked
    against the manufacturability bounds and *recorded* if it violates them,
    rather than assumed feasible.
    """

    if design_cl <= 0.0:
        raise ValueError(f"design_cl must be positive, got {design_cl}")

    r_over_R = np.asarray(r_over_R, dtype=float)
    phi = schmitz_inflow_angle(r_over_R, tip_speed_ratio)
    radius = radius_m * r_over_R

    return (16.0 * math.pi * radius) / (n_blades * design_cl) * np.sin(phi / 2.0) ** 2


def schmitz_twist(r_over_R, tip_speed_ratio, alpha_design_rad):
    """
    theta(r) = phi(r) - alpha_design, radians.

    This is absolute local twist in the convention `bem/station.py` uses
    (alpha = phi - twist), so it already includes whatever pitch the design
    angle of attack implies. No separate tip-pitch term is added.
    """

    return schmitz_inflow_angle(r_over_R, tip_speed_ratio) - alpha_design_rad


def max_lift_to_drag_point(interpolant, reynolds, alpha_range_deg=(-5.0, 15.0),
                           n_samples=2001):
    """
    The airfoil's maximum-L/D operating point at one Reynolds number.

    Returns `(alpha_rad, cl, lift_to_drag)`.

    Sampled on a fine grid rather than solved for, deliberately: L/D near its
    maximum is flat, so a root-find on d(L/D)/dalpha is poorly conditioned
    exactly where the answer matters, while a 0.01 deg grid resolves the peak
    far more tightly than the design decision needs. The search is confined to
    the attached-flow band -- past stall the ratio is meaningless as a design
    point even where it is numerically defined.

    Reynolds dependence is real and should be stated wherever this is used:
    for SG6043 the peak moves from L/D = 84.8 at alpha = 6.1 deg at Re = 150k
    to L/D = 132.4 at alpha = 3.9 deg at Re = 400k. Plan step 1.7 says "at a
    representative Re", which makes the choice of that Re a documented input,
    not a detail.
    """

    lo, hi = alpha_range_deg
    best = None
    for i in range(n_samples):
        alpha_deg = lo + (hi - lo) * i / (n_samples - 1)
        cl = interpolant.cl(alpha_deg, reynolds)
        cd = interpolant.cd(alpha_deg, reynolds)
        if cd <= 0.0 or cl <= 0.0:
            continue
        ratio = cl / cd
        if best is None or ratio > best[2]:
            best = (math.radians(alpha_deg), cl, ratio)

    if best is None:
        raise ValueError(
            f"no positive-lift, positive-drag point found in "
            f"alpha = {alpha_range_deg} deg at Re = {reynolds:,.0f}")

    return best


def schmitz_distribution(stations_r_over_R, design_cl, alpha_design_rad,
                         tip_speed_ratio=None, radius_m=None, n_blades=None):
    """
    Chord and twist arrays at the given normalised stations.

    Rotor properties default to `config/rotor_design.yaml`, so a caller that
    wants the project's actual baseline supplies only the airfoil operating
    point -- which is the part that comes from the polar data rather than from
    configuration.

    Returns
    -------
    (chord_m, twist_rad) : ndarray, ndarray
    """

    design = load_design_rotor()
    tip_speed_ratio = design.design_tsr if tip_speed_ratio is None else tip_speed_ratio
    radius_m = design.radius_m if radius_m is None else radius_m
    n_blades = design.n_blades if n_blades is None else n_blades

    chord = schmitz_chord(stations_r_over_R, tip_speed_ratio, radius_m,
                          n_blades, design_cl)
    twist = schmitz_twist(stations_r_over_R, tip_speed_ratio, alpha_design_rad)

    return chord, twist
