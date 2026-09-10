"""
The Schmitz baseline `x0`, and the reference numbers measured from it
(plan step 1.7).

The plan is emphatic about what this is: the evaluated performance of this
blade gives "**the reference numbers for the entire results chapter**". Every
later claim about what the optimiser achieved is a comparison against these,
so they are constructed once, committed as a fixed artefact, and not
recomputed casually.

Construction, in order:

  1. The analytic Schmitz distribution (`schmitz.py`) at the SG6043 maximum-L/D
     point, sampled at the BEM strip stations.
  2. Projected onto the *same* spline parameterisation the optimiser uses,
     with the fitting error reported explicitly rather than absorbed. `x0` is
     the fitted control-point vector, not the analytic blade -- the optimiser
     can only start somewhere it can represent.
  3. Feasibility checked against the design-variable bounds, clipped and
     recorded if violated.
  4. Evaluated with the same solver and settings as everything else.

Step 3 is currently a no-op that reports itself as such: the bounds are still
`TODO` in `config/rotor_design.yaml` (see `docs/OUTSTANDING-INPUTS.md`), so
`feasibility` comes back as "not checked" with the reason attached rather than
as a quiet pass. The distinction matters -- a baseline that was never checked
must not read as a baseline that passed.

AEP is likewise absent rather than approximated: it needs the Weibull
parameters, which are `TODO`. `evaluate_baseline` returns everything that does
not depend on the wind resource -- Cp-lambda, spanwise loading, root bending
moment, peak thrust -- and names AEP as outstanding.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
from dataclasses import dataclass, field

import numpy as np

from bem.rotor import solve_rotor
from config import load_design_rotor, load_site
from config.unresolved import UnresolvedConfigError
from design.bounds import DesignBounds
from design.parameterisation import BladeParameterisation
from design.schmitz import (
    DEFAULT_DESIGN_REYNOLDS,
    max_lift_to_drag_point,
    schmitz_distribution,
)
from polars.polar import interpolant_for


@dataclass
class BaselineBlade:
    """The fitted Schmitz baseline and everything needed to reproduce it."""

    design_vector: np.ndarray
    parameterisation: BladeParameterisation
    chord_target_m: np.ndarray
    twist_target_rad: np.ndarray
    chord_rms_error_m: float
    twist_rms_error_rad: float
    design_reynolds: float
    alpha_design_rad: float
    design_cl: float
    max_lift_to_drag: float
    feasibility: dict = field(default_factory=dict)

    @property
    def chord_m(self):
        """Fitted chord at the stations -- what the solver actually sees."""

        return self.parameterisation.chord(self.design_vector)

    @property
    def twist_rad(self):
        return self.parameterisation.twist(self.design_vector)

    def to_geometry(self, polar_cache="sg6043"):
        return self.parameterisation.to_geometry(self.design_vector,
                                                 polar_cache=polar_cache)


def _check_feasibility(design_vector, parameterisation):
    """
    Feasibility against the design-variable bounds, or an explicit "not
    checked" when they are unresolved.

    Never returns a bare pass. Plan step 1.7 asks for the feasibility status to
    be *recorded*, and "the bounds do not exist yet" is a status -- one that
    must not be mistaken for "checked and satisfied" by anyone reading the
    committed artefact later.
    """

    try:
        bounds = DesignBounds.from_config(
            n_chord=parameterisation.n_chord, n_twist=parameterisation.n_twist)
    except UnresolvedConfigError as error:
        return {
            "checked": False,
            "reason": str(error),
            "clipped": False,
            "violations": [],
        }

    clipped, violations = bounds.clip_physical(design_vector)
    return {
        "checked": True,
        "reason": None,
        "clipped": bool(violations),
        "violations": [
            {"index": index, "value": value, "bound": bound}
            for index, value, bound in violations
        ],
        "clipped_vector": [float(value) for value in clipped],
    }


def build_schmitz_baseline(parameterisation=None,
                           design_reynolds=DEFAULT_DESIGN_REYNOLDS,
                           polar_cache="sg6043"):
    """
    Construct the baseline: analytic Schmitz, fitted, feasibility-checked.

    Parameters
    ----------
    parameterisation : BladeParameterisation or None
        Defaults to the configured 5 + 5 control points over 25 strips, i.e.
        the count the representation study selected.
    design_reynolds : float
        Where the airfoil's maximum-L/D point is read. A documented input, not
        a detail: the point moves from L/D = 84.8 at alpha = 6.1 deg at
        Re = 150k to L/D = 132.4 at alpha = 3.9 deg at Re = 400k.

    Returns
    -------
    BaselineBlade
    """

    parameterisation = parameterisation or BladeParameterisation()

    alpha_design, design_cl, lift_to_drag = max_lift_to_drag_point(
        interpolant_for(polar_cache), design_reynolds)

    r_over_R = parameterisation.radii / parameterisation.radius_m
    chord_target, twist_target = schmitz_distribution(
        r_over_R, design_cl, alpha_design)

    design_vector, chord_rms, twist_rms = parameterisation.fit(
        chord_target, twist_target)

    return BaselineBlade(
        design_vector=design_vector,
        parameterisation=parameterisation,
        chord_target_m=chord_target,
        twist_target_rad=twist_target,
        chord_rms_error_m=chord_rms,
        twist_rms_error_rad=twist_rms,
        design_reynolds=design_reynolds,
        alpha_design_rad=alpha_design,
        design_cl=design_cl,
        max_lift_to_drag=lift_to_drag,
        feasibility=_check_feasibility(design_vector, parameterisation),
    )


def root_bending_moment(stations, chords, air_density, n_blades, r_hub):
    """
    Flapwise root bending moment for one blade, N.m.

    M = integral of dT/dr * (r - r_hub) dr, with dT/dr the same normal-force
    integrand `solve_rotor` uses for thrust. Per *blade*, not per rotor: the
    root attachment carries one blade's load, and quoting a rotor-summed figure
    here would overstate it by a factor of B.
    """

    radii = [station["r"] for station in stations]
    lever_arms = [r - r_hub for r in radii]
    load = [
        0.5 * air_density * station["w"] ** 2 * chord * station["Cn"] * arm
        for station, chord, arm in zip(stations, chords, lever_arms)
    ]
    return sum(0.5 * (load[i] + load[i + 1]) * (radii[i + 1] - radii[i])
               for i in range(len(radii) - 1))


def evaluate_baseline(baseline, wind_speeds=None, tsr_values=None,
                      polar_cache="sg6043"):
    """
    The reference numbers, everything that does not need the wind resource.

    The operating line is fixed tip-speed ratio from cut-in to rated, which is
    the design strategy below rated (plan step 1.5). Above rated the strategy
    involves power limiting, which lands with the AEP model; the sweep
    therefore stops at rated and says so rather than extrapolating a strategy
    that has not been specified.

    Returns a dict with `cp_lambda`, `operating_line`, `spanwise` at the design
    point, `root_bending_moment_nm`, `peak_thrust_n`, and `outstanding` naming
    what could not be computed.
    """

    design = load_design_rotor()
    site = load_site()
    air = {"air_density": site.air_density,
           "kinematic_viscosity": site.kinematic_viscosity}

    geometry = baseline.to_geometry(polar_cache=polar_cache)
    chords = geometry.chord

    if wind_speeds is None:
        wind_speeds = [design.cut_in_wind_speed_ms + i * 0.5 for i in range(
            int((design.rated_wind_speed_ms - design.cut_in_wind_speed_ms) / 0.5) + 1)]
    if tsr_values is None:
        tsr_values = [round(2.0 + 0.25 * i, 2) for i in range(33)]  # 2.0 .. 10.0

    # Operating line at the design tip-speed ratio.
    operating_line = []
    peak_thrust = 0.0
    peak_thrust_speed = None
    for v_inf in wind_speeds:
        result = solve_rotor(geometry, tsr=design.design_tsr, v_inf=v_inf, **air)
        area = math.pi * geometry.R ** 2
        thrust = result["Ct"] * 0.5 * site.air_density * v_inf ** 2 * area
        power = result["Cp"] * 0.5 * site.air_density * v_inf ** 3 * area
        if thrust > peak_thrust:
            peak_thrust, peak_thrust_speed = thrust, v_inf
        operating_line.append({
            "v_inf": v_inf,
            "Cp": result["Cp"],
            "Ct": result["Ct"],
            "power_w": power,
            "thrust_n": thrust,
            "converged": result["converged"],
            "failed_stations": result["failed_stations"],
        })

    # Cp-lambda at the rated wind speed.
    cp_lambda = []
    for tsr in tsr_values:
        result = solve_rotor(geometry, tsr=tsr,
                             v_inf=design.rated_wind_speed_ms, **air)
        cp_lambda.append({
            "tsr": tsr,
            "Cp": result["Cp"],
            "Ct": result["Ct"],
            "converged": result["converged"],
        })

    # Spanwise detail and the root bending moment at the design point.
    design_point = solve_rotor(geometry, tsr=design.design_tsr,
                               v_inf=design.rated_wind_speed_ms, **air)
    stations = design_point["stations"]
    moment = root_bending_moment(stations, chords, site.air_density,
                                 geometry.n_blades, geometry.r_hub)

    return {
        "operating_line": operating_line,
        "cp_lambda": cp_lambda,
        "design_point": {
            "v_inf": design.rated_wind_speed_ms,
            "tsr": design.design_tsr,
            "Cp": design_point["Cp"],
            "Ct": design_point["Ct"],
            "converged": design_point["converged"],
        },
        "spanwise": [
            {
                "r": station["r"],
                "r_over_R": station["r"] / geometry.R,
                "chord_m": chord,
                "twist_deg": math.degrees(twist),
                "alpha_deg": math.degrees(station["alpha"]),
                "a": station["a"],
                "a_prime": station["a_prime"],
                "Cl": station["Cl"],
                "Cd": station["Cd"],
                "reynolds": station["reynolds"],
                "converged": station["converged"],
            }
            for station, chord, twist in zip(stations, chords, geometry.twist)
        ],
        "root_bending_moment_nm": moment,
        "peak_thrust_n": peak_thrust,
        "peak_thrust_wind_speed_ms": peak_thrust_speed,
        "outstanding": {
            "aep_mwh_per_year": (
                "BLOCKED: needs site.weibull_k and site.weibull_c_ms, still "
                "TODO. See docs/OUTSTANDING-INPUTS.md section 1."
            ),
            "above_rated_operating_line": (
                "Not swept: power limiting above the rated wind speed is part "
                "of the AEP model (plan step 1.5) and is not yet specified."
            ),
        },
    }
