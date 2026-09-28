"""
Power curve and Cp-lambda sweep generation on top of the rotor solver.

`rotor.solve_rotor` answers a single operating point and reports only the
dimensionless Ct/Cp. Everything downstream -- the Cp-lambda validation
curve and the AEP integration -- needs the same solver
swept over a range of operating points, in dimensional units (watts, newtons,
newton-metres). That sweep is what this module adds; it contains no new
aerodynamics of its own and deliberately calls straight through to
`solve_rotor` so there is exactly one BEM implementation in the project.

Two machine types, because the sweep convention differs and mixing them up is
the single easiest way to produce a meaningless comparison (exactly that
mistake -- wind speed vs TSR -- once faked a 40-degree AoA disagreement
against QBlade):

  fixed-speed    (`rpm=`)  rotor speed is held constant and TSR therefore
                           falls as wind speed rises. This is the NREL Phase
                           VI convention (71.63 RPM synchronous) and how most
                           small stall-regulated turbines actually run.
  variable-speed (`tsr=`)  TSR is held at its design value and rotor speed
                           tracks the wind. Cp is then constant with wind
                           speed apart from Reynolds-number drift, which is a
                           useful check on the polar cache (see
                           tests/test_powercurve.py, check 3).

Power here is aerodynamic rotor power: no generator/gearbox efficiency, no
rated-power cap or pitch regulation above rated. Those belong with the AEP
objective (`objective/power.py`), not in the aerodynamic solver.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math

from bem.rotor import RotorGeometry, solve_rotor

#: Betz limit, 16/27 -- the theoretical maximum Cp for any actuator disc.
#: Used as a sanity bound, never as a correction.
BETZ_LIMIT = 16.0 / 27.0


def swept_area(geometry: RotorGeometry):
    """Rotor swept area, m^2."""

    return math.pi * geometry.R ** 2


def rpm_to_omega(rpm):
    """Rotational speed, RPM -> rad/s."""

    return rpm * 2.0 * math.pi / 60.0


def omega_to_rpm(omega):
    """Rotational speed, rad/s -> RPM."""

    return omega * 60.0 / (2.0 * math.pi)


def operating_point(geometry: RotorGeometry, v_inf, tsr, *, air_density,
                    kinematic_viscosity, airfoil_for_reynolds=None,
                    airfoils=None, keep_stations=False):
    """
    One operating point, with `solve_rotor`'s dimensionless result converted
    to dimensional power/thrust/torque.

    Parameters
    ----------
    geometry : RotorGeometry
    v_inf : float
        Freestream wind speed, m/s.
    tsr : float
        Rotor tip-speed ratio, Omega * R / v_inf.
    air_density : float
        kg/m^3. Required and keyword-only -- see `solve_rotor`'s docstring for
        why there is no default. Supplied from `config/`.
    kinematic_viscosity : float
        m^2/s. Likewise required; passed through to `solve_rotor` for the
        per-station Reynolds estimate.
    airfoil_for_reynolds : callable or None
        Passed straight through to `solve_rotor` -- see its docstring. The
        normal path leaves this None; the blade names its own airfoil.
    airfoils : list or None
        Passed straight through to `solve_rotor` -- see its docstring.
    keep_stations : bool
        Include the per-station dicts in the result. Off by default: a sweep
        of a few hundred points holding every station would be large and is
        rarely what the caller wants.

    Returns
    -------
    dict
        v_inf, tsr, omega (rad/s), rpm, Cp, Ct, power (W), thrust (N),
        torque (N m), and "stations" if `keep_stations`.
    """

    if v_inf <= 0.0:
        raise ValueError("v_inf must be positive")
    if tsr <= 0.0:
        raise ValueError("tsr must be positive")

    result = solve_rotor(geometry, tsr=tsr, v_inf=v_inf,
                         air_density=air_density,
                         kinematic_viscosity=kinematic_viscosity,
                         airfoil_for_reynolds=airfoil_for_reynolds,
                         airfoils=airfoils)

    area = swept_area(geometry)
    omega = tsr * v_inf / geometry.R

    # Cp/Ct are defined by solve_rotor against exactly these references, so
    # inverting them here is an identity, not an independent estimate --
    # tests/test_powercurve.py's check 1 asserts the round-trip holds.
    power = result["Cp"] * 0.5 * air_density * v_inf ** 3 * area
    thrust = result["Ct"] * 0.5 * air_density * v_inf ** 2 * area
    torque = power / omega

    point = {
        "v_inf": v_inf,
        "tsr": tsr,
        "omega": omega,
        "rpm": omega_to_rpm(omega),
        "Cp": result["Cp"],
        "Ct": result["Ct"],
        "power": power,
        "thrust": thrust,
        "torque": torque,
    }
    if keep_stations:
        point["stations"] = result["stations"]
    return point


def power_curve(geometry: RotorGeometry, wind_speeds, rpm=None, tsr=None, *,
                air_density, kinematic_viscosity,
                airfoil_for_reynolds=None, airfoils=None,
                keep_stations=False):
    """
    Sweep the solver across a range of wind speeds.

    Exactly one of `rpm` (fixed-speed machine) or `tsr` (variable-speed
    machine) must be given -- see the module docstring for why the two are
    kept explicitly distinct rather than defaulted.

    Parameters
    ----------
    geometry : RotorGeometry
    wind_speeds : iterable of float
        Freestream wind speeds to evaluate, m/s. All must be positive.
    rpm : float or None
        Constant rotor speed, RPM. TSR is derived per wind speed.
    tsr : float or None
        Constant tip-speed ratio. Rotor speed is derived per wind speed.
    air_density : float
        kg/m^3. Required and keyword-only, from `config/`.
    kinematic_viscosity : float
        m^2/s. Required and keyword-only, from `config/`.
    airfoil_for_reynolds : callable or None
        Passed through to `solve_rotor`.
    airfoils : list or None
        Passed through to `solve_rotor`.
    keep_stations : bool
        Include per-station detail on every point.

    Returns
    -------
    dict
        "mode": "fixed-speed" or "variable-speed"
        "points": list of `operating_point` dicts, in the order given
        plus flat per-quantity lists ("v_inf", "tsr", "Cp", "Ct", "power",
        "thrust", "torque", "rpm") for convenient plotting.
    """

    if (rpm is None) == (tsr is None):
        raise ValueError(
            "give exactly one of rpm= (fixed-speed) or tsr= (variable-speed); "
            "see the module docstring"
        )

    wind_speeds = [float(v) for v in wind_speeds]
    if not wind_speeds:
        raise ValueError("wind_speeds must not be empty")

    mode = "fixed-speed" if rpm is not None else "variable-speed"
    omega_fixed = rpm_to_omega(rpm) if rpm is not None else None

    points = []
    for v in wind_speeds:
        tsr_v = (omega_fixed * geometry.R / v) if omega_fixed is not None else tsr
        points.append(operating_point(
            geometry, v_inf=v, tsr=tsr_v, air_density=air_density,
            kinematic_viscosity=kinematic_viscosity,
            airfoil_for_reynolds=airfoil_for_reynolds,
            airfoils=airfoils, keep_stations=keep_stations,
        ))

    return _with_flat_lists({"mode": mode, "points": points})


def cp_lambda_curve(geometry: RotorGeometry, tsr_values, v_inf=7.0, *,
                    air_density, kinematic_viscosity,
                airfoil_for_reynolds=None, airfoils=None,
                    keep_stations=False):
    """
    Sweep TSR at fixed wind speed -- the Cp-lambda (and Ct-lambda) curve.

    Wind speed is held constant and rotor speed varies, matching the
    convention QBlade's and most textbooks' Cp-lambda sweeps use. `v_inf`
    still matters, despite Cp-lambda being nominally dimensionless, because
    it sets each station's Reynolds number for the polar lookup.

    Returns
    -------
    dict
        Same shape as `power_curve`, with "mode": "cp-lambda".
    """

    tsr_values = [float(t) for t in tsr_values]
    if not tsr_values:
        raise ValueError("tsr_values must not be empty")

    points = [
        operating_point(geometry, v_inf=v_inf, tsr=t, air_density=air_density,
                        kinematic_viscosity=kinematic_viscosity,
                        airfoil_for_reynolds=airfoil_for_reynolds,
                        airfoils=airfoils, keep_stations=keep_stations)
        for t in tsr_values
    ]
    return _with_flat_lists({"mode": "cp-lambda", "points": points})


def peak_cp(curve):
    """
    The maximum-Cp point of a curve, and its neighbours' Cp values.

    Returns
    -------
    dict
        The winning `operating_point` dict, plus "index" -- its position in
        curve["points"]. Useful both for reporting a rotor's design point and
        for the unimodality check in tests/test_powercurve.py.
    """

    points = curve["points"]
    index = max(range(len(points)), key=lambda i: points[i]["Cp"])
    best = dict(points[index])
    best["index"] = index
    return best


def _with_flat_lists(curve):
    """Attach per-quantity lists alongside curve["points"], for plotting."""

    for key in ("v_inf", "tsr", "rpm", "Cp", "Ct", "power", "thrust", "torque"):
        curve[key] = [p[key] for p in curve["points"]]
    return curve
