"""
Multi-station spanwise loop over the single-station solver, using real
cached polar data (polars.polar.CachedPolar) instead of a synthetic linear
polar. Which cache is a property of the blade --
see RotorGeometry.polar_cache.

Each station is solved independently via station.solve_station -- there is
deliberately no spanwise coupling/smoothing between stations. This module only adds: (a) a container for
spanwise geometry, (b) a per-station Reynolds estimate so the real polar
can be queried, (c) the loop itself, and (d) trapezoidal spanwise
integration to rotor-level Ct, Cp.

This module answers one operating point. Sweeping it across wind speeds or
tip-speed ratios -- the power curve and the Cp-lambda curve, in dimensional
units -- lives in powercurve.py, which calls straight through to
`solve_rotor` and adds no aerodynamics of its own.

Not done here (explicitly deferred): validation against NREL Phase VI
*experimental* performance data (still unsourced), AEP integration, and anything adjoint-related. The demo
geometry below is a synthetic, smoothly-tapered/twisted blade chosen only to
exercise the pipeline sensibly -- it is NOT NREL Phase VI's published
chord/twist table, which will be sourced properly in the dedicated
Phase VI validation stage.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math
from dataclasses import dataclass

from bem.station import StationParams, solve_station
from polars.polar import polar_factory_for

# No module-level air properties, deliberately. A sea-level kinematic
# viscosity used to be a module constant here and a sea-level density a
# default argument on solve_rotor below. The design site is at 1800 m, where
# both differ by around a fifth: omitting the density gave an answer 25 % high
# that still looked plausible, and the viscosity would have put every
# design-rotor Reynolds number 23 % high. Both are now required arguments,
# supplied from `config/` -- see src/config/. No literal air property appears
# anywhere under src/; grep is the check.


@dataclass
class RotorGeometry:
    """
    Spanwise blade geometry, one entry per station.

    Parameters
    ----------
    r : list of float
        Station radii, m. Must be strictly increasing, strictly greater
        than r_hub (if set), and strictly less than R.
    chord : list of float
        Station chord lengths, m.
    twist : list of float
        Station twist (+ pitch), radians.
    R : float
        Rotor (blade tip) radius, m.
    polar_cache : str
        Which polar cache this blade's sections are, named by its
        `config/polars_<name>.yaml` file -- "sg6043" for the design rotor,
        "s809" for the NREL Phase VI validation rotor.

        Required, with no default, and carried here rather than passed to
        `solve_rotor`, for one reason: the airfoil is a property of the blade,
        not of the call. Two caches are in use at the same time, and anything that lets the two drift apart -- a module-level
        "active airfoil" (since removed), or a default on the solver --
        produces a plausible number computed from the wrong table. This is the
        same rule that removed the air density and viscosity defaults: a default is
        worse than no default when the wrong answer stays inside the
        believable band.
    n_blades : int
        Number of blades, B.
    r_hub : float or None
        Hub radius, m. Optional -- see station.StationParams.
    """

    r: list
    chord: list
    twist: list
    R: float
    polar_cache: str
    n_blades: int = 3
    r_hub: float = None

    def __post_init__(self):
        n = len(self.r)
        if not (len(self.chord) == n and len(self.twist) == n):
            raise ValueError("r, chord, twist must all have the same length")
        if any(self.r[i] >= self.r[i + 1] for i in range(n - 1)):
            raise ValueError("r must be strictly increasing")


def demo_rotor_geometry(n_stations=15, R=5.5, r_hub=0.5, n_blades=3,
                        polar_cache="s809"):
    """
    A synthetic, smoothly-tapered/twisted blade spanning r_hub to R --
    linear chord taper (0.5 m -> 0.1 m) and linear twist (15 deg -> -2 deg),
    both monotonic so the spanwise loop has no built-in discontinuities to
    confound smoothness checks. This is a placeholder geometry
    for exercising the pipeline, not NREL Phase VI's actual blade (see
    module docstring).

    `polar_cache` defaults to "s809" only because that is what this
    placeholder has always been solved with -- naming it here keeps
    `RotorGeometry`'s required field satisfied at the one place this blade is
    defined, rather than reintroducing a solver-side default.
    """

    stations = n_stations
    tip_margin = 0.1   # keep the last station strictly inside R (see StationParams)
    hub_margin = 0.3   # keep the first station strictly outside r_hub, likewise --
                        # found the hard way: placing a station exactly at r_hub
                        # degenerates F_hub to 0 exactly as r=R does for F_tip
                        # (see station.StationParams.__post_init__).
    r_lo, r_hi = r_hub + hub_margin, R - tip_margin
    r = [r_lo + (r_hi - r_lo) * i / (stations - 1) for i in range(stations)]
    chord_root, chord_tip = 0.5, 0.1
    twist_root, twist_tip = math.radians(15.0), math.radians(-2.0)
    chord = []
    twist = []
    for ri in r:
        frac = (ri - r_lo) / (r_hi - r_lo)
        chord.append(chord_root + (chord_tip - chord_root) * frac)
        twist.append(twist_root + (twist_tip - twist_root) * frac)

    return RotorGeometry(r=r, chord=chord, twist=twist, R=R,
                         polar_cache=polar_cache, n_blades=n_blades, r_hub=r_hub)


#: NREL Phase VI, Table A-1 "Blade chord and twist distributions"
#: (Hand, M.M. et al. (2001), NREL/TP-500-29955, p.64), airfoil-bearing
#: stations only (r >= 1.2575 m -- the report states "Except for the root
#: [0 to 1.257 m], the blade uses the S809 at all span locations"; the
#: cylindrical root/transition stations below that are excluded here, see
#: PHASE_VI_ROOT_EXCLUDED_NOTE). Twist is relative to zero twist at the
#: 3.772 m (75% span) station, per the table's footnote 2 -- see
#: phase_vi_geometry() for how the operating tip-pitch angle is added to
#: get the absolute local twist this solver expects.
PHASE_VI_TABLE_A1_STATIONS = [
    # r (m), chord (m), twist_relative_to_75pct_span (deg)
    (1.2575, 0.737, 20.040),
    (1.343, 0.728, 18.074),
    (1.510, 0.711, 14.292),
    (1.648, 0.697, 11.909),
    (1.952, 0.666, 7.979),
    (2.257, 0.636, 5.308),
    (2.343, 0.627, 4.715),
    (2.562, 0.605, 3.425),
    (2.867, 0.574, 2.083),
    (3.172, 0.543, 1.150),
    (3.185, 0.542, 1.115),
    (3.476, 0.512, 0.494),
    (3.781, 0.482, -0.015),
    (4.023, 0.457, -0.381),
    (4.086, 0.451, -0.475),
    (4.391, 0.420, -0.920),
    (4.696, 0.389, -1.352),
    (4.780, 0.381, -1.469),
    (5.000, 0.358, -1.775),
]

PHASE_VI_ROOT_EXCLUDED_NOTE = (
    "Table A-1 also lists three stations from r=0.5083 m (hub attachment) "
    "to r=0.8835 m that are a cylindrical, non-airfoil root adapter (chord "
    "constant at 0.218/0.183 m, twist 0 deg, no S809 polar applies), plus a "
    "transition region blending into the S809 profile by r=1.2575 m. "
    "PHASE_VI_TABLE_A1_STATIONS omits all of these, consistent with how "
    "the demo geometry's root is handled: no published 2D polar "
    "exists for this non-airfoil section, and it contributes negligible "
    "torque at its small radius. r_hub is still set to the real 0.508 m "
    "hub-attachment point for the Prandtl hub-loss factor."
)

#: Standard rotor radius (10.058 m diameter / 2), matching every published
#: Phase VI reference; the table's two 5.305/5.532 m rows belong to the
#: 11.064 m tip-extension configuration (Sequence W only) and are excluded.
PHASE_VI_R = 5.029
PHASE_VI_R_HUB = 0.508
PHASE_VI_N_BLADES = 2
PHASE_VI_RATED_RPM = 71.63  # "Rotational speed: 71.63 RPM synchronous speed" (Appendix A)
PHASE_VI_SEQUENCE_S_TIP_PITCH_DEG = 3.0  # standard Sequence S operator pitch setting


def phase_vi_geometry(tip_pitch_deg=PHASE_VI_SEQUENCE_S_TIP_PITCH_DEG):
    """
    NREL Phase VI blade geometry (Sequence S configuration), from the
    published Table A-1 station table -- see PHASE_VI_TABLE_A1_STATIONS
    and PHASE_VI_ROOT_EXCLUDED_NOTE above for exactly what is and is not
    included, and the module docstring for how this was sourced.

    This is real published geometry, unlike demo_rotor_geometry()'s
    synthetic blade -- used for the pyBEMT and CCBlade cross-checks, not (yet)
    for validation against Phase VI experimental performance data, which
    remains unavailable.

    Parameters
    ----------
    tip_pitch_deg : float
        Operator-set collective pitch, added uniformly to the table's
        twist-relative-to-75%-span values to get each station's absolute
        local twist (see module docstring). Default 3.0 deg, the standard
        Sequence S setting.

    Returns
    -------
    RotorGeometry
    """

    r = [row[0] for row in PHASE_VI_TABLE_A1_STATIONS]
    chord = [row[1] for row in PHASE_VI_TABLE_A1_STATIONS]
    twist = [math.radians(row[2] + tip_pitch_deg) for row in PHASE_VI_TABLE_A1_STATIONS]

    return RotorGeometry(
        r=r, chord=chord, twist=twist, R=PHASE_VI_R,
        # Phase VI is an S809 blade -- "Except for the root, the blade uses
        # the S809 at all span locations" (Hand et al. 2001, Appendix A).
        # Pinned here, at the blade's definition, so this validation rotor can
        # never inherit the design rotor's SG6043 sections from a caller that
        # did not say which airfoil it meant.
        polar_cache="s809",
        n_blades=PHASE_VI_N_BLADES, r_hub=PHASE_VI_R_HUB,
    )


def solve_rotor(geometry: RotorGeometry, tsr, v_inf=7.0, *,
                air_density, kinematic_viscosity,
                airfoil_for_reynolds=None, airfoils=None):
    """
    Solve every station independently (single-station solver, real cached polar)
    and integrate to rotor-level Ct, Cp.

    Parameters
    ----------
    geometry : RotorGeometry
    tsr : float
        Rotor tip-speed ratio, Omega * R / v_inf.
    v_inf : float
        Freestream wind speed, m/s. Only affects the per-station Reynolds
        estimate (used to query the real polar) -- the induction solve
        itself is dimensionless (local tsr only).
    air_density : float
        kg/m^3, used only for the rotor-level Ct/Cp integration. Required,
        and keyword-only: there is no defensible default. The validation
        rotor's sea-level value is a property of that experiment and lives in
        `config/rotor_phase_vi.yaml`; the design rotor's is a property of the
        site (1800 m, about 20 % less dense) and lives in `config/site.yaml`.
    kinematic_viscosity : float
        m^2/s, used only for the per-station Reynolds estimate that queries
        the polar. Required and keyword-only for the same reason -- the site
        value is 23 % above the sea-level one, and that difference lands
        directly on every Reynolds number.
    airfoil_for_reynolds : callable or None
        Optional `reynolds -> polar object` factory, overriding the
        `polars.polar.CachedPolar` built over `geometry.polar_cache`. The
        normal path needs nothing here: the blade already names its own
        airfoil (see `RotorGeometry.polar_cache`), which is why there is no
        default airfoil on this signature to get wrong.
    airfoils : list of airfoil objects or None
        Optional, one per station, overriding both of the above -- used by the
        pyBEMT cross-check (compare_pybemt.py) and the QBlade
        comparison to force both solvers onto the exact same discretized
        (alpha, Re-bucket) polar table rather than this module's normal
        continuous Re interpolation, so any difference in results comes from
        the solver, not the input data. Reynolds is still estimated and
        reported either way.

    Raises
    ------
    polars.interpolant.PolarDomainError
        If a station's estimated Reynolds number or solved angle of attack
        falls outside what `geometry.polar_cache` covers. This
        raises rather than clamping to the table edge: an out-of-envelope
        station is a cache-coverage fact the caller needs, not something to
        be silently replaced with the nearest value that happens to exist.

    Returns
    -------
    dict
        "stations": list of per-station dicts, each with r, phi, a,
        a_prime, alpha, Cl, Cd, Cn, Ct, Cq, F, reynolds (Ct/Cq here follow
        station.py's convention: the blade-element tangential-force
        coefficient, not the rotor-integrated thrust coefficient below).
        "Ct": rotor-integrated thrust coefficient.
        "Cp": rotor-integrated power coefficient.
    """

    omega = tsr * v_inf / geometry.R
    stations = []

    if airfoils is None and airfoil_for_reynolds is None:
        # Built once per call and shared by every station: the interpolant
        # behind it is immutable and memoised per cache, so this costs one
        # dict lookup, not a rebuild (see polars.polar.interpolant_for).
        airfoil_for_reynolds = polar_factory_for(geometry.polar_cache)

    for i, (r, chord, twist) in enumerate(zip(geometry.r, geometry.chord, geometry.twist)):
        # Zero-induction relative-velocity estimate for the polar's Reynolds
        # number -- a fixed Re per station, not re-solved iteratively;
        # adequate for this stage's sanity check, not for AEP-grade accuracy.
        # This estimate is load-bearing: an estimate outside the cache's Reynolds range raises here
        # instead of clamping to the ceiling, which is what hid the total loss
        # of Reynolds dependence in every Phase VI run before 2026-07-28.
        w_approx = math.hypot(v_inf, omega * r)
        reynolds = w_approx * chord / kinematic_viscosity

        airfoil = airfoils[i] if airfoils is not None else airfoil_for_reynolds(reynolds)

        tsr_local = omega * r / v_inf
        station = StationParams(
            r=r, chord=chord, twist=twist, airfoil=airfoil, tsr=tsr_local,
            R=geometry.R, n_blades=geometry.n_blades, r_hub=geometry.r_hub,
        )
        result = solve_station(station)

        alpha = result["phi"] - twist
        cn = result["Cl"] * math.cos(result["phi"]) + result["Cd"] * math.sin(result["phi"])
        w = v_inf * (1.0 - result["a"]) / math.sin(result["phi"])

        stations.append({
            "r": r,
            "phi": result["phi"],
            "a": result["a"],
            "a_prime": result["a_prime"],
            "alpha": alpha,
            "Cl": result["Cl"],
            "Cd": result["Cd"],
            "Cn": cn,
            "Ct": result["Ct"],
            "Cq": result["Cq"],
            "F": result["F"],
            "reynolds": reynolds,
            "w": w,
            # Solver status, carried per station. A station that
            # did not converge is reported, not raised -- see solve_station.
            "converged": result["converged"],
            "residual": result["residual"],
            "residual_initial": result["residual_initial"],
            "iterations": result["iterations"],
            "failure": result["failure"],
        })

    r_arr = [s["r"] for s in stations]
    dt_dr = [
        0.5 * air_density * s["w"] ** 2 * geometry.n_blades * chord * s["Cn"]
        for s, chord in zip(stations, geometry.chord)
    ]
    dq_dr = [
        0.5 * air_density * s["w"] ** 2 * geometry.n_blades * chord * s["Ct"] * s["r"]
        for s, chord in zip(stations, geometry.chord)
    ]

    thrust = _trapz(dt_dr, r_arr)
    torque = _trapz(dq_dr, r_arr)
    power = torque * omega

    ct_rotor = thrust / (0.5 * air_density * v_inf ** 2 * math.pi * geometry.R ** 2)
    cp_rotor = power / (0.5 * air_density * v_inf ** 3 * math.pi * geometry.R ** 2)

    # Rotor-level status. Ct and Cp are still returned when a station failed,
    # deliberately: the smoothness gate has to be able to plot the point
    # to see where the objective breaks down, and a None there would make it
    # blind exactly where it most needs to see. The flag sits alongside the
    # number rather than replacing it, and `failed_stations` names which.
    failed = [i for i, s in enumerate(stations) if not s["converged"]]

    return {
        "stations": stations,
        "Ct": ct_rotor,
        "Cp": cp_rotor,
        "converged": not failed,
        "failed_stations": failed,
    }


def _trapz(y, x):
    """Trapezoidal integration, no numpy dependency."""

    return sum(0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]) for i in range(len(x) - 1))
