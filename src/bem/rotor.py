"""
Stage 4: multi-station spanwise loop over the Stage 1-3 single-station
solver, using real S809 polar data (airfoil.S809Polar) instead of Stage 1's
synthetic linear polar.

Each station is solved independently via station.solve_station -- there is
deliberately no spanwise coupling/smoothing between stations at this stage
(see PROJECT_PLAN.md Phase 1). This module only adds: (a) a container for
spanwise geometry, (b) a per-station Reynolds estimate so the real polar
can be queried, (c) the loop itself, and (d) trapezoidal spanwise
integration to rotor-level Ct, Cp.

Not done here (explicitly deferred): the Cp-lambda sweep validation against
NREL Phase VI, AEP integration, and anything adjoint-related. The demo
geometry below is a synthetic, smoothly-tapered/twisted blade chosen only to
exercise the pipeline sensibly -- it is NOT NREL Phase VI's published
chord/twist table, which will be sourced properly in the dedicated
Phase VI validation stage.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
from dataclasses import dataclass

from bem.airfoil import S809Polar
from bem.station import StationParams, solve_station

AIR_KINEMATIC_VISCOSITY = 1.5e-5  # m^2/s, standard sea-level air


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
    n_blades : int
        Number of blades, B.
    r_hub : float or None
        Hub radius, m. Optional -- see station.StationParams.
    """

    r: list
    chord: list
    twist: list
    R: float
    n_blades: int = 3
    r_hub: float = None

    def __post_init__(self):
        n = len(self.r)
        if not (len(self.chord) == n and len(self.twist) == n):
            raise ValueError("r, chord, twist must all have the same length")
        if any(self.r[i] >= self.r[i + 1] for i in range(n - 1)):
            raise ValueError("r must be strictly increasing")


def demo_rotor_geometry(n_stations=15, R=5.5, r_hub=0.5, n_blades=3):
    """
    A synthetic, smoothly-tapered/twisted blade spanning r_hub to R --
    linear chord taper (0.5 m -> 0.1 m) and linear twist (15 deg -> -2 deg),
    both monotonic so the spanwise loop has no built-in discontinuities to
    confound the Stage 4 smoothness checks. This is a placeholder geometry
    for exercising the pipeline, not NREL Phase VI's actual blade (see
    module docstring).
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

    return RotorGeometry(r=r, chord=chord, twist=twist, R=R, n_blades=n_blades, r_hub=r_hub)


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
    "Stage 4/5 handled the demo geometry's root: no published 2D polar "
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
    included, and module docstring / the 2026-07-25 journal entry for how
    this was sourced and verified.

    This is real published geometry, unlike demo_rotor_geometry()'s
    synthetic blade -- built for the Stage 6 pyBEMT cross-check, not (yet)
    for validation against Phase VI experimental performance data, which
    remains unavailable (see the Stage 5 journal entry).

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
        n_blades=PHASE_VI_N_BLADES, r_hub=PHASE_VI_R_HUB,
    )


def solve_rotor(geometry: RotorGeometry, tsr, v_inf=7.0, air_density=1.225, airfoils=None):
    """
    Solve every station independently (Stage 1-3 solver, real S809 polar)
    and integrate to rotor-level Ct, Cp.

    Parameters
    ----------
    geometry : RotorGeometry
    tsr : float
        Rotor tip-speed ratio, Omega * R / v_inf.
    v_inf : float
        Freestream wind speed, m/s. Only affects the per-station Reynolds
        estimate (used to query the real polar) -- the induction solve
        itself is dimensionless (local tsr only), as in Stage 1-3.
    air_density : float
        kg/m^3, used only for the rotor-level Ct/Cp integration.
    airfoils : list of airfoil objects or None
        Optional, one per station, overriding the default per-station
        S809Polar(reynolds) construction -- used by the Stage 6 pyBEMT
        cross-check (compare_pybemt.py) to force both solvers onto the
        exact same discretized (alpha, Re-bucket) polar table rather than
        this module's normal continuous Re interpolation, so any
        difference in results comes from the solver, not the input data.
        Reynolds is still estimated and reported either way. Default
        (None) preserves Stage 4's original behaviour exactly.

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

    for i, (r, chord, twist) in enumerate(zip(geometry.r, geometry.chord, geometry.twist)):
        # Zero-induction relative-velocity estimate for the Reynolds lookup
        # (see airfoil.S809Polar) -- a fixed Re per station, not re-solved
        # iteratively; adequate for this stage's sanity check, not for
        # AEP-grade accuracy.
        w_approx = math.hypot(v_inf, omega * r)
        reynolds = w_approx * chord / AIR_KINEMATIC_VISCOSITY

        airfoil = airfoils[i] if airfoils is not None else S809Polar(reynolds)

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

    return {"stations": stations, "Ct": ct_rotor, "Cp": cp_rotor}


def _trapz(y, x):
    """Trapezoidal integration, no numpy dependency."""

    return sum(0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]) for i in range(len(x) - 1))
