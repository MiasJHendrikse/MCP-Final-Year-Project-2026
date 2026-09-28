"""
The golden-file snapshot definition: one function, `build_golden()`, that
computes every quantity the golden regression watches.

Both sides of the regression use this one definition -- `generate_golden.py`
writes its output to `tests/golden/*.json`, and `test_golden_regression.py`
calls it again and compares. That is deliberate: if the snapshot were
recomputed independently in the test, the two definitions could drift and the
regression would start comparing two different things without saying so.

What is captured, and why each piece is here
---------------------------------------------
1. `cp_lambda`     -- the full Cp(lambda) and Ct(lambda) arrays for the NREL
                      Phase VI rotor across the cross-validated lambda range
                      (1.5 to 7.5, the span covered by the Sequence S
                      fixed-RPM sweep in docs/validation/bem-cross-validation.md).
                      This is the rotor-level integrated answer -- the thing a
                      refactor is most likely to be judged by and least likely
                      to visibly break.

2. `spanwise`      -- per-station a, a', phi, alpha, Re, F, Cl, Cd and the
                      local loads dT/dr and dQ/dr at four operating points on
                      the Sequence S fixed-RPM line. A mean Cp barely moves
                      when the spanwise distribution shifts materially, so the
                      rotor-level arrays above are not sufficient on their own.
                      The four points span the envelope:

                        v = 5.0 m/s   alpha  1.9 to  4.3 deg  fully attached
                        v = 7.0 m/s   alpha  3.4 to 10.2 deg  attached, stall approaching inboard
                        v = 10.0 m/s  alpha  5.9 to 19.7 deg  stall onset / near-stall
                        v = 15.0 m/s  alpha 10.6 to 32.4 deg  deep stall, past the S809 table's
                                                              +18 deg XFOIL-converged band, so
                                                              this point reads the Viterna
                                                              extrapolation

                      The last point was included precisely because
                      bem.airfoil.S809Polar's alpha clamp, due to be removed, used
                      to return Cl(18 deg) here for a station sitting at 32
                      deg. It changed, deliberately and visibly (-25 % in Cp);
                      see tests/golden/README.md's change log.

3. `cross_tool`    -- our solver's side of the three cross-tool comparisons,
                      recomputed through the exact airfoil construction each
                      comparison script uses. See the `cross_tool` section
                      below and tests/golden/README.md for the important
                      caveat about the CCBlade and pyBEMT reference values.

Everything here runs against the polar cache in `data/polars/s809/` as
committed. Nothing in this module calls XFOIL, and nothing writes to `data/`.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import math
import os

from bem.rotor import phase_vi_geometry, solve_rotor
from config import load_phase_vi_rotor

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
GOLDEN_DIR = os.path.join(_HERE, "golden")
DOCS_VALIDATION_DIR = os.path.join(REPO_ROOT, "docs", "validation")
QBLADE_DIR = os.path.join(REPO_ROOT, "data", "qblade")

# --------------------------------------------------------------------------
# Operating conditions. These are fixed inputs to the snapshot, not results;
# changing any of them invalidates the golden files and requires a regenerate.
# --------------------------------------------------------------------------

#: Cp-lambda sweep. Fixed wind speed with lambda varying is the convention
#: bem.powercurve.cp_lambda_curve uses; v_inf still matters because it sets
#: each station's Reynolds number for the polar lookup.
CP_LAMBDA_V_INF = 7.0
CP_LAMBDA_VALUES = [1.5 + 0.5 * i for i in range(13)]  # 1.5 .. 7.5

#: Spanwise detail points, on the Sequence S fixed-RPM line (71.63 RPM
#: synchronous), matching the wind speeds used by the cross-tool comparison.
SPANWISE_WIND_SPEEDS = [5.0, 7.0, 10.0, 15.0]

#: The cross-tool sweep's own wind speeds, matching compare_pybemt.py and
#: compare_ccblade.py exactly (NREL Sequence S test points).
CROSS_TOOL_WIND_SPEEDS = [5.0, 7.0, 10.0, 13.0, 15.0, 20.0, 25.0]

#: The Phase VI operating condition, from config/rotor_phase_vi.yaml.
#: Read from config rather than restated here on purpose: if the rewiring of
#: rho and nu through the config layer changed either number, every golden
#: comparison below would move, and the regression would say so. The values are
#: the same sea-level pair the snapshot was originally generated with.
_PHASE_VI = load_phase_vi_rotor()
PHASE_VI_RPM = _PHASE_VI.rated_rpm
AIR_DENSITY = _PHASE_VI.air_density              # kg/m^3, sea level, as per the experiment
KINEMATIC_VISCOSITY = _PHASE_VI.kinematic_viscosity  # m^2/s, sea level

#: The two fixed Reynolds numbers the QBlade export/comparison pair uses
#: (data/qblade/S809.plr and S809_Re500k.plr).
QBLADE_FIXED_REYNOLDS = [100000, 500000]


def _omega(rpm=PHASE_VI_RPM):
    return rpm * 2.0 * math.pi / 60.0


def _local_loads(stations, geometry):
    """
    Per-station dT/dr and dQ/dr, in N/m and N.m/m.

    Reproduces bem.rotor.solve_rotor's own integrand exactly -- solve_rotor
    computes these internally on its way to Ct/Cp but does not return them,
    and the local loads are part of what the snapshot captures. Defined here once so
    both the generator and the test derive them identically.
    """

    dt_dr = [
        0.5 * AIR_DENSITY * s["w"] ** 2 * geometry.n_blades * chord * s["Cn"]
        for s, chord in zip(stations, geometry.chord)
    ]
    dq_dr = [
        0.5 * AIR_DENSITY * s["w"] ** 2 * geometry.n_blades * chord * s["Ct"] * s["r"]
        for s, chord in zip(stations, geometry.chord)
    ]
    return dt_dr, dq_dr


#: Per-station fields carried into the golden file. "Ct" and "Cq" here follow
#: station.py's convention (blade-element tangential-force coefficients), not
#: the rotor-integrated thrust coefficient -- see solve_rotor's docstring.
STATION_FIELDS = ("r", "phi", "a", "a_prime", "alpha", "Cl", "Cd",
                  "Cn", "Ct", "Cq", "F", "reynolds", "w")


def build_cp_lambda():
    """Cp(lambda) and Ct(lambda) for the Phase VI rotor, default polar path."""

    geometry = phase_vi_geometry()
    cp, ct = [], []
    for lam in CP_LAMBDA_VALUES:
        result = solve_rotor(geometry, tsr=lam, v_inf=CP_LAMBDA_V_INF,
                             air_density=AIR_DENSITY,
                             kinematic_viscosity=KINEMATIC_VISCOSITY)
        cp.append(result["Cp"])
        ct.append(result["Ct"])

    return {
        "v_inf": CP_LAMBDA_V_INF,
        "air_density": AIR_DENSITY,
        "tsr": list(CP_LAMBDA_VALUES),
        "Cp": cp,
        "Ct": ct,
    }


def build_spanwise():
    """Per-station state and local loads at the four envelope-spanning points."""

    geometry = phase_vi_geometry()
    omega = _omega()
    points = []

    for v_inf in SPANWISE_WIND_SPEEDS:
        tsr = omega * geometry.R / v_inf
        result = solve_rotor(geometry, tsr=tsr, v_inf=v_inf,
                             air_density=AIR_DENSITY,
                             kinematic_viscosity=KINEMATIC_VISCOSITY)
        stations = result["stations"]
        dt_dr, dq_dr = _local_loads(stations, geometry)

        point = {
            "v_inf": v_inf,
            "rpm": PHASE_VI_RPM,
            "tsr": tsr,
            "air_density": AIR_DENSITY,
            "Cp": result["Cp"],
            "Ct": result["Ct"],
            "dT_dr": dt_dr,
            "dQ_dr": dq_dr,
        }
        for field in STATION_FIELDS:
            point[field] = [s[field] for s in stations]
        points.append(point)

    return {"points": points}


def _cross_tool_table_airfoils(geometry, omega):
    """
    The pinned per-station Reynolds-bucket polar tables that compare_pybemt.py
    and compare_ccblade.py both build, imported from compare_pybemt rather
    than reimplemented so this snapshot tracks the comparison scripts rather
    than a copy of them.
    """

    from validation.compare_pybemt import (
        TableS809Polar, build_dense_table, reynolds_bucket_for_station,
    )

    buckets = [
        reynolds_bucket_for_station(r, c, omega)
        for r, c in zip(geometry.r, geometry.chord)
    ]
    tables = {b: build_dense_table(b) for b in sorted(set(buckets))}
    return buckets, [TableS809Polar(*tables[b]) for b in buckets]


def _load_stored_external(subdir, key):
    """Stored external-solver Cp values from a committed comparison results.json."""

    path = os.path.join(DOCS_VALIDATION_DIR, subdir, "results.json")
    with open(path) as f:
        rows = json.load(f)
    return {row["v_inf"]: row[key] for row in rows}


def build_cross_tool():
    """
    Our solver's side of all three cross-tool comparisons.

    CCBlade / pyBEMT
    -----------------
    Both scripts build the same pinned per-station Reynolds-bucket tables and
    call solve_rotor over the same seven Sequence S wind speeds, so one sweep
    covers both. What is recorded is:

      * `Cp_ours` / `Ct_ours` at each point -- reproducible in-repo, and the
        actual regression quantity.
      * `mean_abs_cp_deviation_pct` against the *stored* external values.

    That second figure is recorded for tracking, NOT as a validation result.
    The committed docs/validation/{ccblade,pybemt}_case/results.json predate
    both the Ncrit=5 S809 cache rebuild and the extension of the Reynolds
    bucket list past 500k, so their external columns were produced from polar
    tables the repo no longer contains. Re-establishing the published 0.51 %
    and 1.83 % figures needs the external tools re-run against the current
    cache. See tests/golden/README.md.

    QBlade
    -------
    The QBlade pair pins every station to a single Reynolds number matching
    the exported .plr, exactly as compare_qblade.py does. This one *is*
    current: the committed data/qblade/bem_solver_phase_vi_result_Re*.csv
    reproduce bitwise from today's tree.
    """

    geometry = phase_vi_geometry()
    omega = _omega()
    buckets, airfoils = _cross_tool_table_airfoils(geometry, omega)

    stored_cc = _load_stored_external("ccblade_case", "Cp_ccblade")
    stored_pb = _load_stored_external("pybemt_case", "Cp_pybemt")

    rows, dev_cc, dev_pb = [], [], []
    for v_inf in CROSS_TOOL_WIND_SPEEDS:
        tsr = omega * geometry.R / v_inf
        result = solve_rotor(geometry, tsr=tsr, v_inf=v_inf,
                             air_density=AIR_DENSITY,
                             kinematic_viscosity=KINEMATIC_VISCOSITY,
                             airfoils=airfoils)
        rows.append({
            "v_inf": v_inf,
            "tsr": tsr,
            "Cp_ours": result["Cp"],
            "Ct_ours": result["Ct"],
        })
        dev_cc.append(abs(result["Cp"] - stored_cc[v_inf]) / stored_cc[v_inf])
        dev_pb.append(abs(result["Cp"] - stored_pb[v_inf]) / stored_pb[v_inf])

    # QBlade side: single pinned Reynolds per run, all stations.
    from polars.polar import CachedPolar, interpolant_for

    qblade = []
    interpolant = interpolant_for(geometry.polar_cache)
    for reynolds in QBLADE_FIXED_REYNOLDS:
        pinned = [CachedPolar(interpolant, reynolds) for _ in geometry.r]
        result = solve_rotor(geometry, tsr=omega * geometry.R / 7.0, v_inf=7.0,
                             air_density=AIR_DENSITY,
                             kinematic_viscosity=KINEMATIC_VISCOSITY,
                             airfoils=pinned)
        stations = result["stations"]
        entry = {
            "fixed_reynolds": reynolds,
            "v_inf": 7.0,
            "rpm": PHASE_VI_RPM,
            "Cp": result["Cp"],
            "Ct": result["Ct"],
        }
        for field in STATION_FIELDS:
            entry[field] = [s[field] for s in stations]
        qblade.append(entry)

    n = len(CROSS_TOOL_WIND_SPEEDS)
    return {
        "reynolds_buckets_per_station": buckets,
        "sweep": rows,
        "mean_abs_cp_deviation_pct": {
            "vs_stored_ccblade": 100.0 * sum(dev_cc) / n,
            "vs_stored_pybemt": 100.0 * sum(dev_pb) / n,
        },
        "qblade_pinned": qblade,
    }


def build_golden():
    """The complete snapshot, as written to / compared against tests/golden/."""

    return {
        "phase_vi_cp_lambda": build_cp_lambda(),
        "phase_vi_spanwise": build_spanwise(),
        "cross_tool_summary": build_cross_tool(),
    }


#: Golden filenames, keyed by the top-level `build_golden()` key they hold.
GOLDEN_FILES = {
    "phase_vi_cp_lambda": "phase_vi_cp_lambda.json",
    "phase_vi_spanwise": "phase_vi_spanwise.json",
    "cross_tool_summary": "cross_tool_summary.json",
}


def golden_path(key):
    return os.path.join(GOLDEN_DIR, GOLDEN_FILES[key])


def load_golden(key):
    with open(golden_path(key)) as f:
        return json.load(f)
