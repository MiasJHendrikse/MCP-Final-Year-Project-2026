"""
Stage 6: cross-check our BEM solver (Stages 1-4) against pyBEMT
(https://github.com/kegiljarhus/pyBEMT, MIT licensed) on the real NREL
Phase VI rotor geometry (bem.rotor.phase_vi_geometry).

THIS IS NOT VALIDATION AGAINST EXPERIMENTAL DATA. Phase VI's own measured
Cp-lambda curve could not be sourced this session (see the Stage 5 journal
entry) -- this script only compares two independent BEM implementations
against each other on the same geometry and (as close as practical) the
same input polar data. Agreement here says "two different pieces of BEM
code, with different corrections and root-finding, land in the same
neighbourhood" -- it is not evidence that either one matches the real
rotor. See docs/validation/pybemt-comparison.md for the full writeup and
how to reproduce this from a clean checkout.

pyBEMT is NOT vendored into this repository. It is expected to already be
cloned and installed as a sibling directory, by default
"../../pybemt-reference" relative to this file's repo root (override with
the PYBEMT_REFERENCE_PATH environment variable), inside its own virtualenv
-- see the "How to reproduce" section of the doc above for exact commands,
including a one-line Python-3.12 compatibility patch pyBEMT's own
SafeConfigParser usage needs (removed from the stdlib in 3.12; this repo
does not carry that patch, it is applied to the external clone only). Run
this script with THAT virtualenv's python, invoked from src/ so `bem` and
`xfoil` resolve the same way every validate_stageN script already relies
on:

    cd src
    "<path-to-pybemt-reference>/.venv/Scripts/python.exe" -m bem.compare_pybemt

Fairness of the comparison
---------------------------
Both solvers are given the exact same discretized (alpha -> Cl, Cd) table
per Reynolds bucket -- generated once here from data/polars/s809/*.csv,
extended flat (constant) outside the XFOIL-converged [-8, 18] deg range out
to +/-180 deg (neither solver's root-search needs that extrapolation to be
aerodynamically accurate, only defined, since a converged solution never
lands out there -- see station._select_bracket's docstring for why our own
solver already assumes this). The same table is written as an AeroDyn-style
.dat file for pyBEMT (bilinear/quadratic-interpolated by pyBEMT's own
loader) and used directly (linearly interpolated, via numpy) for our
solver, bypassing bem.rotor's normal continuous-Reynolds S809Polar for
this comparison only, via solve_rotor's `airfoils=` override. The remaining
known difference is linear (ours) vs. quadratic (pyBEMT's) interpolation
between table points -- judged negligible given the table is generated at
a 0.5 deg resolution, and left as-is rather than reimplementing pyBEMT's
interpolation choice.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

from bem.rotor import phase_vi_geometry, solve_rotor  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data")
DOCS_VALIDATION_DIR = os.path.join(REPO_ROOT, "docs", "validation")

PYBEMT_REFERENCE_PATH = os.environ.get(
    "PYBEMT_REFERENCE_PATH",
    os.path.join(os.path.dirname(REPO_ROOT), "pybemt-reference"),
)

RE_BUCKETS = [100000, 150000, 200000, 300000, 400000, 500000]
AIR_DENSITY = 1.225        # kg/m^3
AIR_KINEMATIC_VISCOSITY = 1.5e-5  # m^2/s -- matches bem.rotor.AIR_KINEMATIC_VISCOSITY
AIR_DYNAMIC_VISCOSITY = AIR_DENSITY * AIR_KINEMATIC_VISCOSITY  # kept consistent with the above

ALPHA_GRID_DEG = np.arange(-180.0, 180.01, 0.5)


def _load_cache_csv(re_bucket):
    path = os.path.join(DATA_DIR, "polars", "s809", f"S809_Re{re_bucket}.csv")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1], data[:, 2]  # alpha, cl, cd


def build_dense_table(re_bucket):
    """
    (alpha_grid_deg, cl, cd) for one Reynolds bucket: linear interpolation
    of the cached XFOIL data over its converged range, flat-extrapolated
    to +/-180 deg outside it. Shared, identical input for both solvers --
    see module docstring.
    """

    alpha_c, cl_c, cd_c = _load_cache_csv(re_bucket)
    alpha_lo, alpha_hi = alpha_c.min(), alpha_c.max()

    alpha_clamped = np.clip(ALPHA_GRID_DEG, alpha_lo, alpha_hi)
    cl = np.interp(alpha_clamped, alpha_c, cl_c)
    cd = np.interp(alpha_clamped, alpha_c, cd_c)
    return ALPHA_GRID_DEG.copy(), cl, cd


class TableS809Polar:
    """
    Our-solver-side airfoil adapter over a pre-built dense (alpha, Cl, Cd)
    table (see build_dense_table) -- linear interpolation via numpy, radians
    in, matching station.py's duck-typed .cl(alpha)/.cd(alpha) interface.
    Used only by this comparison script, to force bit-close agreement with
    the AeroDyn table written for pyBEMT (see module docstring).
    """

    def __init__(self, alpha_deg, cl, cd):
        self._alpha_deg = alpha_deg
        self._cl = cl
        self._cd = cd

    def cl(self, alpha_rad):
        return float(np.interp(math.degrees(alpha_rad), self._alpha_deg, self._cl))

    def cd(self, alpha_rad):
        return float(np.interp(math.degrees(alpha_rad), self._alpha_deg, self._cd))


def write_aerodyn_dat(path, alpha_deg, cl, cd, label):
    """Write an AeroDyn-format single-table airfoil file (pyBEMT's expected format)."""

    lines = [
        f"AeroDyn airfoil file, generated by compare_pybemt.py for {label}",
        f'Polar "{label}" :: generated from data/polars/s809 (Phase 0 XFOIL cache)',
        "1              Number of airfoil tables in this file",
        "0              Table ID parameter",
        "0.0            Stall angle (deg) -- not used by pyBEMT's loader",
        "0              No longer used, enter zero",
        "0              No longer used, enter zero",
        "0              No longer used, enter zero",
        "0.0            Angle of attack for zero Cn for linear Cn curve (deg) -- not used",
        "0.0            Cn slope for zero lift for linear Cn curve (1/rad) -- not used",
        "0.0            Cn at stall value for positive angle of attack -- not used",
        "0.0            Cn at stall value for negative angle of attack -- not used",
        "0.0            Angle of attack for minimum CD (deg) -- not used",
        "0.0            Minimum CD value -- not used",
    ]
    assert len(lines) == 14, "pyBEMT's load_airfoil uses skiprows=14"

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
        for a, l, d in zip(alpha_deg, cl, cd):
            f.write(f"{a:10.2f}{l:10.4f}{d:10.4f}\n")


def reynolds_bucket_for_station(r, chord, omega, v_inf_reference=7.0):
    """Nearest cached Re bucket for a station, using the same zero-induction
    estimate as bem.rotor.solve_rotor, at a fixed reference wind speed (kept
    fixed across the whole sweep -- see module docstring on why Re is not
    re-estimated per sweep point)."""

    w_approx = math.hypot(v_inf_reference, omega * r)
    reynolds = w_approx * chord / AIR_KINEMATIC_VISCOSITY
    return min(RE_BUCKETS, key=lambda b: abs(b - reynolds))


def trapezoidal_dr(r):
    """
    Per-station integration width w_i such that sum(f_i * w_i) exactly
    equals the trapezoidal-rule integral of f(r) sampled at the (possibly
    irregularly spaced) stations r -- i.e. w_i = 0.5*(r[i+1]-r[i-1]) for
    interior stations, half that at each end. Phase VI's real station
    table (bem.rotor.PHASE_VI_TABLE_A1_STATIONS) is far from uniformly
    spaced (e.g. 0.014 m between two stations, 0.304 m between two others),
    so this cannot be left to pyBEMT's own default -- Rotor.__init__ falls
    back to a single dr = r[1]-r[0] applied to *every* station when `dr`
    isn't given in the config, which is silently wrong for a non-uniform
    table (found the hard way: it inflated our solver's Ct/Cp by ~2.5x
    relative to pyBEMT before this was added, entirely from the
    integration setup, not from any solver/correction difference -- see
    the Stage 6 journal entry). Matches bem.rotor._trapz's own weighting
    exactly, so the two solvers integrate identical per-station loads the
    same way.
    """

    n = len(r)
    if n == 1:
        return [1.0]
    w = [0.0] * n
    w[0] = 0.5 * (r[1] - r[0])
    w[-1] = 0.5 * (r[-1] - r[-2])
    for i in range(1, n - 1):
        w[i] = 0.5 * (r[i + 1] - r[i - 1])
    return w


def write_pybemt_ini(path, geometry, rpm, rho, mu, station_buckets):
    sections = " ".join(f"S809_Re{b}" for b in station_buckets)
    radii = " ".join(f"{r:.4f}" for r in geometry.r)
    chords = " ".join(f"{c:.4f}" for c in geometry.chord)
    pitches = " ".join(f"{math.degrees(t):.4f}" for t in geometry.twist)
    dr = " ".join(f"{w:.5f}" for w in trapezoidal_dr(geometry.r))

    content = f"""\
; Auto-generated by compare_pybemt.py -- NREL Phase VI (Sequence S) geometry,
; Table A-1 stations (see bem.rotor.phase_vi_geometry). Regenerate by
; re-running compare_pybemt.py rather than hand-editing.
[case]
rpm = {rpm}
v_inf = 7.0
twist = 0.0

[turbine]
nblades = {geometry.n_blades}
diameter = {2 * geometry.R}
radius_hub = {geometry.r_hub}
section = {sections}
radius = {radii}
chord  = {chords}
pitch = {pitches}
dr = {dr}

[fluid]
rho = {rho}
mu = {mu}

[solver]
solver = bisect
"""
    with open(path, "w") as f:
        f.write(content)


def main():
    if not os.path.isdir(PYBEMT_REFERENCE_PATH):
        raise SystemExit(
            f"pyBEMT reference clone not found at '{PYBEMT_REFERENCE_PATH}'. "
            "See docs/validation/pybemt-comparison.md 'How to reproduce' -- "
            "clone/install it externally first, or set PYBEMT_REFERENCE_PATH."
        )
    try:
        from pybemt.solver import Solver
    except ImportError as e:
        raise SystemExit(
            "Could not import pybemt. This script must be run with the "
            "pyBEMT reference clone's own virtualenv python (it is not a "
            "dependency of this repo) -- see docs/validation/"
            f"pybemt-comparison.md. Original error: {e}"
        )

    geometry = phase_vi_geometry()
    rpm = 71.63
    omega = rpm * 2 * math.pi / 60.0

    station_buckets = [
        reynolds_bucket_for_station(r, c, omega) for r, c in zip(geometry.r, geometry.chord)
    ]
    unique_buckets = sorted(set(station_buckets))
    print(f"Reynolds buckets used across {len(geometry.r)} stations: {unique_buckets}")

    tables = {b: build_dense_table(b) for b in unique_buckets}

    airfoils_dir = os.path.join(PYBEMT_REFERENCE_PATH, "pybemt", "airfoils")
    for b, (alpha_deg, cl, cd) in tables.items():
        write_aerodyn_dat(
            os.path.join(airfoils_dir, f"S809_Re{b}.dat"), alpha_deg, cl, cd, f"S809_Re{b}"
        )

    our_airfoils = [TableS809Polar(*tables[b]) for b in station_buckets]

    os.makedirs(os.path.join(DOCS_VALIDATION_DIR, "pybemt_case"), exist_ok=True)
    ini_path = os.path.join(DOCS_VALIDATION_DIR, "pybemt_case", "phase_vi.ini")
    write_pybemt_ini(ini_path, geometry, rpm, AIR_DENSITY, AIR_DYNAMIC_VISCOSITY, station_buckets)

    wind_speeds = [5.0, 7.0, 10.0, 13.0, 15.0, 20.0, 25.0]  # matches Sequence S test points

    solver = Solver(ini_path)

    rows = []
    for v_inf in wind_speeds:
        tsr = omega * geometry.R / v_inf

        ours = solve_rotor(geometry, tsr=tsr, v_inf=v_inf, air_density=AIR_DENSITY, airfoils=our_airfoils)

        solver.v_inf = v_inf
        solver.rpm = rpm
        T, Q, P, _sec_df = solver.run()
        _tsr_check, cp_pybemt, ct_pybemt = solver.turbine_coeffs(T, Q, P)

        rows.append({
            "v_inf": v_inf,
            "tsr": tsr,
            "Cp_ours": ours["Cp"],
            "Ct_ours": ours["Ct"],
            "Cp_pybemt": cp_pybemt,
            "Ct_pybemt": ct_pybemt,
        })

    return rows, ini_path


def make_plot(rows, out_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tsr = [r["tsr"] for r in rows]
    fig, (ax_cp, ax_ct) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax_cp.plot(tsr, [r["Cp_ours"] for r in rows], "o-", label="Our solver (Stages 1-4)")
    ax_cp.plot(tsr, [r["Cp_pybemt"] for r in rows], "s--", label="pyBEMT")
    ax_cp.set_xlabel("Tip-speed ratio")
    ax_cp.set_ylabel("$C_p$")
    ax_cp.set_title("Power coefficient")
    ax_cp.legend()
    ax_cp.grid(alpha=0.3)

    ax_ct.plot(tsr, [r["Ct_ours"] for r in rows], "o-", label="Our solver (Stages 1-4)")
    ax_ct.plot(tsr, [r["Ct_pybemt"] for r in rows], "s--", label="pyBEMT")
    ax_ct.set_xlabel("Tip-speed ratio")
    ax_ct.set_ylabel("$C_t$")
    ax_ct.set_title("Thrust coefficient")
    ax_ct.legend()
    ax_ct.grid(alpha=0.3)

    fig.suptitle(
        "Our solver vs. pyBEMT on NREL Phase VI geometry (solver cross-check, not vs. experiment)"
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote plot: {out_path}")


def make_deviation_table(rows):
    """Markdown deviation table: peak Cp (within this sweep) + two off-design points."""

    peak = max(rows, key=lambda r: r["Cp_ours"])
    off_design = [r for r in rows if r["v_inf"] in (10.0, 25.0)]
    highlighted = [peak] + off_design

    lines = [
        "| Wind speed (m/s) | TSR | Cp (ours) | Cp (pyBEMT) | Cp deviation | Ct (ours) | Ct (pyBEMT) | Ct deviation |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in highlighted:
        cp_dev = r["Cp_ours"] - r["Cp_pybemt"]
        cp_pct = 100.0 * cp_dev / r["Cp_pybemt"]
        ct_dev = r["Ct_ours"] - r["Ct_pybemt"]
        ct_pct = 100.0 * ct_dev / r["Ct_pybemt"]
        tag = " (peak $C_p$ in this sweep)" if r is peak else ""
        lines.append(
            f"| {r['v_inf']:.1f}{tag} | {r['tsr']:.3f} | {r['Cp_ours']:.4f} | {r['Cp_pybemt']:.4f} | "
            f"{cp_dev:+.4f} ({cp_pct:+.1f}%) | {r['Ct_ours']:.4f} | {r['Ct_pybemt']:.4f} | "
            f"{ct_dev:+.4f} ({ct_pct:+.1f}%) |"
        )
    return "\n".join(lines)


def make_full_table(rows):
    lines = [
        "| Wind speed (m/s) | TSR | Cp (ours) | Cp (pyBEMT) | Ct (ours) | Ct (pyBEMT) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['v_inf']:.1f} | {r['tsr']:.3f} | {r['Cp_ours']:.4f} | {r['Cp_pybemt']:.4f} | "
            f"{r['Ct_ours']:.4f} | {r['Ct_pybemt']:.4f} |"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    rows, ini_path = main()
    print(f"\nWrote pyBEMT config: {ini_path}\n")
    print(f"{'v_inf':>6} {'TSR':>7} {'Cp_ours':>9} {'Cp_pybemt':>10} {'Ct_ours':>9} {'Ct_pybemt':>10}")
    for row in rows:
        print(f"{row['v_inf']:6.1f} {row['tsr']:7.3f} {row['Cp_ours']:9.4f} "
              f"{row['Cp_pybemt']:10.4f} {row['Ct_ours']:9.4f} {row['Ct_pybemt']:10.4f}")

    plot_path = os.path.join(DOCS_VALIDATION_DIR, "pybemt_cp_lambda_comparison.png")
    make_plot(rows, plot_path)

    print("\n--- Full sweep table (markdown) ---\n")
    print(make_full_table(rows))
    print("\n--- Deviation table: peak Cp + off-design points (markdown) ---\n")
    print(make_deviation_table(rows))
