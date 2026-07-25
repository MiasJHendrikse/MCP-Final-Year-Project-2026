"""
Stage 6: cross-check our BEM solver (Stages 1-4) against CCBlade
(https://github.com/WISDEM/CCBlade, part of NREL's WISDEM stack, Apache-2.0
licensed) on the real NREL Phase VI rotor geometry (bem.rotor.phase_vi_geometry).

CCBlade is the reference implementation of S. Andrew Ning's BEM formulation
(Ning 2014) -- the same underlying method this project's solver uses,
unlike pyBEMT's more classic BEM implementation (compare_pybemt.py). Closer
agreement with CCBlade than with pyBEMT would therefore be a meaningful
result, not just "two more digits of agreement" -- see
docs/validation/bem-cross-validation.md for whether that's what was found.

This is the second of three scripts in the Stage 6 cross-check pipeline
(see compare_pybemt.py's module docstring for the split rationale): this
one, run with CCBlade's own virtualenv, writes
docs/validation/ccblade_case/results.json; plot_bem_comparison.py (plain
repo python) reads it together with pyBEMT's results and produces the
combined three-way plot and deviation tables.

THIS IS NOT VALIDATION AGAINST EXPERIMENTAL DATA -- see compare_pybemt.py's
module docstring; the same caveat applies here.

CCBlade is NOT vendored into this repository. It ships as part of WISDEM
and needs a compiled Fortran extension (`wisdem.ccblade._bem`); a prebuilt
wheel exists on PyPI for this platform/Python version (`pip install
wisdem`), so no local Fortran/meson/ninja toolchain was actually needed
here -- see the "How to reproduce" section of the doc above for what to do
if a prebuilt wheel isn't available for your platform (build from CCBlade's
own source, which does need that toolchain). Expected as a sibling
directory, by default "../../ccblade-reference" relative to this file's
repo root (override with CCBLADE_REFERENCE_PATH), inside its own
virtualenv:

    cd src
    "<path-to-ccblade-reference>/.venv/Scripts/python.exe" -m bem.compare_ccblade

Fairness of the comparison
---------------------------
Same dense per-Reynolds-bucket (alpha, Cl, Cd) tables as compare_pybemt.py
(built once here from data/polars/s809/*.csv via the shared
build_dense_table/reynolds_bucket_for_station helpers, imported directly
from compare_pybemt -- not reimplemented) -- given to CCAirfoil as a
single-Reynolds-column table (RectBivariateSpline requires >=2 Re points
internally, so CCAirfoil duplicates a single-Re table itself; see its own
docstring). The remaining known difference is CCBlade's smoothed bivariate
spline interpolation (s=0.01 for Cl, s=0.001 for Cd) vs. our solver's
linear interpolation over the same table -- judged negligible at the same
0.5 deg grid resolution used for pyBEMT, for the same reason.

CCBlade defaults to shearExp=0.2 (wind shear across the rotor disk by
hub height) -- explicitly set to 0.0 here, since neither our solver nor
pyBEMT model shear; both assume uniform inflow equal to the given wind
speed.

Integration domain -- a second fairness issue found and fixed here
----------------------------------------------------------------------
CCBlade's own `evaluate()` integrates each station's load out to the
*physical* Rhub/Rtip boundaries (effectively assuming the load tapers to
zero there), not just between the first and last given station. Since our
geometry's stations start at r=1.2575 m while the real hub is at r=0.508 m
(the root/transition region is deliberately unmodeled -- no S809 polar
applies there, see bem.rotor.PHASE_VI_ROOT_EXCLUDED_NOTE), this let
CCBlade silently integrate over a wider span than our solver or pyBEMT do,
inflating its Cp/Ct by ~4-5% relative to ours for a reason that had
nothing to do with solver physics -- confirmed by comparing per-station
`a`, `alpha`, `Cl` (which matched closely) against the integrated totals
(which didn't), then by re-running with Rhub pulled in to just under the
first station, which closed most of the gap. Rather than move Rhub itself
(which is also, correctly, an input to the Prandtl hub-loss factor --
moving it would change the physics, not just the integration domain), this
script uses CCBlade's `distributedAeroLoads()` for per-station Np/Tp and
integrates them itself with the exact same trapezoidal rule as
bem.rotor._trapz, over exactly the same station range as our solver and
the pyBEMT config's `dr` -- see the Stage 6 journal entry.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

from bem.compare_pybemt import (  # noqa: E402
    AIR_DENSITY,
    AIR_DYNAMIC_VISCOSITY,
    DOCS_VALIDATION_DIR,
    build_dense_table,
    reynolds_bucket_for_station,
)
from bem.rotor import _trapz, phase_vi_geometry, solve_rotor  # noqa: E402

CCBLADE_REFERENCE_PATH = os.environ.get(
    "CCBLADE_REFERENCE_PATH",
    os.path.join(
        os.path.dirname(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))),
        "ccblade-reference",
    ),
)


class TableS809Polar:
    """Duplicated from compare_pybemt.py deliberately (not imported) -- keeping
    each solver-facing script's own-solver-side airfoil adapter local avoids a
    surprising cross-script coupling if one is ever edited independently. See
    compare_pybemt.TableS809Polar for the identical implementation/docstring."""

    def __init__(self, alpha_deg, cl, cd):
        self._alpha_deg = alpha_deg
        self._cl = cl
        self._cd = cd

    def cl(self, alpha_rad):
        return float(np.interp(math.degrees(alpha_rad), self._alpha_deg, self._cl))

    def cd(self, alpha_rad):
        return float(np.interp(math.degrees(alpha_rad), self._alpha_deg, self._cd))


def main():
    if not os.path.isdir(CCBLADE_REFERENCE_PATH):
        raise SystemExit(
            f"CCBlade reference environment not found at '{CCBLADE_REFERENCE_PATH}'. "
            "See docs/validation/bem-cross-validation.md 'How to reproduce' -- "
            "install it externally first, or set CCBLADE_REFERENCE_PATH."
        )
    try:
        from wisdem.ccblade.ccblade import CCAirfoil, CCBlade
    except ImportError as e:
        raise SystemExit(
            "Could not import wisdem.ccblade. This script must be run with "
            "the CCBlade reference environment's own virtualenv python (it "
            "is not a dependency of this repo) -- see docs/validation/"
            f"bem-cross-validation.md. Original error: {e}"
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
    our_airfoils = [TableS809Polar(*tables[b]) for b in station_buckets]

    ccblade_airfoils = []
    for b in station_buckets:
        alpha_deg, cl, cd = tables[b]
        ccblade_airfoils.append(CCAirfoil(alpha_deg, [b], cl, cd, AFName=f"S809_Re{b}"))

    twist_deg = [math.degrees(t) for t in geometry.twist]

    rotor = CCBlade(
        r=geometry.r,
        chord=geometry.chord,
        theta=twist_deg,
        af=ccblade_airfoils,
        Rhub=geometry.r_hub,
        Rtip=geometry.R,
        B=geometry.n_blades,
        rho=AIR_DENSITY,
        mu=AIR_DYNAMIC_VISCOSITY,
        precone=0.0,
        tilt=0.0,
        yaw=0.0,
        shearExp=0.0,  # no wind shear -- match our solver's/pyBEMT's uniform-inflow assumption
        tiploss=True,
        hubloss=True,
        wakerotation=True,
        usecd=True,
        iterRe=1,
    )

    wind_speeds = [5.0, 7.0, 10.0, 13.0, 15.0, 20.0, 25.0]  # matches Sequence S test points, compare_pybemt.py
    area = math.pi * geometry.R ** 2

    rows = []
    for v_inf in wind_speeds:
        tsr = omega * geometry.R / v_inf
        ours = solve_rotor(geometry, tsr=tsr, v_inf=v_inf, air_density=AIR_DENSITY, airfoils=our_airfoils)

        # distributedAeroLoads + our own trapezoidal integration, over exactly
        # the given station range -- not rotor.evaluate()'s own integrated
        # CP/CT, which extrapolates to Rhub/Rtip (see module docstring).
        loads, _derivs = rotor.distributedAeroLoads(v_inf, rpm, 0.0, 0.0)
        Np, Tp = list(loads["Np"]), list(loads["Tp"])
        Tp_r = [tp * r for tp, r in zip(Tp, geometry.r)]
        thrust = geometry.n_blades * _trapz(Np, geometry.r)
        torque = geometry.n_blades * _trapz(Tp_r, geometry.r)
        power = torque * omega
        cp_ccblade = power / (0.5 * AIR_DENSITY * v_inf ** 3 * area)
        ct_ccblade = thrust / (0.5 * AIR_DENSITY * v_inf ** 2 * area)

        rows.append({
            "v_inf": v_inf,
            "tsr": tsr,
            "Cp_ours": ours["Cp"],
            "Ct_ours": ours["Ct"],
            "Cp_ccblade": cp_ccblade,
            "Ct_ccblade": ct_ccblade,
        })

    os.makedirs(os.path.join(DOCS_VALIDATION_DIR, "ccblade_case"), exist_ok=True)
    results_path = os.path.join(DOCS_VALIDATION_DIR, "ccblade_case", "results.json")
    with open(results_path, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"Wrote results: {results_path}")

    return rows


if __name__ == "__main__":
    rows = main()
    print(f"\n{'v_inf':>6} {'TSR':>7} {'Cp_ours':>9} {'Cp_ccblade':>11} {'Ct_ours':>9} {'Ct_ccblade':>11}")
    for row in rows:
        print(f"{row['v_inf']:6.1f} {row['tsr']:7.3f} {row['Cp_ours']:9.4f} "
              f"{row['Cp_ccblade']:11.4f} {row['Ct_ours']:9.4f} {row['Ct_ccblade']:11.4f}")
    print(
        "\nNext: run plot_bem_comparison.py (plain repo python) for the combined "
        "three-way plot and deviation tables -- see docs/validation/bem-cross-validation.md."
    )
