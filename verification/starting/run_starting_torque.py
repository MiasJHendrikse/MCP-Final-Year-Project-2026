"""
Starting torque of the three committed blades, relative to the reference.

The design takes the
3 m/s cut-in as a machine requirement and never asks whether the
minimum-material blade starts as readily as the Schmitz reference. Small
rotors start slowly because at rest every section meets the wind at
alpha = 90 deg - theta(r), deep in stall, and the aerodynamic torque there
comes mostly from the inboard, highly twisted span (Wood 2011, the report's
[3]). This script evaluates that parked torque for `x0`, `x_c` and `x_m`:

    Q_0 = B int_{r_hub}^{R} (1/2) rho V^2 c(r) C_l(90 deg - theta(r)) r dr

(phi = 90 deg with the rotor at rest and no induction; the drag term drops
because cos phi = 0), trapezoid rule on the 25 BEM stations.

What it can and cannot say. The torque is a ratio between blades, not a
prediction: every alpha here lies in 66..90 deg, outside the XFOIL band, so
the coefficients are the cache's Viterna extension, and the station Reynolds
numbers at 3 m/s (about 10 000 to 45 000) are below the cache's 40 000 floor,
so the coefficients are read on the 40 000 row. Viterna is almost independent
of Reynolds number at these angles; both simplifications apply to the three
blades alike. Starting also depends on the generator's cogging and friction
torque and on the rotor inertia, neither of which is specified; the inertia
is compared through the shell mass, and the ratio of torque to that mass is
reported as a rough figure of merit for the acceleration.

Outputs, next to this script: `starting_torque.json`.

    python verification/starting/run_starting_torque.py     # ~2 s

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import datetime
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
MASS_DIR = os.path.join(REPO_ROOT, "verification", "mass_optimisation")
for path in (os.path.join(REPO_ROOT, "src"), MASS_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)

import _common as C  # noqa: E402
from polars.polar import CachedPolar, interpolant_for  # noqa: E402

XM_PATH = os.path.join(MASS_DIR, "result_delta0.json")
OUT_PATH = os.path.join(_HERE, "starting_torque.json")


def trapz(y, x):
    return float(sum(0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i]) for i in range(len(x) - 1)))


def parked_torque(problem, d, v, reynolds_row):
    """Rotor torque at rest [N m] and its spanwise integrand, one wind speed."""

    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    polar = CachedPolar(interpolant_for(problem.polar_cache), reynolds_row)
    rho = float(problem.site.air_density)
    nu = float(problem.site.kinematic_viscosity)
    r = np.array(geometry.r, dtype=float)
    c = np.array(geometry.chord, dtype=float)
    theta = np.array(geometry.twist, dtype=float)
    alpha = 0.5 * math.pi - theta
    cl = np.array([polar.cl(a) for a in alpha])
    dq = 0.5 * rho * v ** 2 * c * cl * r
    torque = geometry.n_blades * trapz(dq, r)
    return {
        "torque_nm": torque,
        "alpha_deg_range": [float(np.degrees(alpha).min()), float(np.degrees(alpha).max())],
        "station_reynolds_range": [float((v * c / nu).min()), float((v * c / nu).max())],
        "inboard_half_share": trapz(dq[r <= 0.5 * (r[0] + r[-1])], r[r <= 0.5 * (r[0] + r[-1])])
        / trapz(dq, r),
    }


def main():
    problem, _ = C.prepared_problem()
    material = problem.material_model()
    v = float(problem.design.cut_in_wind_speed_ms)
    reynolds_row = float(interpolant_for(problem.polar_cache).re_values[0])
    blades = {"x0": C.load_x0(), "x_c": C.load_xc(),
              "x_m": np.array(C.load_json(XM_PATH)["x_m"], dtype=float)}
    records = {}
    for name, d in blades.items():
        rec = parked_torque(problem, d, v, reynolds_row)
        rec["shell_mass_kg_per_blade"] = float(material.mass_kg(d))
        rec["torque_per_rotor_shell_mass_nm_per_kg"] = rec["torque_nm"] / (
            geometry_blades(problem, d) * rec["shell_mass_kg_per_blade"])
        records[name] = rec
    q0 = records["x0"]["torque_nm"]
    m0 = records["x0"]["torque_per_rotor_shell_mass_nm_per_kg"]
    for name, rec in records.items():
        rec["torque_vs_x0_pct"] = 100.0 * (rec["torque_nm"] / q0 - 1.0)
        rec["torque_per_mass_vs_x0_pct"] = 100.0 * (rec["torque_per_rotor_shell_mass_nm_per_kg"] / m0 - 1.0)
        print(f"{name:3s}: parked torque at {v:g} m/s {rec['torque_nm']:.3f} N m "
              f"({rec['torque_vs_x0_pct']:+.2f} % vs x0); per unit shell mass "
              f"{rec['torque_per_mass_vs_x0_pct']:+.2f} %; alpha {np.round(rec['alpha_deg_range'], 1)} deg; "
              f"inboard-half share {rec['inboard_half_share']:.2f}")
    record = {
        "description": __doc__.split("\n\n")[0].strip(),
        "command": "python verification/starting/run_starting_torque.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "wind_speed_ms": v,
        "reynolds_row_read": reynolds_row,
        "caveats": ("ratio between blades only; coefficients from the cache's Viterna "
                    "extension on its 40 000 row (station Re at cut-in is below the cache "
                    "floor); no generator resistance or inertia model"),
        "blades": records,
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"wrote {OUT_PATH}")


def geometry_blades(problem, d):
    return problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache).n_blades


if __name__ == "__main__":
    main()
