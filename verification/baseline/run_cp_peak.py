"""
Where does the Schmitz reference's power coefficient peak, and how exact is it?

The report called `x0` "the exact
maximiser of the power coefficient at the design tip-speed ratio" and the
peak's position at lambda = 6.5 "the strongest available check" on its
construction. Schmitz's closed form leaves out tip loss and treats drag and
the Reynolds number only through the one design polar point, while the BEM
model includes all three, so neither statement can be exact. This script
measures how far from exact they are: the tip-speed ratio at which `x0`'s
C_P peaks along a line of fixed tip-speed ratio, at several wind speeds in
the site atmosphere (the Reynolds number rises with wind speed).

Outputs, next to this script: `cp_peak.json`.

    python verification/baseline/run_cp_peak.py     # ~20 s

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import datetime
import json
import os
import sys

from scipy.optimize import minimize_scalar

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
MASS_DIR = os.path.join(REPO_ROOT, "verification", "mass_optimisation")
for path in (os.path.join(REPO_ROOT, "src"), MASS_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)

import _common as C  # noqa: E402
from bem.rotor import solve_rotor  # noqa: E402

OUT_PATH = os.path.join(_HERE, "cp_peak.json")
WIND_SPEEDS_MS = (5.0, 7.0, 9.0, 11.0)


def main():
    problem, _ = C.prepared_problem()
    geometry = problem.parameterisation.to_geometry(C.load_x0(), polar_cache=problem.polar_cache)
    rho = float(problem.site.air_density)
    nu = float(problem.site.kinematic_viscosity)
    design = float(problem.design.design_tsr)
    rows = []
    for v in WIND_SPEEDS_MS:
        def negative_cp(lam, v=v):
            return -solve_rotor(geometry, tsr=lam, v_inf=v, air_density=rho,
                                kinematic_viscosity=nu)["Cp"]
        best = minimize_scalar(negative_cp, bounds=(5.5, 7.5), method="bounded",
                               options={"xatol": 1e-4})
        rows.append({"v_ms": v, "peak_tsr": float(best.x), "peak_cp": float(-best.fun),
                     "cp_at_design_tsr": float(-negative_cp(design)),
                     "cp_gain_at_peak_pct": 100.0 * (float(-best.fun) / float(-negative_cp(design)) - 1.0)})
        print(f"V = {v:4.1f} m/s: C_P peaks at lambda = {best.x:.3f} ({-best.fun:.4f}); "
              f"C_P(6.5) = {rows[-1]['cp_at_design_tsr']:.4f}", flush=True)
    record = {
        "description": __doc__.split("\n\n")[0].strip(),
        "command": "python verification/baseline/run_cp_peak.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "design_tsr": design,
        "note": "x0 along fixed-lambda lines in the site atmosphere; bounded scalar search, xatol 1e-4",
        "rows": rows,
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
