"""
The thin-shell bending stress along the span, not only at the root.

The stress row of the mass problem
constrains `KS / c0^2`, the thin-shell stress at the root section, and
`run_absolute_material.py` turns that one section into megapascals. This
script evaluates the same thin-shell stress at every BEM station of the three
committed blades at the rated point,

    sigma(r) = M(r) / (k_Z t c(r)^2),

with `M(r)` the flapwise moment about the station from the outboard load
alone (`objective.loads.spanwise_moments`, the quantity of the report's
Table 9.1), and reports where along the span it peaks and by how much the
peak exceeds the root value. Nothing is optimised and no row changes.

Inputs: `verification/baseline/x0.json`, `load_constraint/result_eps0.json`
(`x_c`), `mass_optimisation/result_delta0.json` (`x_m`),
`config/rotor_design.yaml::structure` (`k_Z`, `t`). Output:
`spanwise_stress.json`, next to this script.

    python verification/absolute_material/run_spanwise_stress.py     # ~5 s

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import datetime
import json
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

XM_PATH = os.path.join(MASS_DIR, "result_delta0.json")
OUT_PATH = os.path.join(_HERE, "spanwise_stress.json")


def blade(problem, label, d, kz_t):
    """Rated-point spanwise thin-shell stress of one blade, MPa."""

    loads = C.loads_at(problem, d, with_spanwise=True)
    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    radii = np.array(loads["radii_m"])
    moments = np.array(loads["spanwise_moment_nm"])
    chords = np.array(loads["chord_m"])
    sigma = moments / (kz_t * chords ** 2) / 1e6
    c_root = float(d[0])  # the clamped spline passes through its first control point
    sigma_root = loads["root_bending_moment_nm"] / (kz_t * c_root ** 2) / 1e6
    k = int(np.argmax(sigma))
    return {
        "label": label,
        "rated_point": {"v_ms": loads["v_ms"], "tsr": loads["tsr"], "rpm": loads["rpm"]},
        "root": {"r_m": float(geometry.r_hub),
                 "chord_mm": 1e3 * c_root,
                 "moment_nm": loads["root_bending_moment_nm"],
                 "stress_mpa": float(sigma_root)},
        "stations": {"r_m": radii.tolist(), "chord_mm": (1e3 * chords).tolist(),
                     "moment_nm": moments.tolist(), "stress_mpa": sigma.tolist()},
        "peak": {"station": k, "r_m": float(radii[k]), "r_over_R": float(radii[k] / radii[-1]),
                 "chord_mm": float(1e3 * chords[k]), "moment_nm": float(moments[k]),
                 "stress_mpa": float(sigma[k]),
                 "over_root": float(sigma[k] / sigma_root)},
    }


def main():
    problem, _ = C.prepared_problem()
    material = problem.material_model()
    kz_t = float(material.section_modulus_factor_m())
    blades = {"x0": C.load_x0(), "x_c": C.load_xc(),
              "x_m": np.array(C.load_json(XM_PATH)["x_m"], dtype=float)}
    records = {name: blade(problem, name, d, kz_t) for name, d in blades.items()}
    for name, rec in records.items():
        p = rec["peak"]
        print(f"{name:3s}: root {rec['root']['stress_mpa']:.2f} MPa; peak {p['stress_mpa']:.2f} MPa "
              f"at r = {p['r_m']:.3f} m (r/R {p['r_over_R']:.2f}), {p['over_root']:.2f} x root")
    peak_x0 = records["x0"]["peak"]["stress_mpa"]
    record = {
        "description": __doc__.split("\n\n")[0].strip(),
        "command": "python verification/absolute_material/run_spanwise_stress.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "section_modulus_factor_kZ_t_m": kz_t,
        "note": ("sigma(r) = M(r) / (k_Z t c(r)^2) at the 25 BEM stations, rated point; "
                 "M(r) about the station from the outboard load; the root value is about "
                 "r_hub with the root control-point chord. Station values only: the peak "
                 "between stations is not interpolated."),
        "blades": records,
        "x_m_peak_vs_x0_peak_pct": 100.0 * (records["x_m"]["peak"]["stress_mpa"] / peak_x0 - 1.0),
        "x_c_peak_vs_x0_peak_pct": 100.0 * (records["x_c"]["peak"]["stress_mpa"] / peak_x0 - 1.0),
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"x_m peak vs x0 peak: {record['x_m_peak_vs_x0_peak_pct']:+.2f} %; "
          f"x_c: {record['x_c_peak_vs_x0_peak_pct']:+.2f} %")
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
