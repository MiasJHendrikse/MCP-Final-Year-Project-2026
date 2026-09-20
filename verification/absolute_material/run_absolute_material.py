"""
The absolute material figures (2026-09-20, evening): the committed blades
`x0`, `x_c`, `x_m` in kilograms, megapascals and millimetres, the absolute
rows' Tier 3 at the new constants, and the measured answer to "does the
optimum move when the relative stress row is replaced by the absolute one".

What is run
-----------
1. `x0`, `x_c`, `x_m` evaluated, not optimised: shell and solid mass per
   blade and per rotor (`MaterialModel.mass_kg`, 1920 kg/m^3 x 2 mm skin),
   the thin-shell root stress `KS M_ref / (k_Z c0^2 t)` against the 196.5 MPa
   design allowable, the tip deflection `delta_ref D / (E k_I t)` in mm, and
   every mass-problem slack including the absolute rows.
2. Tier 3 (adjoint vs central FD at the committed `h*`, acceptance 3 eps_j)
   and the Taylor remainder for `stress_absolute` at `x0` and `x_m`, and for
   `deflection_absolute` at both points under an INJECTED 150 mm clearance
   (a test value: `tip_clearance_m` is TODO and the row is not assembled in
   any optimisation here). The Tier 1-2 checks are inherited unchanged: the
   rows read the same `RootMomentSystem` / `DeflectionSystem` results as the
   relative rows and differ by a constant.
3. Two optimisations at `delta = 0`:
   * `relative + stress_absolute`: the committed row set with the absolute
     stress row added. Warm from `x_m`. Expected: the same point (the row
     is slack there).
   * `stress -> stress_absolute`: the relative stress row REPLACED by the
     absolute one, cold from `x0` and warm from `x_m`. Expected: the
     "no stress" ablation optimum of `mass_optimisation/ablation.json`
     returns, because 196.5 MPa is never approached at the operating loads.

Inputs: `verification/baseline/x0.json`, `load_constraint/result_eps0.json`
(`x_c`), `mass_optimisation/result_delta0.json` (`x_m`),
`mass_optimisation/ablation.json` (the comparison), `fd_step_size/sweep.json`
(`h*`), `config/rotor_design.yaml::structure`. Output: `result.json`,
`absolute_material.png`. Shared code is `mass_optimisation/_common.py` and
the Tier 3 / Taylor helpers of `mass_optimisation/run_mass_checks.py`.

    python verification/absolute_material/run_absolute_material.py     # ~2 min

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import dataclasses
import datetime
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
MASS_DIR = os.path.join(REPO_ROOT, "verification", "mass_optimisation")
for path in (os.path.join(REPO_ROOT, "src"), MASS_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)

import _common as C  # noqa: E402
import run_mass_checks as checks  # noqa: E402
from gradients.problem import ABSOLUTE_ROWS, MASS_PROBLEM_ROWS  # noqa: E402

RESULT_PATH = os.path.join(_HERE, "result.json")
FIGURE_PATH = os.path.join(_HERE, "absolute_material.png")
XM_PATH = os.path.join(MASS_DIR, "result_delta0.json")
ABLATION_PATH = os.path.join(MASS_DIR, "ablation.json")

#: The clearance the deflection row's Tier 3 is run under. A TEST VALUE so
#: the row's Jacobian can be checked; not a design assumption, not used in
#: any optimisation, and recorded as such in `result.json`.
INJECTED_CLEARANCE_M = 0.150

#: Rows of the two optimisation runs.
ROWS_ADDED = MASS_PROBLEM_ROWS + ("stress_absolute",)
ROWS_SWAPPED = tuple(r for r in MASS_PROBLEM_ROWS if r != "stress") + ("stress_absolute",)


def blade(problem, label, u, delta=0.0):
    """`C.blade_record` plus the absolute figures in engineering units."""

    record = C.blade_record(problem, u, delta, label)
    material = record["material"]
    n_blades = int(problem.design.n_blades)
    slacks = record["slacks"]
    record["absolute"] = {
        "shell_mass_kg_per_blade": material["shell_mass_kg"],
        "solid_mass_kg_per_blade": material["solid_mass_kg"],
        "shell_mass_kg_per_rotor": material["shell_mass_kg"] * n_blades,
        "solid_mass_kg_per_rotor": material["solid_mass_kg"] * n_blades,
        "root_stress_mpa": slacks["stress_absolute"]["stress_mpa"],
        "design_allowable_mpa": slacks["stress_absolute"]["allowable_mpa"],
        "stress_over_allowable": slacks["stress_absolute"]["stress_mpa"]
        / slacks["stress_absolute"]["allowable_mpa"],
        "tip_deflection_mm": slacks["deflection_absolute"]["tip_deflection_mm"],
        "tip_deflection_over_radius": slacks["deflection_absolute"]["tip_deflection_mm"]
        / (1e3 * float(problem.design.radius_m)),
        "tip_clearance": slacks["deflection_absolute"]["available"],
    }
    return record


def verify_rows(problem, points, h_star, names, clearance_problem):
    """Tier 3 and Taylor for the absolute rows at every point in `points`."""

    out = {}
    stress = problem.stress_constraint_absolute()
    deflection = clearance_problem.deflection_constraint_absolute()
    for label, u in points.items():
        print(f"  Tier 3 / Taylor at {label} ...", flush=True)
        out[label] = {
            "stress_absolute": {
                "tier3": checks.tier3(problem, stress, u, h_star, names),
                "taylor": checks.taylor(problem, stress, u, seed=31),
            },
            "deflection_absolute": {
                "injected_clearance_m": INJECTED_CLEARANCE_M,
                "tier3": checks.tier3(clearance_problem, deflection, u, h_star, names),
                "taylor": checks.taylor(clearance_problem, deflection, u, seed=37),
            },
        }
    return out


def optimise(problem, u0, u_start, include, label, names):
    """One SLSQP solve of the mass problem at `delta = 0` over `include`."""

    print(f"\n== {label}: rows {list(include)} ==", flush=True)
    t0 = time.perf_counter()
    result, recorder, wall = C.run_mass_slsqp(problem, u_start, 0.0, include, verbose=True)
    u_opt = np.asarray(result.x, dtype=float)
    record = blade(problem, label, u_opt)
    kkt = C.kkt_report(problem, u_opt, 0.0, include, names)
    return {
        "label": label,
        "rows": list(include),
        "success": bool(result.success), "status": int(result.status),
        "message": str(result.message), "nit": int(result.nit),
        "wall_time_s": float(time.perf_counter() - t0),
        "u": [float(x) for x in u_opt],
        "x": [float(x) for x in problem.physical(u_opt)],
        "control_points": C.control_points_record(problem, problem.physical(u_opt)),
        "shell_pct_vs_x0": record["shell_pct_vs_x0"],
        "solid_pct_vs_x0": record["solid_pct_vs_x0"],
        "aep_pct_vs_x0": record["aep_pct_vs_x0"],
        "absolute": record["absolute"],
        "stress_ratio": record["stress_ratio"],
        "deflection_ratio": record["deflection_ratio"],
        "KS_over_KS0": record["KS_over_KS0"],
        "slacks": record["slacks"],
        "kkt": kkt,
        "n_evaluation_failures": len(problem.evaluation_failures),
    }


def plot(blades, path):
    """
    Two panels, one axis each: the masses in kg (shell and solid, per blade)
    and the thin-shell root stress against its allowable. The deflection is
    a third panel in mm with no allowable line (there is none yet). House
    colours (the `mass_optimisation` figures'); every bar carries its value.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [b["label"] for b in blades]
    x = np.arange(len(labels))
    blue, light, ink, grid = "#1f5fbf", "#5aa9e6", "#333333", "#dddddd"
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))

    ax = axes[0]
    shell = [b["absolute"]["shell_mass_kg_per_blade"] for b in blades]
    solid = [b["absolute"]["solid_mass_kg_per_blade"] for b in blades]
    bars = ax.bar(x - 0.21, shell, 0.4, color=blue, label="shell skin, 2 mm E-glass/epoxy")
    bars2 = ax.bar(x + 0.21, solid, 0.4, color=light, label="solid section, same laminate")
    for group in (bars, bars2):
        for bar in group:
            ax.annotate(f"{bar.get_height():.2f}", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                        ha="center", va="bottom", fontsize=8, color=ink, xytext=(0, 2),
                        textcoords="offset points")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("mass per blade [kg]")
    ax.set_title("material in kilograms (1920 kg/m^3)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.set_ylim(0.0, max(solid) * 1.25)
    ax.grid(True, axis="y", color=grid, lw=0.6)

    ax = axes[1]
    stress = [b["absolute"]["root_stress_mpa"] for b in blades]
    allowable = blades[0]["absolute"]["design_allowable_mpa"]
    bars = ax.bar(x, stress, 0.5, color=blue)
    for bar in bars:
        ax.annotate(f"{bar.get_height():.1f}", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    ha="center", va="bottom", fontsize=8, color=ink, xytext=(0, 2),
                    textcoords="offset points")
    ax.axhline(allowable, color=ink, lw=1.0, ls="--")
    ax.annotate(f"design allowable {allowable:.1f} MPa (702 / 3.5725)", (0.03, allowable),
                xycoords=("axes fraction", "data"), fontsize=8, color=ink, va="bottom",
                xytext=(0, 2), textcoords="offset points")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("thin-shell root stress [MPa]")
    ax.set_title("KS root moment over k_Z c0^2 t, operating loads", fontsize=10)
    ax.set_ylim(0.0, allowable * 1.15)
    ax.grid(True, axis="y", color=grid, lw=0.6)

    ax = axes[2]
    deflection = [b["absolute"]["tip_deflection_mm"] for b in blades]
    bars = ax.bar(x, deflection, 0.5, color=blue)
    for bar in bars:
        ax.annotate(f"{bar.get_height():.1f}", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    ha="center", va="bottom", fontsize=8, color=ink, xytext=(0, 2),
                    textcoords="offset points")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("tip deflection [mm]")
    ax.set_title("delta_ref D / (E k_I t), E = 41.8 GPa; no clearance specified", fontsize=10)
    ax.set_ylim(0.0, max(deflection) * 1.25)
    ax.grid(True, axis="y", color=grid, lw=0.6)

    fig.suptitle("x0 / x_c / x_m with the recorded laminate -- thin-shell blade, "
                 "E-LT-5500/EP-3 (Griffith & Ashwill 2011, Table 19)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    started = time.perf_counter()
    h_star = float(C.load_json(C.SWEEP_PATH)["h_star_global"])
    problem, u0 = C.prepared_problem()
    names = C.variable_names(problem)
    u_c = problem.scaled(C.load_xc())
    artefact = C.load_json(XM_PATH)
    u_m = np.array(artefact["u_m"], dtype=float)
    material = problem.material_model()

    print("== the committed blades in engineering units ==", flush=True)
    blades = [blade(problem, "x0", u0), blade(problem, "x_c", u_c), blade(problem, "x_m", u_m)]
    for b in blades:
        a = b["absolute"]
        print(f"  {b['label']:3s}  shell {a['shell_mass_kg_per_blade']:.4f} kg  solid "
              f"{a['solid_mass_kg_per_blade']:.4f} kg  sigma {a['root_stress_mpa']:.2f} MPa  "
              f"delta {a['tip_deflection_mm']:.1f} mm", flush=True)

    # A second problem carrying the injected clearance, for the deflection
    # row's Tier 3 only.
    clearance_problem = C.build_problem()
    clearance_problem.design = dataclasses.replace(clearance_problem.design,
                                                   tip_clearance_m=INJECTED_CLEARANCE_M)
    clearance_problem.set_reference(clearance_problem.scaled(clearance_problem.reference_design()))
    print("\n== Tier 3 and Taylor for the absolute rows ==", flush=True)
    verification = verify_rows(problem, {"x0": u0, "x_m": u_m}, h_star, names, clearance_problem)
    for label, rows in verification.items():
        for row, r in rows.items():
            print(f"  {label} {row}: worst |adj - fd| / eps = "
                  f"{r['tier3']['worst_abs_error_over_eps']:.3f} ({r['tier3']['worst_variable']}), "
                  f"Taylor min ratio {r['taylor']['min_ratio']:.1f}", flush=True)

    runs = {
        "added_warm_from_x_m": optimise(problem, u0, u_m, ROWS_ADDED,
                                        "relative + stress_absolute, warm from x_m", names),
        "swapped_cold_from_x0": optimise(problem, u0, u0, ROWS_SWAPPED,
                                         "stress -> stress_absolute, cold from x0", names),
        "swapped_warm_from_x_m": optimise(problem, u0, u_m, ROWS_SWAPPED,
                                          "stress -> stress_absolute, warm from x_m", names),
    }

    ablation = {c["label"]: c for c in C.load_json(ABLATION_PATH)["cases"]}
    no_stress = ablation["no stress"]
    swapped = runs["swapped_cold_from_x0"]
    comparison = {
        "added_row_moves_x_m_by_u_inf": float(np.max(np.abs(np.array(runs["added_warm_from_x_m"]["u"]) - u_m))),
        "added_row_shell_pct": runs["added_warm_from_x_m"]["shell_pct_vs_x0"],
        "committed_x_m_shell_pct": float(artefact["shell_pct_vs_x0"]),
        "swapped_cold_vs_warm_u_inf": float(np.max(np.abs(np.array(swapped["u"])
                                                          - np.array(runs["swapped_warm_from_x_m"]["u"])))),
        "swapped_vs_no_stress_ablation_u_inf": float(np.max(np.abs(np.array(swapped["u"])
                                                                   - np.array(no_stress["u_m"])))),
        "swapped_shell_pct": swapped["shell_pct_vs_x0"],
        "no_stress_ablation_shell_pct": float(no_stress["shell_pct_vs_x0"]),
        "swapped_stress_ratio_vs_x0": swapped["stress_ratio"],
        "swapped_root_stress_mpa": swapped["absolute"]["root_stress_mpa"],
        "swapped_moves_x_m_by_u_inf": float(np.max(np.abs(np.array(swapped["u"]) - u_m))),
        "agreement_criterion_u_inf": float(C.load_json(C.MULTISTART_PATH)["spread_of_optima_u_inf"])
        if "spread_of_optima_u_inf" in C.load_json(C.MULTISTART_PATH) else None,
    }

    out = {
        "problem": "the committed mass problem with the recorded laminate: masses in kg, "
                   "the thin-shell root stress and tip deflection in MPa and mm, the absolute "
                   "rows' Tier 3, and the optimum with the absolute stress row added / swapped in",
        "bounds_label": C.BOUNDS_LABEL,
        "law_label": C.LAW_LABEL,
        "command": "python verification/absolute_material/run_absolute_material.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "structural_inputs": material.structural_inputs(),
        "design_allowable_stress_mpa": problem.design_allowable_stress_pa() / 1e6,
        "stiffness_factor_E_kI_t_pa_m": material.stiffness_factor_pa_m(),
        "section_modulus_factor_kZ_t_m": material.section_modulus_factor_m(),
        "absolute_rows_available": problem.absolute_rows_available(),
        "references": {
            "m_ref_nm": float(problem.m_ref), "KS0": float(problem.KS0),
            "delta_ref_per_unit_stiffness": float(problem.delta_ref), "D0": float(problem.D0),
            "material_ref_shell_m2": float(problem.material_ref),
            "rated_tip_deflection_x0_mm": float(problem.delta_ref) / material.stiffness_factor_pa_m() * 1e3,
        },
        "blades": blades,
        "verification": verification,
        "all_tier3_pass": bool(all(r[row]["tier3"]["passes"] for r in verification.values()
                                   for row in ("stress_absolute", "deflection_absolute"))),
        "all_taylor_pass": bool(all(r[row]["taylor"]["passes"] for r in verification.values()
                                    for row in ("stress_absolute", "deflection_absolute"))),
        "runs": runs,
        "comparison": comparison,
        "rows_added": list(ROWS_ADDED),
        "rows_swapped": list(ROWS_SWAPPED),
        "absolute_rows": list(ABSOLUTE_ROWS),
        "evaluation_failures": problem.evaluation_failures,
        "wall_time_s": float(time.perf_counter() - started),
    }
    with open(RESULT_PATH, "w", encoding="utf-8") as handle:
        json.dump(out, handle, indent=1)
    plot(blades, FIGURE_PATH)
    print(f"\nwrote {RESULT_PATH} and {FIGURE_PATH} in {out['wall_time_s']:.0f} s", flush=True)
    print(json.dumps(comparison, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
