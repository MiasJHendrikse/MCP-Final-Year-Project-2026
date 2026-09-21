"""
Phase 5 (2026-09-20) -- the checks behind the mass optimum, and the
reference-blades table.

At `x0` and at the production optimum `x_m` (`result_delta0.json`):

  * the material objective: its exact gradient against central FD (the
    objective is polynomial in `d`, so the agreement is round-off), and the
    solid proxy's gradient likewise;
  * the manufacturing rows: the difference matrix reproduces the control-
    point differences to the bit, and `x0` is strictly interior;
  * forward twins: the deflection system's per-point deflections against the
    forward path (`objective.loads.tip_deflection_at`), the moment system's
    against `root_moment`;
  * Tier 1 -- `dD/dphi`, `dD/dd` (deflection functional) and `dKS/dx`,
    `dKS/dd` (moment functional) against a complex step of the value code;
  * Tier 2 -- the deflection tangent against the adjoint direction, and the
    adjoint assembly identity;
  * Tier 3 -- the **stress** and **deflection** constraint Jacobians against
    central FD of their `fun` at the committed `h*`, with `eps_j` re-measured
    as the three-step jitter and the round-off floor reported on failure
    (`verification/load_constraint/run_load_checks.py`'s rule, verbatim);
  * Taylor remainders of both rows, second order in three random directions.

Also: the **reference-blades table** `x0 / x_c / x_m` (AEP, both material
proxies, KS, rated-point moment, thrust and deflection, stress and
deflection ratios, cut-in `reynolds_min`, control points) -- where the
"+6.2 % shell for +0.15 % energy" motivation and the "-4.2 % shell for
0 % energy" result both come from -- and the re-evaluation of `x_m` from
its JSON through the forward path only, confirming the recorded AEP, root
moment and tip deflection to 1e-10 (plan, Verification 5).

Outputs, next to this script: `checks.json`, `reference_blades.json`,
`reference_blades.png`, `blades_rendered.png` and `blades_rendered_edge.png`
(the three blades drawn with a cosmetic root cylinder and tip rounding: both
are drawn only, never modelled).

Run from the repo root (about three minutes):

    python verification/mass_optimisation/run_mass_checks.py
    python verification/mass_optimisation/run_mass_checks.py --replot

`--replot` redraws `reference_blades.png` from the committed
`reference_blades.json`; no check is run and no JSON is written.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import datetime
import json
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import _common as C  # noqa: E402
from gradients import central_difference  # noqa: E402
from objective.loads import load_operating_points, root_moment, tip_deflection_at  # noqa: E402
from objective.objective import objective  # noqa: E402
from objective.power import tsr_schedule  # noqa: E402

RESULT_PATH = os.path.join(_HERE, "result_delta0.json")
CHECKS_PATH = os.path.join(_HERE, "checks.json")
TABLE_PATH = os.path.join(_HERE, "reference_blades.json")
FIGURE_PATH = os.path.join(_HERE, "reference_blades.png")

H = 1e-30
SQRT10 = math.sqrt(10.0)
EPS_FACTOR = 3.0
TAYLOR_STEPS = (1e-2, 1e-3, 1e-4)
TAYLOR_SEEDS = (0, 1, 2)
FD_H_MASS = 1e-6


def mixed_error(estimate, partial):
    err = np.abs(np.asarray(estimate) - np.asarray(partial))
    return float(np.max(err / np.maximum(1.0, np.abs(partial))))


# ---------------------------------------------------------------------------
# the objective and the linear rows
# ---------------------------------------------------------------------------

def mass_gradient_check(problem, u):
    model = problem.material_model()
    d = problem.physical(u)
    out = {}
    for name, power in (("shell", 1), ("solid", 2)):
        exact = model.planform_integral_gradient(d, power)
        fd = central_difference(lambda w: model.planform_integral(w, power), d, FD_H_MASS)[0]
        out[name] = {"max_abs_error": float(np.max(np.abs(exact - fd))),
                     "max_rel_error": float(np.max(np.abs(exact - fd))
                                            / max(1.0, np.max(np.abs(exact)))),
                     "twist_block_is_zero": bool(np.all(exact[problem.parameterisation.n_chord:] == 0.0))}
    jac = problem.mass_jac(u)
    fd = central_difference(problem.mass, u, FD_H_MASS)[0]
    out["scaled_objective"] = {"max_abs_error": float(np.max(np.abs(jac - fd))),
                               "constant_for_shell": bool(np.array_equal(jac, problem.mass_jac(problem.u0)))}
    return out


def manufacturing_check(problem, u):
    matrix, labels = problem.manufacturing_rows()
    d = problem.physical(u)
    n_c = problem.parameterisation.n_chord
    expected = np.concatenate([d[:n_c - 1] - d[1:n_c], d[n_c:-1] - d[n_c + 1:]])
    values = problem.manufacturing_constraints()["fun"](u)
    return {"n_rows": int(matrix.shape[0]), "labels": labels,
            "values": [float(x) for x in values],
            "exact_to_the_bit": bool(np.array_equal(matrix @ d, expected)),
            "min_slack": float(values.min()), "strictly_interior": bool(values.min() > 0.0)}


# ---------------------------------------------------------------------------
# forward twins and Tiers 1-2
# ---------------------------------------------------------------------------

def forward_twins(problem, u):
    system = problem._deflection_system()
    d = problem.physical(u)
    state = system.solve(d)
    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    air = (float(problem.site.air_density), float(problem.site.kinematic_viscosity))
    points = load_operating_points()
    forward_delta = np.array([tip_deflection_at(geometry, v, lam, *air)[0] for v, lam in points])
    forward_moment = np.array([root_moment(geometry, v, lam, *air)[0] for v, lam in points])
    system_delta = system.deflections(state.phi, state.d)
    system_moment = system.moments_from_m(system.partials(state.phi, state.d, derivatives=False)["m"])
    return {
        "deflection_max_rel_error": float(np.max(np.abs(system_delta / forward_delta - 1.0))),
        "moment_max_rel_error": float(np.max(np.abs(system_moment / forward_moment - 1.0))),
        "deflections_per_point": [float(x) for x in system_delta],
        "moments_per_point_nm": [float(x) for x in system_moment],
    }, system, state


def tier1(system, state):
    parts = system.partials(state.phi, state.d)
    out = {}
    for name, value, partial_x, partial_d in (
        ("deflection", system.D, system.dD_dx, system.dD_dd),
        ("moment", system.KS, system.dKS_dx, system.dKS_dd),
    ):
        dx = partial_x(state.phi, state.d, parts)
        estimate = np.empty_like(dx)
        for b in range(system.n_points):
            for i in range(system.n_stations):
                phi = state.phi.astype(complex)
                phi[b, i] += 1j * H
                estimate[b, i] = value(phi, state.d).imag / H
        dd = partial_d(state.phi, state.d, parts)
        dd_estimate = np.empty_like(dd)
        for j in range(system.n_design):
            d = state.d.astype(complex)
            d[j] += 1j * H
            dd_estimate[j] = value(state.phi, d).imag / H
        out[name] = {"d_dphi": mixed_error(estimate, dx), "d_dd_explicit": mixed_error(dd_estimate, dd)}
    return out, parts


def tier2(system, state, parts):
    result = system.gradient(state.d, state=state)
    assembled = result.dD_dd_explicit + system.apply_dR_dd_T(parts, result.psi)
    identity = float(np.max(np.abs(assembled - result.dD_dd)) / max(1.0, np.max(np.abs(result.dD_dd))))
    rng = np.random.default_rng(42)
    worst = 0.0
    for _ in range(3):
        v = rng.normal(size=system.n_design)
        tangent = system.tangent(state.phi, state.d, v)
        adjoint = float(result.dD_dd @ v)
        worst = max(worst, abs(tangent - adjoint) / max(1.0, abs(adjoint)))
    return {"tangent_vs_adjoint": float(worst), "assembly_identity": identity,
            "D": float(result.D), "rated_point_weight": float(result.weights[-1])}


# ---------------------------------------------------------------------------
# Tier 3 and Taylor, the committed rule verbatim
# ---------------------------------------------------------------------------

def roundoff_floor(g, u, h_star):
    ts = np.linspace(-4e-12, 4e-12, 9)
    design = np.vstack([ts, np.ones_like(ts)]).T
    floors = np.empty(len(u))
    for j in range(len(u)):
        e = np.zeros(len(u))
        e[j] = 1.0
        values = np.array([g(u + t * e) for t in ts])
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        residual = values - design @ coefficients
        floors[j] = float(np.sqrt(np.mean(residual ** 2))) / h_star
    return floors


def tier3(problem, constraint, u, h_star, names):
    jac = constraint["jac"](u)[0]
    g = lambda w: float(constraint["fun"](w)[0])  # noqa: E731

    def fd(h):
        return central_difference(g, u, h)[0]

    g_mid, g_lo, g_hi = fd(h_star), fd(h_star / SQRT10), fd(h_star * SQRT10)
    eps = np.maximum(np.abs(g_mid - g_lo), np.abs(g_mid - g_hi))
    abs_err = np.abs(jac - g_mid)
    ratio = abs_err / eps
    within = abs_err <= EPS_FACTOR * eps
    out = {
        "h_star": float(h_star),
        "adjoint": [float(x) for x in jac],
        "fd_at_h_star": [float(x) for x in g_mid],
        "eps_per_u": [float(x) for x in eps],
        "abs_error": [float(x) for x in abs_err],
        "abs_error_over_eps": [float(x) for x in ratio],
        "within_3eps": [bool(x) for x in within],
        "worst_variable": names[int(np.argmax(ratio))],
        "worst_abs_error_over_eps": float(np.max(ratio)),
        "passes": bool(np.all(within)),
    }
    if not out["passes"]:
        # The floor is measured only on failure (the committed rule): it
        # describes the failure, it does not pass it.
        fit_floor = roundoff_floor(g, u, h_star)
        eps_below = eps < fit_floor
        out["measured_roundoff_floor"] = [float(x) for x in fit_floor]
        out["eps_below_roundoff_floor"] = [bool(x) for x in eps_below]
        out["floor_limited"] = bool(np.all(within | eps_below))
    else:
        out["floor_limited"] = False
    return out


def taylor(problem, constraint, u, seed):
    jac = constraint["jac"](u)[0]
    g = lambda w: float(constraint["fun"](w)[0])  # noqa: E731
    envelope = problem.envelope_constraint()["fun"]
    rng = np.random.default_rng(seed)
    for _attempt in range(50):
        v = rng.normal(size=problem.n)
        v /= np.linalg.norm(v)
        if all(np.all(envelope(np.clip(u + h * v, 0.0, 1.0)) > 0.0) for h in TAYLOR_STEPS) \
                and np.all(u + max(TAYLOR_STEPS) * v >= 0.0) and np.all(u + max(TAYLOR_STEPS) * v <= 1.0):
            break
    else:
        raise RuntimeError("no feasible direction found for the Taylor test")
    g0 = g(u)
    slope = float(jac @ v)
    remainders = [abs(g(u + h * v) - g0 - h * slope) for h in TAYLOR_STEPS]
    ratios = [remainders[k] / remainders[k + 1] for k in range(len(TAYLOR_STEPS) - 1)]
    return {"seed": seed, "remainder": [float(r) for r in remainders],
            "ratios_per_decade": [float(r) for r in ratios],
            "min_ratio": float(min(ratios)), "passes": bool(min(ratios) >= 30.0)}


def verify_point(problem, label, u, h_star, names):
    u = np.asarray(u, dtype=float)
    print(f"\n== {label} ==", flush=True)
    twins, system, state = forward_twins(problem, u)
    t1, parts = tier1(system, state)
    t2 = tier2(system, state, parts)
    rows = {"stress": problem.stress_constraint(), "deflection": problem.deflection_constraint()}
    t3 = {name: tier3(problem, row, u, h_star, names) for name, row in rows.items()}
    tay = {name: [taylor(problem, row, u, seed) for seed in TAYLOR_SEEDS] for name, row in rows.items()}
    point = {
        "label": label, "u": [float(x) for x in u], "x": [float(x) for x in problem.physical(u)],
        "mass_gradient": mass_gradient_check(problem, u),
        "manufacturing": manufacturing_check(problem, u),
        "forward_twins": twins,
        "tier1": t1, "tier2": t2, "tier3": t3, "taylor": tay,
        "tier3_all_pass": bool(all(t["passes"] for t in t3.values())),
        "taylor_all_pass": bool(all(t["passes"] for ts in tay.values() for t in ts)),
    }
    print(f"  mass gradient: shell {point['mass_gradient']['shell']['max_rel_error']:.2e}, "
          f"solid {point['mass_gradient']['solid']['max_rel_error']:.2e}; manufacturing exact "
          f"{point['manufacturing']['exact_to_the_bit']} (min slack {point['manufacturing']['min_slack']:.4f})")
    print(f"  twins: deflection {twins['deflection_max_rel_error']:.2e}, moment {twins['moment_max_rel_error']:.2e}")
    print(f"  Tier 1: dD/dphi {t1['deflection']['d_dphi']:.2e}, dD/dd {t1['deflection']['d_dd_explicit']:.2e}, "
          f"dKS/dx {t1['moment']['d_dphi']:.2e}, dKS/dd {t1['moment']['d_dd_explicit']:.2e}")
    print(f"  Tier 2: tangent {t2['tangent_vs_adjoint']:.2e}, assembly {t2['assembly_identity']:.2e}")
    for name, t in t3.items():
        print(f"  Tier 3 {name}: worst |diff|/eps {t['worst_abs_error_over_eps']:.3f} at "
              f"{t['worst_variable']} (passes {t['passes']}, floor-limited {t['floor_limited']})")
    print(f"  Taylor: {[round(t['min_ratio'], 1) for ts in tay.values() for t in ts]}", flush=True)
    return point


# ---------------------------------------------------------------------------
# the reference-blades table and the forward re-evaluation of x_m
# ---------------------------------------------------------------------------

def reevaluate_from_json(problem, artefact):
    """`x_m` from its JSON through the forward path only, against the recorded numbers."""

    x_m = np.array(artefact["x_m"], dtype=float)
    recorded = artefact["optimum"]
    geometry = problem.parameterisation.to_geometry(x_m, polar_cache=problem.polar_cache)
    air = (float(problem.site.air_density), float(problem.site.kinematic_viscosity))
    v_rated = float(problem.design.rated_wind_speed_ms)
    tsr_rated = tsr_schedule(v_rated, problem.design.design_tsr, problem.design.max_tip_speed_ms)
    aep = -float(objective(x_m, resource=problem.resource, parameterisation=problem.parameterisation,
                           polar_cache=problem.polar_cache))
    moment = float(root_moment(geometry, v_rated, tsr_rated, *air)[0])
    deflection = float(tip_deflection_at(geometry, v_rated, tsr_rated, *air)[0])
    shell = float(problem.material_model().shell_area(x_m))
    comparisons = {
        "aep_mwh_per_yr": (aep, recorded["aep_mwh_per_yr"]),
        "root_bending_moment_nm": (moment, recorded["loads_rated"]["root_bending_moment_nm"]),
        "tip_deflection_per_unit_stiffness": (deflection, recorded["loads_rated"]["tip_deflection_per_unit_stiffness"]),
        "shell_area_m2": (shell, recorded["material"]["shell_area_m2"]),
    }
    out = {name: {"forward": f, "recorded": r, "rel_error": abs(f / r - 1.0)}
           for name, (f, r) in comparisons.items()}
    out["all_within_1e-10"] = bool(all(v["rel_error"] <= 1e-10 for v in out.values()))
    return out


def plot_reference_blades(table, path):
    """
    The three committed blades as two panels side by side, written through
    `figstyle.save`:

    * "Material proxies" -- the shell and the solid material of each blade
      relative to `x0`, in per cent (grouped bars);
    * "Constrained quantities" -- the AEP, KS-moment, stress-proxy and
      deflection-proxy ratios to `x0`, with a line at 1.0.

    The blade labels and tick labels come from `figstyle.BLADES`, so the
    bars read the same as every other figure. `table` is the list of blade
    records (`reference_blades.json`'s `blades`). Returns the path.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from plotting import figstyle

    figstyle.apply()

    keys = [b["label"] for b in table]
    tick_labels = [figstyle.BLADES[key]["label"] for key in keys]
    x = np.arange(len(table))
    blue, light, ink = "#1f3b8b", "#5aa9e6", "#888888"

    fig, axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE)

    width = 0.38
    axes[0].bar(x - width / 2.0, [b["shell_pct_vs_x0"] for b in table], width,
                color=blue, label=r"shell $k_P\int c\,dr$")
    axes[0].bar(x + width / 2.0, [b["solid_pct_vs_x0"] for b in table], width,
                color=light, label=r"solid $k_A\int c^2\,dr$")
    axes[0].axhline(0.0, color=ink, lw=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(tick_labels, rotation=35, ha="right",
                            rotation_mode="anchor", fontsize=7)
    axes[0].set_ylabel(figstyle.LABELS["shell_material"])
    figstyle.title(axes[0], "Material proxies")
    axes[0].legend(loc="upper left", fontsize=7)

    keys = (("aep_over_x0", "energy"), ("KS_over_KS0", "KS moment"),
            ("stress_ratio", "stress proxy"), ("deflection_ratio", "deflection proxy"))
    width = 0.2
    for k, (key, name) in enumerate(keys):
        axes[1].bar(x + (k - 1.5) * width, [b[key] for b in table], width,
                    label=name)
    axes[1].axhline(1.0, color=ink, lw=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(tick_labels, rotation=35, ha="right",
                            rotation_mode="anchor", fontsize=7)
    axes[1].set_ylabel(figstyle.LABELS["ratio_to_reference"])
    axes[1].set_ylim(0.72, 1.16)
    figstyle.title(axes[1], "Constrained quantities")
    axes[1].legend(loc="upper left", ncol=2, fontsize=7)

    fig.tight_layout()
    figstyle.save(fig, path)
    plt.close(fig)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--replot", action="store_true",
                        help="redraw reference_blades.png from the committed "
                             "reference_blades.json; no check runs, no JSON "
                             "is written")
    args = parser.parse_args(argv)

    if args.replot:
        table = C.load_json(TABLE_PATH)["blades"]
        plot_reference_blades(table, FIGURE_PATH)
        print(f"redrew {FIGURE_PATH} from {TABLE_PATH}; no check run",
              flush=True)
        return table

    artefact = C.load_json(RESULT_PATH)
    h_star = float(C.load_json(C.SWEEP_PATH)["h_star_global"])
    problem, u0 = C.prepared_problem()
    names = C.variable_names(problem)
    u_c = problem.scaled(C.load_xc())
    u_m = np.array(artefact["u_m"], dtype=float)
    print(f"J(x0) = {problem.J0:.6f}; KS0 = {problem.KS0!r}; D0 = {problem.D0!r}; "
          f"delta_ref = {problem.delta_ref!r}; material_ref = {problem.material_ref!r}", flush=True)

    started = time.perf_counter()
    points = [verify_point(problem, "x0", u0, h_star, names),
              verify_point(problem, "mass optimum x_m", u_m, h_star, names)]
    table = [C.blade_record(problem, u0, 0.0, "x0", u0=u0),
             C.blade_record(problem, u_c, 0.0, "x_c", u0=u0),
             C.blade_record(problem, u_m, 0.0, "x_m", u0=u0)]
    reevaluation = reevaluate_from_json(problem, artefact)
    wall = time.perf_counter() - started
    plot_reference_blades(table, FIGURE_PATH)
    render_paths = C.render_blades(problem, [("x0", table[0]),
                                             ("x_c", table[1]),
                                             ("x_m", table[2])], _HERE)

    checks = {
        "problem": C.PROBLEM_LABEL,
        "bounds_label": C.BOUNDS_LABEL, "law_label": C.LAW_LABEL,
        "command": "python verification/mass_optimisation/run_mass_checks.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "h_star": h_star,
        "acceptance": {"tier1": "complex step, mixed tolerance",
                       "tier2": "tangent vs adjoint, assembly identity",
                       "tier3": f"|adj - fd| <= {EPS_FACTOR:g} eps_j; floor measured on failure only",
                       "taylor": "remainder ratio >= 30 per decade",
                       "mass_gradient": f"central FD at h = {FD_H_MASS:g}, round-off"},
        "references": {"J0": float(problem.J0), "KS0": float(problem.KS0), "D0": float(problem.D0),
                       "delta_ref": float(problem.delta_ref), "material_ref": float(problem.material_ref)},
        "points": points,
        "all_tier3_pass": bool(all(p["tier3_all_pass"] for p in points)),
        "tier3_failures": [{"point": p["label"], "row": name, "variable": t["worst_variable"],
                            "floor_limited": t["floor_limited"]}
                           for p in points for name, t in p["tier3"].items() if not t["passes"]],
        "all_taylor_pass": bool(all(p["taylor_all_pass"] for p in points)),
        "reevaluation_of_x_m_from_json": reevaluation,
        "wall_time_s": wall,
    }
    with open(CHECKS_PATH, "w", encoding="utf-8") as handle:
        json.dump(checks, handle, indent=1)
    with open(TABLE_PATH, "w", encoding="utf-8") as handle:
        json.dump({
            "problem": C.PROBLEM_LABEL,
            "bounds_label": C.BOUNDS_LABEL, "law_label": C.LAW_LABEL,
            "command": checks["command"], "generated": checks["generated"],
            "src_commit": checks["src_commit"],
            "sources": {"x0": "verification/baseline/x0.json",
                        "x_c": "verification/load_constraint/result_eps0.json (x_c)",
                        "x_m": "verification/mass_optimisation/result_delta0.json (x_m)"},
            "blades": table,
            "headline": {
                "x_c_shell_pct_vs_x0": table[1]["shell_pct_vs_x0"],
                "x_c_aep_pct_vs_x0": table[1]["aep_pct_vs_x0"],
                "x_m_shell_pct_vs_x0": table[2]["shell_pct_vs_x0"],
                "x_m_solid_pct_vs_x0": table[2]["solid_pct_vs_x0"],
                "x_m_aep_pct_vs_x0": table[2]["aep_pct_vs_x0"],
            },
        }, handle, indent=1)

    print("\nreference blades:")
    for b in table:
        print(f"  {b['label']:>4}: AEP {b['aep_pct_vs_x0']:+.4f} %  shell {b['shell_pct_vs_x0']:+.3f} %  "
              f"solid {b['solid_pct_vs_x0']:+.3f} %  KS/KS0 {b['KS_over_KS0']:.4f}  stress {b['stress_ratio']:.4f}  "
              f"defl {b['deflection_ratio']:.4f}  M_rated {b['loads_rated']['root_bending_moment_nm']:.2f} N m  "
              f"thrust {b['loads_rated']['rotor_thrust_n']:.1f} N  Re_min {b['reynolds_min']:.0f}")
    print(f"re-evaluation of x_m from JSON within 1e-10: {reevaluation['all_within_1e-10']} "
          f"({ {k: v['rel_error'] for k, v in reevaluation.items() if isinstance(v, dict)} })")
    print(f"Tier 3 all pass: {checks['all_tier3_pass']}; failures {checks['tier3_failures']}; "
          f"Taylor all pass: {checks['all_taylor_pass']}  ({wall:.0f} s)")
    print(f"wrote {CHECKS_PATH}\n      {TABLE_PATH}\n      {FIGURE_PATH}\n      "
          + "\n      ".join(render_paths))
    return checks


if __name__ == "__main__":
    main()
