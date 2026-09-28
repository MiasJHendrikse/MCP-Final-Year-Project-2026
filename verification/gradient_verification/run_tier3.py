"""
Tier 3: the adjoint gradient against the FD noise floor.

At three points -- `x0`, the mid-run iterate `k = 17` of the FD-driven SLSQP
run (`verification/fd_optimisation/iterates.json`, `nit/2`), and that run's
optimum `u*` -- for all 10 scaled design variables:

    adjoint_j      `ScaledProblem.jac_adjoint`, undone to MWh/yr per unit u
    FD_j           central difference at `h*_j` from the step-size study
    eps_j          the FD noise floor: `max(|g(h*) - g(h*/sqrt10)|, |g(h*) - g(h* sqrt10)|)`.
                   At `x0` this is the committed A3 value; at the other two
                   points it is re-measured by the same three-step local
                   sweep, because the floor is a property of the point.

Acceptance, derived not invented:

    |adjoint_j - FD_j| <= 3 eps_j          and        |adjoint_j - FD_j| / |FD_j| <= 1e-3

The first is the reference's own noise; the second is a gross-error
tripwire. A violation is a bug until Tier 4 (`run_tier4.py`) attributes it
to the polar interpolation specifically.

Whole-chain Taylor remainder at every point: random unit `v`,
`r(eps) = |fun(u + eps v) - fun(u) - eps g.v|` for `eps = 1e-2, 1e-3, 1e-4`
with `g` the adjoint; `r` must fall ~100x per decade, asserted >= 30 (a
wrong gradient gives ~10; a C2 break inside the step can pull it below 100
without the gradient being wrong).

Also: wall time of `J`, the FD gradient and the adjoint gradient at `x0`,
and A3's V-curves re-plotted with the adjoint as the reference
(`--reference` machinery of `run_sweep.py`, output kept here so A3's
committed figure is untouched).

Outputs, next to this script: `tier3.json`, `adjoint_x0.json` (the
reference file for `run_sweep.py --reference`), `v_curve_vs_adjoint.png`,
`tier3_agreement.png`, and `README.md` (by hand, from the JSON).

Bounds: the configured set (`DesignBounds.from_config()`), grounded
2026-09-19 -- `chord_max_m = 0.30 m`. The scaling and every gradient here
are stated under those bounds; the 0.45 m placeholder is retired.

Run from the repo root (about 1.5 min):

    python verification/gradient_verification/run_tier3.py
    python verification/gradient_verification/run_tier3.py --replot   # figure only

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import argparse
import datetime
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from design import BladeParameterisation, DesignBounds  # noqa: E402
from gradients import ScaledProblem, central_difference  # noqa: E402
from objective import WeibullResource  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
SWEEP_PATH = os.path.join(REPO_ROOT, "verification", "fd_step_size", "sweep.json")
SWEEP_SCRIPT = os.path.join(REPO_ROOT, "verification", "fd_step_size", "run_sweep.py")
ITERATES_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation", "iterates.json")
RESULT_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation", "result.json")

TIER3_PATH = os.path.join(_HERE, "tier3.json")
ADJOINT_X0_PATH = os.path.join(_HERE, "adjoint_x0.json")
VCURVE_PATH = os.path.join(_HERE, "v_curve_vs_adjoint.png")
AGREEMENT_PATH = os.path.join(_HERE, "tier3_agreement.png")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")

SQRT10 = math.sqrt(10.0)
EPS_FACTOR = 3.0           # |adj - FD| <= EPS_FACTOR * eps_j
REL_TRIPWIRE = 1e-3        # gross-error tripwire on |adj - FD| / |FD|
TAYLOR_STEPS = (1e-2, 1e-3, 1e-4)
TAYLOR_MIN_RATIO = 30.0
TAYLOR_DRAWS = 3
#: The mid-run iterate Tier 3 checks is `nit // 2` of the committed A4 run,
#: read from `fd_optimisation/result.json` (35 iterations -> 17 today; the
#: 27-iteration run before the 300 rpm law gave 13). `MID_ITERATE_TODAY` pins
#: what the committed `tier3.json` was produced at.
MID_ITERATE_TODAY = 17


def mid_iterate_from_result():
    return int(load_json(RESULT_PATH)["nit"]) // 2


def load_x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def build_problem():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config())


def variable_names(parameterisation):
    return ([f"chord_{j}" for j in range(parameterisation.n_chord)]
            + [f"twist_{j}" for j in range(parameterisation.n_twist)])


# ---------------------------------------------------------------------------
# one point
# ---------------------------------------------------------------------------

def fd_at(problem, u, steps):
    """Central FD of `fun` at per-variable steps, in MWh/yr per unit u."""

    grad, _n = central_difference(problem.fun, u, steps)
    return problem.unscale(grad)


def local_noise_floor(problem, u, h_star):
    """`eps_j` re-measured at `u`: the three-step local sweep at `h*_j`."""

    g_mid = fd_at(problem, u, h_star)
    g_lo = fd_at(problem, u, h_star / SQRT10)
    g_hi = fd_at(problem, u, h_star * SQRT10)
    eps = np.maximum(np.abs(g_mid - g_lo), np.abs(g_mid - g_hi))
    return g_mid, eps, {"h_star_over_sqrt10": g_lo, "h_star": g_mid, "h_star_times_sqrt10": g_hi}


def taylor_remainder(problem, u, g_scaled, seed):
    """`|fun(u + eps v) - fun(u) - eps g.v|` for the three steps, one unit `v`."""

    rng = np.random.default_rng(seed)
    envelope = problem.envelope_constraint()["fun"]
    for _attempt in range(50):
        v = rng.normal(size=problem.n)
        v /= np.linalg.norm(v)
        if all(np.all(envelope(u + eps * v) > 0.0) for eps in TAYLOR_STEPS):
            break
    else:
        raise RuntimeError("no envelope-feasible direction found for the Taylor test")

    f0 = problem.fun(u)
    slope = float(g_scaled @ v)
    remainders = []
    for eps in TAYLOR_STEPS:
        remainders.append(abs(problem.fun(u + eps * v) - f0 - eps * slope))
    ratios = [remainders[k] / remainders[k + 1] if remainders[k + 1] > 0.0 else float("inf")
              for k in range(len(TAYLOR_STEPS) - 1)]
    return {
        "seed": seed,
        "v": [float(x) for x in v],
        "g_dot_v_scaled": slope,
        "steps": list(TAYLOR_STEPS),
        "remainder_fun_units": [float(r) for r in remainders],
        "ratios_per_decade": [float(r) for r in ratios],
        "min_ratio": float(min(ratios)),
        "passes": bool(min(ratios) >= TAYLOR_MIN_RATIO),
    }


def verify_point(problem, label, u, h_star, eps_reference, names, time_it=False):
    print(f"\n== {label} ==", flush=True)
    u = np.asarray(u, dtype=float)

    timings = {}
    if time_it:
        started = time.perf_counter()
        problem.fun(u)
        timings["J_s"] = time.perf_counter() - started
        started = time.perf_counter()
        problem.jac_fd(u, h_star)
        timings["fd_gradient_s"] = time.perf_counter() - started

    started = time.perf_counter()
    adjoint_scaled = problem.jac_adjoint(u)
    timings["adjoint_gradient_s"] = time.perf_counter() - started
    adjoint = problem.unscale(adjoint_scaled)

    fd, eps_local, local = local_noise_floor(problem, u, h_star)
    eps = eps_local if eps_reference is None else np.asarray(eps_reference, dtype=float)

    abs_err = np.abs(adjoint - fd)
    rel_err = abs_err / np.abs(fd)
    within_eps = abs_err <= EPS_FACTOR * eps
    within_rel = rel_err <= REL_TRIPWIRE

    print(f"{'variable':<9} {'adjoint':>16} {'FD(h*)':>16} {'|diff|':>10} {'rel':>9} "
          f"{'h*':>8} {'eps':>10} {'|diff|/eps':>10}")
    for j, name in enumerate(names):
        print(f"{name:<9} {adjoint[j]:>16.10f} {fd[j]:>16.10f} {abs_err[j]:>10.2e} "
              f"{rel_err[j]:>9.2e} {h_star[j]:>8.0e} {eps[j]:>10.2e} {abs_err[j] / eps[j]:>10.2f}"
              + ("" if within_eps[j] and within_rel[j] else "   <-- FAIL"))

    taylor = [taylor_remainder(problem, u, adjoint_scaled, seed) for seed in range(TAYLOR_DRAWS)]
    for t in taylor:
        print(f"Taylor seed {t['seed']}: remainders {['%.3e' % r for r in t['remainder_fun_units']]} "
              f"ratios {['%.1f' % r for r in t['ratios_per_decade']]} "
              f"{'ok' if t['passes'] else '<-- FAIL'}")

    return {
        "label": label,
        "u": [float(x) for x in u],
        "x": [float(x) for x in problem.physical(u)],
        "J_mwh_per_yr": float(problem.J(u)),
        "adjoint_mwh_per_u": [float(x) for x in adjoint],
        "adjoint_scaled": [float(x) for x in adjoint_scaled],
        "fd_at_h_star_mwh_per_u": [float(x) for x in fd],
        "fd_local_sweep_mwh_per_u": {k: [float(x) for x in v] for k, v in local.items()},
        "h_star_per_variable": [float(x) for x in h_star],
        "eps_source": "A3 sweep.json (epsilon_mwh_per_u)" if eps_reference is not None
                      else "re-measured here (three-step local sweep at h*_j)",
        "eps_mwh_per_u": [float(x) for x in eps],
        "eps_remeasured_mwh_per_u": [float(x) for x in eps_local],
        "abs_error_mwh_per_u": [float(x) for x in abs_err],
        "rel_error": [float(x) for x in rel_err],
        "abs_error_over_eps": [float(x) for x in abs_err / eps],
        "within_3eps": [bool(x) for x in within_eps],
        "within_rel_tripwire": [bool(x) for x in within_rel],
        "passes": bool(np.all(within_eps) and np.all(within_rel)),
        "worst_variable": names[int(np.argmax(abs_err / eps))],
        "worst_abs_error_over_eps": float(np.max(abs_err / eps)),
        "worst_rel_error": float(np.max(rel_err)),
        "taylor": taylor,
        "taylor_passes": bool(all(t["passes"] for t in taylor)),
        "timings_s": timings,
    }


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def replot_v_curves(sweep, adjoint_x0):
    spec = importlib.util.spec_from_file_location("run_sweep", SWEEP_SCRIPT)
    run_sweep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_sweep)

    steps = np.array([r["h"] for r in sweep["sweep"]])
    grads = run_sweep.gradient_matrix(sweep["sweep"])
    reference = run_sweep.load_reference(ADJOINT_X0_PATH, sweep["J0_mwh_per_yr"])[0]
    run_sweep.plot(steps, grads, reference, sweep["variables"], sweep["n_chord"],
                   sweep["h_star_global"], VCURVE_PATH,
                   run_sweep.spans_from_bounds(sweep["bounds"], sweep["variables"]),
                   sweep["J0_mwh_per_yr"])


#: Display names for the three points, in the order they are measured.
POINT_NAMES = ("reference", "mid-run iterate", "optimum")


def variable_symbol(name):
    """`chord_3` -> `$c_3$`, `twist_2` -> `$\\theta_2$`."""

    family, index = name.rsplit("_", 1)
    symbol = "c" if family.startswith("chord") else r"\theta"
    return rf"${symbol}_{{{index}}}$"


def plot_agreement(data, out_dir):
    """
    `tier3_agreement.png`: `|adjoint - FD(h*_j)| / eps_j` and the relative
    error, per variable at the three points, against their acceptance lines.
    Writes through `figstyle.save` and returns the path written.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.lines
    import matplotlib.pyplot as plt
    from plotting import figstyle

    figstyle.apply()

    points, names = data["points"], data["variables"]
    path = os.path.join(out_dir, "tier3_agreement.png")

    colours = figstyle.PALETTE[:3]
    fig, axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE)
    x = np.arange(len(names))
    width = 0.26

    for k, point in enumerate(points):
        label = POINT_NAMES[k] if k < len(POINT_NAMES) else point["label"]
        axes[0].bar(x + (k - 1) * width, point["abs_error_over_eps"], width,
                    color=colours[k], label=label)
        axes[1].bar(x + (k - 1) * width, point["rel_error"], width,
                    color=colours[k], label=label)

    axes[0].axhline(EPS_FACTOR, color="#444444", lw=0.9, ls="--",
                    label=rf"${EPS_FACTOR:g}\,\varepsilon_j$")
    axes[0].set_yscale("log")
    axes[0].set_ylabel(figstyle.label(
        "Disagreement",
        r"\|g_j^{\mathrm{FD}} - g_j^{\mathrm{adj}}\|/\varepsilon_j", "--"))
    figstyle.title(axes[0], "Disagreement over acceptance scale")


    axes[1].axhline(REL_TRIPWIRE, color="#444444", lw=0.9, ls="--",
                    label=r"$10^{-3}$ tripwire")
    axes[1].set_yscale("log")
    axes[1].set_ylabel(figstyle.label(
        "Relative error",
        r"\|g_j^{\mathrm{FD}} - g_j^{\mathrm{adj}}\|/\|g_j^{\mathrm{FD}}\|", "--"))
    figstyle.title(axes[1], "Relative error")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels([variable_symbol(name) for name in names])
        # A log axis labels its decades in math text already; the
        # ScalarFormatter behind figstyle.sci would label mantissas only.

    # The bars and the acceptance lines leave no free corner, and the three
    # points are the same in both panels: one legend, under the figure.
    # Each acceptance line is named where it is drawn, at the right-hand end,
    # clear of the bars; the legend carries the three points only.
    for ax, value, text, above in ((axes[0], EPS_FACTOR,
                                    rf"${EPS_FACTOR:g}\,\varepsilon_j$", True),
                                   (axes[1], REL_TRIPWIRE,
                                    r"$10^{-3}$ tripwire", False)):
        ax.annotate(text, (1.0, value), xycoords=ax.get_yaxis_transform(),
                    xytext=(-3, 3 if above else -3), textcoords="offset points",
                    ha="right", va="bottom" if above else "top", fontsize=8,
                    color="#444444")
    bars = [h for h in axes[0].get_legend_handles_labels()[0]
            if not isinstance(h, matplotlib.lines.Line2D)]
    fig.legend(bars, [b.get_label() for b in bars], loc="lower center",
               ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    figstyle.save(fig, path)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mid-iterate", type=int, default=None,
                        help="default: nit // 2 of fd_optimisation/result.json")
    parser.add_argument("--replot", action="store_true",
                        help="redraw tier3_agreement.png from the committed tier3.json; "
                             "no re-run")
    args = parser.parse_args(argv)
    if args.replot:
        print(f"redrew {plot_agreement(load_json(TIER3_PATH), _HERE)} "
              f"from {TIER3_PATH}")
        return
    if args.mid_iterate is None:
        args.mid_iterate = mid_iterate_from_result()
        assert args.mid_iterate == MID_ITERATE_TODAY, (
            f"nit // 2 is {args.mid_iterate}, the committed tier3.json was produced "
            f"at {MID_ITERATE_TODAY}: fd_optimisation/result.json moved")

    problem = build_problem()
    names = variable_names(problem.parameterisation)
    x0 = load_x0()
    u0 = problem.scaled(x0)
    J0 = problem.set_reference(u0)
    print(f"J(x0) = {J0:.6f} MWh/yr")

    sweep = load_json(SWEEP_PATH)
    h_star = np.array(sweep["h_star_per_variable"], dtype=float)
    eps_x0 = np.array(sweep["epsilon_mwh_per_u"], dtype=float)
    if abs(sweep["J0_mwh_per_yr"] - J0) > 1e-12 * abs(J0):
        raise RuntimeError("J0 differs from the step-size study's: the objective changed")

    iterates = load_json(ITERATES_PATH)["iterates"]
    mid = iterates[args.mid_iterate]
    assert mid["k"] == args.mid_iterate
    result = load_json(RESULT_PATH)
    u_star = np.array(result["u_star"], dtype=float)

    started = time.perf_counter()
    points = [
        verify_point(problem, "x0", u0, h_star, eps_x0, names, time_it=True),
        verify_point(problem, f"FD-SLSQP iterate k={args.mid_iterate}", mid["u"], h_star, None, names),
        verify_point(problem, "FD-SLSQP optimum u*", u_star, h_star, None, names),
    ]
    wall = time.perf_counter() - started

    x0_point = points[0]
    with open(ADJOINT_X0_PATH, "w", encoding="utf-8") as handle:
        json.dump({
            "label": "adjoint",
            "description": "Discrete-adjoint gradient of fun(u) = J(u)/|J(u0)| at x0, "
                           "for run_sweep.py --reference.",
            "bounds_label": BOUNDS_LABEL,
            "gradient_mwh_per_u": x0_point["adjoint_mwh_per_u"],
            "gradient_scaled": x0_point["adjoint_scaled"],
            "J0_mwh_per_yr": J0,
        }, handle, indent=1)
    replot_v_curves(sweep, x0_point)

    summary = {
        "description": "Tier 3: discrete-adjoint gradient vs central FD at h*_j at three "
                       "points, with the FD noise floor as the acceptance scale, plus the "
                       "whole-chain Taylor-remainder test.",
        "bounds_label": BOUNDS_LABEL,
        "bounds": problem.bounds.as_record(),
        "command": "python verification/gradient_verification/run_tier3.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "variables": names,
        "J0_mwh_per_yr": J0,
        "acceptance": {
            "abs": f"|adjoint_j - FD_j| <= {EPS_FACTOR:g} eps_j",
            "rel": f"|adjoint_j - FD_j| / |FD_j| <= {REL_TRIPWIRE:g}",
            "taylor": f"remainder ratio per decade >= {TAYLOR_MIN_RATIO:g} for eps = {TAYLOR_STEPS}",
        },
        "mid_iterate": args.mid_iterate,
        "sources": {"h_star_and_eps": os.path.relpath(SWEEP_PATH, REPO_ROOT),
                    "iterates": os.path.relpath(ITERATES_PATH, REPO_ROOT),
                    "optimum": os.path.relpath(RESULT_PATH, REPO_ROOT)},
        "points": points,
        "all_pass": bool(all(p["passes"] and p["taylor_passes"] for p in points)),
        "wall_time_s": wall,
        "objective_evaluations_total": int(problem.n_fun_evals),
        "adjoint_evaluations_total": int(problem.n_adjoint_evals),
    }
    with open(TIER3_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    plot_agreement(summary, _HERE)

    print(f"\nall pass: {summary['all_pass']}   ({wall:.1f} s, "
          f"{problem.n_fun_evals} objective evaluations, {problem.n_adjoint_evals} adjoint)")
    t = x0_point["timings_s"]
    print(f"timings at x0: J {t['J_s']:.3f} s, FD gradient {t['fd_gradient_s']:.2f} s, "
          f"adjoint gradient {t['adjoint_gradient_s']:.3f} s")
    print(f"wrote {TIER3_PATH}\n      {ADJOINT_X0_PATH}\n      {VCURVE_PATH}\n      {AGREEMENT_PATH}")
    return summary


if __name__ == "__main__":
    main()
