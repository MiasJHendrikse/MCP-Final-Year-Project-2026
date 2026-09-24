"""
Phase 3, B4 -- Tier 4: attribute what Tier 3 left, to the polar interpolation
or to round-off, and to nothing else.

Tier 3 measured `|adjoint_j - FD_j|` at `h*_j` against the FD noise floor
`eps_j` and found every pair inside `3 eps_j`. This script explains the
residual disagreement -- and, in particular, the three pairs above `1 eps_j`
-- by asking what a central-difference stencil `[u - h, u + h]` actually
crosses:

  (i)   alpha knots: the polar interpolant is C1 with C2 breaks at every
        0.5 deg alpha node. For every (b, i), `|d alpha_{b,i} / d u_j| h`
        against the distance to the nearest node, with
        `d alpha/d u_j = (d phi/d u_j - N_theta[i, j]) span_j` and
        `d phi/d u_j = -(dR/dd_j)/(dR/dphi)` from the adjoint system itself.
  (i')  Reynolds rows: the same C2 break in the Reynolds direction (PCHIP in
        log Re), `|d Re_{b,i} / d u_j| h` against the distance to the nearest
        cached row; `d Re/d u_j = (Re/c_i) N_c[i, j] span_j`.
  (iii) Buhl `a = 0.4`: `a(Y, F)` is C1 there, so `dR/dphi` is continuous
        and the adjoint is exact, but the C2 break degrades FD by O(h) delta
        f''. `|d a_{b,i} / d u_j| h` against `|a - 0.4|`, reported on its own
        line, never folded into the polar figure.
  (ii)  shrink (or, where the stencil is already clean, widen to the largest
        clean step on the A3 grid), recompute FD, and show the disagreement
        against the adjoint sits on the round-off floor of `J`, which is
        measured directly: `J` along a line of 1e-12 steps, linear fit,
        residual std = `delta J`; central FD round-off is `delta J / h`.

All of this at the three Tier 3 points, for all ten variables (the worst
Tier 3 variable at each point is the headline). For the worst variable a
15-step sweep on A3's grid (1e-2 .. 1e-9) is run and each step's
`|FD(h) - adjoint|` is tabulated next to its crossing counts and the
round-off floor, so the whole V-curve is attributed: truncation slope on
the large-`h` side, crossings counted there, round-off on the small-`h`
side. Nothing left unexplained.

Outputs, next to this script: `tier4.json`, `tier4_attribution.png` (reference
| optimum), `tier4_attribution_iterate.png` (the mid-run iterate alone), and
the Tier 4 section of `README.md` (by hand, from the JSON).

Bounds: the configured set (`DesignBounds.from_config()`), grounded
2026-09-19 -- `chord_max_m = 0.30 m`; the 0.45 m placeholder is retired.

Run from the repo root (about 1 min):

    python verification/gradient_verification/run_tier4.py
    python verification/gradient_verification/run_tier4.py --replot   # figures only

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
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from bem.corrections import BUHL_AC  # noqa: E402
from design import BladeParameterisation, DesignBounds  # noqa: E402
from gradients import ScaledProblem  # noqa: E402
from gradients.finite_difference import _central_component  # noqa: E402
from objective import WeibullResource  # noqa: E402

SWEEP_PATH = os.path.join(REPO_ROOT, "verification", "fd_step_size", "sweep.json")
TIER3_PATH = os.path.join(_HERE, "tier3.json")
TIER4_PATH = os.path.join(_HERE, "tier4.json")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")

NOISE_STEPS = np.arange(-4, 5) * 1e-12   # the line J is sampled along for delta J


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def build_problem():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config())


# ---------------------------------------------------------------------------
# crossings inside a stencil, from the linearised state
# ---------------------------------------------------------------------------

def _margin(distance, movement):
    """`min(distance / movement)` over stations that move; `inf` if none does."""

    moving = movement > 0.0
    if not np.any(moving):
        return float("inf")
    return float(np.min(distance[moving] / movement[moving]))


class Stencil:
    """Per-station sensitivities to `u_j` at one point, and what a step `h` crosses."""

    def __init__(self, system, state, parts, alpha_nodes_deg, re_rows, span):
        self.system = system
        self.state = state
        self.parts = parts
        self.span = span
        alpha_deg = np.degrees(parts["alpha"])
        self.alpha_deg = alpha_deg
        # distance to the nearest interpolant node, on the interpolant's own grid
        idx = np.searchsorted(alpha_nodes_deg, alpha_deg)
        lo = alpha_nodes_deg[np.clip(idx - 1, 0, len(alpha_nodes_deg) - 1)]
        hi = alpha_nodes_deg[np.clip(idx, 0, len(alpha_nodes_deg) - 1)]
        self.dist_alpha = np.minimum(np.abs(alpha_deg - lo), np.abs(alpha_deg - hi))
        self.dist_re = np.min(np.abs(state.reynolds[:, :, None] - re_rows[None, None, :]), axis=2)
        self.dist_a = np.abs(state.a - BUHL_AC)
        self.chord = np.array(system.chord_twist(state.d)[0])

    def sensitivities(self, j):
        s = self.span[j]
        dphi = self.system.state_sensitivity(self.parts, j) * s
        dalpha_deg = np.degrees(dphi - (self.system.N_theta[:, j] * s)[None, :])
        dre = (self.state.reynolds / self.chord[None, :]) * (self.system.N_c[:, j] * s)[None, :]
        da = (self.parts["da_dphi"] * dphi
              + self.parts["da_dc"] * (self.system.N_c[:, j] * s)[None, :]
              + self.parts["da_dtheta"] * (self.system.N_theta[:, j] * s)[None, :])
        return dalpha_deg, dre, da

    def crossings(self, j, h):
        dalpha, dre, da = self.sensitivities(j)
        move_alpha, move_re, move_a = np.abs(dalpha) * h, np.abs(dre) * h, np.abs(da) * h
        return {
            "h": float(h),
            "alpha_knot_crossings": int(np.sum(move_alpha > self.dist_alpha)),
            "reynolds_row_crossings": int(np.sum(move_re > self.dist_re)),
            "buhl_crossings": int(np.sum(move_a > self.dist_a)),
            # margin = distance / movement: how many stencil half-widths to the
            # nearest break, over the stations this variable actually moves
            "alpha_margin_min": _margin(self.dist_alpha, move_alpha),
            "reynolds_margin_min": _margin(self.dist_re, move_re),
            "buhl_margin_min": _margin(self.dist_a, move_a),
            "max_alpha_move_deg": float(move_alpha.max()),
            "max_reynolds_move": float(move_re.max()),
            "max_a_move": float(move_a.max()),
        }

    def clean(self, j, h):
        c = self.crossings(j, h)
        return (c["alpha_knot_crossings"] + c["reynolds_row_crossings"] + c["buhl_crossings"]) == 0


# ---------------------------------------------------------------------------
# the round-off floor of J
# ---------------------------------------------------------------------------

def evaluation_noise(problem, u, j):
    """
    `delta J`: `J` sampled along `u + t e_j`, `t` in +-4e-12, linear fit,
    residual standard deviation, in MWh/yr. Steps this small move `J` by
    ~1e-12 (its slope), so the residual is the evaluation's own round-off
    (root-finder termination, 450-term sums), not curvature.
    """

    values = []
    for t in NOISE_STEPS:
        v = np.array(u, dtype=float)
        v[j] += t
        values.append(problem.J(v))
    values = np.array(values)
    coefficients = np.polyfit(NOISE_STEPS, values, 1)
    residual = values - np.polyval(coefficients, NOISE_STEPS)
    return float(np.std(residual, ddof=2)), [float(x) for x in values]


# ---------------------------------------------------------------------------
# one point
# ---------------------------------------------------------------------------

def analyse_point(problem, point, sweep_steps, alpha_nodes, re_rows, names, a3_sweep=None):
    label = point["label"]
    print(f"\n== {label} ==", flush=True)
    u = np.array(point["u"], dtype=float)
    h_star = np.array(point["h_star_per_variable"], dtype=float)
    adjoint = np.array(point["adjoint_mwh_per_u"], dtype=float)
    fd = np.array(point["fd_at_h_star_mwh_per_u"], dtype=float)
    eps = np.array(point["eps_mwh_per_u"], dtype=float)
    diff = np.abs(adjoint - fd)

    system = problem.adjoint_system()
    state = system.solve(problem.physical(u))
    parts = system.partials(state.phi, state.d)
    stencil = Stencil(system, state, parts, alpha_nodes, re_rows, problem.bounds.span())

    # (i), (i'), (iii) at h*_j for every variable
    per_variable = []
    for j, name in enumerate(names):
        c = stencil.crossings(j, h_star[j])
        c.update({"variable": name, "h_star": float(h_star[j]),
                  "abs_error_mwh_per_u": float(diff[j]), "eps_mwh_per_u": float(eps[j]),
                  "abs_error_over_eps": float(diff[j] / eps[j])})
        per_variable.append(c)
        print(f"  {name:<8} h*={h_star[j]:.0e} |diff|/eps={diff[j] / eps[j]:5.2f}  "
              f"knot crossings {c['alpha_knot_crossings']:3d} (margin {c['alpha_margin_min']:8.1f}x)  "
              f"Re rows {c['reynolds_row_crossings']:3d} (margin {c['reynolds_margin_min']:8.1f}x)  "
              f"Buhl {c['buhl_crossings']:3d} (margin {c['buhl_margin_min']:8.1f}x)")

    # the headline variable: worst Tier 3 disagreement at this point
    j = int(np.argmax(diff / eps))
    name = names[j]

    # (ii) the round-off floor, measured
    delta_J, samples = evaluation_noise(problem, u, j)
    floor_at_h_star = delta_J / h_star[j]
    print(f"  worst: {name}; delta J = {delta_J:.2e} MWh/yr -> round-off floor at h* "
          f"{floor_at_h_star:.2e} MWh/yr per u; observed |diff| {diff[j]:.2e} "
          f"({diff[j] / floor_at_h_star:.2f} x floor)")

    # (ii) the largest clean step on the grid, FD there
    clean_steps = [h for h in sweep_steps if stencil.clean(j, h)]
    h_clean = max(clean_steps)
    fd_clean = problem.unscale(_central_component(problem.fun, u, j, h_clean))
    diff_clean = abs(fd_clean - adjoint[j])
    print(f"  largest crossing-free step {h_clean:.0e}: |FD - adjoint| = {diff_clean:.2e} "
          f"(round-off floor there {delta_J / h_clean:.2e})")

    # the whole V for the worst variable: sweep, crossings, floor
    sweep_rows = []
    for h in sweep_steps:
        if a3_sweep is not None:
            g = a3_sweep[h][j]
            g = None if g is None else float(problem.unscale(g))
        else:
            g = float(problem.unscale(_central_component(problem.fun, u, j, h)))
        c = stencil.crossings(j, h)
        c["fd_mwh_per_u"] = g
        c["abs_error_mwh_per_u"] = None if g is None else abs(g - adjoint[j])
        c["roundoff_floor_mwh_per_u"] = delta_J / h
        sweep_rows.append(c)
    # slope of log|err| vs log h over the crossing-populated side (h >= 1e-4)
    big = [(r["h"], r["abs_error_mwh_per_u"]) for r in sweep_rows
           if r["h"] >= 1e-4 and r["abs_error_mwh_per_u"]]
    slope = float(np.polyfit(np.log10([h for h, _ in big]), np.log10([e for _, e in big]), 1)[0])

    # The smooth truncation law C h^2, fitted on steps whose stencil crosses
    # nothing and whose error is well clear of the round-off floor; the
    # excess of every crossing-populated step over that law is the polar /
    # Buhl contribution, with its sign.
    def _is_clean(r):
        return (r["alpha_knot_crossings"] + r["reynolds_row_crossings"] + r["buhl_crossings"]) == 0

    # Window: crossing-free steps whose error is well clear of the round-off
    # floor, i.e. steps that are genuinely truncation-dominated. If no step
    # clears that bar the truncation regime is not observable at this point:
    # no law is reported rather than a fit through pure round-off, and the
    # reason is recorded in place of the coefficient.
    fit_rows = [r for r in sweep_rows if _is_clean(r) and r["abs_error_mwh_per_u"]
                and r["abs_error_mwh_per_u"] > 30.0 * r["roundoff_floor_mwh_per_u"]]
    fit_rule = "crossing-free and above 30x the round-off floor"
    C = (float(np.exp(np.mean([math.log(r["abs_error_mwh_per_u"] / r["h"] ** 2)
                               for r in fit_rows]))) if fit_rows else None)
    if C is None:
        fit_rule = ("no law fitted: no crossing-free step rises clear of its round-off floor, "
                    "so the O(h^2) truncation regime is not observable at this point")
    for r in sweep_rows:
        r["smooth_h2_prediction_mwh_per_u"] = None if C is None else C * r["h"] ** 2
        r["crossing_free"] = _is_clean(r)
        r["excess_over_h2_mwh_per_u"] = (None if C is None or r["abs_error_mwh_per_u"] is None
                                         else r["abs_error_mwh_per_u"] - C * r["h"] ** 2)
    crossing_rows = [r for r in sweep_rows
                     if not r["crossing_free"] and r["abs_error_mwh_per_u"] is not None]
    buhl_rows = [r for r in crossing_rows if r["buhl_crossings"] > 0]
    excesses = [r["excess_over_h2_mwh_per_u"] for r in crossing_rows
                if r["excess_over_h2_mwh_per_u"] is not None]
    smallest_crossing_h = min(r["h"] for r in crossing_rows) if crossing_rows else None
    smallest_buhl_h = min(r["h"] for r in buhl_rows) if buhl_rows else None

    law = ("no law fitted" if C is None
           else f"C = {C:.3e} fitted on {len(fit_rows)} steps")
    print(f"  {name}: truncation-side slope d log|err| / d log h = {slope:.2f} over h >= 1e-4 "
          f"(2 = smooth O(h^2); 1 = crossing-dominated O(h)); smooth law C h^2: {law} "
          f"({fit_rule})")
    for r in sweep_rows:
        e = r["abs_error_mwh_per_u"]
        x = r["excess_over_h2_mwh_per_u"]
        p = r["smooth_h2_prediction_mwh_per_u"]
        print(f"    h={r['h']:.0e}  |FD-adj|={'   n/a  ' if e is None else f'{e:.2e}'}  "
              f"h^2 law={'  n/a ' if p is None else f'{p:.1e}'}  "
              f"excess={'   n/a  ' if x is None else f'{x:+.1e}'}  "
              f"floor={r['roundoff_floor_mwh_per_u']:.1e}  knots {r['alpha_knot_crossings']:3d}  "
              f"Re {r['reynolds_row_crossings']:3d}  Buhl {r['buhl_crossings']:3d}")
    print(f"  smallest step with any crossing: "
          f"{'none' if smallest_crossing_h is None else f'{smallest_crossing_h:.0e}'}; "
          f"with a Buhl crossing: "
          f"{'none' if smallest_buhl_h is None else f'{smallest_buhl_h:.0e}'}; h* = {h_star[j]:.0e}")

    return {
        "label": label,
        "u": point["u"],
        "worst_variable": name,
        "worst_index": j,
        "per_variable_at_h_star": per_variable,
        "any_crossing_at_h_star": bool(any(
            c["alpha_knot_crossings"] + c["reynolds_row_crossings"] + c["buhl_crossings"]
            for c in per_variable)),
        "delta_J_mwh_per_yr": delta_J,
        "delta_J_samples_mwh_per_yr": samples,
        "roundoff_floor_at_h_star_mwh_per_u": floor_at_h_star,
        "abs_error_at_h_star_mwh_per_u": float(diff[j]),
        "abs_error_over_roundoff_floor": float(diff[j] / floor_at_h_star),
        "largest_clean_step": float(h_clean),
        "fd_at_clean_step_mwh_per_u": float(fd_clean),
        "abs_error_at_clean_step_mwh_per_u": float(diff_clean),
        "roundoff_floor_at_clean_step_mwh_per_u": float(delta_J / h_clean),
        "truncation_side_slope": slope,
        "smooth_h2_coefficient": C,
        "smooth_h2_fit_rule": fit_rule,
        "smooth_h2_fit_steps": [r["h"] for r in fit_rows],
        "smallest_step_with_any_crossing": (None if smallest_crossing_h is None
                                            else float(smallest_crossing_h)),
        "smallest_step_with_buhl_crossing": smallest_buhl_h,
        "max_abs_excess_at_crossing_steps_mwh_per_u": (None if not excesses
                                                       else float(max(abs(x) for x in excesses))),
        "excess_at_buhl_steps_mwh_per_u": [[r["h"], r["excess_over_h2_mwh_per_u"]] for r in buhl_rows],
        "sweep_worst_variable": sweep_rows,
        "n_stations_above_buhl": int(np.sum(state.a > BUHL_AC)),
        "min_distance_to_alpha_knot_deg": float(stencil.dist_alpha.min()),
        "min_distance_to_reynolds_row": float(stencil.dist_re.min()),
        "min_distance_to_buhl": float(stencil.dist_a.min()),
    }


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

#: Display names for the three points, in the order they are measured.
POINT_NAMES = ("reference", "mid-run iterate", "optimum")


def variable_symbol(name):
    """`chord_3` -> `$c_3$`, `twist_2` -> `$\\theta_2$`."""

    family, index = name.rsplit("_", 1)
    symbol = "c" if family.startswith("chord") else r"\theta"
    return rf"${symbol}_{{{index}}}$"


def crossing_counts(row):
    """The crossings of one step as a short annotation, e.g. `37$\\alpha$+4Re`."""

    parts = []
    if row["alpha_knot_crossings"]:
        parts.append(f"{row['alpha_knot_crossings']}$\\alpha$")
    if row["reynolds_row_crossings"]:
        parts.append(f"{row['reynolds_row_crossings']}Re")
    if row["buhl_crossings"]:
        parts.append(f"{row['buhl_crossings']}B")
    return "+".join(parts)


def draw_attribution_panel(ax, point, panel_title=None):
    """
    One point's `|FD(h) - adjoint|` V over the 15-step grid, with the fitted
    `C h^2` law, the measured round-off floor and the crossing counts.
    """

    from plotting import figstyle

    rows = point["sweep_worst_variable"]
    h = np.array([r["h"] for r in rows])
    err = np.array([np.nan if r["abs_error_mwh_per_u"] is None
                    else r["abs_error_mwh_per_u"] for r in rows])
    floor = np.array([r["roundoff_floor_mwh_per_u"] for r in rows])
    law = np.array([np.nan if r["smooth_h2_prediction_mwh_per_u"] is None
                    else r["smooth_h2_prediction_mwh_per_u"] for r in rows])

    ax.plot(h, err, "-o", color=figstyle.PALETTE[0], lw=1.4, ms=4,
            label="FD $-$ adjoint")
    if np.isfinite(law).any():
        ax.plot(h, law, "--", color="#888888", lw=1.0, label="$C h^2$ fit")
    ax.plot(h, floor, ":", color="#888888", lw=1.2, label="round-off floor")
    # A step whose stencil crosses a polar node, a Reynolds row or the
    # high-induction blend is drawn hollow. The per-step counts are in
    # tier4.json (`crossing_counts` formats them); as text on the curve they
    # would sit on the data.
    crossed = np.array([bool(r["alpha_knot_crossings"]
                             + r["reynolds_row_crossings"] + r["buhl_crossings"])
                        for r in rows])
    if crossed.any():
        ax.plot(h[crossed], err[crossed], "o", ms=4, color=figstyle.PALETTE[0],
                markerfacecolor="white", markeredgewidth=1.1,
                label="stencil crosses a node or blend")

    h_star = next(r["h_star"] for r in point["per_variable_at_h_star"]
                  if r["variable"] == point["worst_variable"])
    ax.axvline(h_star, color="#444444", lw=0.8, alpha=0.6)
    ax.annotate(r"$h^*$", (h_star, 0.04), xycoords=ax.get_xaxis_transform(),
                xytext=(3, 0), textcoords="offset points", fontsize=8,
                color="#444444")

    ax.set_xscale("log")
    ax.set_yscale("log")
    # A log axis already labels its decades in math text (`10^{-3}`); the
    # ScalarFormatter behind figstyle.sci would label mantissas only.
    ax.invert_xaxis()
    ax.set_xlabel(figstyle.LABELS["fd_step"])
    if panel_title is not None:
        figstyle.title(ax, panel_title)


def plot_attribution(points, out_dir):
    """
    The two Tier 4 figures, both written through `figstyle.save`:
    `tier4_attribution.png` = reference | optimum (two panels);
    `tier4_attribution_iterate.png` = the mid-run iterate alone.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plotting import figstyle

    figstyle.apply()

    ylabel = figstyle.label(r"$\|g_j^{\mathrm{FD}}(h) - g_j^{\mathrm{adj}}\|$",
                            None, r"MWh yr$^{-1}$")
    # The points are stored in the measured order: reference, mid-run iterate,
    # optimum (the order of `tier3.json`'s points).
    reference, iterate, optimum = points[0], points[1], points[2]

    attribution = os.path.join(out_dir, "tier4_attribution.png")
    figure, axes = plt.subplots(1, 2, figsize=figstyle.DOUBLE, sharey=True)
    for ax, point, name in ((axes[0], reference, POINT_NAMES[0]),
                            (axes[1], optimum, POINT_NAMES[2])):
        draw_attribution_panel(
            ax, point, f"{name.capitalize()}: {variable_symbol(point['worst_variable'])}")
    axes[0].set_ylabel(ylabel)
    # One legend under the figure, carrying every series either panel draws.
    entries = {}
    for ax in axes:
        for handle, text in zip(*ax.get_legend_handles_labels()):
            entries.setdefault(text, handle)
    figure.legend(list(entries.values()), list(entries), loc="lower center",
                  ncol=len(entries), frameon=False, bbox_to_anchor=(0.5, 0.0))
    figure.tight_layout(rect=(0, 0.1, 1, 1))
    figstyle.save(figure, attribution)
    plt.close(figure)

    iterate_path = os.path.join(out_dir, "tier4_attribution_iterate.png")
    iterate_figure, iterate_ax = plt.subplots(figsize=figstyle.SINGLE)
    draw_attribution_panel(iterate_ax, iterate)
    iterate_ax.set_ylabel(ylabel)
    figstyle.legend_below(iterate_ax, ncol=2, offset=0.28)
    iterate_figure.tight_layout()
    figstyle.save(iterate_figure, iterate_path)
    plt.close(iterate_figure)

    return [attribution, iterate_path]


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--replot", action="store_true",
                        help="redraw the two Tier 4 figures from the committed "
                             "tier4.json; no solve")
    args = parser.parse_args(argv)
    if args.replot:
        for written in plot_attribution(load_json(TIER4_PATH)["points"], _HERE):
            print(f"redrew {written} from {TIER4_PATH}")
        return

    tier3 = load_json(TIER3_PATH)
    sweep = load_json(SWEEP_PATH)
    names = tier3["variables"]
    sweep_steps = [float(h) for h in sweep["steps"]]
    a3 = {r["h"]: [None if g == "out of cache" else g for g in r["gradient"]] for r in sweep["sweep"]}

    problem = build_problem()
    system = problem.adjoint_system()
    alpha_nodes = np.asarray(system.interpolant.alpha_values, dtype=float)
    re_rows = np.asarray(system.interpolant.re_values, dtype=float)
    u0 = np.array(tier3["points"][0]["u"])
    J0 = problem.set_reference(u0)
    if abs(J0 - tier3["J0_mwh_per_yr"]) > 1e-12 * abs(J0):
        raise RuntimeError("J0 differs from tier3.json's: the objective changed")

    started = time.perf_counter()
    points = []
    for k, point in enumerate(tier3["points"]):
        points.append(analyse_point(problem, point, sweep_steps, alpha_nodes, re_rows, names,
                                    a3_sweep=a3 if k == 0 else None))
    wall = time.perf_counter() - started

    for written in plot_attribution(points, _HERE):
        print(f"wrote {written}")

    summary = {
        "description": "Tier 4: what the FD stencil crosses (alpha knots, Reynolds rows, "
                       "Buhl a = 0.4) at h*_j and across the step grid, the measured "
                       "round-off floor of J, and the attribution of |FD - adjoint|.",
        "bounds_label": BOUNDS_LABEL,
        "command": "python verification/gradient_verification/run_tier4.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "variables": names,
        "alpha_node_spacing_deg": float(np.median(np.diff(alpha_nodes))),
        "reynolds_rows": [float(r) for r in re_rows],
        "buhl_ac": BUHL_AC,
        "steps": sweep_steps,
        "points": points,
        "any_crossing_at_h_star_anywhere": bool(any(p["any_crossing_at_h_star"] for p in points)),
        "wall_time_s": wall,
        "objective_evaluations_total": int(problem.n_fun_evals),
    }
    with open(TIER4_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    print(f"\ncrossings at h* anywhere: {summary['any_crossing_at_h_star_anywhere']}   "
          f"({wall:.1f} s, {problem.n_fun_evals} objective evaluations)")
    print(f"wrote {TIER4_PATH}")
    return summary


if __name__ == "__main__":
    main()
