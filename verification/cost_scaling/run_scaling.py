"""
Phase 4, Step 3 -- the cost-scaling study: gradient wall time against the
number of design variables.

Plan section 8.6: the parameterisation is run at 5, 10, 20, 40 and 80 control
points per block (`n = 2k = 10 ... 160` design variables, 25 BEM strips
throughout), and at each count the cost of one objective evaluation and of one
gradient by each of the three methods is measured:

  * `J`               -- `ScaledProblem.J(u)`: 17 forward solves, one per bin.
  * central FD        -- `ScaledProblem.jac_fd(u, h*)`: `2n` objective
                         evaluations at the committed global step
                         `h* = 3.162277660168379e-06` in `u`.
  * tangent (direct)  -- one forward solve, one `partials`, then `n`
                         directional products `BEMSystem.tangent_from_parts`
                         along `e_j`. Timed twice: the total, and the `n`
                         products alone (the part that grows with `n`).
  * adjoint           -- `BEMSystem.gradient(d)`: one solve, the partials over
                         17 x 25 stations, `psi`, one assembly.
  * moment adjoint    -- `RootMomentSystem.gradient(d)`: the Phase 4
                         constraint's gradient, the same chain over the
                         9-point load set `L` (9 x 25 stations).

Forward solves per gradient are counted, not timed: FD `2n`, tangent 1,
adjoint 1, moment adjoint 1.

The point at each `k` is the least-squares projection of the analytic Schmitz
blade onto that block count (`build_schmitz_baseline(parameterisation)`, the
same construction as `x0`; `k = 5` *is* `x0`). Beyond 25 control points the
25-station fit is underdetermined and `lstsq` returns the minimum-norm
solution: the station chords and twists are reproduced exactly, some control
points fall outside the configured box, and neither matters for timing. The
AEP at each projected point is recorded so the reader can see it is the same
blade (within 0.02 %).

Wall times are best-of-5 (`time.perf_counter`): a wall-clock measurement on a
shared machine has a hard floor and an unbounded tail, so the minimum is the
robust estimator (the same rule as `tests/test_cost.py`). Each method is
fitted with `t = a n^p` on the five points; the tangent total is also fitted
affinely, `t = c + b n`, because its constant part (one solve + partials)
dominates over this range and a power law would misreport it.

Outputs, next to this script: `scaling.json`, `scaling.png`, and a
hand-written `README.md`.

Run from the repo root (about 12 minutes; the FD gradient at n = 160 is 320
objective evaluations, five times):

    python verification/cost_scaling/run_scaling.py

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import datetime
import json
import os
import platform
import subprocess
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from adjoint.loads import RootMomentSystem  # noqa: E402
from adjoint.system import BEMSystem  # noqa: E402
from design import BladeParameterisation, DesignBounds, build_schmitz_baseline  # noqa: E402
from gradients import ScaledProblem  # noqa: E402
from objective import WeibullResource  # noqa: E402

BASELINE_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "baseline_reference.json")
SWEEP_PATH = os.path.join(REPO_ROOT, "verification", "fd_step_size", "sweep.json")
JSON_PATH = os.path.join(_HERE, "scaling.json")
FIGURE_PATH = os.path.join(_HERE, "scaling.png")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")
LAW_LABEL = ("lambda(V) = min(6.5, Omega_max R / V), Omega_max = 300 rpm "
             "(V_c = 9.67 m/s), fixed rating 3822.189755449124 W")

#: Control points per block. Plan section 8.6.
CONTROL_POINTS = (5, 10, 20, 40, 80)

#: Best-of repeats, the same rule as `tests/test_cost.py`.
REPEATS = 5

#: `M(x0)` at the design condition, the Phase 4 normalising moment. One number
#: in one place: `ScaledProblem.load_system()` reads it from the committed
#: baseline, but that path only works for the 5 + 5 parameterisation (it scales
#: `x0` with the problem's own bounds), so here it is passed explicitly and
#: cross-checked once against the `k = 5` problem.
M_REF_NM = 177.3755406092969

#: The projected blades must be the same blade to this, in AEP.
SAME_BLADE_PCT = 0.02


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def src_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       cwd=REPO_ROOT, text=True).strip()
    except Exception:  # pragma: no cover - reporting only
        return None


def best_of(fn, repeats=REPEATS):
    """Minimum wall time of `fn()` over `repeats` calls, and the last result."""

    best = np.inf
    result = None
    for _ in range(repeats):
        started = time.perf_counter()
        result = fn()
        best = min(best, time.perf_counter() - started)
    return best, result


class SolveCounter:
    """
    Count `BEMSystem.solve` calls on one system without touching `src`.

    The adjoint's structural claim is that its forward-solve count does not
    depend on `n`; the count is read here rather than inferred from timing.
    """

    def __init__(self, system):
        self.system = system
        self.count = 0
        self._solve = system.solve

    def __enter__(self):
        def counted(d):
            self.count += 1
            return self._solve(d)
        self.system.solve = counted
        return self

    def __exit__(self, *exc):
        self.system.solve = self._solve
        return False


def fit_power(n, t):
    """`t = a n^p` by least squares in log space; returns `(a, p)`."""

    p, log_a = np.polyfit(np.log(n), np.log(t), 1)
    return float(np.exp(log_a)), float(p)


def fit_affine(n, t):
    """`t = c + b n` by least squares; returns `(c, b)`."""

    b, c = np.polyfit(n, t, 1)
    return float(c), float(b)


# ---------------------------------------------------------------------------
# one control-point count
# ---------------------------------------------------------------------------

def measure(k, resource, h_star, aep_x0):
    parameterisation = BladeParameterisation(n_chord=k, n_twist=k)
    bounds = DesignBounds.from_config(n_chord=k, n_twist=k)
    n = parameterisation.n_design_variables
    baseline = build_schmitz_baseline(parameterisation=parameterisation)
    d = baseline.design_vector

    problem = ScaledProblem(parameterisation, bounds, resource)
    u = problem.scaled(d)
    J0 = problem.set_reference(u)
    aep = -J0
    aep_delta_pct = 100.0 * (aep - aep_x0) / aep_x0
    print(f"\nk = {k:3d}  n = {n:3d}  degree {parameterisation.degree}  "
          f"AEP {aep:.6f} MWh/yr ({aep_delta_pct:+.4f} % vs x0)  "
          f"fit rms chord {baseline.chord_rms_error_m:.2e} m, twist "
          f"{baseline.twist_rms_error_rad:.2e} rad, control points outside the box: "
          f"{len(baseline.feasibility['violations'])}")
    if abs(aep_delta_pct) > SAME_BLADE_PCT:
        raise RuntimeError(f"k = {k}: the projected blade differs from x0 by "
                           f"{aep_delta_pct:+.4f} % in AEP, over {SAME_BLADE_PCT} %")

    system = BEMSystem(parameterisation, bounds, resource)
    moment_system = RootMomentSystem(parameterisation, bounds, resource, m_ref_nm=M_REF_NM)
    if k == 5:
        m_ref_problem = float(problem.load_system().m_ref)
        if m_ref_problem != M_REF_NM:
            raise RuntimeError(f"M_REF_NM {M_REF_NM!r} != ScaledProblem's {m_ref_problem!r}")

    # Warm every path once: the memoised polar interpolant, the solver's first
    # call, the load system's first solve. Not timed.
    problem.J(u)
    system.gradient(d)
    moment_system.gradient(d)

    # J ----------------------------------------------------------------------
    t_J, _ = best_of(lambda: problem.J(u))

    # central FD ---------------------------------------------------------------
    before = problem.n_fun_evals
    t_fd, grad_fd = best_of(lambda: problem.jac_fd(u, h_star))
    fd_solves = (problem.n_fun_evals - before) // REPEATS

    # tangent (direct) ---------------------------------------------------------
    def tangent_products(parts, state):
        grad = np.empty(n)
        for j in range(n):
            e = np.zeros(n)
            e[j] = 1.0
            grad[j] = system.tangent_from_parts(parts, state.phi, state.d, e, state.limited)
        return grad

    def tangent_total():
        with SolveCounter(system) as counter:
            state = system.solve(d)
            parts = system.partials(state.phi, state.d)
            grad = tangent_products(parts, state)
        return grad, counter.count

    t_tangent, (grad_tangent, tangent_solves) = best_of(tangent_total)
    state = system.solve(d)
    parts = system.partials(state.phi, state.d)
    t_tangent_products, _ = best_of(lambda: tangent_products(parts, state))

    # adjoint ------------------------------------------------------------------
    def adjoint():
        with SolveCounter(system) as counter:
            result = system.gradient(d)
        return result, counter.count

    t_adjoint, (adjoint_result, adjoint_solves) = best_of(adjoint)

    # moment adjoint -----------------------------------------------------------
    def moment_adjoint():
        with SolveCounter(moment_system) as counter:
            result = moment_system.gradient(d)
        return result, counter.count

    t_moment, (moment_result, moment_solves) = best_of(moment_adjoint)

    # agreement, so the timed gradients are the same gradient -------------------
    span = bounds.span()
    grad_adjoint_fun = adjoint_result.dJ_dd * span / abs(J0)
    grad_tangent_fun = grad_tangent * span / abs(J0)
    scale = max(1.0, float(np.max(np.abs(grad_adjoint_fun))))
    tangent_vs_adjoint = float(np.max(np.abs(grad_tangent_fun - grad_adjoint_fun)) / scale)
    fd_vs_adjoint = float(np.max(np.abs(grad_fd - grad_adjoint_fun)) / scale)

    print(f"  J {t_J:.4f} s | FD {t_fd:.3f} s ({fd_solves} solves) | tangent "
          f"{t_tangent:.4f} s total, {t_tangent_products:.5f} s products "
          f"({tangent_solves} solve) | adjoint {t_adjoint:.4f} s ({adjoint_solves} solve) | "
          f"moment adjoint {t_moment:.4f} s ({moment_solves} solve)")
    print(f"  tangent vs adjoint {tangent_vs_adjoint:.2e}, FD vs adjoint {fd_vs_adjoint:.2e} "
          f"(relative to max |dfun/du|)")

    return {
        "k": int(k),
        "n": int(n),
        "spline_degree": int(parameterisation.degree),
        "n_stations": int(parameterisation.n_stations),
        "aep_mwh_per_yr": float(aep),
        "aep_delta_pct_vs_x0": float(aep_delta_pct),
        "fit_rms_error_chord_m": float(baseline.chord_rms_error_m),
        "fit_rms_error_twist_rad": float(baseline.twist_rms_error_rad),
        "control_points_outside_box": int(len(baseline.feasibility["violations"])),
        "wall_s": {
            "J": float(t_J),
            "fd": float(t_fd),
            "tangent_total": float(t_tangent),
            "tangent_products": float(t_tangent_products),
            "adjoint": float(t_adjoint),
            "moment_adjoint": float(t_moment),
        },
        "forward_solves_per_gradient": {
            "fd": int(fd_solves),
            "tangent": int(tangent_solves),
            "adjoint": int(adjoint_solves),
            "moment_adjoint": int(moment_solves),
        },
        "ratios": {
            "fd_over_adjoint": float(t_fd / t_adjoint),
            "fd_over_J": float(t_fd / t_J),
            "adjoint_over_J": float(t_adjoint / t_J),
            "moment_adjoint_over_adjoint": float(t_moment / t_adjoint),
        },
        "agreement_max_abs_diff_relative": {
            "tangent_vs_adjoint": tangent_vs_adjoint,
            "fd_vs_adjoint": fd_vs_adjoint,
        },
        "KS_at_point": float(moment_result.KS),
    }


# ---------------------------------------------------------------------------
# fits and figure
# ---------------------------------------------------------------------------

METHODS = (
    ("fd", "central FD (2n evaluations)"),
    ("tangent_total", "tangent, total (1 solve + partials + n products)"),
    ("tangent_products", "tangent, the n products alone"),
    ("adjoint", "adjoint (1 solve + partials + assembly)"),
    ("moment_adjoint", "moment adjoint (9-point load set)"),
    ("J", "J (one objective evaluation)"),
)


def fits(rows):
    n = np.array([row["n"] for row in rows], dtype=float)
    out = {}
    for key, label in METHODS:
        t = np.array([row["wall_s"][key] for row in rows], dtype=float)
        a, p = fit_power(n, t)
        out[key] = {"label": label, "power_law": {"a": a, "p": p}}
        if key == "tangent_total":
            c, b = fit_affine(n, t)
            out[key]["affine"] = {"c": c, "b": b}
    return out


#: The methods drawn in the figure, in draw order, with the legend name the
#: report uses; the fitted exponent is appended from the JSON's `fits` block.
PLOTTED = (
    ("J", "objective"),
    ("fd", "central FD"),
    ("tangent_total", "tangent"),
    ("adjoint", "adjoint"),
    ("moment_adjoint", "moment adjoint"),
)

#: Marker and colour per method.
STYLES = {
    "J": ("x", "0.4"),
    "fd": ("o", "C3"),
    "tangent_total": ("s", "C1"),
    "adjoint": ("^", "C0"),
    "moment_adjoint": ("v", "C2"),
}


def power_text(p):
    """`p` to two decimals, without a signed zero."""

    return f"{0.0 if abs(p) < 0.005 else p:.2f}"


def plot(rows, fitted, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plotting import figstyle

    figstyle.apply()

    n = np.array([row["n"] for row in rows], dtype=float)
    grid = np.geomspace(n[0], n[-1], 100)

    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    for key, name in PLOTTED:
        t = np.array([row["wall_s"][key] for row in rows])
        marker, colour = STYLES[key]
        entry = fitted.get(key, {}).get("power_law")
        label = name if entry is None else f"{name}, $p = {power_text(entry['p'])}$"
        ax.plot(n, t, marker, color=colour, ms=5, label=label)
        if entry is not None:
            ax.plot(grid, entry["a"] * grid ** entry["p"], "-", color=colour,
                    lw=1, alpha=0.7)

    ax.set_xscale("log")
    ax.set_yscale("log")
    # A log axis already labels its decades in math text (`10^{-3}`); the
    # ScalarFormatter behind figstyle.sci would label mantissas only.
    ax.set_xlabel(figstyle.LABELS["design_variables"])
    ax.set_ylabel(figstyle.LABELS["wall_time"])
    ax.grid(True, which="both")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2)
    fig.tight_layout(rect=(0.0, 0.22, 1.0, 1.0))
    figstyle.save(fig, path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--control-points", type=int, nargs="+", default=list(CONTROL_POINTS))
    parser.add_argument("--replot", action="store_true",
                        help="redraw scaling.png from the committed scaling.json; no timing")
    args = parser.parse_args(argv)
    if args.replot:
        summary = load_json(JSON_PATH)
        plot(summary["rows"], summary["fits"], FIGURE_PATH)
        print(f"redrew {FIGURE_PATH} from {JSON_PATH}")
        return

    resource = WeibullResource.from_config()
    sweep = load_json(SWEEP_PATH)
    h_star = float(sweep["h_star_global"])
    aep_x0 = float(load_json(BASELINE_PATH)["aep_mwh_per_year"])
    print(f"h* (global, in u) = {h_star:.6e}; AEP(x0) = {aep_x0:.9f} MWh/yr")

    started = time.perf_counter()
    rows = [measure(k, resource, h_star, aep_x0) for k in args.control_points]
    wall = time.perf_counter() - started

    fitted = fits(rows)
    print("\nfits t = a n^p:")
    for key, entry in fitted.items():
        line = f"  {entry['label']:52s} p = {entry['power_law']['p']:+.3f}"
        if "affine" in entry:
            line += (f"   (affine: c = {entry['affine']['c']:.4f} s, "
                     f"b = {entry['affine']['b']:.2e} s per variable)")
        print(line)

    plot(rows, fitted, FIGURE_PATH)

    summary = {
        "description": ("Phase 4 Step 3 cost-scaling study: wall time of one objective "
                        "evaluation and of one gradient by central FD, tangent (direct) "
                        "mode, the discrete adjoint and the root-moment adjoint, against "
                        "the number of design variables n = 2k, at the least-squares "
                        "projection of the analytic Schmitz blade onto each block count."),
        "bounds_label": BOUNDS_LABEL,
        "law_label": LAW_LABEL,
        "command": "python verification/cost_scaling/run_scaling.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": src_commit(),
        "machine": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "repeats_best_of": REPEATS,
        "h_star_global_u": h_star,
        "m_ref_nm": M_REF_NM,
        "aep_x0_mwh_per_yr": aep_x0,
        "same_blade_tolerance_pct": SAME_BLADE_PCT,
        "note_minimum_norm": ("For k > 25 the 25-station least-squares fit is underdetermined "
                              "and lstsq returns the minimum-norm control points: the station "
                              "chords and twists are reproduced to round-off, some control "
                              "points fall outside the configured box, and neither matters "
                              "for timing. The AEP column shows the same blade throughout."),
        "note_honest_cost": ("With a diagonal dR/dx and one scalar state per station the "
                             "adjoint's advantage over central FD is the 2n factor and "
                             "accuracy, not a linear-solve saving (derivation section 8). "
                             "Its own cost is the partials over 17 x 25 stations, independent "
                             "of n; the N_c^T assembly is 25 x n and negligible."),
        "wall_time_s": float(wall),
        "rows": rows,
        "fits": fitted,
    }
    with open(JSON_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    print(f"\nwrote {JSON_PATH} and {FIGURE_PATH} in {wall / 60:.1f} min")


if __name__ == "__main__":
    main()
