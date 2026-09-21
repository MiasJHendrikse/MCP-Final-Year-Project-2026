"""
Phase 2, Stage A4: FD-driven SLSQP, end to end.

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10
    subject to the polar-cache Reynolds envelope (linear, 50 rows, margin 5 %)
               and the configured local-solidity cap (25 rows, inactive while
               chord_max_m = 0.30 m binds first) -- `problem.constraints()`
    gradient   central finite differences at the global h* from
               verification/fd_step_size/sweep.json

`scipy.optimize.minimize(method="SLSQP", ftol=1e-8, maxiter=200)`. Every
accepted iterate is saved to `iterates.json` -- `(k, u_k, J_k, g_k,
wall_time)` -- because the adjoint's Tier 3 check needs a mid-run point, not
just the endpoints. `result.json` carries the optimum, the counters, the
active set and the post-checks; `optimised_blade.png` overlays chord and
twist of `x0` and the optimum.

Sanity gates (the implementation plan §5 A4): the expected AEP improvement is
2-6 %. Above 15 %, or negative, the run is flagged `suspected_defect` in
`result.json` and must be written up as such, not reported as a result.

`PolarDomainError` during the run: the offending `u` is already logged by
`ScaledProblem.fun`; it is written to `result.json`, the margin is raised
from 0.05 to 0.10 once, recorded, and the run restarted. A second failure
propagates.

Bounds: the configured set (`DesignBounds.from_config()`), grounded
2026-09-19 -- `chord_max_m = 0.30 m` binds, `chord_min_m`, `twist_min_rad`
and `twist_max_rad` as in `config/`. `BOUNDS_LABEL` is stamped into every
artefact this script writes; the 0.45 m placeholder is retired.

Run from the repo root:

    python verification/fd_optimisation/run_fd_slsqp.py [--step H] [--maxiter N]

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
from scipy.optimize import Bounds, minimize

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from design import BladeParameterisation, DesignBounds, clamped_knots  # noqa: E402
from gradients import ScaledProblem  # noqa: E402
from gradients.problem import DEFAULT_ENVELOPE_MARGIN  # noqa: E402
from objective import WeibullResource  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
SWEEP_PATH = os.path.join(REPO_ROOT, "verification", "fd_step_size", "sweep.json")
RESULT_PATH = os.path.join(_HERE, "result.json")
ITERATES_PATH = os.path.join(_HERE, "iterates.json")
FIGURE_PATH = os.path.join(_HERE, "optimised_blade.png")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")

ACTIVE_TOL = 1e-6          # |u| or |1 - u| below this: bound active
ENVELOPE_ACTIVE_TOL = 1e-6  # constraint row (metres of chord) below this: active
EXPECTED_GAIN_PCT = (2.0, 6.0)
DEFECT_GAIN_PCT = 15.0


def load_x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"]
                    + artefact["twist_control_points_rad"])


def load_h_star():
    with open(SWEEP_PATH, encoding="utf-8") as handle:
        return float(json.load(handle)["h_star_global"])


def build_problem(margin):
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config(),
                         margin=margin)


class Recorder:
    """
    Caches the last `fun`/`jac` evaluation so the callback can record the
    iterate's objective and gradient without spending extra evaluations
    (SLSQP evaluates both at the accepted point before calling back).
    """

    def __init__(self, problem, h):
        self.problem = problem
        self.h = h
        self.started = time.perf_counter()
        self.iterates = []
        self._last_fun = None
        self._last_jac = None
        self.extra_evals = 0

    def fun(self, u):
        value = self.problem.fun(u)
        self._last_fun = (np.array(u, dtype=float), value)
        return value

    def jac(self, u):
        grad = self.problem.jac_fd(u, self.h)
        self._last_jac = (np.array(u, dtype=float), grad)
        return grad

    def _cached(self, cache, u, compute):
        if cache is not None and np.array_equal(cache[0], u):
            return cache[1]
        self.extra_evals += 1
        return compute(u)

    def callback(self, u):
        u = np.array(u, dtype=float)
        value = self._cached(self._last_fun, u, self.problem.fun)
        grad = self._cached(self._last_jac, u, self.jac)
        self.iterates.append({
            "k": len(self.iterates),
            "u": [float(x) for x in u],
            "J_mwh_per_yr": float(self.problem.unscale(value)),
            "fun": float(value),
            "g_scaled": [float(x) for x in grad],
            "g_mwh_per_u": [float(x) for x in self.problem.unscale(grad)],
            "wall_time_s": time.perf_counter() - self.started,
        })
        record = self.iterates[-1]
        print(f"  k={record['k']:3d}  J={record['J_mwh_per_yr']:.6f} MWh/yr  "
              f"|g|={np.linalg.norm(grad):.3e}  t={record['wall_time_s']:.1f}s",
              flush=True)

    def record_start(self, u0, value, grad):
        self.iterates.append({
            "k": 0,
            "u": [float(x) for x in u0],
            "J_mwh_per_yr": float(self.problem.unscale(value)),
            "fun": float(value),
            "g_scaled": [float(x) for x in grad],
            "g_mwh_per_u": [float(x) for x in self.problem.unscale(grad)],
            "wall_time_s": time.perf_counter() - self.started,
        })


def run(problem, u0, h, maxiter, ftol):
    constraints = problem.constraints()
    recorder = Recorder(problem, h)

    value0 = recorder.fun(u0)
    grad0 = recorder.jac(u0)
    recorder.record_start(u0, value0, grad0)
    print(f"J(x0) = {problem.J0:.6f} MWh/yr; h* = {h:.0e}; margin = {problem.margin}")

    started = time.perf_counter()
    result = minimize(
        recorder.fun, u0, jac=recorder.jac, method="SLSQP",
        bounds=Bounds(np.zeros(problem.n), np.ones(problem.n)),
        constraints=constraints,
        options=dict(ftol=ftol, maxiter=maxiter, disp=True),
        callback=recorder.callback,
    )
    wall = time.perf_counter() - started
    return result, recorder, wall


def plot(problem, x0, x_star, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = problem.parameterisation
    r = p.radii
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    for ax, (label, unit, get, block) in zip(axes, [
        ("chord", "m", p.chord, slice(0, p.n_chord)),
        ("twist", "deg", lambda d: np.degrees(p.twist(d)), slice(p.n_chord, None)),
    ]):
        ax.plot(r, get(x0), "-", color="#888888", lw=1.6, label=r"$\mathbf{x}_0$ Schmitz reference")
        ax.plot(r, get(x_star), "-", color="#1f5fbf", lw=1.6, label="FD-SLSQP optimum")
        scale = np.degrees(1.0) if label == "twist" else 1.0
        # control points drawn at their Greville abscissae (mean of `degree`
        # consecutive interior knots), mapped onto the station span
        n_pts = p.n_chord if label == "chord" else p.n_twist
        knots = clamped_knots(n_pts, p.degree)
        greville = np.array([knots[j + 1:j + p.degree + 1].mean() for j in range(n_pts)])
        s = r[0] + greville * (r[-1] - r[0])
        ax.plot(s, x0[block] * scale, "o", color="#888888", ms=5, mfc="white")
        ax.plot(s, x_star[block] * scale, "o", color="#1f5fbf", ms=5, mfc="white")
        ax.set_xlabel("radius r [m]")
        ax.set_ylabel(f"{label} [{unit}]")
        ax.set_title(f"{label} distribution", fontsize=11)
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--step", type=float, default=None,
                        help="FD step (default: h*_global from the step-size study)")
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--ftol", type=float, default=1e-8)
    args = parser.parse_args(argv)

    h = args.step if args.step is not None else load_h_star()
    x0 = load_x0()

    margin = DEFAULT_ENVELOPE_MARGIN
    margin_raised = False
    domain_errors = []
    while True:
        problem = build_problem(margin)
        u0 = problem.scaled(x0)
        try:
            result, recorder, wall = run(problem, u0, h, args.maxiter, args.ftol)
            break
        except PolarDomainError as error:
            domain_errors.append({"margin": margin, "u": problem.domain_errors[-1],
                                  "error": str(error)})
            if margin_raised:
                _write_partial(domain_errors)
                raise
            print(f"PolarDomainError at margin {margin}: raising margin to 0.10 once "
                  f"and restarting", flush=True)
            margin, margin_raised = 0.10, True

    u_star = np.array(result.x, dtype=float)
    x_star = problem.physical(u_star)
    aep0 = -problem.J0
    aep_star = -problem.J(u_star)
    gain_pct = 100.0 * (aep_star - aep0) / aep0

    envelope = problem.envelope_constraint()
    g_env = envelope["fun"](u_star)
    labels = problem.envelope_row_labels()
    active_env = [labels[i] for i in np.flatnonzero(g_env < ENVELOPE_ACTIVE_TOL)]
    active_lower = [i for i in range(problem.n) if abs(u_star[i]) < ACTIVE_TOL]
    active_upper = [i for i in range(problem.n) if abs(1.0 - u_star[i]) < ACTIVE_TOL]
    names = ([f"chord_{j}" for j in range(problem.parameterisation.n_chord)]
             + [f"twist_{j}" for j in range(problem.parameterisation.n_twist)])

    post = problem.alpha_check(u_star)
    post0 = problem.alpha_check(u0)

    suspected = None
    if gain_pct > DEFECT_GAIN_PCT:
        suspected = (f"AEP gain {gain_pct:.2f} % exceeds {DEFECT_GAIN_PCT} %: suspected "
                     "defect (unfair baseline, units, or gradient) -- not a result")
    elif gain_pct < 0.0:
        suspected = f"AEP gain {gain_pct:.2f} % is negative: suspected defect -- not a result"

    summary = {
        "description": "FD-driven SLSQP on fun(u) = J(u)/|J(u0)| with the polar-cache "
                       "envelope constraint, from x0.",
        "bounds_label": BOUNDS_LABEL,
        "bounds": problem.bounds.as_record(),
        "command": "python verification/fd_optimisation/run_fd_slsqp.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "variables": names,
        "fd_step": h,
        "envelope_margin": problem.margin,
        "envelope_margin_raised": margin_raised,
        "domain_errors": domain_errors,
        "slsqp_options": {"ftol": args.ftol, "maxiter": args.maxiter},
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nit": int(result.nit),
        "nfev": int(result.nfev),
        "njev": int(result.njev),
        "objective_evaluations_total": int(problem.n_fun_evals),
        "callback_extra_evaluations": recorder.extra_evals,
        "wall_time_s": wall,
        "J0_mwh_per_yr": float(problem.J0),
        "aep_x0_mwh_per_yr": float(aep0),
        "aep_optimum_mwh_per_yr": float(aep_star),
        "delta_aep_pct_vs_x0": float(gain_pct),
        "expected_gain_pct": list(EXPECTED_GAIN_PCT),
        "suspected_defect": suspected,
        "u_star": [float(v) for v in u_star],
        "x_star": [float(v) for v in x_star],
        "x0": [float(v) for v in x0],
        "u0": [float(v) for v in u0],
        "fun_star": float(result.fun),
        "gradient_star_scaled": [float(v) for v in result.jac],
        "gradient_star_mwh_per_u": [float(v) for v in problem.unscale(result.jac)],
        "active_bounds": {"lower": [names[i] for i in active_lower],
                          "upper": [names[i] for i in active_upper]},
        "active_envelope_rows": active_env,
        "envelope_min_row_value_m": float(g_env.min()),
        "envelope_min_row": labels[int(np.argmin(g_env))],
        "active_set": problem.active_set(u_star),
        "post_check_x0": post0,
        "post_check_optimum": post,
    }
    with open(RESULT_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    with open(ITERATES_PATH, "w", encoding="utf-8") as handle:
        json.dump({"bounds_label": BOUNDS_LABEL, "fd_step": h,
                   "iterates": recorder.iterates}, handle, indent=1)

    plot(problem, x0, x_star, FIGURE_PATH)

    print(f"\n{result.message}")
    print(f"nit={result.nit} nfev={result.nfev} njev={result.njev} "
          f"objective evaluations={problem.n_fun_evals} wall={wall:.1f}s")
    print(f"AEP: {aep0:.6f} -> {aep_star:.6f} MWh/yr  ({gain_pct:+.3f} %)")
    print(f"active bounds: {summary['active_bounds']}")
    print(f"active envelope rows: {active_env}  (min row {summary['envelope_min_row']} "
          f"= {g_env.min():.4e} m)")
    print(f"alpha at optimum: {post['alpha_min_deg']:.2f}..{post['alpha_max_deg']:.2f} deg, "
          f"within={post['within']}; Re {post['reynolds_min']:.0f}..{post['reynolds_max']:.0f}")
    if suspected:
        print(f"\n*** {suspected} ***")
    print(f"wrote {RESULT_PATH}\n      {ITERATES_PATH}\n      {FIGURE_PATH}")


def _write_partial(domain_errors):
    with open(RESULT_PATH, "w", encoding="utf-8") as handle:
        json.dump({"bounds_label": BOUNDS_LABEL, "aborted": True,
                   "domain_errors": domain_errors}, handle, indent=1)


if __name__ == "__main__":
    main()
