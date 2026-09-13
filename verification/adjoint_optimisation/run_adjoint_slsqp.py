"""
Phase 3, B5: adjoint-driven SLSQP, end to end -- A4's run with the adjoint
gradient in place of central differences.

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10
    subject to the polar-cache Reynolds envelope (linear, 50 rows, margin 5 %)
    gradient   ScaledProblem.jac_adjoint  (discrete adjoint, Phase 3)

Everything else is A4's (`verification/fd_optimisation/run_fd_slsqp.py`):
`scipy.optimize.minimize(method="SLSQP", ftol=1e-8, maxiter=200)`, the same
starting point `x0`, the same envelope constraint and margin escalation,
the same iterate recording, the same sanity gates on the AEP gain, the same
provisional-bounds label. The one change is the `jac` argument.

The comparison with A4 is the point. the implementation plan §6 B5: the two
optima "must agree to SLSQP's tolerance; a different optimum is a suspected
silent gradient error -- investigate before reporting". `result.json`
records `||u*_adj - u*_fd||_inf`, the AEP difference, and the counters side
by side, and `comparison_agrees` is decided against the tolerances stated
below -- SLSQP's `ftol` on `fun`, and the spread the multi-start study
(`verification/fd_optimisation_multistart/`) measured between equally
converged optima of this flat objective in `u`.

Provisional bounds: `chord_max_m = 0.45 m` is a placeholder pending the
hub-radius / root-attachment decision; the optimum is "under provisional
bounds" in every artefact here.

Run from the repo root:

    python verification/adjoint_optimisation/run_adjoint_slsqp.py [--maxiter N]

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import datetime
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import Bounds, minimize

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))
sys.path.insert(0, os.path.join(REPO_ROOT, "tests"))

from design import BladeParameterisation, DesignBounds, clamped_knots  # noqa: E402
from gradients import ScaledProblem  # noqa: E402
from gradients.problem import DEFAULT_ENVELOPE_MARGIN  # noqa: E402
from objective import WeibullResource  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402
from test_parameterisation import PROVISIONAL_BOUNDS  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
FD_RESULT_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation", "result.json")
MULTISTART_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation_multistart", "results.json")
RESULT_PATH = os.path.join(_HERE, "result.json")
ITERATES_PATH = os.path.join(_HERE, "iterates.json")
FIGURE_PATH = os.path.join(_HERE, "optimised_blade.png")

PROVISIONAL_LABEL = ("under provisional bounds (chord_max_m = 0.45 m provisional; "
                     "chord_min_m, twist_min, twist_max grounded)")

ACTIVE_TOL = 1e-6
ENVELOPE_ACTIVE_TOL = 1e-6
EXPECTED_GAIN_PCT = (2.0, 6.0)
DEFECT_GAIN_PCT = 15.0

#: Agreement with A4. `ftol = 1e-8` is SLSQP's absolute tolerance on `fun`;
#: two runs that both stopped on it can differ in `fun` by a few `ftol`,
#: i.e. a few 1e-7 MWh/yr. In `u` the objective is flat at the optimum
#: (gradient ~1e-3 MWh/yr per u), and the multi-start study measured a
#: spread of 8.6e-4 in `u` across nine equally converged FD optima; the
#: `u` tolerance is that measured spread with a factor of ten, stated here
#: and checked against the multi-start JSON at run time.
FUN_AGREEMENT_TOL = 1e-7          # |fun*_adj - fun*_fd|, i.e. ~1e-6 MWh/yr
U_AGREEMENT_FACTOR = 10.0         # times the multi-start spread in u


def load_x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"]
                    + artefact["twist_control_points_rad"])


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def build_problem(margin):
    parameterisation = BladeParameterisation()
    bounds = DesignBounds(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist, **PROVISIONAL_BOUNDS)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config(),
                         margin=margin)


class Recorder:
    """A4's recorder with the adjoint in place of `jac_fd`."""

    def __init__(self, problem):
        self.problem = problem
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
        grad = self.problem.jac_adjoint(u)
        self._last_jac = (np.array(u, dtype=float), grad)
        return grad

    def _cached(self, cache, u, compute):
        if cache is not None and np.array_equal(cache[0], u):
            return cache[1]
        self.extra_evals += 1
        return compute(u)

    def _record(self, u, value, grad):
        self.iterates.append({
            "k": len(self.iterates),
            "u": [float(x) for x in u],
            "J_mwh_per_yr": float(self.problem.unscale(value)),
            "fun": float(value),
            "g_scaled": [float(x) for x in grad],
            "g_mwh_per_u": [float(x) for x in self.problem.unscale(grad)],
            "wall_time_s": time.perf_counter() - self.started,
        })
        return self.iterates[-1]

    def callback(self, u):
        u = np.array(u, dtype=float)
        value = self._cached(self._last_fun, u, self.problem.fun)
        grad = self._cached(self._last_jac, u, self.jac)
        record = self._record(u, value, grad)
        print(f"  k={record['k']:3d}  J={record['J_mwh_per_yr']:.6f} MWh/yr  "
              f"|g|={np.linalg.norm(grad):.3e}  t={record['wall_time_s']:.1f}s", flush=True)

    def record_start(self, u0, value, grad):
        self._record(u0, value, grad)


def run(problem, u0, maxiter, ftol):
    envelope = problem.envelope_constraint()
    recorder = Recorder(problem)

    value0 = recorder.fun(u0)
    grad0 = recorder.jac(u0)
    recorder.record_start(u0, value0, grad0)
    print(f"J(x0) = {problem.J0:.6f} MWh/yr; gradient = adjoint; margin = {problem.margin}")

    started = time.perf_counter()
    result = minimize(
        recorder.fun, u0, jac=recorder.jac, method="SLSQP",
        bounds=Bounds(np.zeros(problem.n), np.ones(problem.n)),
        constraints=[envelope],
        options=dict(ftol=ftol, maxiter=maxiter, disp=True),
        callback=recorder.callback,
    )
    wall = time.perf_counter() - started
    return result, recorder, wall


def plot(problem, x0, x_star, x_fd, path):
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
        ax.plot(r, get(x0), "-", color="#888888", lw=1.6, label="x0 (fitted Schmitz)")
        ax.plot(r, get(x_fd), "-", color="#d1495b", lw=2.6, alpha=0.45, label="FD-SLSQP optimum (A4)")
        ax.plot(r, get(x_star), "-", color="#1f5fbf", lw=1.4, label="adjoint-SLSQP optimum")
        scale = np.degrees(1.0) if label == "twist" else 1.0
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

    fig.suptitle("Adjoint-driven SLSQP optimum vs x0 and the FD optimum -- " + PROVISIONAL_LABEL,
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--ftol", type=float, default=1e-8)
    args = parser.parse_args(argv)

    x0 = load_x0()
    fd = load_json(FD_RESULT_PATH)
    multistart = load_json(MULTISTART_PATH)

    margin = DEFAULT_ENVELOPE_MARGIN
    margin_raised = False
    domain_errors = []
    while True:
        problem = build_problem(margin)
        u0 = problem.scaled(x0)
        try:
            result, recorder, wall = run(problem, u0, args.maxiter, args.ftol)
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

    suspected = None
    if gain_pct > DEFECT_GAIN_PCT:
        suspected = (f"AEP gain {gain_pct:.2f} % exceeds {DEFECT_GAIN_PCT} %: suspected "
                     "defect (unfair baseline, units, or gradient) -- not a result")
    elif gain_pct < 0.0:
        suspected = f"AEP gain {gain_pct:.2f} % is negative: suspected defect -- not a result"

    # -- comparison with A4 -------------------------------------------------
    u_fd = np.array(fd["u_star"], dtype=float)
    x_fd = problem.physical(u_fd)
    du_inf = float(np.max(np.abs(u_star - u_fd)))
    d_aep = float(aep_star - fd["aep_optimum_mwh_per_yr"])
    d_fun = float(result.fun - fd["fun_star"])
    # the FD gradient's own value of J at the adjoint optimum, and vice versa,
    # so the comparison does not depend on either run's last evaluation
    aep_fd_here = -problem.J(u_fd)
    u_spread = float(multistart["spread_of_optima_u_inf"])
    u_tol = U_AGREEMENT_FACTOR * u_spread
    agrees_fun = abs(d_fun) <= FUN_AGREEMENT_TOL
    agrees_u = du_inf <= u_tol
    comparison = {
        "fd_result": os.path.relpath(FD_RESULT_PATH, REPO_ROOT),
        "u_star_fd": [float(v) for v in u_fd],
        "x_star_fd": [float(v) for v in x_fd],
        "du_inf": du_inf,
        "du_inf_per_variable": {names[j]: float(u_star[j] - u_fd[j]) for j in range(problem.n)},
        "aep_fd_mwh_per_yr": float(fd["aep_optimum_mwh_per_yr"]),
        "aep_adjoint_mwh_per_yr": float(aep_star),
        "delta_aep_mwh_per_yr": d_aep,
        "delta_fun": d_fun,
        "aep_at_fd_optimum_evaluated_here_mwh_per_yr": float(aep_fd_here),
        "delta_aep_pct_fd": float(fd["delta_aep_pct_vs_x0"]),
        "delta_aep_pct_adjoint": float(gain_pct),
        "nit_fd": int(fd["nit"]), "nfev_fd": int(fd["nfev"]), "njev_fd": int(fd["njev"]),
        "objective_evaluations_total_fd": int(fd["objective_evaluations_total"]),
        "wall_time_s_fd": float(fd["wall_time_s"]),
        "fun_agreement_tol": FUN_AGREEMENT_TOL,
        "multistart_u_spread_inf": u_spread,
        "u_agreement_tol": u_tol,
        "agrees_fun": bool(agrees_fun),
        "agrees_u": bool(agrees_u),
        "comparison_agrees": bool(agrees_fun and agrees_u),
    }

    summary = {
        "description": "Adjoint-driven SLSQP on fun(u) = J(u)/|J(u0)| with the polar-cache "
                       "envelope constraint, from x0 -- A4's run with jac = jac_adjoint.",
        "provisional_bounds": PROVISIONAL_LABEL,
        "bounds": {k: float(v) for k, v in PROVISIONAL_BOUNDS.items()},
        "command": "python verification/adjoint_optimisation/run_adjoint_slsqp.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "variables": names,
        "gradient": "discrete adjoint (ScaledProblem.jac_adjoint)",
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
        "adjoint_evaluations_total": int(problem.n_adjoint_evals),
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
        "post_check_optimum": post,
        "comparison_with_fd": comparison,
    }
    with open(RESULT_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    with open(ITERATES_PATH, "w", encoding="utf-8") as handle:
        json.dump({"provisional_bounds": PROVISIONAL_LABEL, "gradient": "adjoint",
                   "iterates": recorder.iterates}, handle, indent=1)

    plot(problem, x0, x_star, x_fd, FIGURE_PATH)

    print(f"\n{result.message}")
    print(f"nit={result.nit} nfev={result.nfev} njev={result.njev} "
          f"objective evaluations={problem.n_fun_evals} adjoint evaluations="
          f"{problem.n_adjoint_evals} wall={wall:.1f}s")
    print(f"AEP: {aep0:.6f} -> {aep_star:.6f} MWh/yr  ({gain_pct:+.3f} %)")
    print(f"active bounds: {summary['active_bounds']}; active envelope rows: {active_env}")
    print(f"alpha at optimum: {post['alpha_min_deg']:.2f}..{post['alpha_max_deg']:.2f} deg, "
          f"within={post['within']}; Re {post['reynolds_min']:.0f}..{post['reynolds_max']:.0f}")
    print(f"\nvs A4 (FD): ||du||_inf = {du_inf:.2e} (tol {u_tol}), "
          f"dAEP = {d_aep:+.3e} MWh/yr, dfun = {d_fun:+.3e} (tol {FUN_AGREEMENT_TOL:g}); "
          f"A4 nit/nfev/wall = {fd['nit']}/{fd['nfev']}/{fd['wall_time_s']:.0f}s; "
          f"agrees = {comparison['comparison_agrees']}")
    if suspected:
        print(f"\n*** {suspected} ***")
    if not comparison["comparison_agrees"]:
        print("\n*** optimum differs from A4 beyond tolerance: suspected silent gradient "
              "error -- investigate before reporting ***")
    print(f"wrote {RESULT_PATH}\n      {ITERATES_PATH}\n      {FIGURE_PATH}")


def _write_partial(domain_errors):
    with open(RESULT_PATH, "w", encoding="utf-8") as handle:
        json.dump({"provisional_bounds": PROVISIONAL_LABEL, "aborted": True,
                   "domain_errors": domain_errors}, handle, indent=1)


if __name__ == "__main__":
    main()
