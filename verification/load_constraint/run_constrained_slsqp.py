"""
Phase 4, Step 2d -- SLSQP on the AEP objective with the root-moment KS
constraint at reduction fraction `eps`.

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10
    subject to the polar-cache Reynolds envelope (50 rows) and the local
               solidity cap (25 rows)       -- `problem.constraints()`, unchanged
               g_eps(u) = (1 - eps) KS0 - KS(u) >= 0
                                             -- `problem.moment_constraint(eps)`,
                                                adjoint Jacobian
    gradient   discrete adjoint (`ScaledProblem.jac_adjoint`)
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-8, maxiter=200,
               from x0 (cold)

The recorder is B5's, with `KS_k` added. The summary reports the optimum, the
AEP and moment at `x0`, at the unconstrained optimum and at this constrained
one, the full active set (including the moment row's slack), a KKT check on the
slope, the post-checks, the SLSQP counters, and the moment-solve count. The
`--eps 0` run is the Phase 5 production optimum for now.

Outputs, next to this script: `result_eps{E}.json`, `iterates_eps{E}.json`,
`constrained_blade_eps{E}.png` (`E = 0, 0.02, 0.05, 0.1`).

Run from the repo root (about one minute with adjoint gradients):

    python verification/load_constraint/run_constrained_slsqp.py --eps 0
    python verification/load_constraint/run_constrained_slsqp.py --eps 0.05

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import datetime
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import Bounds, minimize

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from bem.rotor import solve_rotor  # noqa: E402
from design import BladeParameterisation, DesignBounds, clamped_knots  # noqa: E402
from gradients import ScaledProblem  # noqa: E402
from gradients.problem import DEFAULT_ENVELOPE_MARGIN  # noqa: E402
from objective import WeibullResource  # noqa: E402
from objective.power import tsr_schedule  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
ADJOINT_RESULT_PATH = os.path.join(REPO_ROOT, "verification", "adjoint_optimisation", "result.json")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")
LAW_LABEL = ("lambda(V) = min(6.5, Omega_max R / V), Omega_max = 300 rpm "
             "(V_c = 9.67 m/s), fixed rating 3822.189755449124 W")

ACTIVE_TOL = 1e-6

#: The stated sanity band for the eps = 0 AEP gain, from HANDOFF-2026-09-19
#: section 3.1: it cannot exceed the unconstrained +0.14666 % and losing more
#: than a third of it to a 0.31 % moment cap is a suspected defect.
EXPECTED_GAIN_PCT = (0.10, 0.14666)
#: The eps = 0.10 AEP-cost band; above 5 % is a suspected defect.
EXPECTED_COST_AT_10_PCT = (0.3, 3.0)
DEFECT_COST_PCT = 5.0


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def load_x0():
    artefact = load_json(X0_PATH)
    return np.array(artefact["chord_control_points_m"]
                    + artefact["twist_control_points_rad"])


def src_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       cwd=REPO_ROOT, text=True).strip()
    except Exception:  # pragma: no cover - reporting only
        return None


def build_problem(margin=DEFAULT_ENVELOPE_MARGIN):
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config(),
                         margin=margin)


def variable_names(problem):
    return ([f"chord_{j}" for j in range(problem.parameterisation.n_chord)]
            + [f"twist_{j}" for j in range(problem.parameterisation.n_twist)])


def loads_at(problem, d):
    """Forward-path loads at the design condition (11 m/s on the ceiling)."""

    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    v_inf = float(problem.design.rated_wind_speed_ms)
    tsr = tsr_schedule(v_inf, problem.design.design_tsr, problem.design.max_tip_speed_ms)
    result = solve_rotor(geometry, tsr=tsr, v_inf=v_inf,
                         air_density=problem.site.air_density,
                         kinematic_viscosity=problem.site.kinematic_viscosity)
    area = math.pi * geometry.R ** 2
    from objective.loads import root_bending_moment
    moment = root_bending_moment(result["stations"], geometry.chord,
                                 problem.site.air_density, geometry.n_blades,
                                 geometry.r_hub)
    return {
        "v_ms": v_inf, "tsr": float(tsr),
        "rpm": float(tsr * v_inf / geometry.R * 60.0 / (2.0 * math.pi)),
        "Cp": float(result["Cp"]), "Ct": float(result["Ct"]),
        "rotor_thrust_n": float(result["Ct"] * 0.5 * problem.site.air_density
                                * v_inf ** 2 * area),
        "root_bending_moment_nm": float(moment),
        "converged": bool(result["converged"]),
    }


class Recorder:
    """B5's recorder with `KS_k` added."""

    def __init__(self, problem, eps):
        self.problem = problem
        self.eps = eps
        self.started = time.perf_counter()
        self.iterates = []
        self._last_fun = None
        self._last_jac = None
        self._last_ks = None
        self.extra_evals = 0

    def fun(self, u):
        value = self.problem.fun(u)
        self._last_fun = (np.array(u, dtype=float), value)
        return value

    def jac(self, u):
        grad = self.problem.jac_adjoint(u)
        self._last_jac = (np.array(u, dtype=float), grad)
        return grad

    def ks(self, u):
        value = self.problem.moment_ks(u)
        self._last_ks = (np.array(u, dtype=float), value)
        return value

    def _cached(self, cache, u, compute):
        if cache is not None and np.array_equal(cache[0], u):
            return cache[1]
        self.extra_evals += 1
        return compute(u)

    def _record(self, u, value, grad, ks):
        self.iterates.append({
            "k": len(self.iterates),
            "u": [float(x) for x in u],
            "J_mwh_per_yr": float(self.problem.unscale(value)),
            "fun": float(value),
            "g_scaled": [float(x) for x in grad],
            "g_mwh_per_u": [float(x) for x in self.problem.unscale(grad)],
            "KS": float(ks),
            "moment_active_slack": float((self.problem.KS0 - ks) - self.eps * self.problem.KS0),
            "wall_time_s": time.perf_counter() - self.started,
        })
        return self.iterates[-1]

    def callback(self, u):
        u = np.array(u, dtype=float)
        value = self._cached(self._last_fun, u, self.problem.fun)
        grad = self._cached(self._last_jac, u, self.jac)
        ks = self._cached(self._last_ks, u, self.ks)
        record = self._record(u, value, grad, ks)
        print(f"  k={record['k']:3d}  J={record['J_mwh_per_yr']:.6f} MWh/yr  "
              f"KS={record['KS']:.8f}  |g|={np.linalg.norm(grad):.3e}  "
              f"t={record['wall_time_s']:.1f}s", flush=True)

    def record_start(self, u0, value, grad, ks):
        self._record(u0, value, grad, ks)


def run_constrained(problem, u0, eps, maxiter=200, ftol=1e-8, x_start=None):
    """One SLSQP run for reduction fraction `eps`. Returns `(result, recorder, wall)`."""

    constraints = problem.constraints_with_moment(eps)
    recorder = Recorder(problem, eps)

    # Iterate 0 is the point SLSQP actually starts from: x0 cold, or the
    # previous eps's optimum when warm-started by run_pareto.py.
    start = u0 if x_start is None else np.asarray(x_start, dtype=float)
    value0 = recorder.fun(start)
    grad0 = recorder.jac(start)
    ks0 = recorder.ks(start)
    recorder.record_start(start, value0, grad0, ks0)
    print(f"J(x0) = {problem.J0:.6f} MWh/yr; KS0 = {problem.KS0!r}; "
          f"eps = {eps:g}; margin = {problem.margin}; "
          f"start = {'warm' if x_start is not None else 'x0'}", flush=True)

    started = time.perf_counter()
    result = minimize(
        recorder.fun, start, jac=recorder.jac, method="SLSQP",
        bounds=Bounds(np.zeros(problem.n), np.ones(problem.n)),
        constraints=constraints,
        options=dict(ftol=ftol, maxiter=maxiter, disp=True),
        callback=recorder.callback,
    )
    wall = time.perf_counter() - started
    return result, recorder, wall


def kkt_report(problem, u, eps, g_scaled, names):
    """
    SciPy's SLSQP does not expose multipliers, so estimate them: the
    least-squares coefficients of `grad fun` on the active constraint/bound
    normals are the multipliers (all must be >= 0 at a minimiser), the
    projection is the part they explain, and the residual after removing it is
    zero at a KKT point up to the solver tolerance.
    """

    u = np.asarray(u, dtype=float)
    rows, labels = [], []

    envelope = problem.envelope_constraint()
    g_env = envelope["fun"](u)
    jac_env = envelope["jac"](u)
    for k in np.flatnonzero(g_env < ACTIVE_TOL):
        labels.append(problem.envelope_row_labels()[k])
        rows.append(jac_env[k])

    solidity = problem.solidity_constraint()
    g_sol = solidity["fun"](u)
    jac_sol = solidity["jac"](u)
    for k in np.flatnonzero(g_sol < ACTIVE_TOL):
        labels.append(problem.solidity_row_labels()[k])
        rows.append(jac_sol[k])

    moment = problem.moment_active(u, eps)
    if moment["active"]:
        labels.append("moment")
        rows.append(problem.moment_constraint(eps)["jac"](u)[0])

    for j in range(problem.n):
        if abs(u[j]) < ACTIVE_TOL:
            labels.append(f"lower:{names[j]}")
            row = np.zeros(problem.n)
            row[j] = 1.0
            rows.append(row)
        elif abs(1.0 - u[j]) < ACTIVE_TOL:
            labels.append(f"upper:{names[j]}")
            row = np.zeros(problem.n)
            row[j] = -1.0
            rows.append(row)

    # Every row is written so that its KKT multiplier is non-negative at a
    # minimiser: constraints as `g >= 0` (SciPy's convention, `grad fun = sum
    # lambda_i grad g_i`), lower bounds as `+e_j`, upper bounds as `-e_j`. The
    # least-squares coefficients ARE the multiplier estimates; a negative one
    # would mean the row is pushing the wrong way and the point is not KKT.
    g = np.asarray(g_scaled, dtype=float)
    if rows:
        A = np.vstack(rows)
        multipliers, *_ = np.linalg.lstsq(A.T, g, rcond=None)
        projection = A.T @ multipliers
    else:
        A = np.zeros((0, problem.n))
        multipliers = np.zeros(0)
        projection = np.zeros(problem.n)
    residual = g - projection
    free = np.array([abs(u[j]) > ACTIVE_TOL and abs(1.0 - u[j]) > ACTIVE_TOL
                     for j in range(problem.n)])
    multiplier_by_row = {label: float(value) for label, value in zip(labels, multipliers)}
    moment_multiplier = multiplier_by_row.get("moment")
    # The moment row's multiplier is the shadow price of the cap: in scaled
    # units it is d(fun)/d(g), and `|J0|` times it is the AEP given up per unit
    # of KS -- the number that says whether a reduction "costs nothing".
    shadow_price = (None if moment_multiplier is None
                    else float(moment_multiplier * abs(problem.J0)))

    return {
        "active_rows": labels,
        "n_active_rows": len(labels),
        "multipliers_scaled": multiplier_by_row,
        "all_multipliers_nonnegative": bool(np.all(multipliers >= -1e-12)),
        "moment_multiplier_scaled": moment_multiplier,
        "moment_shadow_price_mwh_per_yr_per_unit_ks": shadow_price,
        "gradient_norm_scaled": float(np.linalg.norm(g)),
        "projection_norm_onto_active_normals": float(np.linalg.norm(projection)),
        "residual_norm_after_projection": float(np.linalg.norm(residual)),
        "free_variables": [names[j] for j in range(problem.n) if free[j]],
        "gradient_norm_on_free_variables": float(np.linalg.norm(g[free])) if free.any() else 0.0,
    }


def plot(problem, x0, x_star, x_c, eps, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = problem.parameterisation
    radii = p.radii
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    for ax, (label, unit, get, block) in zip(axes, [
        ("chord", "m", p.chord, slice(0, p.n_chord)),
        ("twist", "deg", lambda d: np.degrees(p.twist(d)), slice(p.n_chord, None)),
    ]):
        ax.plot(radii, get(x0), "-", color="#888888", lw=1.6, label="x0 (fitted Schmitz)")
        ax.plot(radii, get(x_star), "--", color="#d1495b", lw=1.8,
                label="unconstrained optimum u*")
        ax.plot(radii, get(x_c), "-", color="#1f5fbf", lw=2.4,
                label=f"constrained eps = {eps:g}")
        scale = np.degrees(1.0) if label == "twist" else 1.0
        n_pts = p.n_chord if label == "chord" else p.n_twist
        knots = clamped_knots(n_pts, p.degree)
        greville = np.array([knots[j + 1:j + p.degree + 1].mean() for j in range(n_pts)])
        s = radii[0] + greville * (radii[-1] - radii[0])
        ax.plot(s, x0[block] * scale, "o", color="#888888", ms=5, mfc="white")
        ax.plot(s, x_c[block] * scale, "o", color="#1f5fbf", ms=5, mfc="white")
        ax.set_xlabel("radius r [m]")
        ax.set_ylabel(f"{label} [{unit}]")
        ax.set_title(f"{label} distribution", fontsize=11)
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    fig.suptitle(f"Phase 4 constrained optimum (eps = {eps:g}, moment cap) vs x0 and u* -- "
                 + BOUNDS_LABEL, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def summarise(problem, u0, eps, result, recorder, wall, args, margin_raised, domain_errors):
    x0 = load_x0()
    unconstrained = load_json(ADJOINT_RESULT_PATH)
    u_star = np.array(unconstrained["u_star"], dtype=float)
    x_star = problem.physical(u_star)
    aep_star = float(unconstrained["aep_optimum_mwh_per_yr"])
    names = variable_names(problem)

    u_c = np.array(result.x, dtype=float)
    x_c = problem.physical(u_c)
    aep0 = -problem.J0
    aep_c = -problem.J(u_c)
    gain_pct = 100.0 * (aep_c - aep0) / aep0
    # The constraint's cost is the AEP given up relative to the UNCONSTRAINED
    # optimum (section 3.1's definition), not the change relative to x0.
    cost_vs_unconstrained = 100.0 * (aep_star - aep_c) / aep_star

    ks0 = float(problem.KS0)
    ks_c = problem.moment_ks(u_c)
    moment_row = problem.moment_active(u_c, eps)

    loads = {
        "x0": loads_at(problem, x0),
        "x_star_unconstrained": loads_at(problem, x_star),
        "x_c": loads_at(problem, x_c),
    }
    base_moment = loads["x0"]["root_bending_moment_nm"]
    base_thrust = loads["x0"]["rotor_thrust_n"]
    for entry in loads.values():
        entry["root_bending_moment_pct_vs_x0"] = 100.0 * (
            entry["root_bending_moment_nm"] / base_moment - 1.0)
        entry["rotor_thrust_pct_vs_x0"] = 100.0 * (entry["rotor_thrust_n"] / base_thrust - 1.0)

    suspected = None
    if eps == 0.0:
        if not (EXPECTED_GAIN_PCT[0] <= gain_pct <= EXPECTED_GAIN_PCT[1]):
            suspected = (f"eps = 0 AEP gain {gain_pct:+.3f} % is outside the stated "
                         f"band {EXPECTED_GAIN_PCT} -- a finding, not a defect gate")
        if not moment_row["active"] and moment_row["slack"] > ACTIVE_TOL:
            suspected = (suspected + "; " if suspected else "") + (
                f"eps = 0 constraint is not active at the optimum (slack "
                f"{moment_row['slack']:.3e})")
    within_eps10_band = None
    if eps >= 0.10:
        within_eps10_band = bool(EXPECTED_COST_AT_10_PCT[0] <= cost_vs_unconstrained
                                 <= EXPECTED_COST_AT_10_PCT[1])
        if cost_vs_unconstrained > DEFECT_COST_PCT:
            suspected = (f"eps = {eps:g} AEP cost {cost_vs_unconstrained:.2f} % exceeds "
                         f"{DEFECT_COST_PCT} %: suspected defect")

    return {
        "description": f"Constrained SLSQP with the root-moment KS constraint at eps = {eps:g}.",
        "bounds_label": BOUNDS_LABEL,
        "law_label": LAW_LABEL,
        "command": f"python verification/load_constraint/run_constrained_slsqp.py --eps {eps:g}",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": src_commit(),
        "bounds": problem.bounds.as_record(),
        "eps": float(eps),
        "rho": float(problem.load_system().rho),
        "m_ref_nm": float(problem.m_ref),
        "KS0": ks0,
        "KS_c": ks_c,
        "KS_reduction": ks0 - ks_c,
        "KS_reduction_over_eps_ks0": (ks0 - ks_c) / (eps * ks0) if eps > 0 else None,
        "envelope_margin": problem.margin,
        "envelope_margin_raised": bool(margin_raised),
        "domain_errors": domain_errors,
        "slsqp_options": {"ftol": args.ftol, "maxiter": args.maxiter},
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nit": int(result.nit), "nfev": int(result.nfev), "njev": int(result.njev),
        "objective_evaluations_total": int(problem.n_fun_evals),
        "adjoint_evaluations_total": int(problem.n_adjoint_evals),
        "moment_solves_total": int(problem.n_moment_solves),
        "callback_extra_evaluations": recorder.extra_evals,
        "wall_time_s": wall,
        "J0_mwh_per_yr": float(problem.J0),
        "aep_x0_mwh_per_yr": float(aep0),
        "aep_optimum_mwh_per_yr": float(aep_c),
        "aep_star_unconstrained_mwh_per_yr": float(aep_star),
        "delta_aep_pct_vs_x0": float(gain_pct),
        "aep_cost_pct_vs_unconstrained": float(cost_vs_unconstrained),
        "expected_gain_pct_eps0": list(EXPECTED_GAIN_PCT),
        "expected_cost_pct_eps10": list(EXPECTED_COST_AT_10_PCT),
        "within_expected_cost_band_eps10": within_eps10_band,
        "suspected_defect": suspected,
        "u_c": [float(x) for x in u_c],
        "x_c": [float(x) for x in x_c],
        "x0": [float(x) for x in x0],
        "u0": [float(x) for x in u0],
        "u_star_unconstrained": [float(x) for x in u_star],
        "x_star_unconstrained": [float(x) for x in x_star],
        "fun_c": float(result.fun),
        "gradient_c_scaled": [float(x) for x in result.jac],
        "gradient_c_mwh_per_u": [float(x) for x in problem.unscale(result.jac)],
        "moment_row": moment_row,
        "active_set": problem.active_set(u_c, moment_eps=eps),
        "kkt": kkt_report(problem, u_c, eps, result.jac, names),
        "loads": loads,
        "post_check_optimum": problem.alpha_check(u_c),
    }


def optimise_with_escalation(eps, args):
    """
    Run `eps`, raising the envelope margin to 0.10 once and restarting if a
    `PolarDomainError` escapes -- B5's rule, recorded. Returns
    `(problem, u0, result, recorder, wall, margin_raised, domain_errors)`.
    """

    x0 = load_x0()
    margin = DEFAULT_ENVELOPE_MARGIN
    margin_raised = False
    domain_errors = []
    while True:
        problem = build_problem(margin)
        u0 = problem.scaled(x0)
        problem.set_reference(u0)
        problem.load_system()
        try:
            result, recorder, wall = run_constrained(problem, u0, eps,
                                                     args.maxiter, args.ftol)
            return problem, u0, result, recorder, wall, margin_raised, domain_errors
        except PolarDomainError as error:
            domain_errors.append({"margin": margin, "u": problem.domain_errors[-1],
                                  "error": str(error)})
            if margin_raised:
                raise
            print(f"PolarDomainError at margin {margin}: raising margin to 0.10 once "
                  f"and restarting", flush=True)
            margin, margin_raised = 0.10, True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--eps", type=float, default=0.0,
                        help="required KS reduction fraction (0, 0.02, 0.05, 0.10)")
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--ftol", type=float, default=1e-8)
    args = parser.parse_args(argv)

    problem, u0, result, recorder, wall, margin_raised, domain_errors = \
        optimise_with_escalation(args.eps, args)
    summary = summarise(problem, u0, args.eps, result, recorder, wall, args,
                        margin_raised, domain_errors)

    tag = f"{args.eps:g}"
    result_path = os.path.join(_HERE, f"result_eps{tag}.json")
    iterates_path = os.path.join(_HERE, f"iterates_eps{tag}.json")
    figure_path = os.path.join(_HERE, f"constrained_blade_eps{tag}.png")

    with open(result_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    with open(iterates_path, "w", encoding="utf-8") as handle:
        json.dump({"bounds_label": BOUNDS_LABEL, "eps": args.eps,
                   "gradient": "adjoint", "iterates": recorder.iterates},
                  handle, indent=1)
    plot(problem, np.array(summary["x0"], dtype=float),
         np.array(summary["x_star_unconstrained"], dtype=float),
         np.array(summary["x_c"], dtype=float), args.eps, figure_path)

    print(f"\n{result.message}")
    print(f"eps={args.eps:g}: nit={result.nit} nfev={result.nfev} njev={result.njev} "
          f"moment solves={problem.n_moment_solves} wall={wall:.1f}s")
    print(f"AEP: {summary['aep_x0_mwh_per_yr']:.6f} -> "
          f"{summary['aep_optimum_mwh_per_yr']:.6f} MWh/yr "
          f"({summary['delta_aep_pct_vs_x0']:+.3f} % vs x0; cost "
          f"{summary['aep_cost_pct_vs_unconstrained']:.4f} % vs the unconstrained optimum)")
    print(f"KS: {summary['KS0']:.8f} -> {summary['KS_c']:.8f} "
          f"(reduction {summary['KS_reduction']:.3e}, "
          f"{100.0 * summary['KS_reduction'] / summary['KS0']:+.3f} %)")
    print(f"active set: {summary['active_set']['bounds_lower']} lower, "
          f"{summary['active_set']['bounds_upper']} upper; moment slack "
          f"{summary['moment_row']['slack']:.3e}")
    kkt = summary["kkt"]
    print(f"KKT: multipliers {kkt['multipliers_scaled']} (all >= 0: "
          f"{kkt['all_multipliers_nonnegative']}); moment shadow price "
          f"{kkt['moment_shadow_price_mwh_per_yr_per_unit_ks']} MWh/yr per unit KS; "
          f"residual {kkt['residual_norm_after_projection']:.2e}")
    if summary["suspected_defect"]:
        print(f"\n*** {summary['suspected_defect']} ***")
    print(f"wrote {result_path}\n      {iterates_path}\n      {figure_path}")
    return summary


if __name__ == "__main__":
    main()
