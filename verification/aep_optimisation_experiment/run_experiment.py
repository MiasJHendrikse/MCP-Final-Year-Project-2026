"""
EXPERIMENT (2026-09-19, not a project result): adjoint-driven SLSQP on the
site-specific AEP objective under a rotor-speed ceiling.

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10             (the project's 5 + 5 control points)
    subject to the polar-cache Reynolds envelope (linear, 50 rows, margin 5 %)
    gradient   CeilingProblem.jac_adjoint  (discrete adjoint, per-bin lambda)
    schedule   lambda(V) = min(6.5, V_tip,max / V)

One case per ceiling in {none, 60, 55, 50} m/s. `none` is the control: it
is the project's own objective and must reproduce
`verification/adjoint_optimisation/result.json`. Everything about the
optimiser is B5's (`run_adjoint_slsqp.py`): SLSQP, `ftol = 1e-8`,
`maxiter = 200`, the envelope constraint, `x0` as the start, the same
iterate recorder. Only the objective and its adjoint differ, and they
differ only by the per-bin TSR (`ceiling_objective.py`).

After each optimisation the result is checked, independently of the
optimiser's own evaluations:

  1. AEP at x_AEP* re-evaluated through a fresh parameterisation and
     resource, and re-assembled from the station integrand `q` by the
     adjoint system (a second code path); for the control, also through
     `objective.annual_energy_mwh` in `src/`.
  2. Adjoint vs central FD at x_AEP* (two steps, the project's 3e-6 and
     1e-5 in `u`), and again at x0.
  3. SLSQP status, iteration and evaluation counts, gradient at the optimum.
  4. Every `u` in [0, 1]; active bounds and active envelope rows listed.
  5. BEM convergence at every bin at x_AEP*; alpha inside the polar cache.
  6. The bin masses sum to CDF(20) - CDF(3); the per-bin energies sum to the
     total; the adjoint's weights equal -(T/1e6) m_b on the uncapped bins.
  7. Against the FD-driven optimum the audit found for the same ceiling
     (`verification/aep_gain_audit/opt_<v>_fixed.json`, rating 3822.2 W
     rounded): ||du||_inf and the AEP difference, both evaluated here.

Provisional bounds (`chord_max_m = 0.45 m`), as every optimisation
artefact in this repository. Run from the repo root:

    python verification/aep_optimisation_experiment/run_experiment.py [--cases none,60,55,50]

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
sys.path.insert(0, os.path.join(REPO_ROOT, "tests"))
sys.path.insert(0, _HERE)

from config import load_design_rotor, load_site  # noqa: E402
from design import BladeParameterisation, DesignBounds  # noqa: E402
from objective import WeibullResource  # noqa: E402
from objective.objective import HOURS_PER_YEAR, annual_energy_mwh  # noqa: E402
from objective.power import aerodynamic_power, wind_speed_bins  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402
from test_parameterisation import PROVISIONAL_BOUNDS  # noqa: E402

from ceiling_objective import CeilingProblem, evaluate_aep, tsr_schedule  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
XSTAR_PATH = os.path.join(REPO_ROOT, "verification", "adjoint_optimisation", "result.json")
AUDIT_DIR = os.path.join(REPO_ROOT, "verification", "aep_gain_audit")
RESULTS_DIR = os.path.join(_HERE, "results")

PROVISIONAL_LABEL = ("under provisional bounds (chord_max_m = 0.45 m provisional; "
                     "chord_min_m, twist_min, twist_max grounded)")
CASES = {"none": None, "60": 60.0, "55": 55.0, "50": 50.0}
FD_STEPS = (3e-6, 1e-5)
ACTIVE_TOL = 1e-6
ENVELOPE_ACTIVE_TOL = 1e-6
#: The wind speed the Cp reporting is done at: the spanwise
#: reference speed, near the energy-weighted centre of the below-rated bins.
CP_REPORT_SPEED_MS = 8.5
CP_SWEEP_TSR = [3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5]


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def load_x0():
    artefact = load_json(X0_PATH)
    return np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"],
                    dtype=float)


def build_problem(vtip_max, margin=0.05):
    parameterisation = BladeParameterisation()
    bounds = DesignBounds(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist, **PROVISIONAL_BOUNDS)
    return CeilingProblem(parameterisation, bounds, WeibullResource.from_config(),
                          vtip_max=vtip_max, margin=margin)


class Recorder:
    """B5's iterate recorder, unchanged in substance."""

    def __init__(self, problem):
        self.problem = problem
        self.started = time.perf_counter()
        self.iterates = []
        self._last_fun = None
        self._last_jac = None

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
        return compute(u)

    def record(self, u, value, grad):
        self.iterates.append({
            "k": len(self.iterates),
            "u": [float(x) for x in u],
            "aep_mwh_per_yr": float(-self.problem.unscale(value)),
            "fun": float(value),
            "g_mwh_per_u": [float(x) for x in self.problem.unscale(grad)],
            "wall_time_s": time.perf_counter() - self.started,
        })
        return self.iterates[-1]

    def callback(self, u):
        u = np.array(u, dtype=float)
        value = self._cached(self._last_fun, u, self.problem.fun)
        grad = self._cached(self._last_jac, u, self.jac)
        rec = self.record(u, value, grad)
        print(f"    k={rec['k']:3d}  AEP={rec['aep_mwh_per_yr']:.6f} MWh/yr  "
              f"|g|={np.linalg.norm(grad):.3e}  t={rec['wall_time_s']:.1f}s", flush=True)


def optimise(problem, u0, maxiter=200, ftol=1e-8):
    envelope = problem.envelope_constraint()
    recorder = Recorder(problem)
    value0 = recorder.fun(u0)
    grad0 = recorder.jac(u0)
    recorder.record(u0, value0, grad0)
    started = time.perf_counter()
    result = minimize(
        recorder.fun, u0, jac=recorder.jac, method="SLSQP",
        bounds=Bounds(np.zeros(problem.n), np.ones(problem.n)),
        constraints=[envelope],
        options=dict(ftol=ftol, maxiter=maxiter, disp=False),
        callback=recorder.callback,
    )
    return result, recorder, time.perf_counter() - started


# -- verification ------------------------------------------------------------

def fd_gradient(problem, u, h):
    """Central differences of `fun` in `u`, the project's own convention."""

    u = np.asarray(u, dtype=float)
    grad = np.zeros(problem.n)
    for j in range(problem.n):
        up, um = u.copy(), u.copy()
        up[j] += h
        um[j] -= h
        grad[j] = (problem.fun(up) - problem.fun(um)) / (2.0 * h)
    return grad


def gradient_check(problem, u, label):
    adjoint = problem.jac_adjoint(u)
    out = {"point": label, "adjoint_scaled": [float(v) for v in adjoint],
           "adjoint_mwh_per_u": [float(v) for v in problem.unscale(adjoint)], "fd": {}}
    for h in FD_STEPS:
        fd = fd_gradient(problem, u, h)
        diff = np.abs(adjoint - fd)
        out["fd"][f"{h:g}"] = {
            "fd_scaled": [float(v) for v in fd],
            "abs_diff_max": float(diff.max()),
            "rel_diff_max_to_gradient_norm": float(diff.max() / np.linalg.norm(fd)),
            "rel_diff_per_component": [float(a / max(abs(f), 1e-300)) for a, f in zip(diff, fd)],
        }
    return out


def alpha_by_cap(problem, u):
    """
    Alpha extremes split by whether the bin is capped. A capped bin adds a
    constant to J and nothing to the gradient, so alpha outside the XFOIL
    band there cannot steer the optimiser; on an uncapped bin it can, and
    it also means the blade is being scored on the Viterna extrapolation.
    """

    state = problem.operating_state(u)
    p_rated = float(problem.design.rated_power_w)
    lo, hi = problem.alpha_min_deg, problem.alpha_max_deg
    unc = [i for i, p in enumerate(state["power_w"]) if p <= p_rated]
    cap = [i for i, p in enumerate(state["power_w"]) if p > p_rated]
    a_unc = np.array([state["alpha_deg"][i] for i in unc])
    a_cap = np.array([state["alpha_deg"][i] for i in cap]) if cap else np.zeros((0, 1))
    return {
        "uncapped_bins_ms": [state["speeds"][i] for i in unc],
        "alpha_uncapped_min_deg": float(a_unc.min()), "alpha_uncapped_max_deg": float(a_unc.max()),
        "uncapped_within_cache_band": bool(a_unc.min() >= lo and a_unc.max() <= hi),
        "uncapped_station_bins_outside_band": int(((a_unc < lo) | (a_unc > hi)).sum()),
        "alpha_capped_max_deg": float(a_cap.max()) if cap else None,
        "all_converged": bool(all(state["converged"])),
    }


def independent_aep(problem, x, vtip_max):
    """
    AEP at `x` by routes that share nothing with the optimiser's evaluation
    object: a fresh parameterisation and resource through `evaluate_aep`;
    the adjoint system's `J(phi, d)` assembled from `q`; and, when there is
    no ceiling, `src/objective.annual_energy_mwh` itself.
    """

    fresh_par = BladeParameterisation()
    fresh_res = WeibullResource.from_config()
    design, site = load_design_rotor(), load_site()
    aep_fresh, rows, converged = evaluate_aep(x, fresh_par, fresh_res, design, site, vtip_max)

    system = problem.adjoint_system()
    state = system.solve(x)
    aep_from_q = -float(np.real(system.J(state.phi, x, state.limited)))

    out = {"aep_fresh_forward_chain_mwh": aep_fresh,
           "aep_assembled_from_q_mwh": aep_from_q,
           "all_bins_converged": bool(converged),
           "bins": rows}
    if vtip_max is None:
        out["aep_src_objective_mwh"] = float(annual_energy_mwh(x, resource=fresh_res,
                                                               parameterisation=fresh_par))
    return out


def weighting_check(problem):
    """The bin masses and the adjoint's weights, against the resource's CDF."""

    edges, mids, _ = wind_speed_bins()
    res = problem.resource
    mass = np.asarray(res.probability_between(edges[:-1], edges[1:]), dtype=float)
    cdf_span = float(res.cdf(edges[-1]) - res.cdf(edges[0])) if hasattr(res, "cdf") else None
    system = problem.adjoint_system()
    omega_uncapped = system.weights(np.zeros(len(mids), dtype=bool))
    return {
        "n_bins": int(len(mids)),
        "edges_ms": [float(v) for v in edges],
        "midpoints_ms": [float(v) for v in mids],
        "mass": [float(m) for m in mass],
        "mass_sum": float(mass.sum()),
        "cdf_cut_out_minus_cut_in": cdf_span,
        "adjoint_weight_equals_minus_T_over_1e6_times_mass": bool(
            np.allclose(omega_uncapped, -HOURS_PER_YEAR / 1e6 * mass, rtol=0, atol=1e-15)),
        "hours_per_year": HOURS_PER_YEAR,
    }


def cp_metrics(problem, x, vtip_max):
    """Cp at lambda = 6.5 and a Cp-lambda sweep at the report speed; Ct too."""

    par = problem.parameterisation
    rho, nu = float(problem.site.air_density), float(problem.site.kinematic_viscosity)
    geometry = par.to_geometry(x)
    lam_opt = float(problem.design.design_tsr)

    _p, at_design = aerodynamic_power(geometry, CP_REPORT_SPEED_MS, lam_opt, rho, nu)
    sweep = []
    for lam in CP_SWEEP_TSR:
        _p, res = aerodynamic_power(geometry, CP_REPORT_SPEED_MS, lam, rho, nu)
        sweep.append({"tsr": lam, "Cp": float(res["Cp"]), "Ct": float(res["Ct"]),
                      "converged": bool(res["converged"])})
    best = max(sweep, key=lambda r: r["Cp"])

    # Energy-weighted mean Cp over the bins, on this case's schedule, weighted by
    # the *unlimited* energy each bin would carry (m_b V_b^3): the Cp the site sees.
    aep, rows, _ = evaluate_aep(x, par, problem.resource, problem.design, problem.site, vtip_max)
    w = np.array([r["mass"] * r["v_ms"] ** 3 for r in rows])
    cp = np.array([r["Cp"] for r in rows])
    return {
        "report_speed_ms": CP_REPORT_SPEED_MS,
        "Cp_at_design_tsr": float(at_design["Cp"]),
        "Ct_at_design_tsr": float(at_design["Ct"]),
        "Cp_lambda_sweep": sweep,
        "Cp_max_in_sweep": best["Cp"], "tsr_at_Cp_max": best["tsr"],
        "energy_weighted_mean_Cp_on_schedule": float((w * cp).sum() / w.sum()),
        "per_bin_Cp_on_schedule": [{"v_ms": r["v_ms"], "tsr": r["tsr"], "Cp": r["Cp"],
                                    "limited": r["limited"]} for r in rows],
    }


def geometry_comparison(par, names, x_new, x0, x_star):
    """Design variables side by side; chord in % and twist in degrees."""

    n_c = par.n_chord

    def delta(a, b):
        out = {}
        for j, name in enumerate(names):
            if j < n_c:
                out[name] = {"abs_mm": 1000.0 * (a[j] - b[j]),
                             "pct": 100.0 * (a[j] - b[j]) / b[j]}
            else:
                out[name] = {"abs_deg": math.degrees(a[j] - b[j])}
        return out

    radii = par.radii
    return {
        "variables": names,
        "x0": [float(v) for v in x0],
        "x_star": [float(v) for v in x_star],
        "x_aep": [float(v) for v in x_new],
        "x0_twist_deg": [math.degrees(v) for v in x0[n_c:]],
        "x_star_twist_deg": [math.degrees(v) for v in x_star[n_c:]],
        "x_aep_twist_deg": [math.degrees(v) for v in x_new[n_c:]],
        "x_aep_vs_x0": delta(x_new, x0),
        "x_aep_vs_x_star": delta(x_new, x_star),
        "x_star_vs_x0": delta(x_star, x0),
        "max_abs_diff_u_x_aep_vs_x0": None,  # filled by caller in u
        "spanwise": {
            "radii_m": [float(r) for r in radii],
            "chord_x0_m": [float(v) for v in par.chord(x0)],
            "chord_x_star_m": [float(v) for v in par.chord(x_star)],
            "chord_x_aep_m": [float(v) for v in par.chord(x_new)],
            "twist_x0_deg": [float(v) for v in np.degrees(par.twist(x0))],
            "twist_x_star_deg": [float(v) for v in np.degrees(par.twist(x_star))],
            "twist_x_aep_deg": [float(v) for v in np.degrees(par.twist(x_new))],
        },
    }


# -- one case ---------------------------------------------------------------

def run_case(label, vtip_max, x0, x_star, args):
    print(f"\n=== case {label}: V_tip,max = {vtip_max} ===", flush=True)
    problem = build_problem(vtip_max)
    par, bounds = problem.parameterisation, problem.bounds
    names = ([f"chord_{j}" for j in range(par.n_chord)] + [f"twist_{j}" for j in range(par.n_twist)])
    u0 = problem.scaled(x0)
    u_star_existing = problem.scaled(x_star)

    # The three blades under this objective before any optimisation.
    aep_x0 = -problem.J(u0)
    aep_xstar = -problem.J(u_star_existing)
    print(f"  AEP(x0) = {aep_x0:.6f}   AEP(x*) = {aep_xstar:.6f}   "
          f"x* vs x0 = {100 * (aep_xstar / aep_x0 - 1):+.3f} %", flush=True)

    # Gradient at the start, before the optimiser has touched anything.
    problem.set_reference(u0)
    grad_check_x0 = gradient_check(problem, u0, "x0")
    n_fun_before = problem.n_fun_evals
    n_adj_before = problem.n_adjoint_evals

    result, recorder, wall = optimise(problem, u0, args.maxiter, args.ftol)
    n_fun_opt = problem.n_fun_evals - n_fun_before
    n_adj_opt = problem.n_adjoint_evals - n_adj_before

    u_new = np.array(result.x, dtype=float)
    x_new = problem.physical(u_new)
    aep_new = -problem.J(u_new)
    print(f"  {result.message}: nit={result.nit} nfev={result.nfev} njev={result.njev}  "
          f"wall={wall:.1f}s", flush=True)
    print(f"  AEP: {aep_x0:.6f} -> {aep_new:.6f} MWh/yr  "
          f"({100 * (aep_new / aep_x0 - 1):+.3f} % vs x0, "
          f"{100 * (aep_new / aep_xstar - 1):+.3f} % vs x*)", flush=True)

    # -- checks ----------------------------------------------------------
    indep = independent_aep(problem, x_new, vtip_max)
    grad_check_opt = gradient_check(problem, u_new, "x_aep")
    weights = weighting_check(problem)
    post = problem.alpha_check(u_new)
    alpha_split = {"x0": alpha_by_cap(problem, u0), "x_star": alpha_by_cap(problem, u_star_existing),
                   "x_aep": alpha_by_cap(problem, u_new)}

    envelope = problem.envelope_constraint()
    g_env = envelope["fun"](u_new)
    labels = problem.envelope_row_labels()
    active_env = [labels[i] for i in np.flatnonzero(g_env < ENVELOPE_ACTIVE_TOL)]
    active_lower = [names[i] for i in range(problem.n) if abs(u_new[i]) < ACTIVE_TOL]
    active_upper = [names[i] for i in range(problem.n) if abs(1.0 - u_new[i]) < ACTIVE_TOL]
    within_box = bool(np.all(u_new >= -1e-12) and np.all(u_new <= 1.0 + 1e-12))

    # Gradient on the free variables at the optimum (KKT, ignoring the
    # envelope, which is checked separately for activity).
    g_opt = np.asarray(grad_check_opt["adjoint_mwh_per_u"])
    free = np.array([not (abs(u) < ACTIVE_TOL and g > 0) and not (abs(1 - u) < ACTIVE_TOL and g < 0)
                     for u, g in zip(u_new, g_opt)])
    projected_norm = float(np.linalg.norm(g_opt[free])) if free.any() else 0.0

    # The energies in the report's own rows sum to the total.
    energy_sum = float(sum(r["energy_mwh"] for r in indep["bins"]))

    # Against the independently computed FD-driven optimum for the same ceiling.
    audit_path = os.path.join(AUDIT_DIR, f"opt_{label}_fixed.json")
    fd_comparison = None
    if os.path.exists(audit_path):
        audit = load_json(audit_path)
        u_fd = np.array(audit["u_opt"], dtype=float)
        aep_fd_here = -problem.J(u_fd)
        fd_comparison = {
            "source": os.path.relpath(audit_path, REPO_ROOT),
            "note": "the independent run used the rating rounded to 3822.2 W and central FD "
                    "at h = 3e-6; this run uses the configured 3822.189755 W and the adjoint",
            "du_inf": float(np.max(np.abs(u_new - u_fd))),
            "aep_fd_optimum_evaluated_here_mwh": float(aep_fd_here),
            "aep_adjoint_optimum_mwh": float(aep_new),
            "delta_aep_mwh": float(aep_new - aep_fd_here),
            "delta_aep_pct": float(100.0 * (aep_new / aep_fd_here - 1.0)),
            "audit_gain_pct": float(audit["gain_pct"]),
            "audit_nit": int(audit["nit"]), "audit_nfev_total": int(audit["nfev_total"]),
            "audit_wall_time_s": float(audit["wall_time_s"]),
        }

    geometry = geometry_comparison(par, names, x_new, x0, x_star)
    geometry["max_abs_diff_u_x_aep_vs_x0"] = float(np.max(np.abs(u_new - u0)))
    geometry["max_abs_diff_u_x_aep_vs_x_star"] = float(np.max(np.abs(u_new - u_star_existing)))

    _a0, rows0, _ = evaluate_aep(x0, par, problem.resource, problem.design, problem.site, vtip_max)
    attribution = [{"v_ms": q0["v_ms"], "tsr": q0["tsr"], "Cp_x0": q0["Cp"], "Cp_x_aep": q1["Cp"],
                    "capped_x0": q0["limited"], "capped_x_aep": q1["limited"],
                    "delta_energy_mwh": q1["energy_mwh"] - q0["energy_mwh"],
                    "share_of_gain_pct": 100.0 * (q1["energy_mwh"] - q0["energy_mwh"]) / (aep_new - aep_x0)}
                   for q0, q1 in zip(rows0, indep["bins"])]

    cp = {"x0": cp_metrics(problem, x0, vtip_max),
          "x_star": cp_metrics(problem, x_star, vtip_max),
          "x_aep": cp_metrics(problem, x_new, vtip_max)}

    record = {
        "experiment": "site-specific AEP optimisation under a rotor-speed ceiling, "
                      "discrete-adjoint gradient -- EXPERIMENT, not a project result",
        "case": label,
        "vtip_max_ms": vtip_max,
        "omega_max_rpm": None if vtip_max is None else vtip_max / float(problem.design.radius_m) * 60 / (2 * math.pi),
        "v_c_ms": None if vtip_max is None else vtip_max / float(problem.design.design_tsr),
        "schedule": "lambda(V) = min(6.5, V_tip,max / V); P = min(P_aero, P_rated)",
        "provisional_bounds": PROVISIONAL_LABEL,
        "bounds": {k: float(v) for k, v in PROVISIONAL_BOUNDS.items()},
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "variables": names,
        "gradient": "discrete adjoint (CeilingProblem.jac_adjoint)",
        "envelope_margin": problem.margin,
        "slsqp_options": {"ftol": args.ftol, "maxiter": args.maxiter},
        "convergence": {
            "success": bool(result.success), "status": int(result.status),
            "message": str(result.message),
            "nit": int(result.nit), "nfev": int(result.nfev), "njev": int(result.njev),
            "objective_evaluations_in_optimisation": int(n_fun_opt),
            "adjoint_evaluations_in_optimisation": int(n_adj_opt),
            "wall_time_s": wall,
            "aep_initial_mwh": float(aep_x0), "aep_final_mwh": float(aep_new),
            "fun_star": float(result.fun),
            "gradient_star_mwh_per_u": [float(v) for v in g_opt],
            "projected_gradient_norm_mwh_per_u": projected_norm,
            "envelope_max_violation": float(max(0.0, -g_env.min())),
        },
        "aep": {
            "x0_mwh": float(aep_x0), "x_star_mwh": float(aep_xstar), "x_aep_mwh": float(aep_new),
            "gain_x_star_vs_x0_pct": float(100 * (aep_xstar / aep_x0 - 1)),
            "gain_x_aep_vs_x0_pct": float(100 * (aep_new / aep_x0 - 1)),
            "gain_x_aep_vs_x_star_pct": float(100 * (aep_new / aep_xstar - 1)),
        },
        "checks": {
            "independent_aep": {k: v for k, v in indep.items() if k != "bins"},
            "independent_aep_minus_optimiser_mwh": {
                "fresh_forward_chain": float(indep["aep_fresh_forward_chain_mwh"] - aep_new),
                "assembled_from_q": float(indep["aep_assembled_from_q_mwh"] - aep_new),
                **({"src_objective": float(indep["aep_src_objective_mwh"] - aep_new)}
                   if "aep_src_objective_mwh" in indep else {}),
            },
            "per_bin_energy_sum_minus_total_mwh": float(energy_sum - indep["aep_fresh_forward_chain_mwh"]),
            "gradient_at_x0": grad_check_x0,
            "gradient_at_x_aep": grad_check_opt,
            "within_box": within_box,
            "active_bounds": {"lower": active_lower, "upper": active_upper},
            "active_envelope_rows": active_env,
            "envelope_min_row_value_m": float(g_env.min()),
            "envelope_min_row": labels[int(np.argmin(g_env))],
            "bem_all_bins_converged_at_x_aep": bool(indep["all_bins_converged"]),
            "alpha_check_at_x_aep": post,
            "alpha_uncapped_vs_capped_bins": alpha_split,
            "weighting": weights,
            "domain_errors": problem.domain_errors,
        },
        "comparison_with_audit_fd_optimum": fd_comparison,
        "cp": cp,
        "geometry": geometry,
        "bins_at_x_aep": indep["bins"],
        "gain_attribution_by_bin": attribution,
        "u0": [float(v) for v in u0], "u_aep": [float(v) for v in u_new],
        "iterates": recorder.iterates,
    }
    with open(os.path.join(RESULTS_DIR, f"case_{label}.json"), "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    return record


# -- figures -----------------------------------------------------------------

COLOURS = {"x0": "#888888", "x_star": "#d1495b", "x_aep": "#1f5fbf"}
LABELS = {"x0": "Schmitz x0", "x_star": "existing x* (Cp @ 6.5)", "x_aep": "AEP optimum (this case)"}


def plot_geometry(records, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(records)
    fig, axes = plt.subplots(2, n, figsize=(4.2 * n, 7.2), sharex=True)
    axes = np.atleast_2d(axes).reshape(2, n)
    for col, rec in enumerate(records):
        span = rec["geometry"]["spanwise"]
        r = span["radii_m"]
        title = "no ceiling (control)" if rec["vtip_max_ms"] is None else \
            f"V_tip,max = {rec['vtip_max_ms']:g} m/s ({rec['omega_max_rpm']:.0f} rpm)"
        for row, (key, unit, scale) in enumerate([("chord", "m", 1.0), ("twist", "deg", 1.0)]):
            ax = axes[row, col]
            for name in ("x0", "x_star", "x_aep"):
                y = span[f"{key}_{name}_{'m' if key == 'chord' else 'deg'}"]
                ax.plot(r, y, "-", color=COLOURS[name], lw=2.0, label=LABELS[name])
            ax.grid(True, color="#e6e6e6", lw=0.6)
            ax.set_ylabel(f"{key} [{unit}]")
            if row == 0:
                ax.set_title(f"{title}\nAEP gain vs x0: {rec['aep']['gain_x_aep_vs_x0_pct']:+.2f} %",
                             fontsize=10)
            if row == 1:
                ax.set_xlabel("radius r [m]")
            if row == 0 and col == 0:
                ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_cp_lambda(records, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(records)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, rec in zip(axes, records):
        for name in ("x0", "x_star", "x_aep"):
            sweep = rec["cp"][name]["Cp_lambda_sweep"]
            ax.plot([s["tsr"] for s in sweep], [s["Cp"] for s in sweep], "-",
                    color=COLOURS[name], lw=2.0, label=LABELS[name])
        if rec["vtip_max_ms"] is not None:
            lam_rated = tsr_schedule(float(rec["bins_at_x_aep"][7]["v_ms"]), 6.5, rec["vtip_max_ms"])
            ax.axvspan(lam_rated, 6.5, color="#1f5fbf", alpha=0.08, lw=0)
        title = "no ceiling (control)" if rec["vtip_max_ms"] is None else \
            f"V_tip,max = {rec['vtip_max_ms']:g} m/s"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("tip-speed ratio lambda")
        ax.grid(True, color="#e6e6e6", lw=0.6)
    axes[0].set_ylabel(f"Cp at V = {CP_REPORT_SPEED_MS} m/s")
    axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# -- main --------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--cases", default="none,60,55,50")
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--ftol", type=float, default=1e-8)
    args = parser.parse_args(argv)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    x0 = load_x0()
    x_star = np.array(load_json(XSTAR_PATH)["x_star"], dtype=float)

    records = []
    for label in args.cases.split(","):
        label = label.strip()
        try:
            records.append(run_case(label, CASES[label], x0, x_star, args))
        except PolarDomainError as error:
            print(f"  case {label} aborted: PolarDomainError {error}", flush=True)
            records.append({"case": label, "vtip_max_ms": CASES[label], "aborted": str(error)})

    complete = [r for r in records if "aborted" not in r]

    # Robustness: every case's optimum evaluated under every schedule.
    blades = {"x0": x0, "x_star": x_star}
    blades.update({f"x_aep_{r['case']}": np.array(r["geometry"]["x_aep"]) for r in complete})
    par, res = BladeParameterisation(), WeibullResource.from_config()
    design, site = load_design_rotor(), load_site()
    cross = {}
    for label, vtip in CASES.items():
        row = {}
        for name, x in blades.items():
            aep, _rows, ok = evaluate_aep(x, par, res, design, site, vtip)
            row[name] = {"aep_mwh": aep, "all_converged": ok}
        for name in blades:
            row[name]["vs_x0_pct"] = 100.0 * (row[name]["aep_mwh"] / row["x0"]["aep_mwh"] - 1.0)
        cross[label] = row
    with open(os.path.join(RESULTS_DIR, "cross_evaluation.json"), "w", encoding="utf-8") as handle:
        json.dump({"description": "AEP [MWh/yr] of every blade under every schedule; rows are the "
                                  "schedule (ceiling), columns the blade", "cases": cross}, handle, indent=1)
    print("\n=== cross-evaluation (AEP, % vs x0 under the same schedule) ===")
    for label, row in cross.items():
        print(f"{label:>5}: " + "  ".join(f"{n} {v['aep_mwh']:.4f} ({v['vs_x0_pct']:+.2f}%)" for n, v in row.items()))

    summary = {
        "experiment": "site-specific AEP optimisation under a rotor-speed ceiling, "
                      "discrete-adjoint gradient -- EXPERIMENT, not a project result",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "command": "python verification/aep_optimisation_experiment/run_experiment.py",
        "provisional_bounds": PROVISIONAL_LABEL,
        "cases": {
            r["case"]: {
                "vtip_max_ms": r["vtip_max_ms"], "omega_max_rpm": r["omega_max_rpm"],
                **r["aep"],
                "nit": r["convergence"]["nit"], "nfev": r["convergence"]["nfev"],
                "njev": r["convergence"]["njev"], "status": r["convergence"]["status"],
                "message": r["convergence"]["message"], "wall_time_s": r["convergence"]["wall_time_s"],
                "Cp_at_6p5_x0": r["cp"]["x0"]["Cp_at_design_tsr"],
                "Cp_at_6p5_x_star": r["cp"]["x_star"]["Cp_at_design_tsr"],
                "Cp_at_6p5_x_aep": r["cp"]["x_aep"]["Cp_at_design_tsr"],
                "adjoint_vs_fd_at_optimum_rel_max": r["checks"]["gradient_at_x_aep"]["fd"]["3e-06"]["rel_diff_max_to_gradient_norm"],
                "active_bounds": r["checks"]["active_bounds"],
                "active_envelope_rows": r["checks"]["active_envelope_rows"],
                "all_converged": r["checks"]["bem_all_bins_converged_at_x_aep"],
                "alpha_uncapped_within_band_x_aep": r["checks"]["alpha_uncapped_vs_capped_bins"]["x_aep"]["uncapped_within_cache_band"],
                "alpha_uncapped_within_band_x0": r["checks"]["alpha_uncapped_vs_capped_bins"]["x0"]["uncapped_within_cache_band"],
                "x_aep": r["geometry"]["x_aep"],
                "du_inf_vs_audit_fd_optimum": None if r["comparison_with_audit_fd_optimum"] is None
                else r["comparison_with_audit_fd_optimum"]["du_inf"],
            } for r in complete
        },
        "aborted": [r for r in records if "aborted" in r],
    }
    with open(os.path.join(RESULTS_DIR, "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)

    if complete:
        plot_geometry(complete, os.path.join(RESULTS_DIR, "geometry.png"))
        plot_cp_lambda(complete, os.path.join(RESULTS_DIR, "cp_lambda.png"))

    print("\n=== summary ===")
    for label, c in summary["cases"].items():
        print(f"{label:>5}: x0 {c['x0_mwh']:.4f}  x* {c['x_star_mwh']:.4f}  x_AEP {c['x_aep_mwh']:.4f}  "
              f"gain vs x0 {c['gain_x_aep_vs_x0_pct']:+.3f} %  vs x* {c['gain_x_aep_vs_x_star_pct']:+.3f} %  "
              f"nit {c['nit']} status {c['status']}  adj/FD {c['adjoint_vs_fd_at_optimum_rel_max']:.1e}")
    print(f"wrote {RESULTS_DIR}")


if __name__ == "__main__":
    main()
