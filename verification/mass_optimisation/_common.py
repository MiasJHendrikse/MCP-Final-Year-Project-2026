"""
Phase 5 (2026-09-20) -- shared pieces of the mass-optimisation artefact.

The four runners in this directory (`run_mass_slsqp.py`, `run_multistart.py`,
`run_sweep.py`, `run_mass_checks.py`, `run_cross_evaluation.py`) solve or
evaluate ONE problem, stated in `docs/PLAN-mass-objective-2026-09-20.md`:

    minimise    f(u) = m_shell(d(u)) / m_shell(x0)          geometric, no BEM
    subject to  g_AEP = AEP(u)/AEP(x0) - (1 - delta)  >= 0   the energy floor
                g_M   = KS0 - KS(u)                   >= 0   the Phase 4 load cap
                g_s   = KS0/c00^2 - KS(u)/c0(u)^2     >= 0   root-stress proxy
                g_d   = D0 - D(u)                     >= 0   tip-deflection proxy
                envelope (50), solidity (25), monotone chord/twist (4 + 4)
                u in [0, 1]^10

The rows are `ScaledProblem.constraints_for_mass_problem(delta)`; the
Jacobians are the objective adjoint (AEP floor), the moment adjoint (cap and
stress) and the deflection adjoint. Everything here that the Phase 4 scripts
already had (`build_problem`, `loads_at`, the recorder, the KKT estimate) is
lifted from `verification/load_constraint/run_constrained_slsqp.py` rather
than copied a fourth time, and extended to the mass problem's rows.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import datetime
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import Bounds, minimize, nnls

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
if os.path.join(REPO_ROOT, "src") not in sys.path:
    sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from bem.rotor import solve_rotor  # noqa: E402
from design import BladeParameterisation, DesignBounds, clamped_knots  # noqa: E402
from gradients import ScaledProblem  # noqa: E402
from gradients.problem import (ACTIVE_TOL, DEFAULT_ENVELOPE_MARGIN,  # noqa: E402
                               MASS_PROBLEM_ROWS)
from objective import WeibullResource  # noqa: E402
from objective.loads import root_bending_moment, spanwise_moments, tip_deflection  # noqa: E402
from objective.power import tsr_schedule  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
XC_PATH = os.path.join(REPO_ROOT, "verification", "load_constraint", "result_eps0.json")
STARTS_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation_multistart", "starts.json")
MULTISTART_PATH = os.path.join(REPO_ROOT, "verification", "fd_optimisation_multistart",
                               "results.json")
SWEEP_PATH = os.path.join(REPO_ROOT, "verification", "fd_step_size", "sweep.json")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved 2026-09-19; "
                "chord_min_m, twist_min, twist_max grounded 2026-09-13) and the 60 mm "
                "buildable-tip floor of the mass problem (manufacturing.min_chord_m, 2026-09-20)")
LAW_LABEL = ("lambda(V) = min(6.5, Omega_max R / V), Omega_max = 300 rpm "
             "(V_c = 9.67 m/s), fixed rating 3822.189755449124 W")
PROBLEM_LABEL = ("minimise shell material k_P int c dr subject to AEP >= (1 - delta) AEP(x0), "
                 "the Phase 4 moment cap, the root-stress and tip-deflection proxies, "
                 "monotone chord and twist control points, the 60 mm min-chord floor, "
                 "the envelope and solidity")

#: The energy floors of the sweep; `0.0` is the production optimum.
DELTAS = (0.0, 0.0025, 0.005, 0.01, 0.02)
SLSQP_FTOL = 1e-10
SLSQP_MAXITER = 400

#: The sanity band stated in the plan BEFORE the run (plan, Verification 3):
#: at delta = 0 with the full row set, the shell saving lies between 2 % and
#: 10 %, the solid proxy between 5 % and 20 %; a saving beyond the scratch's
#: deflection-free 7.5 % (shell) is a suspected defect.
EXPECTED_SHELL_SAVING_PCT = (2.0, 10.0)
EXPECTED_SOLID_SAVING_PCT = (5.0, 20.0)
DEFECT_SHELL_SAVING_PCT = 7.5


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def load_x0():
    artefact = load_json(X0_PATH)
    return np.array(artefact["chord_control_points_m"]
                    + artefact["twist_control_points_rad"])


def load_xc():
    """The Phase 4 `eps = 0` optimum `x_c` -- the energy-optimal blade. Its tip
    (47.7 mm) predates the 60 mm min-chord row of 2026-09-20: inside the box,
    infeasible for that row; evaluated as it is, and a legitimate (infeasible)
    start."""

    return np.array(load_json(XC_PATH)["x_c"], dtype=float)


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


def prepared_problem(margin=DEFAULT_ENVELOPE_MARGIN):
    """
    A `ScaledProblem` with every reference measured from the committed `x0`:
    `J0`, `KS0`/`m_ref`, `material_ref`, `delta_ref`/`D0`. Returns
    `(problem, u0)`.
    """

    problem = build_problem(margin)
    u0 = problem.scaled(load_x0())
    problem.set_reference(u0)
    problem.load_system()
    problem.material_model()
    problem._deflection_system()
    return problem, u0


def variable_names(problem):
    return ([f"chord_{j}" for j in range(problem.parameterisation.n_chord)]
            + [f"twist_{j}" for j in range(problem.parameterisation.n_twist)])


def control_points_record(problem, d):
    """Chord in mm and twist in degrees, the planform table's units."""

    n_c = problem.parameterisation.n_chord
    return {"chord_mm": [float(1e3 * x) for x in d[:n_c]],
            "twist_deg": [float(math.degrees(x)) for x in d[n_c:]]}


# ---------------------------------------------------------------------------
# forward-path evaluation of one blade
# ---------------------------------------------------------------------------

def loads_at(problem, d, with_spanwise=False):
    """Forward-path loads at the design condition (11 m/s on the ceiling)."""

    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    v_inf = float(problem.design.rated_wind_speed_ms)
    tsr = tsr_schedule(v_inf, problem.design.design_tsr, problem.design.max_tip_speed_ms)
    rho_air = float(problem.site.air_density)
    result = solve_rotor(geometry, tsr=tsr, v_inf=v_inf, air_density=rho_air,
                         kinematic_viscosity=problem.site.kinematic_viscosity)
    area = math.pi * geometry.R ** 2
    moment = root_bending_moment(result["stations"], geometry.chord, rho_air,
                                 geometry.n_blades, geometry.r_hub)
    deflection = tip_deflection(result["stations"], geometry.chord, rho_air, geometry.R)
    record = {
        "v_ms": v_inf, "tsr": float(tsr),
        "rpm": float(tsr * v_inf / geometry.R * 60.0 / (2.0 * math.pi)),
        "Cp": float(result["Cp"]), "Ct": float(result["Ct"]),
        "rotor_thrust_n": float(result["Ct"] * 0.5 * rho_air * v_inf ** 2 * area),
        "root_bending_moment_nm": float(moment),
        "tip_deflection_per_unit_stiffness": float(deflection),
        "converged": bool(result["converged"]),
    }
    if with_spanwise:
        record["radii_m"] = [float(s["r"]) for s in result["stations"]]
        record["spanwise_moment_nm"] = [float(m) for m in
                                        spanwise_moments(result["stations"], geometry.chord, rho_air)]
        record["chord_m"] = [float(c) for c in geometry.chord]
    return record


def blade_record(problem, u, delta, label, u0=None):
    """
    Everything the artefacts quote for one blade: AEP and its ratio to
    `x0`, both material proxies, the load rows' values and ratios, every
    mass-problem slack at `delta`, the rated-point loads, the post-check
    and the control points. Unguarded: a blade that cannot be evaluated
    raises.
    """

    u = np.asarray(u, dtype=float)
    d = problem.physical(u)
    aep = problem.aep_mwh(u)
    aep0 = -float(problem.J0)
    material = problem.material_report(u)
    slacks = problem.mass_problem_slacks(u, delta)
    loads = loads_at(problem, d, with_spanwise=True)
    check = problem.alpha_check(u)
    return {
        "label": label,
        "u": [float(x) for x in u],
        "x": [float(x) for x in d],
        "control_points": control_points_record(problem, d),
        "aep_mwh_per_yr": float(aep),
        "aep_over_x0": float(aep / aep0),
        "aep_pct_vs_x0": float(100.0 * (aep / aep0 - 1.0)),
        "material": material,
        "shell_pct_vs_x0": float(100.0 * (material["shell_over_reference"] - 1.0)),
        "solid_pct_vs_x0": float(100.0 * (material["solid_over_reference"] - 1.0)),
        "mass_objective": float(problem.mass(u)),
        "KS": float(problem.moment_ks(u)),
        "KS_over_KS0": float(problem.moment_ks(u) / problem.KS0),
        "stress_ratio": float(problem.stress_ratio(u)),
        "deflection_ratio": float(problem.deflection_ks(u) / problem.D0),
        "root_chord_m": float(d[0]),
        "slacks": slacks,
        "loads_rated": loads,
        "post_check": check,
        "reynolds_min": float(check["reynolds_min"]),
        "du_inf_from_x0": None if u0 is None else float(np.max(np.abs(u - np.asarray(u0)))),
    }


# ---------------------------------------------------------------------------
# the recorder and the solve
# ---------------------------------------------------------------------------

class Recorder:
    """
    The Phase 4 recorder for the mass problem. The objective is geometric,
    so the record per accepted iterate is the row values SLSQP has just
    evaluated there (each row's `fun` is wrapped to remember its last
    point), never a second solve. `extra_evals` counts the times the
    callback's point was not the last evaluated one and a row had to be
    re-evaluated.
    """

    def __init__(self, problem, delta, include=MASS_PROBLEM_ROWS):
        self.problem = problem
        self.delta = float(delta)
        self.include = tuple(include)
        self.started = time.perf_counter()
        self.iterates = []
        self.extra_evals = 0
        self._last = {}
        self.rows = []
        self.labels = []
        for label, constraint in problem.mass_problem_rows(delta, include):
            self.labels.append(label)
            self.rows.append(self._remembering(label, constraint))

    def _remembering(self, label, constraint):
        fun = constraint["fun"]

        def remembered_fun(u):
            value = fun(u)
            self._last[label] = (np.array(u, dtype=float), np.array(value, dtype=float))
            return value

        return {"type": "ineq", "fun": remembered_fun, "jac": constraint["jac"]}

    def constraints(self):
        return self.rows

    def row_values(self, u):
        values = {}
        for label, row in zip(self.labels, self.rows):
            cached = self._last.get(label)
            if cached is not None and np.array_equal(cached[0], u):
                values[label] = cached[1]
            else:
                self.extra_evals += 1
                values[label] = np.array(row["fun"](u), dtype=float)
        return values

    def _record(self, u, values):
        problem = self.problem
        record = {
            "k": len(self.iterates),
            "u": [float(x) for x in u],
            "mass": float(problem.mass(u)),
            "shell_pct_vs_x0": float(100.0 * (problem.mass(u) - 1.0)),
            "row_slack_min": {label: float(np.min(v)) for label, v in values.items()},
            "wall_time_s": time.perf_counter() - self.started,
        }
        if "aep_floor" in values:
            record["aep_over_x0"] = float(values["aep_floor"][0] + 1.0 - self.delta)
        # These read the caches the rows just filled at `u`; no new solve.
        if "moment" in values or "stress" in values:
            record["KS_over_KS0"] = float(problem.moment_ks(u) / problem.KS0)
        if "stress" in values:
            record["stress_ratio"] = float(problem.stress_ratio(u))
        if "deflection" in values:
            record["deflection_ratio"] = float(problem.deflection_ks(u) / problem.D0)
        self.iterates.append(record)
        return record

    def callback(self, u):
        u = np.array(u, dtype=float)
        record = self._record(u, self.row_values(u))
        slack = record["row_slack_min"]
        print(f"  k={record['k']:3d}  mass={record['mass']:.6f}  "
              f"AEP/x0={record.get('aep_over_x0', float('nan')):.6f}  "
              f"stress={record.get('stress_ratio', float('nan')):.4f}  "
              f"defl={record.get('deflection_ratio', float('nan')):.4f}  "
              f"minslack={min(slack.values()):.2e}  t={record['wall_time_s']:.1f}s", flush=True)

    def record_start(self, u):
        u = np.array(u, dtype=float)
        return self._record(u, self.row_values(u))


def run_mass_slsqp(problem, u_start, delta, include=MASS_PROBLEM_ROWS,
                   maxiter=SLSQP_MAXITER, ftol=SLSQP_FTOL, verbose=True):
    """
    One SLSQP run of the mass problem at energy floor `delta` from
    `u_start`, with the rows in `include`. Returns `(result, recorder, wall)`.
    A `PolarDomainError` escaping from a row's `jac` (an accepted iterate
    that cannot be evaluated) propagates: the start is a failed start.
    """

    recorder = Recorder(problem, delta, include)
    start = np.asarray(u_start, dtype=float)
    failures_before = len(problem.evaluation_failures)
    recorder.record_start(start)
    if verbose:
        print(f"delta = {delta:g}; rows = {list(include)}; mass(start) = "
              f"{recorder.iterates[0]['mass']:.6f}; margin = {problem.margin}", flush=True)

    started = time.perf_counter()
    result = minimize(
        problem.mass, start, jac=problem.mass_jac, method="SLSQP",
        bounds=Bounds(np.zeros(problem.n), np.ones(problem.n)),
        constraints=recorder.constraints(),
        options=dict(ftol=ftol, maxiter=maxiter, disp=verbose),
        callback=recorder.callback,
    )
    wall = time.perf_counter() - started
    recorder.failures = problem.evaluation_failures[failures_before:]
    return result, recorder, wall


# ---------------------------------------------------------------------------
# KKT estimate at a candidate optimum
# ---------------------------------------------------------------------------

def kkt_report(problem, u, delta, include, names, tol=ACTIVE_TOL):
    """
    SciPy's SLSQP does not expose multipliers, so estimate them: the
    NON-NEGATIVE least-squares coefficients of `grad mass` on the active row
    normals (every row written `g >= 0`, bounds as `+e_j` / `-e_j`) are the
    multipliers; the residual after projection is zero at a KKT point up to
    the solver tolerance. NNLS rather than plain least squares because the
    active set can be degenerate -- with the 60 mm floor two tip control
    points sit on their min-chord rows AND their monotone row is tight, so
    `e_3`, `e_4` and `e_3 - e_4` are dependent and a plain fit can return a
    spurious negative coefficient. With NNLS the sign condition holds by construction
    and the residual alone is the KKT test; `all_multipliers_nonnegative` is
    kept for the record and is always True.

    Both the objective and the AEP row are fractions of `x0`'s, so the AEP
    row's multiplier IS the exchange rate at the margin: material fraction
    per energy fraction, i.e. % material per % energy.
    """

    u = np.asarray(u, dtype=float)
    rows, labels = [], []

    def add_vector_row(name, constraint, row_labels):
        values = constraint["fun"](u)
        jac = constraint["jac"](u)
        for k in np.flatnonzero(values < tol):
            labels.append(row_labels[k])
            rows.append(np.asarray(jac[k], dtype=float))

    if "envelope" in include:
        add_vector_row("envelope", problem.envelope_constraint(), problem.envelope_row_labels())
    if "solidity" in include:
        add_vector_row("solidity", problem.solidity_constraint(), problem.solidity_row_labels())
    if "manufacturing" in include:
        add_vector_row("manufacturing", problem.manufacturing_constraints(),
                       problem.manufacturing_row_labels())
    scalar = {
        "aep_floor": lambda: problem.aep_floor_constraint(delta),
        "moment": lambda: problem.moment_constraint(0.0),
        "stress": lambda: problem.stress_constraint(),
        "deflection": lambda: problem.deflection_constraint(),
    }
    slacks = {}
    for name, builder in scalar.items():
        if name not in include:
            continue
        constraint = builder()
        slack = float(constraint["fun"](u)[0])
        slacks[name] = slack
        if slack < tol:
            labels.append(name)
            rows.append(np.asarray(constraint["jac"](u)[0], dtype=float))

    for j in range(problem.n):
        if abs(u[j]) < tol:
            labels.append(f"lower:{names[j]}")
            row = np.zeros(problem.n)
            row[j] = 1.0
            rows.append(row)
        elif abs(1.0 - u[j]) < tol:
            labels.append(f"upper:{names[j]}")
            row = np.zeros(problem.n)
            row[j] = -1.0
            rows.append(row)

    g = np.asarray(problem.mass_jac(u), dtype=float)
    if rows:
        A = np.vstack(rows)
        multipliers, _rnorm = nnls(A.T, g)
        projection = A.T @ multipliers
    else:
        multipliers = np.zeros(0)
        projection = np.zeros(problem.n)
    residual = g - projection
    free = np.array([abs(u[j]) > tol and abs(1.0 - u[j]) > tol for j in range(problem.n)])
    by_row = {label: float(value) for label, value in zip(labels, multipliers)}
    exchange = by_row.get("aep_floor")

    return {
        "active_rows": labels,
        "n_active_rows": len(labels),
        "scalar_row_slacks": slacks,
        "multipliers": by_row,
        "multiplier_method": "nnls",
        "all_multipliers_nonnegative": bool(np.all(multipliers >= -1e-12)),
        "aep_floor_multiplier": exchange,
        "exchange_rate_pct_material_per_pct_energy": exchange,
        "gradient_norm": float(np.linalg.norm(g)),
        "projection_norm_onto_active_normals": float(np.linalg.norm(projection)),
        "residual_norm_after_projection": float(np.linalg.norm(residual)),
        "free_variables": [names[j] for j in range(problem.n) if free[j]],
        "gradient_norm_on_free_variables": float(np.linalg.norm(g[free])) if free.any() else 0.0,
        "n_rows_vs_n_free": [len(labels), int(free.sum())],
    }


# ---------------------------------------------------------------------------
# summary of one solve, and the sanity band
# ---------------------------------------------------------------------------

def sanity(delta, include, record, kkt, result, failures_at_optimum):
    """The plan's stated band, applied to a `delta = 0` full-set optimum."""

    notes = []
    full = tuple(include) == tuple(MASS_PROBLEM_ROWS)
    shell_saving = -record["shell_pct_vs_x0"]
    solid_saving = -record["solid_pct_vs_x0"]
    if delta == 0.0 and full:
        lo, hi = EXPECTED_SHELL_SAVING_PCT
        if not lo <= shell_saving <= hi:
            notes.append(f"shell saving {shell_saving:.2f} % outside the stated band "
                         f"{EXPECTED_SHELL_SAVING_PCT} -- a finding, not a defect gate")
        lo, hi = EXPECTED_SOLID_SAVING_PCT
        if not lo <= solid_saving <= hi:
            notes.append(f"solid saving {solid_saving:.2f} % outside the stated band "
                         f"{EXPECTED_SOLID_SAVING_PCT}")
        if shell_saving > DEFECT_SHELL_SAVING_PCT:
            notes.append(f"shell saving {shell_saving:.2f} % exceeds the scratch's "
                         f"deflection-free {DEFECT_SHELL_SAVING_PCT} %: suspected defect")
        if not record["slacks"]["aep_floor"]["active"]:
            notes.append("AEP floor not active at the optimum")
        if record["slacks"]["moment"]["active"]:
            notes.append("moment cap active at the optimum (the plan expected it inactive)")
    if kkt["aep_floor_multiplier"] is not None and kkt["aep_floor_multiplier"] < -1e-12:
        notes.append("AEP-floor multiplier negative")
    if not kkt["all_multipliers_nonnegative"]:
        notes.append("a KKT multiplier is negative")
    if not result.success:
        notes.append(f"SLSQP exit {result.status}: {result.message}")
    if failures_at_optimum:
        notes.append("a row could not be evaluated at the accepted optimum")
    if record["post_check"]["uncapped"] is not None and not record["post_check"]["uncapped"]["within"]:
        notes.append("alpha outside the XFOIL band on an uncapped point")
    return notes


def summarise_run(problem, u0, delta, include, result, recorder, wall, start_label,
                  options, names):
    u_m = np.array(result.x, dtype=float)
    record = blade_record(problem, u_m, delta, "x_m", u0=u0)
    kkt = kkt_report(problem, u_m, delta, include, names)
    failures_at_optimum = [f for f in recorder.failures
                           if np.array_equal(np.array(f["u"]), u_m)]
    notes = sanity(delta, include, record, kkt, result, failures_at_optimum)
    return {
        "problem": PROBLEM_LABEL,
        "bounds_label": BOUNDS_LABEL,
        "law_label": LAW_LABEL,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": src_commit(),
        "bounds": problem.bounds.as_record(),
        "delta": float(delta),
        "rows": list(include),
        "start": start_label,
        "slsqp_options": dict(options),
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nit": int(result.nit), "nfev": int(result.nfev), "njev": int(result.njev),
        "wall_time_s": float(wall),
        "callback_extra_evaluations": recorder.extra_evals,
        "references": {
            "J0_mwh_per_yr": float(problem.J0),
            "aep_x0_mwh_per_yr": float(-problem.J0),
            "KS0": float(problem.KS0), "m_ref_nm": float(problem.m_ref),
            "material_ref_shell_m2": float(problem.material_ref),
            "delta_ref": float(problem.delta_ref), "D0": float(problem.D0),
            "rho": float(problem.load_system().rho),
        },
        "optimum": record,
        "u_m": record["u"], "x_m": record["x"],
        "shell_pct_vs_x0": record["shell_pct_vs_x0"],
        "solid_pct_vs_x0": record["solid_pct_vs_x0"],
        "aep_pct_vs_x0": record["aep_pct_vs_x0"],
        "kkt": kkt,
        "active_set": problem.active_set(u_m, mass_delta=delta),
        "evaluation_failures": recorder.failures,
        "n_evaluation_failures": len(recorder.failures),
        "failures_at_optimum": failures_at_optimum,
        "counters": {
            "objective_evaluations_total": int(problem.n_fun_evals),
            "adjoint_evaluations_total": int(problem.n_adjoint_evals),
            "moment_solves_total": int(problem.n_moment_solves),
            "deflection_solves_total": int(problem.n_deflection_solves),
        },
        "expected_band": {"shell_saving_pct": list(EXPECTED_SHELL_SAVING_PCT),
                          "solid_saving_pct": list(EXPECTED_SOLID_SAVING_PCT),
                          "defect_above_shell_saving_pct": DEFECT_SHELL_SAVING_PCT},
        "sanity_notes": notes,
        "suspected_defect": any("defect" in n for n in notes),
    }


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def greville_abscissae(parameterisation, n_pts):
    radii = parameterisation.radii
    knots = clamped_knots(n_pts, parameterisation.degree)
    greville = np.array([knots[j + 1:j + parameterisation.degree + 1].mean()
                         for j in range(n_pts)])
    return radii[0] + greville * (radii[-1] - radii[0])


def plot_blade(problem, blades, delta, path, title_extra=""):
    """
    Planform and twist of `x0`, `x_c` and `x_m` with their control points,
    and the rated-point spanwise moment `M(r)` and the flexibility `1/c^3`
    that the deflection row weighs it with. `blades` is a list of
    `(label, style, blade_record)`.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = problem.parameterisation
    radii = p.radii
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4))
    axes = axes.ravel()

    for ax, (field, unit, get, block) in zip(axes[:2], [
        ("chord", "m", p.chord, slice(0, p.n_chord)),
        ("twist", "deg", lambda d: np.degrees(p.twist(d)), slice(p.n_chord, None)),
    ]):
        n_pts = p.n_chord if field == "chord" else p.n_twist
        s = greville_abscissae(p, n_pts)
        scale = np.degrees(1.0) if field == "twist" else 1.0
        for label, style, record in blades:
            d = np.array(record["x"], dtype=float)
            ax.plot(radii, get(d), style["ls"], color=style["color"], lw=style["lw"], label=label)
            ax.plot(s, d[block] * scale, "o", color=style["color"], ms=4.5, mfc="white")
        ax.set_xlabel("radius r [m]")
        ax.set_ylabel(f"{field} [{unit}]")
        ax.set_title(f"{field} distribution (markers: control points)", fontsize=11)
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    for label, style, record in blades:
        loads = record["loads_rated"]
        r = np.array(loads["radii_m"])
        axes[2].plot(r, loads["spanwise_moment_nm"], style["ls"], color=style["color"],
                     lw=style["lw"], label=label)
        axes[3].plot(r, 1.0 / np.array(loads["chord_m"]) ** 3, style["ls"],
                     color=style["color"], lw=style["lw"], label=label)
    axes[2].set_xlabel("radius r [m]")
    axes[2].set_ylabel("flapwise moment M(r) [N m], one blade")
    axes[2].set_title("spanwise moment at the rated point (11 m/s, 300 rpm)", fontsize=11)
    axes[3].set_xlabel("radius r [m]")
    axes[3].set_ylabel("1 / c(r)^3 [m^-3]")
    axes[3].set_title("flexibility weight of the deflection row", fontsize=11)
    axes[3].set_yscale("log")
    for ax in axes[2:]:
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)

    fig.suptitle(f"Mass problem, delta = {delta:g}{title_extra} -- " + BOUNDS_LABEL, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


STYLE_X0 = {"ls": "-", "color": "#888888", "lw": 1.6}
STYLE_XC = {"ls": "--", "color": "#d1495b", "lw": 1.6}
STYLE_XM = {"ls": "-", "color": "#1f5fbf", "lw": 2.4}


# ---------------------------------------------------------------------------
# rendered blades (2026-09-20, MJ): the surface the BEM sees, plus a
# cylindrical root and a rounded tip that are DRAWN ONLY
# ---------------------------------------------------------------------------

AIRFOIL_PATH = os.path.join(REPO_ROOT, "data", "airfoils", "sg6043.dat")
PITCH_AXIS_FRACTION = 0.30
#: The cosmetic root: a cylinder from the flange radius, blending into the
#: SG6043 section by r_hub (0.30 m, the 15 % cut-out). Outside the station
#: grid, outside the material proxy, outside every constraint.
RENDER_ROOT = {"r_flange_m": 0.10, "r_cylinder_end_m": 0.16, "cylinder_diameter_m": 0.09,
               "tip_rounding_m": 0.06}


def _smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def blade_surface(problem, d, n_span=90):
    """
    `(X, Y, Z)` arrays (span x section points) of the blade `d`: the SG6043
    section scaled by the station chord and rotated by the station twist
    about the pitch axis (30 % chord), with the cosmetic root cylinder and
    tip rounding of `RENDER_ROOT` applied. Metres.
    """

    p = problem.parameterisation
    coords = np.loadtxt(AIRFOIL_PATH, skiprows=1)
    airfoil = np.column_stack([coords[:, 0] - PITCH_AXIS_FRACTION, coords[:, 1]])
    theta = np.linspace(0.0, 2.0 * np.pi, len(coords))
    circle = np.column_stack([0.5 * np.cos(theta), 0.5 * np.sin(theta)])
    r_hub = float(p.radius_m * problem.design.root_fraction)
    root = RENDER_ROOT
    r = np.concatenate([np.linspace(root["r_flange_m"], r_hub, 25, endpoint=False),
                        np.linspace(r_hub, p.radii[-1], n_span)])
    chord_st, twist_st = p.chord(d), p.twist(d)
    chord = np.interp(r, p.radii, chord_st, left=chord_st[0])
    twist = np.interp(r, p.radii, twist_st, left=twist_st[0])
    X, Y, Z = [], [], []
    tip_start = p.radii[-1] - root["tip_rounding_m"]
    for ri, ci, ti in zip(r, chord, twist):
        s = _smoothstep((ri - root["r_cylinder_end_m"]) / (r_hub - root["r_cylinder_end_m"]))
        section = (1.0 - s) * circle * root["cylinder_diameter_m"] + s * airfoil * ci
        if ri > tip_start:
            section = section * (0.15 + 0.85 * math.sqrt(max(0.0, 1.0 - ((ri - tip_start) / root["tip_rounding_m"]) ** 2)))
        tw = s * ti
        X.append(np.full(len(section), ri))
        Y.append(section[:, 0] * np.cos(tw) - section[:, 1] * np.sin(tw))
        Z.append(section[:, 0] * np.sin(tw) + section[:, 1] * np.cos(tw))
    return np.array(X), np.array(Y), np.array(Z)


def render_blades(problem, blades, path):
    """
    Three views (isometric, edge-on, top) of each blade in `blades`
    (`(label, style, blade_record)`), one row per blade, with the cosmetic
    root and tip rounding. The caption says they are cosmetic.
    """

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    R = float(problem.parameterisation.radius_m)
    fig = plt.figure(figsize=(16, 2.7 * len(blades)))
    for k, (label, style, record) in enumerate(blades):
        X, Y, Z = blade_surface(problem, np.array(record["x"], dtype=float))
        for j, (elev, azim) in enumerate([(25, -60), (0, -90), (90, -90)]):
            ax = fig.add_subplot(len(blades), 3, 3 * k + j + 1, projection="3d")
            ax.plot_surface(X, Y, Z, color=style["color"], alpha=0.75, lw=0, rstride=1, cstride=1, shade=True)
            for i in range(0, X.shape[0], 5):
                ax.plot(X[i], Y[i], Z[i], color="k", lw=0.25, alpha=0.5)
            ax.set_box_aspect((2.0, 0.55, 0.55))
            ax.view_init(elev=elev, azim=azim)
            ax.set_xlim(0, R); ax.set_ylim(-0.12, 0.25); ax.set_zlim(-0.12, 0.12)
            ax.set_axis_off()
            if j == 0:
                ax.set_title(label, fontsize=10, loc="left")
    fig.suptitle("Rendered blades -- isometric / edge-on (twist) / top view. The root cylinder "
                 f"(r < {RENDER_ROOT['r_cylinder_end_m']:.2f} m), its transition to the section by "
                 f"r_hub, and the tip rounding are DRAWN ONLY: outside the BEM grid, the material "
                 "proxy and every constraint.", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
