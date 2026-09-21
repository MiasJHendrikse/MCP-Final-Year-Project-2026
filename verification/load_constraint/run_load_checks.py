"""
Step 2d -- the root-moment KS functional and its constraint: the
four-tier checks at `x0` and at the unconstrained optimum `u*`.

At each point:

  * Tier 1 -- every partial of the kernel's `m` (`dm_dphi`, `dm_dd`) and of the
    KS functional (`dKS_dx`, `dKS_dd`) against a complex step of the code's own
    *value* path (`station_partials(..., derivatives=False)`, `KS`), reported as
    the worst mixed error `|estimate - partial| / max(1, |partial|)`.
  * Tier 2 -- the forward-mode tangent against the adjoint direction, and the
    adjoint assembly `dKS_dd == explicit + (dR/dd)^T psi` by construction.
  * Tier 3 -- `moment_constraint(0)["jac"](u)` against central FD of its `fun`
    at the committed global `h* = 3.162277660168379e-06` and at `h* / sqrt(10)`,
    `h* sqrt(10)`, with `eps_j` re-measured as the three-step local jitter (as
    `verification/gradient_verification/run_tier3.py` does) and the measured
    round-off floor `delta_g / h*` at `h*/1000`. Acceptance `|adj - fd| <= 3 eps_j`.

Also recorded, per point: the per-point root moments and softmax weights over
the nine-point load set `L`; `KS` and its conservatism `KS - max` at
`rho in {30, 100, 300}`; the design-condition `Ct` and thrust; the
**B3-dependent** cut-out post-check at 20 m/s; and the wall time of one
objective evaluation, one moment evaluation and one moment-constraint gradient.

Outputs, next to this script: `checks.json`, `moments_x0_xstar.png`, and a
hand-written `README.md`.

Run from the repo root (about two minutes):

    python verification/load_constraint/run_load_checks.py

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

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from design import BladeParameterisation, DesignBounds  # noqa: E402
from gradients import ScaledProblem, central_difference  # noqa: E402
from objective import WeibullResource  # noqa: E402
from objective.loads import load_operating_points  # noqa: E402

X0_PATH = os.path.join(REPO_ROOT, "verification", "baseline", "x0.json")
SWEEP_PATH = os.path.join(REPO_ROOT, "verification", "fd_step_size", "sweep.json")
ADJOINT_RESULT_PATH = os.path.join(REPO_ROOT, "verification", "adjoint_optimisation", "result.json")
CHECKS_PATH = os.path.join(_HERE, "checks.json")
FIGURE_PATH = os.path.join(_HERE, "moments_x0_xstar.png")

BOUNDS_LABEL = ("under the configured bounds (chord_max_m = 0.30 m, resolved "
                "2026-09-19; chord_min_m, twist_min, twist_max grounded 2026-09-13)")
LAW_LABEL = ("lambda(V) = min(6.5, Omega_max R / V), Omega_max = 300 rpm "
             "(V_c = 9.67 m/s), fixed rating 3822.189755449124 W")

H = 1e-30
SQRT10 = math.sqrt(10.0)
EPS_FACTOR = 3.0
TAYLOR_STEPS = (1e-2, 1e-3, 1e-4)
TAYLOR_SEEDS = (0, 1, 2)


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


def build_problem():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config())


def variable_names(parameterisation):
    return ([f"chord_{j}" for j in range(parameterisation.n_chord)]
            + [f"twist_{j}" for j in range(parameterisation.n_twist)])


def mixed_error(estimate, partial):
    """Worst `|estimate - partial| / max(1, |partial|)` (the mixed tolerance)."""

    err = np.abs(np.asarray(estimate) - np.asarray(partial))
    denominator = np.maximum(1.0, np.abs(partial))
    return float(np.max(err / denominator))


# ---------------------------------------------------------------------------
# tiers
# ---------------------------------------------------------------------------

def tier1_m_partials(system, state, parts):
    """Worst mixed error of `dm_dphi` and `dm_dd` against a complex step of `m`."""

    dm_dphi = parts["dm_dphi"]
    estimate = np.empty_like(dm_dphi)
    for b in range(system.n_points):
        for i in range(system.n_stations):
            phi = state.phi.astype(complex)
            phi[b, i] += 1j * H
            estimate[b, i] = system.partials(phi, state.d,
                                             derivatives=False)["m"].imag[b, i] / H

    dm_dd = (parts["dm_dc"][:, :, None] * system.N_c[None, :, :]
             + parts["dm_dtheta"][:, :, None] * system.N_theta[None, :, :])
    dd_estimate = np.empty(dm_dd.shape)
    for j in range(system.n_design):
        d = state.d.astype(complex)
        d[j] += 1j * H
        dd_estimate[:, :, j] = system.partials(state.phi, d,
                                               derivatives=False)["m"].imag / H

    return {"dm_dphi": mixed_error(estimate, dm_dphi),
            "dm_dd": mixed_error(dd_estimate, dm_dd)}


def tier1_ks_partials(system, state, parts):
    """Worst mixed error of `dKS_dx` and `dKS_dd` against a complex step of `KS`."""

    dKS_dx = system.dKS_dx(state.phi, state.d, parts)
    estimate = np.empty_like(dKS_dx)
    for b in range(system.n_points):
        for i in range(system.n_stations):
            phi = state.phi.astype(complex)
            phi[b, i] += 1j * H
            estimate[b, i] = system.KS(phi, state.d).imag / H

    dKS_dd = system.dKS_dd(state.phi, state.d, parts)
    dd_estimate = np.empty_like(dKS_dd)
    for j in range(system.n_design):
        d = state.d.astype(complex)
        d[j] += 1j * H
        dd_estimate[j] = system.KS(state.phi, d).imag / H

    return {"dKS_dx": mixed_error(estimate, dKS_dx),
            "dKS_dd": mixed_error(dd_estimate, dKS_dd)}


def tier2(system, state, parts):
    """Tangent-vs-adjoint directions and the adjoint assembly identity."""

    result = system.gradient(state.d, state=state)
    explicit = result.dKS_dd_explicit
    assembled = explicit + system.apply_dR_dd_T(parts, result.psi)
    identity = float(np.max(np.abs(assembled - result.dKS_dd))
                     / max(1.0, np.max(np.abs(result.dKS_dd))))

    rng = np.random.default_rng(42)
    worst = 0.0
    for _ in range(3):
        v = rng.normal(size=system.n_design)
        tangent = system.tangent(state.phi, state.d, v)
        adjoint = float(result.dKS_dd @ v)
        worst = max(worst, abs(tangent - adjoint) / max(1.0, abs(adjoint)))
    return {"tangent_vs_adjoint": float(worst), "assembly_identity": identity}


def roundoff_floor(g, u, h_star):
    """
    The measured round-off floor of the constraint gradient, Tier 4's way:
    sample `g` along `u + t e_j` for `t = -4e-12 .. 4e-12` (nine points), fit a
    line, take the residual std `delta_g`, and report `delta_g / h_star`. This
    is what `eps_j` is *supposed* to estimate; where it lands below it the
    three-step jitter has hit a degenerate pair (the section 2 fragility).
    """

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


def tier3(problem, u, h_star):
    """The moment constraint Jacobian against central FD of its value."""

    constraint = problem.moment_constraint(0.0)
    jac = constraint["jac"](u)[0]
    g = lambda w: float(constraint["fun"](w)[0])  # noqa: E731

    def fd(h):
        return central_difference(g, u, h)[0]

    g_mid = fd(h_star)
    g_lo = fd(h_star / SQRT10)
    g_hi = fd(h_star * SQRT10)
    eps = np.maximum(np.abs(g_mid - g_lo), np.abs(g_mid - g_hi))
    abs_err = np.abs(jac - g_mid)

    small = h_star / 1000.0
    e = np.eye(problem.n)
    small_step = np.array([abs(g(u + small * e[j]) - g(u - small * e[j])) / h_star
                           for j in range(problem.n)])
    fit_floor = roundoff_floor(g, u, h_star)

    ratio = abs_err / eps
    within = abs_err <= EPS_FACTOR * eps
    eps_below_floor = eps < fit_floor
    # A failure is "floor-limited" when every failing variable's eps_j sits
    # below the measured round-off floor (the section 2 fragility of the
    # acceptance scale). That describes the failure; it does not pass it:
    # `passes` stays False and the reader decides.
    floor_limited = bool((not np.all(within)) and np.all(within | eps_below_floor))
    return {
        "h_star": float(h_star),
        "adjoint": [float(x) for x in jac],
        "fd_at_h_star": [float(x) for x in g_mid],
        "fd_local_sweep": {
            "h_star_over_sqrt10": [float(x) for x in g_lo],
            "h_star_times_sqrt10": [float(x) for x in g_hi],
        },
        "eps_moment_per_u": [float(x) for x in eps],
        "abs_error": [float(x) for x in abs_err],
        "abs_error_over_eps": [float(x) for x in ratio],
        "measured_floor_tiny_step": [float(x) for x in small_step],
        "measured_roundoff_floor": [float(x) for x in fit_floor],
        "eps_below_roundoff_floor": [bool(x) for x in eps_below_floor],
        "within_3eps": [bool(x) for x in within],
        "worst_variable": int(np.argmax(ratio)),
        "worst_abs_error_over_eps": float(np.max(ratio)),
        "passes": bool(np.all(within)),
        "floor_limited": floor_limited,
    }


def taylor(problem, u, seed):
    constraint = problem.moment_constraint(0.0)
    jac = constraint["jac"](u)[0]
    g = lambda w: float(constraint["fun"](w)[0])  # noqa: E731
    envelope = problem.envelope_constraint()["fun"]

    rng = np.random.default_rng(seed)
    for _attempt in range(50):
        v = rng.normal(size=problem.n)
        v /= np.linalg.norm(v)
        if all(np.all(envelope(u + h * v) > 0.0) for h in TAYLOR_STEPS):
            break
    else:
        raise RuntimeError("no envelope-feasible direction found for the moment Taylor test")

    g0 = g(u)
    slope = float(jac @ v)
    remainders = [abs(g(u + h * v) - g0 - h * slope) for h in TAYLOR_STEPS]
    ratios = [remainders[k] / remainders[k + 1] for k in range(len(TAYLOR_STEPS) - 1)]
    return {"seed": seed, "v": [float(x) for x in v],
            "remainder": [float(r) for r in remainders],
            "ratios_per_decade": [float(r) for r in ratios],
            "min_ratio": float(min(ratios)), "passes": bool(min(ratios) >= 30.0)}


def verify_point(problem, label, u, h_star, names, time_it=False):
    u = np.asarray(u, dtype=float)
    system = problem.load_system()
    state = system.solve(problem.physical(u))
    parts = system.partials(state.phi, state.d)

    print(f"\n== {label} ==", flush=True)
    point = {
        "label": label,
        "u": [float(x) for x in u],
        "x": [float(x) for x in problem.physical(u)],
        "J_mwh_per_yr": float(problem.J(u)),
        "KS": problem.moment_ks(u),
        "tier1_m": tier1_m_partials(system, state, parts),
        "tier1_ks": tier1_ks_partials(system, state, parts),
        "tier2": tier2(system, state, parts),
        "tier3": tier3(problem, u, h_star),
        "taylor": [taylor(problem, u, seed) for seed in TAYLOR_SEEDS],
        "moment_report": problem.moment_report(u),
    }
    print(f"  Tier 1 worst mixed: dm_dphi {point['tier1_m']['dm_dphi']:.3e}, "
          f"dm_dd {point['tier1_m']['dm_dd']:.3e}, "
          f"dKS_dx {point['tier1_ks']['dKS_dx']:.3e}, "
          f"dKS_dd {point['tier1_ks']['dKS_dd']:.3e}", flush=True)
    print(f"  Tier 2 tangent {point['tier2']['tangent_vs_adjoint']:.3e}, "
          f"assembly {point['tier2']['assembly_identity']:.3e}", flush=True)
    t3 = point["tier3"]
    print(f"  Tier 3 worst |diff|/eps {t3['worst_abs_error_over_eps']:.3f} "
          f"at {names[t3['worst_variable']]} (passes {t3['passes']}, "
          f"floor-limited {t3['floor_limited']})", flush=True)

    if time_it:
        # Fresh probes, clear of the reference cached at construction, so the
        # timing is one solve and not a cache hit.
        probe = np.clip(u + 1e-4, 0.0, 1.0)
        probe_jac = np.clip(u - 1e-4, 0.0, 1.0)
        started = time.perf_counter()
        problem.J(probe)
        objective_s = time.perf_counter() - started
        started = time.perf_counter()
        problem.moment_state(probe)
        moment_s = time.perf_counter() - started
        started = time.perf_counter()
        problem.moment_constraint(0.0)["jac"](probe_jac)
        moment_jac_s = time.perf_counter() - started
        point["timings_s"] = {"objective": objective_s, "moment_evaluation": moment_s,
                              "moment_constraint_gradient": moment_jac_s}
        print(f"  timings: J {objective_s:.3f}s, moment {moment_s:.3f}s, "
              f"moment jac {moment_jac_s:.3f}s", flush=True)

    return point


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------

def plot(points, names, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    colours = {"x0": "#888888", "unconstrained optimum u*": "#1f5fbf"}

    for k, point in enumerate(points):
        report = point["moment_report"]
        speeds = np.array([p["v_ms"] for p in report["points"]])
        moments = [p["moment_nm"] for p in report["points"]]
        weights = [p["softmax_weight"] for p in report["points"]]
        colour = colours.get(point["label"], None)
        axes[0].plot(speeds, moments, "o-", color=colour, lw=1.8, ms=5,
                     label=point["label"])
        axes[1].bar(speeds + (k - 0.5) * 0.38, weights, width=0.38,
                    color=colour, alpha=0.7, label=point["label"])

    axes[0].set_xlabel("V [m/s]")
    axes[0].set_ylabel("root bending moment M_b [N m]")
    axes[0].set_title("per-blade root moment over the load set L", fontsize=11)
    axes[1].set_xlabel("V [m/s]")
    axes[1].set_ylabel("KS softmax weight")
    axes[1].set_title("softmax weight on each load point (rho = 100)", fontsize=11)
    for ax in axes:
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    args = parser.parse_args(argv)

    x0 = load_x0()
    sweep = load_json(SWEEP_PATH)
    h_star = float(sweep["h_star_global"])
    adjoint_result = load_json(ADJOINT_RESULT_PATH)
    u_star = np.array(adjoint_result["u_star"], dtype=float)

    problem = build_problem()
    names = variable_names(problem.parameterisation)
    u0 = problem.scaled(x0)
    J0 = problem.set_reference(u0)
    problem.load_system()
    print(f"J(x0) = {J0:.6f} MWh/yr; M_ref = {problem.m_ref!r}; "
          f"KS0 = {problem.KS0!r}; rho = {problem.load_system().rho}", flush=True)
    print(f"load set L = {load_operating_points()}", flush=True)

    started = time.perf_counter()
    points = [
        verify_point(problem, "x0", u0, h_star, names, time_it=True),
        verify_point(problem, "unconstrained optimum u*", u_star, h_star, names),
    ]
    wall = time.perf_counter() - started

    plot(points, names, FIGURE_PATH)

    summary = {
        "description": ("Root-moment KS functional and constraint: Tiers 1-3, "
                        "the per-point moments and weights, and the B3-dependent cut-out "
                        "post-check at x0 and the unconstrained optimum u*."),
        "bounds_label": BOUNDS_LABEL,
        "law_label": LAW_LABEL,
        "command": "python verification/load_constraint/run_load_checks.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": src_commit(),
        "bounds": problem.bounds.as_record(),
        "rho": float(problem.load_system().rho),
        "m_ref_nm": float(problem.m_ref),
        "KS0": float(problem.KS0),
        "load_operating_points": [[float(v), float(lam)] for v, lam in load_operating_points()],
        "variables": names,
        "h_star": h_star,
        "acceptance": {"tier1": "complex step, mixed tolerance",
                       "tier2": "tangent vs adjoint, assembly identity",
                       "tier3": f"|adj - fd| <= {EPS_FACTOR:g} eps_j"},
        "points": points,
        "all_tier3_pass": bool(all(p["tier3"]["passes"] for p in points)),
        "tier3_failures": [{"point": p["label"],
                            "variable": names[p["tier3"]["worst_variable"]],
                            "floor_limited": p["tier3"]["floor_limited"]}
                           for p in points if not p["tier3"]["passes"]],
        "all_taylor_pass": bool(all(t["passes"] for p in points for t in p["taylor"])),
        "wall_time_s": wall,
        "n_moment_solves": int(problem.n_moment_solves),
        "objective_evaluations_total": int(problem.n_fun_evals),
        "adjoint_evaluations_total": int(problem.n_adjoint_evals),
    }
    with open(CHECKS_PATH, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)

    print(f"\nTier 3 all pass: {summary['all_tier3_pass']}; "
          f"failures: {summary['tier3_failures']}; "
          f"Taylor all pass: {summary['all_taylor_pass']}  ({wall:.1f} s)")
    print(f"wrote {CHECKS_PATH}\n      {FIGURE_PATH}")
    return summary


if __name__ == "__main__":
    main()
