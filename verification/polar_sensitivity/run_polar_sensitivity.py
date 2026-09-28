"""
Does the minimum-material blade's energy parity survive the polar uncertainty?

The mass optimum `x_m` holds AEP(x0) to ~1e-12 on the production SG6043 polar
cache (XFOIL, n_crit = 9), whose lift differs from the UIUC measurements by
0.11-0.15 RMS at the design Reynolds numbers and which has no measurement at
all below Re = 100 000. The report argues that a polar error cancels in the
energy ratio. That holds for an error the two blades see alike; the saving,
though, moves chord from the tip to mid-span, across different Reynolds
numbers, so an error that varies with Reynolds number need not cancel.

This script re-evaluates the three FIXED blades (`x0`, `x_c`, `x_m`; nothing
is re-optimised) under other polar surfaces and reports AEP(x_m) / AEP(x0)
and AEP(x_c) / AEP(x0) under each:

  * (**XFOIL at n_crit = 7 and 11 -- attempted and abandoned.** Both
    caches were built on 2026-09-26, but XFOIL's gap-closure retries failed
    to converge at Re = 100 000, so the run was stopped and reverted; no
    partial cache is committed. The perturbation cases
    below bound the lift error directly instead.)
  * **uniform lift offsets, +/- s(Re)** -- `s(Re)` is the measured lift RMS
    difference against UIUC at each validated Reynolds number
    (results/sg6043_uiuc_validation/), interpolated in log Re and held at
    its end values outside 100 k - 500 k. An error of this kind should
    largely cancel.
  * **tilted lift offsets, +/- s(Re) w(Re)** -- `w` falls linearly in log Re
    from +1 at Re = 62 000 to -1 at Re = 339 000 (the Reynolds range the
    three blades occupy), so low-Re stations gain lift and high-Re stations
    lose it, or the reverse. This is the error shape that does NOT cancel,
    sized by the measured error; it is an adversarial test, not a model of
    the true error.
  * **drag scaled by 0.8 and 1.2** -- an assumed +/- 20 % on Cd, not derived
    from the UIUC comparison (whose drag RMS includes post-stall points).

The offsets are applied to the production cache's grid over the
XFOIL-converged band (-8 .. 18 deg) only, then a fresh C1 interpolant is
built over the perturbed grid.

Re-evaluation answers "does x_m still match x0?"; it cannot say how much of the
saving survives, because a blade that misses the floor can buy the energy back
with material. So each surface is also **re-optimised**: the delta = 0 mass
problem solved afresh from x0 with every reference re-measured on that
surface. The re-optimised shell saving is the headline under that surface. The rating stays at the configured 3822 W
(it stands for a nameplate, which a polar error does not change), so the
capped bins are whatever each blade reaches under each surface.

Outputs, next to this script: `polar_sensitivity.json`.

    python verification/polar_sensitivity/run_polar_sensitivity.py    # ~3 min

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import datetime
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
MASS_DIR = os.path.join(REPO_ROOT, "verification", "mass_optimisation")
for path in (os.path.join(REPO_ROOT, "src"), MASS_DIR):
    if path not in sys.path:
        sys.path.insert(0, path)

import _common as C  # noqa: E402
import adjoint.system as adjoint_system_module  # noqa: E402
import gradients.problem as problem_module  # noqa: E402
import polars.polar as polar_module  # noqa: E402
from config import load_polar_cache  # noqa: E402
from objective.objective import HOURS_PER_YEAR  # noqa: E402
from objective.power import aerodynamic_power, wind_speed_bins  # noqa: E402
from polars.cache import PolarGrid  # noqa: E402
from polars.interpolant import PolarInterpolant  # noqa: E402

RESULT_PATH = os.path.join(MASS_DIR, "result_delta0.json")
UIUC_PATH = os.path.join(REPO_ROOT, "results", "sg6043_uiuc_validation",
                         "sg6043_uiuc_comparison.json")
OUT_PATH = os.path.join(_HERE, "polar_sensitivity.json")

TILT_RE = (62_000.0, 339_000.0)
DRAG_SCALES = (0.8, 1.2)
#: The production case is the committed optimum; re-solving it is a check that
#: the re-optimisation path reproduces -3.379 %, so it is run as well.
REOPTIMISE_PRODUCTION = True

#: Every module that looks the polar surface up by name. The forward solver
#: goes through `polars.polar`; the adjoint system and the problem imported
#: `interpolant_for` into their own namespaces, so each is redirected too, or
#: a re-optimisation would take forward values from one surface and
#: gradients from another.
_LOOKUP_MODULES = (polar_module, adjoint_system_module, problem_module)
_ORIGINAL_INTERPOLANT_FOR = polar_module.interpolant_for


def cache_dir(name):
    return os.path.join(REPO_ROOT, load_polar_cache(name).directory)


def lift_error_scale():
    """`s(Re)`: the measured lift RMS difference, as a function of Re."""

    rows = sorted(json.load(open(UIUC_PATH, encoding="utf-8")), key=lambda r: r["re"])
    log_re = np.log([r["re"] for r in rows])
    rmse = np.array([r["cl_rmse"] for r in rows])
    return (lambda re: np.interp(np.log(re), log_re, rmse)), {int(r["re"]): r["cl_rmse"] for r in rows}


def tilt(re):
    lo, hi = np.log(TILT_RE[0]), np.log(TILT_RE[1])
    return np.clip(1.0 - 2.0 * (np.log(re) - lo) / (hi - lo), -1.0, 1.0)


def perturbed_interpolant(base_dir, dcl=None, cd_scale=1.0):
    """C1 interpolant over the production grid with a lift offset / drag scale in the XFOIL band."""

    grid = PolarGrid(base_dir)
    band = (grid.alpha_values >= grid.xfoil_alpha_min) & (grid.alpha_values <= grid.xfoil_alpha_max)
    cl = np.array(grid.cl_grid, dtype=float)
    cd = np.array(grid.cd_grid, dtype=float)
    for m, re in enumerate(grid.re_values):
        if dcl is not None:
            cl[m, band] += dcl(re)
        cd[m, band] *= cd_scale
    grid.cl_grid, grid.cd_grid = cl, cd
    return PolarInterpolant(grid)


def use_surface(interpolant):
    """Route every `sg6043` polar lookup to `interpolant` (None restores the cache)."""

    lookup = _ORIGINAL_INTERPOLANT_FOR if interpolant is None else (lambda name: interpolant)
    for module in _LOOKUP_MODULES:
        module.interpolant_for = lookup


def aep(problem, d):
    """AEP [MWh/yr], number of capped bins and convergence of one blade on the active surface."""

    geometry = problem.parameterisation.to_geometry(d, polar_cache=problem.polar_cache)
    lam = float(problem.design.design_tsr)
    vtip = float(problem.design.max_tip_speed_ms)
    p_rated = float(problem.design.rated_power_w)
    air = (float(problem.site.air_density), float(problem.site.kinematic_viscosity))
    edges, mids, _ = wind_speed_bins()
    mass = problem.resource.probability_between(edges[:-1], edges[1:])
    total, n_capped, converged = 0.0, 0, True
    for v, m in zip(mids, mass):
        v = float(v)
        p, res = aerodynamic_power(geometry, v, min(lam, vtip / v), *air)
        converged = converged and bool(res["converged"])
        n_capped += int(p >= p_rated)
        total += min(p, p_rated) * float(m) * HOURS_PER_YEAR / 1e6
    return float(total), n_capped, converged


def reoptimise():
    """
    The mass problem at delta = 0 solved afresh on the active surface, from
    x0: every reference (AEP(x0), KS0, D0, the stress reference) re-measured
    on that surface, the same rows, the same SLSQP settings. Returns the
    saving and the state of the rows.
    """

    problem, u0 = C.prepared_problem()
    result, recorder, wall = C.run_mass_slsqp(problem, u0, 0.0, verbose=False)
    u = np.asarray(result.x, dtype=float)
    slacks = problem.mass_problem_slacks(u, 0.0)
    return {
        "success": bool(result.success), "message": str(result.message), "nit": int(result.nit),
        "wall_s": float(wall), "n_evaluation_failures": len(recorder.failures),
        "shell_pct_vs_x0": 100.0 * (float(problem.mass(u)) - 1.0),
        "aep_pct_vs_x0": 100.0 * (problem.aep_mwh(u) / -float(problem.J0) - 1.0),
        "stress_ratio": float(slacks["stress"]["stress_ratio"]),
        "deflection_ratio": float(slacks["deflection"]["deflection_ratio"]),
        "control_points": C.control_points_record(problem, problem.physical(u)),
    }


def main():
    problem, _ = C.prepared_problem()
    blades = {"x0": C.load_x0(), "x_c": C.load_xc(),
              "x_m": np.array(C.load_json(RESULT_PATH)["x_m"], dtype=float)}
    base_dir = cache_dir(problem.polar_cache)
    s, s_table = lift_error_scale()

    cases = [("production (XFOIL, n_crit = 9)", "production", lambda: None)]
    cases += [
        ("lift + s(Re)", "uniform", lambda: perturbed_interpolant(base_dir, dcl=lambda re: s(re))),
        ("lift - s(Re)", "uniform", lambda: perturbed_interpolant(base_dir, dcl=lambda re: -s(re))),
        ("lift + s(Re) w(Re): low Re up, high Re down", "tilted",
         lambda: perturbed_interpolant(base_dir, dcl=lambda re: s(re) * tilt(re))),
        ("lift - s(Re) w(Re): low Re down, high Re up", "tilted",
         lambda: perturbed_interpolant(base_dir, dcl=lambda re: -s(re) * tilt(re))),
    ]
    cases += [(f"drag x {k:g}", "drag", lambda k=k: perturbed_interpolant(base_dir, cd_scale=k))
              for k in DRAG_SCALES]

    records = []
    for label, kind, build in cases:
        if build is None:
            records.append({"case": label, "kind": kind, "skipped": "cache directory incomplete"})
            print(f"{label:48s} skipped (cache incomplete)", flush=True)
            continue
        use_surface(build())
        try:
            rec = {"case": label, "kind": kind}
            for name, d in blades.items():
                total, n_capped, ok = aep(problem, d)
                rec[name] = {"aep_mwh": total, "n_capped_bins": n_capped, "all_converged": ok}
            a0 = rec["x0"]["aep_mwh"]
            rec["x_m_pct_vs_x0"] = 100.0 * (rec["x_m"]["aep_mwh"] / a0 - 1.0)
            rec["x_c_pct_vs_x0"] = 100.0 * (rec["x_c"]["aep_mwh"] / a0 - 1.0)
            if kind != "production" or REOPTIMISE_PRODUCTION:
                rec["reoptimised"] = reoptimise()
        finally:
            use_surface(None)
        records.append(rec)
        ro = rec.get("reoptimised")
        print(f"{label:48s} AEP x0 {a0:8.4f}  x_m {rec['x_m_pct_vs_x0']:+.4f} %  "
              f"x_c {rec['x_c_pct_vs_x0']:+.4f} %"
              + (f"  | re-optimised: shell {ro['shell_pct_vs_x0']:+.3f} % "
                 f"({'ok' if ro['success'] else ro['message']})" if ro else ""), flush=True)

    run = [r for r in records if "skipped" not in r]
    savings = [r["reoptimised"]["shell_pct_vs_x0"] for r in run if "reoptimised" in r]
    xm = [r["x_m_pct_vs_x0"] for r in run]
    kkt = C.load_json(RESULT_PATH)["kkt"]
    rate = float(kkt["aep_floor_multiplier"])
    worst = min(xm)
    record = {
        "problem": C.PROBLEM_LABEL,
        "description": __doc__.split("\n\n")[0].strip(),
        "command": "python verification/polar_sensitivity/run_polar_sensitivity.py",
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "src_commit": C.src_commit(),
        "blades": {k: [float(v) for v in x] for k, x in blades.items()},
        "lift_error_scale_uiuc_rmse": s_table,
        "tilt_reynolds": list(TILT_RE),
        "drag_scales_assumed": list(DRAG_SCALES),
        "rated_power_w": float(problem.design.rated_power_w),
        "xfoil_ncrit_cases": ("n_crit = 7 and 11 caches attempted 2026-09-26; XFOIL gap-closure "
                              "retries failed at Re = 100 000, run stopped and reverted; not used"),
        "cases": records,
        "x_m_aep_pct_vs_x0_range": [min(xm), max(xm)],
        "reoptimised_shell_pct_range": [min(savings), max(savings)],
        "exchange_rate_pct_per_pct": rate,
        "worst_case_in_material_terms_pct": {
            "note": ("the worst x_m energy shortfall converted to shell material at the "
                     "reference-energy exchange rate: the saving that shortfall is worth "
                     "at the margin, an upper-end reading because the front is concave"),
            "value": -rate * worst if worst < 0 else 0.0},
    }
    with open(OUT_PATH, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=1)
    print(f"\nx_m AEP vs x0 over the cases run: [{min(xm):+.4f}, {max(xm):+.4f}] %")
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
