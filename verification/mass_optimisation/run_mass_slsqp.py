"""
Phase 5 (2026-09-20) -- SLSQP on the mass problem at energy floor `delta`.

    minimise   mass(u) = m_shell(d(u)) / m_shell(x0)
    over       u in [0, 1]^10
    subject to `ScaledProblem.constraints_for_mass_problem(delta)`:
               envelope (50), solidity (25), monotone chord/twist (4 + 4),
               AEP floor, moment cap, root-stress proxy, tip-deflection proxy
    gradients  the objective adjoint (AEP floor), the moment adjoint (cap and
               stress), the deflection adjoint; the objective's gradient is
               exact and constant
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-10,
               maxiter=400, cold from x0 (default) or from x_c

The summary records the optimum `x_m`, both material proxies as % of `x0`,
AEP, KS, the rated-point moment, thrust and deflection, the stress and
deflection ratios, every row's slack, the KKT multipliers (the AEP row's is
the exchange rate, % material per % energy at the margin), the post-check,
`reynolds_min`, the failed trial points, the counters and the wall time.
The `--delta 0` run from `x0` is the production optimum.

Outputs, next to this script: `result_delta{D}.json`, `iterates_delta{D}.json`,
`blade_delta{D}_geometry.png`, `blade_delta{D}_loads.png`.

Run from the repo root (one to two minutes):

    python verification/mass_optimisation/run_mass_slsqp.py --delta 0
    python verification/mass_optimisation/run_mass_slsqp.py --delta 0.01 --start xc

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import _common as C  # noqa: E402
from gradients.problem import (ALL_MASS_PROBLEM_ROWS, DEFAULT_ENVELOPE_MARGIN,  # noqa: E402
                               MASS_PROBLEM_ROWS)
from polars.interpolant import PolarDomainError  # noqa: E402


def start_point(problem, name):
    if name == "x0":
        return np.asarray(problem.u0, dtype=float)
    if name == "xc":
        return problem.scaled(C.load_xc())
    raise ValueError(f"unknown start {name!r}")


def optimise_with_escalation(delta, start, include, args):
    """
    Solve at `delta`, raising the envelope margin to 0.10 once and restarting
    if a `PolarDomainError` escapes from an accepted iterate's Jacobian --
    the margin-escalation rule of the earlier scripts, recorded. Returns
    `(problem, u0, result, recorder, wall, margin_raised, domain_errors)`.
    """

    margin = DEFAULT_ENVELOPE_MARGIN
    margin_raised = False
    domain_errors = []
    while True:
        problem, u0 = C.prepared_problem(margin)
        try:
            result, recorder, wall = C.run_mass_slsqp(
                problem, start_point(problem, start), delta, include,
                maxiter=args.maxiter, ftol=args.ftol)
            return problem, u0, result, recorder, wall, margin_raised, domain_errors
        except PolarDomainError as error:
            domain_errors.append({"margin": margin, "error": str(error),
                                  "u": problem.domain_errors[-1] if problem.domain_errors else None})
            if margin_raised:
                raise
            print(f"PolarDomainError at margin {margin}: raising margin to 0.10 once "
                  f"and restarting", flush=True)
            margin, margin_raised = 0.10, True


def parse_rows(text):
    """
    The `--rows` subset in `MASS_PROBLEM_ROWS` order. Every token must be a
    known row: a mistyped one used to be dropped silently, and the run then
    solved a looser problem under the full problem's label (review 2026-09-20).
    """

    tokens = [t.strip() for t in text.split(",") if t.strip()]
    unknown = [t for t in tokens if t not in ALL_MASS_PROBLEM_ROWS]
    if unknown:
        raise ValueError(f"unknown --rows {unknown}; known rows are {list(ALL_MASS_PROBLEM_ROWS)}")
    return tuple(r for r in ALL_MASS_PROBLEM_ROWS if r in tokens)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--delta", type=float, default=0.0,
                        help="energy floor: AEP >= (1 - delta) AEP(x0)")
    parser.add_argument("--start", choices=("x0", "xc"), default="x0")
    parser.add_argument("--rows", default=",".join(MASS_PROBLEM_ROWS),
                        help="comma-separated subset of the rows (ablation)")
    parser.add_argument("--maxiter", type=int, default=C.SLSQP_MAXITER)
    parser.add_argument("--ftol", type=float, default=C.SLSQP_FTOL)
    parser.add_argument("--tag", default=None, help="output tag (default: the delta)")
    args = parser.parse_args(argv)
    include = parse_rows(args.rows)

    problem, u0, result, recorder, wall, margin_raised, domain_errors = \
        optimise_with_escalation(args.delta, args.start, include, args)
    names = C.variable_names(problem)
    summary = C.summarise_run(problem, u0, args.delta, include, result, recorder, wall,
                              args.start, {"ftol": args.ftol, "maxiter": args.maxiter}, names)
    summary["command"] = (f"python verification/mass_optimisation/run_mass_slsqp.py "
                          f"--delta {args.delta:g} --start {args.start}"
                          + ("" if include == MASS_PROBLEM_ROWS else f" --rows {args.rows}"))
    summary["envelope_margin"] = problem.margin
    summary["envelope_margin_raised"] = bool(margin_raised)
    summary["domain_errors"] = domain_errors
    summary["blades"] = {
        "x0": C.blade_record(problem, u0, args.delta, "x0", u0=u0),
        "x_c": C.blade_record(problem, problem.scaled(C.load_xc()), args.delta, "x_c", u0=u0),
        "x_m": summary["optimum"],
    }

    tag = args.tag if args.tag is not None else f"{args.delta:g}"
    result_path = os.path.join(_HERE, f"result_delta{tag}.json")
    iterates_path = os.path.join(_HERE, f"iterates_delta{tag}.json")
    with open(result_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1)
    with open(iterates_path, "w", encoding="utf-8") as handle:
        json.dump({"bounds_label": C.BOUNDS_LABEL, "delta": args.delta, "rows": list(include),
                   "start": args.start, "iterates": recorder.iterates}, handle, indent=1)
    figure_paths = C.plot_blade(problem, [("x0", summary["blades"]["x0"]),
                                          ("x_c", summary["blades"]["x_c"]),
                                          ("x_m", summary["blades"]["x_m"])],
                                _HERE, stem=f"blade_delta{tag}")

    opt = summary["optimum"]
    print(f"\n{result.message}")
    print(f"delta={args.delta:g}: nit={result.nit} nfev={result.nfev} njev={result.njev} "
          f"wall={wall:.1f}s; solves: J {problem.n_fun_evals}, adjoint {problem.n_adjoint_evals}, "
          f"moment {problem.n_moment_solves}, deflection {problem.n_deflection_solves}; "
          f"failed trials {summary['n_evaluation_failures']}")
    print(f"material: shell {opt['shell_pct_vs_x0']:+.3f} %, solid {opt['solid_pct_vs_x0']:+.3f} % "
          f"vs x0; AEP {opt['aep_pct_vs_x0']:+.4f} % vs x0")
    s = opt["slacks"]
    print(f"rows: AEP slack {s['aep_floor']['slack']:.2e} (active {s['aep_floor']['active']}); "
          f"moment KS/KS0 {s['moment']['ks_over_ks0']:.4f} (active {s['moment']['active']}); "
          f"stress ratio {s['stress']['stress_ratio']:.4f} (active {s['stress']['active']}); "
          f"deflection ratio {s['deflection']['deflection_ratio']:.4f} "
          f"(active {s['deflection']['active']}); monotone active {s['manufacturing']['rows']}")
    print(f"control points: chord {np.round(opt['control_points']['chord_mm'], 1)} mm, "
          f"twist {np.round(opt['control_points']['twist_deg'], 2)} deg")
    kkt = summary["kkt"]
    print(f"KKT: active {kkt['active_rows']}; multipliers {kkt['multipliers']}; "
          f"exchange rate {kkt['exchange_rate_pct_material_per_pct_energy']} % material per % energy; "
          f"residual {kkt['residual_norm_after_projection']:.2e}")
    print(f"post-check: reynolds_min {opt['reynolds_min']:.0f}; alpha within (uncapped) "
          f"{opt['post_check']['uncapped']['within'] if opt['post_check']['uncapped'] else None}")
    for note in summary["sanity_notes"]:
        print(f"*** {note} ***")
    print(f"wrote {result_path}\n      {iterates_path}\n      "
          + "\n      ".join(figure_paths))
    return summary


if __name__ == "__main__":
    main()
