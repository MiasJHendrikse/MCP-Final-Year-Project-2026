# `verification/`

Versioned report figures and data — the artefacts that back a claim in the
writeup. Work order Task 8, and the plan's working conventions.

Distinct from `results/`, which is generated output and partly gitignored.
The rule is: if regenerating it would change what the dissertation says, it
lives here and is committed. If it is scratch, a sweep you ran once, or an
intermediate you would not cite, it belongs in `results/`.

Each subdirectory carries its own README stating what the artefact shows, how
it was produced, and what would make it wrong.

| Directory | Status | Backs |
|---|---|---|
| `polar_interpolant/` | populated (Task 3) | that the C¹ interpolant removed the bilinear derivative staircase |
| `phase_vi/` | populated (Task 5) | that the solver converges across the envelope with checked residuals |
| `wind_resource/` | populated | the 20 m hub-height Weibull fit and the site extrapolation behind AEP |
| `baseline/` | **re-run 2026-09-19** | the Schmitz `x0` and the reference numbers every result is compared to |
| `spline_fit_error/` | **re-run 2026-09-19** | that the baseline's spline projection is not what A4 "gained" |
| `fd_step_size/` | **re-run 2026-09-19** | the per-variable FD step `h*_j` and the disagreement scale `ε_j` the gradient tests use |
| `fd_optimisation/` | **re-run 2026-09-19** | A4: FD-driven SLSQP from `x0`, end to end |
| `adjoint_optimisation/` | **re-run 2026-09-19** | B5: adjoint-driven SLSQP, and its agreement with A4 |
| `gradient_verification/` | **re-run 2026-09-19** | Tier 3 (adjoint vs the FD noise floor) and Tier 4 (attribution of what is left); **one red variable, reported and accepted** |
| `fd_optimisation_multistart/` | **re-run 2026-09-19** | that A4's gain is not an artefact of the starting point |
| `aep_gain_audit/` | **history — one file re-run** | the measurements behind `docs/AEP_GAIN_AUDIT.md`, under the retired provisional bounds; `cross_evaluate_xstar.json` is the exception, re-run 2026-09-19 with the current `x*` |
| `smoothness_gate/` | populated 2026-09-13, **not re-run** | that the objective is smooth enough to differentiate; its numbers predate the fixed rating and the ceiling, its verdict (C¹ not C², from the Buhl blend) does not depend on them |
| `representation_study/` | populated 2026-09-10, **not re-run** | the choice of blade parameterisation and control-point count; a geometric comparison, unaffected by the bounds or the power model |
| `aep_optimisation_experiment/` | **port, not a re-run** | the 2026-09-19 experiment behind inputs B1 and O4 (the 300 rpm ceiling and the 0.30 m chord cap), frozen as a record |

Statuses are current as of 2026-09-19.

## Re-run order

The 2026-09-19 re-run had to respect a dependency chain — several artefacts
consume a file an earlier artefact writes — so this is the order, not an
alphabetical list:

    baseline               x0.json, baseline_reference.json    →
    fd_step_size           sweep.json (h*_j, ε_j)             → consumed by the gradient tests
    fd_optimisation        result.json, iterates.json         → x* consumed by Tier 3/4, fit_error, multistart
    spline_fit_error       fit_error.json                     → consumes result.json's x*
    gradient_verification  tier3.json, tier4.json             → consumes sweep.json + iterates.json + result.json
    fd_optimisation_multistart  starts.json, results.json     → consumes the feasible region and A4's result.json
    adjoint_optimisation   result.json                        → checked against fd_optimisation's; reads the multistart's
                                                                 `spread_of_optima_u_inf` for its u-agreement tolerance,
                                                                 so it runs AFTER the multistart (corrected 2026-09-19, audit)
    cross_evaluate_xstar   cross_evaluate_xstar.json          → consumes the adjoint x* and evaluates it under 5 ceilings
    aep_optimisation_experiment  (port only)

`cross_evaluate_xstar.py` is the one script inside the frozen `aep_gain_audit/`
directory that was re-run: it reads the *current* `x*` from
`adjoint_optimisation/result.json` and evaluates blades rather than optimising
them, so it has no stale input. Everything else in that directory stays as
written, under the retired provisional bounds.

## Regenerating

Every subdirectory has a generator script next to its output, and the output is
committed alongside it. Regeneration is deliberate, not routine — the same rule
`tests/golden/README.md` states for the golden files. A figure that changes
without an explanation in the commit message is a defect, not an update.
