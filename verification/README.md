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
| `load_constraint/` | **new 2026-09-19** | Phase 4: the root-moment KS constraint -- Tiers 1-3, the `eps = 0` constrained optimum (`x_c`, the energy optimum -- the Phase 5 production optimum until the 2026-09-20 re-pitch; now the reference blade the mass optimum is set against), and the 0/2/5/10 % Pareto |
| `mass_optimisation/` | **new 2026-09-20** | Phase 5: the AEP-constrained minimum-material blade -- `x_m` at `delta = 0` (-4.24 % shell material at equal energy), the 10-start agreement, the 0/0.25/0.5/1/2 % energy-floor front with KKT exchange rates, the ablation of the stress, deflection, monotone and moment rows, Tiers 1-3 of the stress and deflection Jacobians at `x0` and `x_m`, the `x0 / x_c / x_m` table, and the cross-evaluation under other resources and ceilings |
| `cost_scaling/` | **new 2026-09-19** | Phase 4 Step 3: wall time of `J`, the FD, tangent, adjoint and moment-adjoint gradients against `n = 10 ... 160`, with the `t = a n^p` fits and the forward-solve counts |
| `gradient_verification/` | **re-run 2026-09-19** | Tier 3 (adjoint vs the FD noise floor) and Tier 4 (attribution of what is left); **one red variable, reported and accepted** |
| `fd_optimisation_multistart/` | **re-run 2026-09-19** | that A4's gain is not an artefact of the starting point |
| `aep_gain_audit/` | **history — one file re-run** | the measurements behind `docs/AEP_GAIN_AUDIT.md`, under the retired provisional bounds; `cross_evaluate_xstar.json` is the exception, re-run 2026-09-19 with the current `x*` |
| `smoothness_gate/` | populated 2026-09-13, **not re-run** | that the objective is smooth enough to differentiate; its numbers predate the fixed rating and the ceiling, its verdict (C¹ not C², from the Buhl blend) does not depend on them |
| `representation_study/` | populated 2026-09-10, **not re-run** | the choice of blade parameterisation and control-point count; a geometric comparison, unaffected by the bounds or the power model |
| `aep_optimisation_experiment/` | **port, not a re-run** | the 2026-09-19 experiment behind inputs B1 and O4 (the 300 rpm ceiling and the 0.30 m chord cap), frozen as a record |

Statuses are current as of 2026-09-20.

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
    load_constraint       checks.json, result_eps*.json,      → consumes adjoint_optimisation/result.json's u* and the
                          pareto.json                            baseline; runs AFTER adjoint_optimisation
    cross_evaluate_xstar   cross_evaluate_xstar.json          → consumes the adjoint x* and evaluates it under 5 ceilings
    mass_optimisation     result_delta*.json, multistart_delta0.json, → consumes baseline/x0.json, load_constraint/result_eps0.json (x_c),
                          pareto.json, ablation.json, checks.json,     the multistart spread and fd_step_size/sweep.json's h*; runs AFTER
                          reference_blades.json, cross_evaluation.json load_constraint (2026-09-20)
    cost_scaling           scaling.json, scaling.png          → reads only h* (fd_step_size) and AEP(x0) (baseline);
                                                                 independent of every optimisation artefact
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
