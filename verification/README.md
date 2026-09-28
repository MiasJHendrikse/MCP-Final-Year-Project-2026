# Verification

This folder holds the evidence behind every number and figure in the project.
Each subfolder contains the script that produced its results, the results
themselves (usually JSON), the figures, and a README explaining what was
checked, what came out, how to reproduce it and what would show it to be
wrong.

It is separate from `results/`, which holds generated output that isn't cited
anywhere. The rule I followed: if regenerating something would change what the
report says, it lives here and is committed.

## What's here

**Setting up the problem**

| Folder | What it shows |
|---|---|
| `wind_resource/` | the site's Weibull wind distribution at the 20 m hub height, extrapolated from 50 m |
| `polar_interpolant/` | the smooth (C¹) polar interpolant removes the derivative "staircase" of plain bilinear interpolation |
| `phase_vi/` | the BEM solver converges, with checked residuals, across the NREL Phase VI operating envelope |
| `representation_study/` | how many spline control points the blade needs (5 chord + 5 twist) |
| `baseline/` | the Schmitz reference blade x₀ that every result is compared with |
| `spline_fit_error/` | fitting the Schmitz blade to the spline costs almost nothing, so it isn't where the optimiser's gain comes from |
| `smoothness_gate/` | the energy function is smooth enough to differentiate (C¹, with explained curvature breaks) |

**Gradients**

| Folder | What it shows |
|---|---|
| `fd_step_size/` | the best finite-difference step for each variable, and the noise floor used to judge agreement |
| `gradient_verification/` | adjoint gradients agree with finite differences to the noise floor (Tier 3), and what's left over is round-off (Tier 4) |
| `cost_scaling/` | an adjoint gradient costs about one objective evaluation, whatever the number of design variables |

**Optimisation**

| Folder | What it shows |
|---|---|
| `fd_optimisation/` | SLSQP driven by finite-difference gradients, run from x₀ to convergence |
| `fd_optimisation_multistart/` | the energy optimum doesn't depend on the starting point |
| `adjoint_optimisation/` | SLSQP driven by adjoint gradients reaches the same optimum, far faster |
| `load_constraint/` | the root-moment constraint, the moment-capped energy optimum x_c, and the energy–moment trade-off |
| `mass_optimisation/` | **the main result**: the minimum-material blade x_m at the reference energy (−3.38 % shell material; −5.27 % within 0.5 % of the energy), multi-start, the energy–material trade-off, which constraints matter, and robustness checks |
| `absolute_material/` | the three blades in kilograms, megapascals and millimetres for an assumed laminate, and the stress along the span |

**Further checks**

| Folder | What it shows |
|---|---|
| `polar_sensitivity/` | how the results change if the SG6043 polars are wrong by a plausible amount |
| `starting/` | the parked-rotor starting torque of the three blades |
| `aep_gain_audit/` | an earlier investigation into why maximising energy gained so little. It used an earlier model configuration and is kept as a record |
| `aep_optimisation_experiment/` | an earlier experiment with rotor-speed ceilings that informed the 300 rpm choice. Also kept as a record |

## Re-running

Several results feed into later ones, so they must be re-run in this order:

    baseline                     x0.json, baseline_reference.json
    fd_step_size                 sweep.json (h*_j, ε_j)                 used by the gradient tests
    fd_optimisation              result.json, iterates.json             x* used by Tiers 3/4, spline_fit_error, multi-start
    spline_fit_error             fit_error.json
    gradient_verification        tier3.json, tier4.json
    fd_optimisation_multistart   starts.json, results.json
    adjoint_optimisation         result.json                            needs the multi-start spread for its tolerance
    load_constraint              checks.json, result_eps*.json, pareto.json
    aep_gain_audit/cross_evaluate_xstar.py                              evaluates the current x* under five ceilings
    mass_optimisation            result_delta*.json, multistart_delta0.json, pareto.json,
                                 ablation.json, checks.json, reference_blades.json,
                                 cross_evaluation.json
    absolute_material            result.json, absolute_material.png
    polar_sensitivity            polar_sensitivity.json
    cost_scaling                 scaling.json                           independent of the optimisation results

`cross_evaluate_xstar.py` is the only script in `aep_gain_audit/` that is
kept up to date. It only evaluates the current energy optimum, so it has no
out-of-date inputs; the rest of that folder stays as it was run.

Re-running is deliberate, not routine. A result or figure that changes
without a reason given in the commit message should be treated as a bug.

## Redrawing the figures

Every figure is drawn through `src/plotting/figstyle.py`, which sets one
consistent style (serif text, labelled axes with units, at most two panels
side by side). It refuses to save a figure with more than two panels, a
super-title, or any project-history text in it, such as a date or a phase
number.

The figures are redrawn from the committed JSON, without re-running the
physics:

    python verification/baseline/generate_baseline.py --replot
    python verification/phase_vi/generate_residual_histories.py --replot
    python verification/polar_interpolant/generate_plots.py
    python verification/representation_study/run_study.py --replot
    python verification/fd_step_size/run_sweep.py --replot
    python verification/gradient_verification/run_tier3.py --replot
    python verification/gradient_verification/run_tier4.py --replot
    python verification/cost_scaling/run_scaling.py --replot
    python verification/mass_optimisation/replot_figures.py
    python verification/absolute_material/run_absolute_material.py --replot
    python src/validation/plot_bem_comparison.py
    python src/validation/compare_sg6043_uiuc.py --replot

The exception is `smoothness_gate/run_gate.py`, which has no `--replot`
option: its figure *is* the run (3,000 objective evaluations, about 13
minutes), so it has to be re-run whenever the objective changes.

A few older figures are no longer produced by any script (`baseline.png`,
`smoothness_gate.png`, `blade_delta0.png`, `pareto.png`, and some figures in
`load_constraint/`, `fd_optimisation/`, `adjoint_optimisation/`,
`aep_optimisation_experiment/` and `spline_fit_error/`). They predate the
style module and are kept as the output of the runs that made them.
