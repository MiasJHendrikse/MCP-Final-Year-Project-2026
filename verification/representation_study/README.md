# Blade representation study

**Empty. Lands with plan step 1.6.** Created by work order Task 8 so the
location is fixed before the work starts.

## What goes here

Plan step 1.6 builds `src/design/parameterisation.py`: the map from a design
vector `d` to a spanwise chord and twist distribution, with analytic `∂c/∂d`
and `∂θ/∂d`. This directory holds the evidence for the *choices* in that map,
which are not obvious and which the writeup has to defend:

- **Control-point count.** Too few and the optimum is not reachable; too many
  and the optimiser has freedom to produce a wavy blade that scores well and
  cannot be manufactured. `config/rotor_design.yaml` carries
  `n_control_points_chord` and `n_control_points_twist` as separate fields
  precisely so the study can vary them independently.
- **Spline family and knot placement**, and what each does to the achievable
  chord and twist distributions.
- **Strip count against control-point count.** `ParameterisationConfig` keeps
  `n_bem_strips` separate from the design-variable count so the two cannot be
  conflated; this is where the separation is justified rather than asserted.

Expected artefacts: recovered-shape comparisons against the Schmitz baseline
(step 1.7), convergence of the achievable optimum with control-point count,
and the complex-step verification plots for `∂c/∂d` and `∂θ/∂d` to ~1e-14 that
`tests/test_parameterisation.py` will assert numerically.

## Note on the derivative tests

`tests/test_parameterisation.py` is named in work order Task 7's list but does
not exist yet — it has nothing to test until this module lands. When it does,
the dtype trap that caught the polar interpolant applies here too: a NumPy
array pre-allocated with a real dtype silently discards the imaginary part of
a complex-step perturbation, capping observed accuracy near 1e-6 while looking
entirely plausible. If verification plateaus there, suspect dtype before
suspecting the mathematics.
