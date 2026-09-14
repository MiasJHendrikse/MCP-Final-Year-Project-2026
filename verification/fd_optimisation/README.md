# FD-driven SLSQP, end to end (Phase 2, Stage A4)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. The optimum below is an
optimum *of this provisional box*; it is not a design. Bounds come from
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` and are not in `config/`.

## What was run

    python verification/fd_optimisation/run_fd_slsqp.py

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10             (5 chord + 5 twist control points, scaled)
    subject to the polar-cache Reynolds envelope: 50 linear rows, margin 5 %
    gradient   central finite differences, h* = 3e-6 (global, from
               verification/fd_step_size/sweep.json)
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-8, maxiter=200

Starting point `x0` = `verification/baseline/x0.json` (fitted Schmitz). No
solidity constraint (no cap has been decided; the constraint exists in
`ScaledProblem.solidity_constraint(cap)` and was not used).

## Result

| quantity | value |
|---|---|
| termination | `Optimization terminated successfully` (exit mode 0) |
| `nit` / `nfev` / `njev` | 34 / 37 / 34 |
| objective evaluations in total | 1419 (37 + 34 × 20 for the FD gradients + the initial gradient + post-checks) |
| wall time | 360.9 s |
| AEP at `x0` | 10.270157 MWh/yr |
| AEP at optimum | 10.292482 MWh/yr |
| **ΔAEP vs x0** | **+0.217 %** |
| active bounds (`|u|` or `|1-u|` < 1e-6) | none |
| active envelope rows (row value < 1e-6 m) | none; tightest row is the floor at `r = 1.966 m` with 17.0 mm of slack |
| α range at optimum (post-check) | -0.44° … 6.13°, inside the cache band [-8°, 18°] |
| Re range at optimum | 62.6 k … 496 k, inside [40 k, 1 M] — the §4.5 cache extension is not triggered |
| `PolarDomainError` during the run | none; margin stayed at 0.05 |

Optimum (physical, `x_star` in `result.json`):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x0 | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x* | 0.317 | 0.167 | 0.121 | 0.077 | 0.048 m | 25.55 | 9.86 | 4.96 | 2.18 | -1.43 ° |

The optimiser thickened the root and mid-span, thinned the tip (chord_4 to
0.048 m, 3 mm above the grounded 0.045 m floor, `u = 0.008`, not active) and
twisted the tip to -1.4° (0.6° above the -2° bound, `u = 0.016`, not active).
`optimised_blade.png` overlays the two distributions with their control
points at the Greville abscissae.

Gradient at the optimum (MWh/yr per unit `u`): largest component 1.6e-3
(`chord_3`), down from 1.06 at `x0`. SLSQP stopped on `ftol`: the last
accepted steps changed `fun` by less than 1e-8, i.e. less than 1e-7 MWh/yr.
The optimum is interior (no active bound, no active envelope row), so the
first-order condition is `g = 0` and the residual 1e-3 is the size of the
gradient SLSQP's step could no longer reduce `fun` against.

## Iterates (`iterates.json`)

35 records, `k = 0 … 34`, each with `u_k`, `J_k` (MWh/yr), `fun_k`, `g_k`
(both units) and wall time. `g_k` is the gradient SLSQP evaluated at the
accepted iterate (cached from the solver's own call; `callback_extra_evaluations`
= 0). **The mid-run point for the adjoint's Tier 3 check is `k = 17`
(`nit / 2`).**

The first iterate overshoots: SLSQP's initial step (identity Hessian
against a gradient of norm 1.46 in `fun` units) takes AEP *down* to
9.419 MWh/yr at `k = 1`, and the inexact line search accepts it. It recovers
at `k = 2` (10.258) and passes `x0` at `k = 3` (10.288). Everything after
`k = 6` is refinement in the fourth decimal.

## The gain is below the expected band — recorded, not hidden

the implementation plan §5 expected 2–6 %. The result is **0.22 %**. That is not
one of the two defect gates (> 15 % or negative), so the run is reported as
a result, but the discrepancy with the expectation needs saying:

- `x0` is the analytic Schmitz optimum (design `Cl` at `(L/D)_max`, design
  `λ = 6.5`) least-squares fitted onto the spline. The objective operates at
  fixed `λ` in every bin, so below rated it is a Weibull-weighted `Cp(λ = 6.5)`
  and above rated it is `P_aero(11 m/s)` — the same `Cp` again. Maximising
  AEP under this operating strategy is therefore very nearly maximising
  `Cp` at the design `λ`, which is exactly what Schmitz already does; what
  is left to gain is the Reynolds dependence across bins, the tip/hub-loss
  and Buhl corrections Schmitz ignores, and the spline fit error.
- Every check on the run itself passes: converged, interior, no domain
  error, α and Re well inside the cache, the gradient three decades down.
- Whether a different starting point (e.g. a deliberately poor blade) reaches
  the same optimum is the natural next check and is not in Stage A's scope.

This is flagged for MJ: the 2–6 % expectation, or the choice of baseline
against which improvement is measured, should be revisited before the
figure is quoted.

## Reproduce

    python verification/fd_optimisation/run_fd_slsqp.py            # ~6 min
    python verification/fd_optimisation/run_fd_slsqp.py --step 1e-5 --maxiter 50

`--step` overrides `h*`; the default reads `h_star_global` from the
step-size study.

## Files

- `run_fd_slsqp.py` — the script (margin escalation 0.05 → 0.10 on a
  `PolarDomainError`, once, recorded; not triggered).
- `result.json` — everything in the table above, plus `u*`, `x*`, the
  gradient at the optimum in both units, the post-checks at `x0` and `x*`,
  the SLSQP counters and options, the provisional label.
- `iterates.json` — the 35 iterates.
- `optimised_blade.png` — chord and twist, `x0` vs optimum.
