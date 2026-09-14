# Adjoint-driven SLSQP, end to end (Phase 3, B5)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. The optimum below is an
optimum *of this provisional box*; it is not a design. Bounds come from
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` and are not in `config/`.

## What was run

    python verification/adjoint_optimisation/run_adjoint_slsqp.py

A4's run (`verification/fd_optimisation/`) with one change: the gradient
handed to SLSQP is the discrete adjoint, `ScaledProblem.jac_adjoint`,
instead of central differences at `h* = 3e-6`.

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10             (5 chord + 5 twist control points, scaled)
    subject to the polar-cache Reynolds envelope: 50 linear rows, margin 5 %
    gradient   discrete adjoint (one forward solve + partials per gradient)
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-8, maxiter=200

Same starting point (`x0`, the fitted Schmitz blade), same envelope
constraint and margin escalation (not triggered), same iterate recording,
same sanity gates on the gain.

## Result, and the comparison with A4

| quantity | adjoint-driven (this run) | FD-driven (A4) |
|---|---|---|
| termination | `Optimization terminated successfully` (exit mode 0) | same |
| `nit` / `nfev` / `njev` | **34 / 37 / 34** | 34 / 37 / 34 |
| objective evaluations in total | 40 (37 + 3 post-checks) | 1419 |
| gradient evaluations in total | 69 adjoint (34 + initial + 34 callback re-evaluations, as A4) | 34 × 20 FD evaluations |
| wall time | **24.7 s** | 360.9 s |
| AEP at `x0` | 10.270157 MWh/yr | same |
| AEP at optimum | 10.292482 MWh/yr | 10.292482 MWh/yr |
| ΔAEP vs x0 | +0.217 % | +0.217 % |
| active bounds / active envelope rows | none / none (tightest: floor at `r = 1.966 m`, 17.0 mm slack) | none / none |
| α at optimum (post-check) | −0.44° … 6.13°, inside [−8°, 18°] | same |
| Re at optimum | 62.6 k … 496 k, inside [40 k, 1 M] | same |
| `PolarDomainError` during the run | none; margin stayed at 0.05 | none |

**The two optima coincide.** the implementation plan §6 B5 requires agreement
to SLSQP's tolerance; the measured differences are

| | value | tolerance (stated in the script) |
|---|---|---|
| `‖u*_adj − u*_fd‖∞` | **8.2e-9** | 8.6e-3 (10 × the spread the multi-start study measured between nine equally converged FD optima; the objective is flat at the optimum, gradient ≈ 1.6e-3 MWh/yr per `u`) |
| `AEP*_adj − AEP*_fd` | **+1.4e-13 MWh/yr** | 1e-7 in `fun` ≈ 1e-6 MWh/yr (`ftol = 1e-8`, a few times over) |
| `fun*_adj − fun*_fd` | −1.3e-14 | 1e-7 |

Not only the endpoints: SLSQP followed the *same trajectory*. Iterate by
iterate against A4's `iterates.json`, `max_k ‖u_k^adj − u_k^fd‖∞ = 6.8e-7`
(at `k = 30`), `max_k |J_k^adj − J_k^fd| = 7.8e-10 MWh/yr`, and the two
gradients at the accepted iterates differ by at most 5.8e-7 MWh/yr per `u`
— the Tier 3 agreement carried through 34 quasi-Newton updates without
the paths separating. The adjoint gradient and the FD gradient at `h*`
are the same gradient to the optimiser.

Optimum (physical, `x_star` in `result.json`; A4's to the same digits):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x0 | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x* | 0.317 | 0.167 | 0.121 | 0.077 | 0.048 m | 25.55 | 9.86 | 4.96 | 2.18 | −1.43 ° |

`optimised_blade.png` overlays `x0`, A4's optimum (wide, pale) and this
run's optimum; the two optima are indistinguishable at plot resolution.

## What this costs, honestly

The run is **14.6 × faster** than A4 for the same 34 iterations: each
gradient is one forward solve plus the partials (≈ 0.25 s) instead of 20
forward solves (≈ 4.3 s). With `n = 10` design variables and a diagonal
state Jacobian that ratio is roughly `2n × J / (1.15 × J) ≈ 17`, and it
would grow linearly with `n`; the adjoint's cost is independent of `n`.
But the brief's point stands: at `n = 10` the adjoint's *decisive*
advantage is not the wall time, which was already tolerable, but a
gradient free of the FD noise floor and of the O(h)·Δf″ contamination at
the polar interpolant's C² breaks (Tier 4). Here that accuracy did not
change the answer — FD at `h*` was already inside its own noise floor —
which is itself the finding: A4's optimum is confirmed by an independent
gradient.

The 34 callback re-evaluations of the gradient (`callback_extra_evaluations`)
are an artefact of recording `g_k` at the accepted iterate when SLSQP's
last gradient call was elsewhere — the same overhead A4 carried, at 0.25 s
each here instead of 4.3 s.

## The gain is still 0.22 %, not 2–6 %

As in A4: below the 2–6 % the brief expected, above neither defect gate
(> 15 % or negative), so reported as a result. The multi-start study
(`verification/fd_optimisation_multistart/`) showed +0.217 % to be the
global optimum of the problem as posed, and this run shows it is not an
artefact of the FD gradient either. The expectation, or the baseline
against which improvement is measured, is the thing to revisit (flagged for
MJ in A4's README).

## Reproduce

    python verification/adjoint_optimisation/run_adjoint_slsqp.py            # ~25 s
    python verification/adjoint_optimisation/run_adjoint_slsqp.py --maxiter 50

## Files

- `run_adjoint_slsqp.py` — the script (A4's with `jac = jac_adjoint` and the comparison block).
- `result.json` — everything in the tables above, plus `u*`, `x*`, the
  gradient at the optimum in both units, the post-checks, the SLSQP
  counters and options, the provisional label, and `comparison_with_fd`
  with the tolerances and the pass flags.
- `iterates.json` — the 35 iterates `(k, u_k, J_k, fun_k, g_k, wall_time)`.
- `optimised_blade.png` — chord and twist: `x0`, A4's optimum, this optimum.
