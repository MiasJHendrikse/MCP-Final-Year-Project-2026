# Adjoint-driven SLSQP, end to end (Phase 3, B5)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. The optimum below is an
optimum *of this provisional box*; it is not a design. Bounds come from
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` and are not in `config/`.

**Re-run 2026-09-19 under the fixed generator rating** (see
`verification/fd_optimisation/README.md` for what changed and why; the
2026-09-13 floating-cap run — +0.217 %, 34 iterations, 24.7 s, interior
optimum — is in the git history). Same script, same A4 comparison, against
the re-run A4.

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
| `nit` / `nfev` / `njev` | **27 / 29 / 27** | 27 / 29 / 27 |
| objective evaluations in total | 32 (29 + 3 post-checks) | 1111 |
| gradient evaluations in total | 54 adjoint (27 + initial + 26 callback re-evaluations, as A4) | 27 × 20 + 26 × 20 FD evaluations |
| wall time | **17.3 s** | 230.3 s |
| AEP at `x0` | 10.270157 MWh/yr | same |
| AEP at optimum | 10.282564 MWh/yr | 10.282564 MWh/yr |
| ΔAEP vs x0 | +0.121 % | +0.121 % |
| active bounds / active envelope rows | **`twist_4` at −2°** / none (tightest: floor at `r = 1.966 m`, 14.9 mm slack) | same / none |
| α at optimum (post-check) | 0.64° … 6.52°, inside [−8°, 18°] | same |
| Re at optimum | 60.1 k … 474 k, inside [40 k, 1 M] | same |
| `PolarDomainError` during the run | none; margin stayed at 0.05 | none |

**The two optima coincide.** the implementation plan §6 B5 requires agreement
to SLSQP's tolerance; the measured differences are

| | value | tolerance (stated in the script) |
|---|---|---|
| `‖u*_adj − u*_fd‖∞` | **1.1e-9** | 0.17 (10 × the spread the multi-start study measured between nine equally converged FD optima, 0.017 — which lies along the near-inert root direction, see that README) |
| `AEP*_adj − AEP*_fd` | **+2.3e-14 MWh/yr** | 1e-7 in `fun` ≈ 1e-6 MWh/yr (`ftol = 1e-8`, a few times over) |
| `fun*_adj − fun*_fd` | −2.2e-15 | 1e-7 |

Not only the endpoints: SLSQP followed the *same trajectory*. Iterate by
iterate against A4's `iterates.json`, `max_k ‖u_k^adj − u_k^fd‖∞ = 3.4e-9`
(at `k = 24`), `max_k |J_k^adj − J_k^fd| = 3.7e-10 MWh/yr`, and the two
gradients at the accepted iterates differ by at most 6.6e-9 MWh/yr per `u`
— the Tier 3 agreement carried through 27 quasi-Newton updates and onto an
active bound without the paths separating. The adjoint gradient and the FD
gradient at `h*` are the same gradient to the optimiser. (The floating-cap
run agreed to 6.8e-7 in `u`; the tighter figure here is not a change in the
gradients but in the trajectory — fewer, shorter terminal steps.)

Optimum (physical, `x_star` in `result.json`; A4's to the same digits):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x0 | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x* | 0.285 | 0.173 | 0.107 | 0.078 | 0.046 m | 23.50 | 10.36 | 3.52 | 2.01 | **−2.00** ° |

`optimised_blade.png` overlays `x0`, A4's optimum (wide, pale) and this
run's optimum; the two optima are indistinguishable at plot resolution.

## What this costs, honestly

The run is **13.3 × faster** than A4 for the same 27 iterations: each
gradient is one forward solve plus the partials (≈ 0.23 s) instead of 20
forward solves (≈ 4.0 s). With `n = 10` design variables and a diagonal
state Jacobian that ratio is roughly `2n × J / (1.15 × J) ≈ 17`, and it
would grow linearly with `n`; the adjoint's cost is independent of `n`.
But the brief's point stands: at `n = 10` the adjoint's *decisive*
advantage is not the wall time, which was already tolerable, but a
gradient free of the FD noise floor and of the O(h)·Δf″ contamination at
the polar interpolant's C² breaks (Tier 4). Here that accuracy did not
change the answer — FD at `h*` was already inside its own noise floor —
which is itself the finding: A4's optimum is confirmed by an independent
gradient.

The 26 callback re-evaluations of the gradient (`callback_extra_evaluations`)
are an artefact of recording `g_k` at the accepted iterate when SLSQP's
last gradient call was elsewhere — the same overhead A4 carried, at 0.23 s
each here instead of 4.0 s.

## The gain is 0.12 %

Below the 2–6 % the brief expected, above neither defect gate (> 15 % or
negative), so reported as a result. `docs/AEP_GAIN_AUDIT.md` is the
explanation (the objective is degenerate in the blade under fixed TSR + an
ideal cap, and Schmitz is its analytic maximiser); the fixed rating removed
the part of the old 0.217 % that was the generator growing with the blade.
The multi-start study (`verification/fd_optimisation_multistart/`, re-run
the same day) shows +0.121 % to be the global optimum of the problem as
posed, and this run shows it is not an artefact of the FD gradient either.
What would make the number larger is a machine fact, not a modelling
choice: inputs B1–B3 in `docs/OUTSTANDING-INPUTS.md` §9.

## Reproduce

    python verification/adjoint_optimisation/run_adjoint_slsqp.py            # ~20 s
    python verification/adjoint_optimisation/run_adjoint_slsqp.py --maxiter 50

## Files

- `run_adjoint_slsqp.py` — the script (A4's with `jac = jac_adjoint` and the comparison block).
- `result.json` — everything in the tables above, plus `u*`, `x*`, the
  gradient at the optimum in both units, the post-checks, the SLSQP
  counters and options, the provisional label, and `comparison_with_fd`
  with the tolerances and the pass flags.
- `iterates.json` — the 28 iterates `(k, u_k, J_k, fun_k, g_k, wall_time)`.
- `optimised_blade.png` — chord and twist: `x0`, A4's optimum, this optimum.
