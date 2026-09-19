# Adjoint-driven SLSQP, end to end (Phase 3, B5)

**Re-run 2026-09-19 under the 300 rpm operating law and the project bounds in
`config/rotor_design.yaml`** (`chord_max_m = 0.30 m`, `chord_min_m = 0.045 m`,
`twist_min_deg = -2`, `twist_max_deg = 35`, local solidity cap 0.5). The
optimum below is an optimum of *that* box under *that* law; it is still not a
design. The retirement of `PROVISIONAL_BOUNDS` and the ordering of this
re-run are in `verification/README.md`.

Two earlier runs are in the git history. The 2026-09-13 floating-cap run:
+0.217 %, 34 iterations, 24.7 s, interior optimum. The 2026-09-19
fixed-rating run: +0.121 %, 27 iterations, 17.3 s, `twist_4` on the −2°
floor. Both were made under the retired `PROVISIONAL_BOUNDS` box
(`chord_max_m = 0.45 m`) and the fixed-λ = 6.5 law. Same script, same A4
comparison, against the re-run A4.

## What was run

    python verification/adjoint_optimisation/run_adjoint_slsqp.py

A4's run (`verification/fd_optimisation/`) with one change: the gradient
handed to SLSQP is the discrete adjoint, `ScaledProblem.jac_adjoint`,
instead of central differences at `h* = 3e-6`.

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10             (5 chord + 5 twist control points, scaled)
    subject to the polar-cache Reynolds envelope: 50 linear rows, margin 5 %
               the root/tip chord envelope and the local solidity cap
               (sigma_i <= 0.5), both from ScaledProblem.constraints()
    gradient   discrete adjoint (one forward solve + partials per gradient)
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-8, maxiter=200

Same starting point (`x0`, the fitted Schmitz blade), same constraint set as
the re-run A4 (both now call `problem.constraints()` rather than the
envelope alone), same margin escalation (not triggered), same iterate
recording, same sanity gates on the gain.

## Result, and the comparison with A4

| quantity | adjoint-driven (this run) | FD-driven (A4) |
|---|---|---|
| termination | `Optimization terminated successfully` (exit mode 0) | same |
| `nit` / `nfev` / `njev` | **35 / 36 / 35** | 35 / 36 / 35 |
| objective evaluations in total | 39 (36 + 3 post-checks) | 1438 |
| gradient evaluations in total | 70 adjoint (35 + initial + 34 callback re-evaluations, as A4) | 35 × 20 + 34 × 20 FD evaluations |
| wall time | **54.4 s** | 367.1 s |
| AEP at `x0` | 10.247707 MWh/yr | same |
| AEP at optimum | 10.262736 MWh/yr | 10.262736 MWh/yr |
| ΔAEP vs x0 | +0.147 % | +0.147 % |
| active bounds / active envelope rows | **`chord_0` at its 0.30 m upper bound** / none (tightest: floor at `r = 1.966 m`, 17.7 mm slack) | same / none |
| active solidity rows | none (σ_max 0.409, cap 0.5; tightest row `r = 0.334 m`, 0.091 slack) | same |
| α at optimum (post-check) | 0.58° … 28.39° over all 17 schedule bins; the 8 bins at or below the rating (up to 10.5 m/s) are 0.58° … 6.54°, inside [−8°, 18°]. The 10 bins above that are on the 300 rpm ceiling, run λ ≈ 3–6 and reach α ≈ 28° in the Viterna extrapolation; they hold a constant power and contribute nothing to the gradient | same |
| Re at optimum | 63.5 k … 339 k, inside [40 k, 1 M] | same |
| `PolarDomainError` during the run | none; margin stayed at 0.05 | none |

The α spread is new and is a consequence of the law, not of the optimum:
the design point sits at λ = 5.712 rather than 6.5, so the bins above
`V_c = 9.67 m/s` see the rotor turning at a lower tip-speed ratio and α
rises with wind speed. `post_check_optimum.uncapped` carries the subset
below the rating; it is the figure that says whether the *optimised* part
of the objective is inside the validated polar band. See
`docs/journal/Session Notes/2026-09-19.md`.

**The two optima coincide.** the implementation plan §6 B5 requires agreement
to SLSQP's tolerance; the measured differences are

| | value | tolerance (stated in the script) |
|---|---|---|
| `‖u*_adj − u*_fd‖∞` | **1.0e-8** | 0.174 (10 × the spread the multi-start study measured between eight equally converged FD optima, 0.0174 — which lies along the near-inert root direction, see that README) |
| `AEP*_adj − AEP*_fd` | **−5.9e-12 MWh/yr** | 1e-7 in `fun` ≈ 1e-6 MWh/yr (`ftol = 1e-8`, a few times over) |
| `fun*_adj − fun*_fd` | +5.7e-13 | 1e-7 |

Not only the endpoints: SLSQP followed the *same trajectory*. Iterate by
iterate against A4's `iterates.json`, `max_k ‖u_k^adj − u_k^fd‖∞ = 1.4e-6`
(at `k = 31`), `max_k |J_k^adj − J_k^fd| = 6.7e-10 MWh/yr`, and the two
gradients at the accepted iterates differ by at most 2.1e-7 MWh/yr per `u`
— the Tier 3 agreement carried through 35 quasi-Newton updates and onto an
active bound without the paths separating. The largest iterate difference
is in the terminal phase, at `k = 29` and `k = 31` of 35, where `chord_0`
is already on its bound and SLSQP's line search is choosing between
near-identical steps; by `k = 35` the two paths are back to 1.0e-8. The
adjoint gradient and the FD gradient at `h*` are the same gradient to the
optimiser. (The floating-cap run agreed to 6.8e-7 in `u`; the tighter
figure here is not a change in the gradients but in the trajectory.)

Optimum (physical, `x_star` in `result.json`; A4's to the same digits):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x0 | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x* | **0.300** | 0.180 | 0.130 | 0.077 | 0.049 m | 24.35 | 10.42 | 5.60 | 2.10 | −1.85 ° |

The first chord control point is on its upper bound, and the tip-twist
point is no longer on the −2° floor — the optimum moved from the twist
floor to the chord ceiling between the two runs, which is the visible
consequence of shortening the box to 0.30 m while the ceiling dropped the
design tip-speed ratio to 5.712.

`optimised_blade.png` overlays `x0`, A4's optimum (wide, pale) and this
run's optimum; the two optima are indistinguishable at plot resolution.

## Loads at the design condition

New in this re-run: `loads_at_rated_ceiling` evaluates both `x0` and `x*`
at the rated wind speed (11 m/s) *on the 300 rpm ceiling*, which is Phase
4's design condition. It is a check on the optimum, not part of the
objective.

| | TSR | rpm | Cp | Ct | rotor thrust | root bending moment |
|---|---|---|---|---|---|---|
| `x0` | 5.712 | 300 | 0.4516 | 0.7029 | 517.5 N | 177.38 N·m |
| `x*` | 5.712 | 300 | 0.4593 | 0.7153 | 526.7 N | 177.92 N·m |

The gain is nearly load-neutral: +0.31 % root bending moment and +1.77 %
thrust for +0.147 % AEP. `x_star_fd` is carried alongside and agrees to
1e-6 N·m, which is the FD optimum re-evaluated (not re-optimised) here.

## What this costs, honestly

The run is **6.8 × faster** than A4 for the same 35 iterations (54.4 s
against 367.1 s) and spends **39 objective evaluations against 1438** —
37 × fewer. Each FD gradient is 20 forward solves at ≈ 0.25 s each
(≈ 5.1 s per gradient); each adjoint gradient is one forward solve plus
the partials, ≈ 0.78 s per adjoint evaluation. With `n = 10` design
variables and a diagonal state Jacobian that ratio is roughly
`2n × J / (1.15 × J) ≈ 17` in evaluations, and it would grow linearly with
`n`; the adjoint's cost is independent of `n`. But the brief's point
stands: at `n = 10` the adjoint's *decisive* advantage is not the wall
time, which was already tolerable, but a gradient free of the FD noise
floor and of the O(h)·Δf″ contamination at the polar interpolant's C²
breaks (Tier 4). Here that accuracy did not change the answer — FD at
`h*` was already inside its own noise floor — which is itself the finding:
A4's optimum is confirmed by an independent gradient, and it is confirmed
*moving to a bound*, where an inaccurate gradient would most easily have
stopped short.

The 34 callback re-evaluations of the gradient (`callback_extra_evaluations`)
are an artefact of recording `g_k` at the accepted iterate when SLSQP's
last gradient call was elsewhere — the same overhead A4 carried, at
≈ 0.78 s each here instead of ≈ 5.1 s.

## The gain is 0.147 %

Below the 2–6 % the brief expected, above neither defect gate (> 15 % or
negative), so reported as a result. `docs/AEP_GAIN_AUDIT.md` is the
explanation (the objective is degenerate in the blade under a fixed TSR
and an ideal cap, and Schmitz is its analytic maximiser); the two
successive weakenings of that degeneracy moved the number from 0.121 % to
0.147 % — the 300 rpm ceiling lowered the design TSR to 5.712, away from
Schmitz's own optimum, so the blade now has something to win and the
optimum has left the interior for the chord ceiling. The fixed rating
removed the part of the old 0.217 % that was the generator growing with
the blade (docs/AEP_GAIN_AUDIT.md section 1.2, section 3.2). The
multi-start study (`verification/fd_optimisation_multistart/`, re-run the
same day) shows +0.147 % to be the global optimum of the problem as
posed, and this run shows it is not an artefact of the FD gradient either.
What would make the number larger is a machine fact, not a modelling
choice: inputs B1–B3 in `docs/OUTSTANDING-INPUTS.md` §9.

## Reproduce

    python verification/adjoint_optimisation/run_adjoint_slsqp.py            # ~55 s
    python verification/adjoint_optimisation/run_adjoint_slsqp.py --maxiter 50

## Files

- `run_adjoint_slsqp.py` — the script (A4's with `jac = jac_adjoint` and the comparison block).
- `result.json` — everything in the tables above, plus `u*`, `x*`, the
  gradient at the optimum in both units, the active set (bounds, envelope
  rows, solidity rows, `sigma_max`), `loads_at_rated_ceiling`, the
  post-checks, the SLSQP counters and options, the bounds label, and
  `comparison_with_fd` with the tolerances and the pass flags.
- `iterates.json` — the 36 iterates `(k, u_k, J_k, fun_k, g_k, wall_time)`.
- `optimised_blade.png` — chord and twist: `x0`, A4's optimum, this optimum.
