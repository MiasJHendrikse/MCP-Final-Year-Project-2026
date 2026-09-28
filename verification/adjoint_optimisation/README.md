# Energy optimisation driven by adjoint gradients

This repeats the finite-difference-driven run in `verification/fd_optimisation/`
with one change: the gradient passed to SLSQP is the discrete adjoint
(`ScaledProblem.jac_adjoint`) instead of central differences at `h* = 3e-6`.

**Result: the same optimum (+0.147 % AEP), reached along the same path, 15
times faster and with 37 times fewer objective evaluations.**

The run uses the 300 rpm operating law and the design bounds in
`config/rotor_design.yaml` (chord 0.045–0.30 m, twist −2° to 35°, local
solidity at most 0.5).

## The problem

    python verification/adjoint_optimisation/run_adjoint_slsqp.py

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10             (5 chord + 5 twist control points, scaled)
    subject to the polar-table Reynolds envelope: 50 linear rows, margin 5 %
               the root/tip chord envelope and the local solidity limit
               (sigma_i <= 0.5), both from ScaledProblem.constraints()
    gradient   discrete adjoint (one forward solve plus the partials per gradient)
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-8, maxiter=200

The starting point (x₀), the constraints, the iterate recording and the
sanity checks on the gain are the same as in the finite-difference run.

This script reads `verification/fd_optimisation_multistart/results.json` to
set its tolerance for agreement in `u` (10 times the spread between the
multi-start optima), so it has to run after the multi-start study.

## Result, compared with the finite-difference run

| quantity | adjoint-driven | finite-difference-driven |
|---|---|---|
| termination | `Optimization terminated successfully` (exit mode 0) | same |
| `nit` / `nfev` / `njev` | **35 / 36 / 35** | 35 / 36 / 35 |
| objective evaluations in total | 39 (36 plus 3 post-checks) | 1438 |
| gradient evaluations in total | 70 adjoint (35, plus the initial one, plus 34 re-evaluations in the callback) | 35 × 20 plus 34 × 20 function evaluations |
| wall time | **24.2 s** | 367.1 s |
| AEP at x₀ | 10.247707 MWh/yr | same |
| AEP at the optimum | 10.262736 MWh/yr | 10.262736 MWh/yr |
| AEP change vs x₀ | +0.147 % | +0.147 % |
| active bounds / envelope rows | **`chord_0` at its 0.30 m upper bound** / none (tightest: the floor at `r = 1.966 m`, 17.7 mm slack) | same / none |
| active solidity rows | none (largest solidity 0.409 against 0.5; tightest at `r = 0.334 m`, 0.091 slack) | same |
| angle of attack at the optimum | 0.58° to 28.39° over all 17 bins. The 8 bins at or below rated (up to 10.5 m/s) are 0.58° to 6.54°, inside −8° to 18°. The bins above that are on the 300 rpm ceiling, run at λ ≈ 3–6 and reach α ≈ 28° in the Viterna extrapolation; they hold constant power and contribute nothing to the gradient | same |
| Reynolds number at the optimum | 63.5 k to 339 k, inside 40 k to 1 M | same |
| polar-range errors | none; the margin stayed at 0.05 | none |

The wide spread in angle of attack comes from the operating law, not the
optimum. Above the ceiling wind speed (9.67 m/s) the rotor runs at a lower
tip-speed ratio, so α rises with wind speed. `post_check_optimum.uncapped` in
`result.json` holds the bins below rated, which is the part of the objective
the optimiser can actually change; those are inside the validated polar
range.

**The two optima are the same.** The run was required to match the
finite-difference optimum to within SLSQP's tolerance:

| | value | tolerance (set in the script) |
|---|---|---|
| `‖u*_adj − u*_fd‖∞` | **1.0e-8** | 0.166 (10 times the 0.0166 spread the multi-start study measured between eight equally converged optima, which lies along the `chord_0` root direction) |
| `AEP*_adj − AEP*_fd` | **−5.9e-12 MWh/yr** | 1e-7 in `fun`, about 1e-6 MWh/yr (a few times `ftol = 1e-8`) |
| `fun*_adj − fun*_fd` | +5.7e-13 | 1e-7 |

It isn't only the end points. SLSQP followed the *same path*. Comparing iterate
by iterate with the finite-difference run's `iterates.json`, the largest
difference in `u` is 1.4e-6 (at `k = 31`), the largest difference in `J` is
6.7e-10 MWh/yr, and the two gradients at the accepted iterates differ by at
most 2.1e-7 MWh/yr per `u`. The Tier 3 agreement carried through 35
quasi-Newton updates, and onto an active bound, without the paths separating.
The largest difference is near the end (`k = 29` and 31 of 35), where
`chord_0` is already on its bound and the line search is choosing between
nearly identical steps; by `k = 35` the paths are back within 1.0e-8. To the
optimiser, the adjoint gradient and the finite-difference gradient at `h*`
are the same gradient.

The optimum in physical units (`x_star` in `result.json`, matching the
finite-difference run to the digits shown):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x₀ | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x\* | **0.300** | 0.180 | 0.130 | 0.077 | 0.049 m | 24.35 | 10.42 | 5.60 | 2.10 | −1.85 ° |

The first chord control point is on its upper bound. `optimised_blade.png`
overlays x₀, the finite-difference optimum (wide, pale line) and this one;
the two optima can't be told apart at plot resolution.

## Loads at the rated point

`loads_at_rated_ceiling` evaluates x₀ and x\* at the rated wind speed (11 m/s)
on the 300 rpm ceiling. This is a check on the optimum, not part of the
objective.

| | TSR | rpm | Cp | Ct | rotor thrust | root bending moment |
|---|---|---|---|---|---|---|
| x₀ | 5.712 | 300 | 0.4516 | 0.7029 | 517.5 N | 177.38 N·m |
| x\* | 5.712 | 300 | 0.4593 | 0.7153 | 526.7 N | 177.92 N·m |

The gain costs little in load: +0.31 % root bending moment and +1.77 % thrust
for +0.147 % AEP. `x_star_fd` (the finite-difference optimum, re-evaluated)
agrees to 1e-6 N·m.

## What the adjoint buys

The run is **15 times faster** than the finite-difference run for the same
35 iterations (24.2 s against 367.1 s) and uses **39 objective evaluations
instead of 1438**, 37 times fewer. Each finite-difference gradient is 20
forward solves at about 0.25 s each (about 5.1 s per gradient); each adjoint
gradient is one forward solve plus the partials, about 0.35 s. With ten
design variables that ratio is roughly `2n × J / (1.15 × J) ≈ 17`, and it
grows linearly with `n`, because the adjoint's cost doesn't depend on the
number of variables (`verification/cost_scaling/`).

At `n = 10`, though, the bigger advantage isn't the time, which was already
tolerable. It's a gradient free of the finite-difference noise floor and of
contamination at the polar interpolant's C² breaks (Tier 4). Here that extra
accuracy didn't change the answer, because the finite difference at `h*` was
already within its own noise floor. That is itself the finding: an
independent gradient confirms the finite-difference optimum, including the
move onto a bound, which is exactly where an inaccurate gradient would most
likely have stopped short.

The 34 extra gradient evaluations in the callback
(`callback_extra_evaluations`) come from recording `g_k` at each accepted
iterate when SLSQP's last gradient call was elsewhere. The finite-difference
run has the same overhead, at about 0.78 s each here instead of 5.1 s.

## Why the gain is 0.147 %

That's below the 2–6 % I originally expected, but the result isn't suspicious:
it isn't negative and it's nowhere near 15 %. `docs/DESIGN-BASIS.md` §2
explains why. In short, with a fixed tip-speed ratio and an ideal power cap,
the Schmitz blade is the analytic optimum. The 300 rpm ceiling moves the
rated point to λ = 5.712, away from Schmitz's own design point, which gives
the blade a little to gain and moves the optimum onto the chord bound. The
multi-start study (`verification/fd_optimisation_multistart/`) shows +0.147 %
is the global optimum of this problem, and this run shows it isn't an
artefact of the finite-difference gradient either. A larger gain would need a
different machine, not a different model choice (`docs/OUTSTANDING-INPUTS.md`
§9).

## Reproducing it

    python verification/adjoint_optimisation/run_adjoint_slsqp.py            # ~55 s
    python verification/adjoint_optimisation/run_adjoint_slsqp.py --maxiter 50

## Files

- `run_adjoint_slsqp.py`: the finite-difference script with `jac = jac_adjoint`
  and a comparison block.
- `result.json`: everything in the tables above, plus `u*`, `x*`, the gradient
  at the optimum in both units, the active set, `loads_at_rated_ceiling`, the
  post-checks, the SLSQP counters and options, the bounds used, and
  `comparison_with_fd` with the tolerances and pass flags.
- `iterates.json`: the 36 iterates `(k, u_k, J_k, fun_k, g_k, wall_time)`.
- `optimised_blade.png`: chord and twist of x₀ and both optima. This figure
  predates the shared figure style and is kept as the run produced it.
