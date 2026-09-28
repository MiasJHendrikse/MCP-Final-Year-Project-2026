# Energy optimisation driven by finite-difference gradients

This is the first end-to-end optimisation: SLSQP maximising annual energy
from the Schmitz blade x₀, with gradients from central finite differences. It
was built first as a safe baseline, so that there would be a working result
before the adjoint existed. `verification/adjoint_optimisation/` repeats the
same run with adjoint gradients.

**Result: +0.147 % AEP**, with the root chord pushed onto its 0.30 m upper
bound.

The run uses the 300 rpm operating law, the design bounds in
`config/rotor_design.yaml` (chord 0.045–0.30 m, twist −2° to 35°, local
solidity at most 0.5) and the fixed generator rating of 3822.19 W (x₀'s own
aerodynamic power at 11 m/s). This is an optimum for that box and that
operating law, not a finished design.

## The problem

    python verification/fd_optimisation/run_fd_slsqp.py

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10             (5 chord + 5 twist control points, scaled)
    subject to the polar-table Reynolds envelope: 50 linear rows, margin 5 %
               the root/tip chord envelope and the local solidity limit
               (sigma_i <= 0.5), both from ScaledProblem.constraints()
    gradient   central finite differences, h* = 3e-6 (global, from
               verification/fd_step_size/sweep.json)
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-8, maxiter=200

It starts from `verification/baseline/x0.json` (the fitted Schmitz blade).
The constraints come from `problem.constraints()`, the same call the adjoint
run makes, so the two runs can't end up solving different problems.

## Result

| quantity | value |
|---|---|
| termination | `Optimization terminated successfully` (exit mode 0) |
| `nit` / `nfev` / `njev` | 35 / 36 / 35 |
| objective evaluations in total | 1438 (36, plus 35 × 20 for the gradients, plus 34 × 20 for gradient re-evaluations in the callback, plus post-checks) |
| wall time | 367.1 s |
| AEP at x₀ | 10.247707 MWh/yr |
| AEP at the optimum | 10.262736 MWh/yr |
| **AEP change vs x₀** | **+0.147 %** |
| active bounds (`u` within 1e-6 of a bound) | **`chord_0` at its 0.30 m upper bound** (`u = 1.0`) |
| active envelope rows | none; the tightest is the floor at `r = 1.966 m`, with 17.7 mm of slack |
| active solidity rows | none; the largest solidity is 0.4089 against the 0.5 limit, at `r = 0.334 m` with 0.091 of slack |
| angle of attack at the optimum | 0.58° to 28.39° over the 17 bins. The 8 at or below rated are 0.58° to 6.54°, inside the polar table's −8° to 18°. The 10 above that are on the rotor-speed ceiling, hold constant power and contribute nothing to the gradient |
| Reynolds number at the optimum | 63.5 k to 338.9 k, inside the table's 40 k to 1 M |
| polar-range errors during the run | none; the envelope margin stayed at 0.05 |

The optimum in physical units (`x_star` in `result.json`):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x₀ | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x\* | **0.300** | 0.180 | 0.130 | 0.077 | 0.049 m | 24.35 | 10.42 | 5.60 | 2.10 | −1.85 ° |

The optimiser pushes the root chord onto its 0.30 m upper bound and moves
chord inboard (c2 from 0.097 to 0.130 m), while leaving the tip twist just
above its −2° floor. The ceiling moves the rated operating point to
λ = 5.712, away from the λ = 6.5 the Schmitz blade is designed for, so there
is something to gain and the optimum lies on a bound rather than inside the
box. `optimised_blade.png` overlays the two blades, with their control points
at the Greville abscissae.

At the optimum, the gradient (MWh/yr per unit `u`) is −7.0e-4 on `chord_0`
(negative, as an upper bound requires: AEP would rise if the chord could
exceed 0.30 m) and at most 3.7e-4 on the other nine variables, down from 0.21
at x₀. SLSQP stopped on `ftol`: the last accepted step changed `fun` by less
than 1e-8, or less than 1e-7 MWh/yr.

## Iterates

`iterates.json` has 36 records (`k = 0 … 35`), each with `u_k`, `J_k`
(MWh/yr), `fun_k`, `g_k` (in both units) and the wall time. `g_k` is the
gradient at the accepted iterate. Where SLSQP's last gradient call was at a
different point, it is re-evaluated (`callback_extra_evaluations = 34`).
**The mid-run point used for the adjoint's Tier 3 check is `k = 17`**, half
the iteration count rounded down.

The first iterate already improves on x₀ (10.2535 MWh/yr). Iterate 3 dips by
7.9e-4, and then the run climbs. Everything after `k = 6` (10.2623) is
refinement in the fourth decimal place and beyond, ending at 10.262736 at
`k = 35`.

## Why the gain is only 0.147 %

With λ fixed at 6.5 in every bin and an ideal power cap, AEP is effectively
Cp at λ = 6.5, and the Schmitz blade built from the same polar is the analytic
maximiser of exactly that. The 300 rpm ceiling breaks that match slightly,
because the rated point now runs at λ = 5.712, and that is why this run finds
+0.147 % and ends on a bound. An earlier version of the model let the rated
power rise with the design, which inflated the gain to 0.217 %; with a fixed
rating that part is gone. `docs/DESIGN-BASIS.md` §2 explains all of this.

The result isn't suspicious: it isn't negative, and it's nowhere near the
15 % that would suggest a bug. The multi-start study
(`verification/fd_optimisation_multistart/`) shows it is the global optimum
of this problem. A larger gain would need the machine model to include
something the real machine has, such as a lower rotor-speed ceiling or a real
above-rated mechanism (`docs/OUTSTANDING-INPUTS.md` §9).

## Reproducing it

    python verification/fd_optimisation/run_fd_slsqp.py            # ~6 min
    python verification/fd_optimisation/run_fd_slsqp.py --step 1e-5 --maxiter 50

`--step` overrides `h*`; by default it reads `h_star_global` from the
step-size study.

## Files

- `run_fd_slsqp.py`: the script. If a polar-range error occurs, it widens the
  envelope margin from 0.05 to 0.10 once and records it; that didn't happen.
- `result.json`: everything in the tables above, plus `u*`, `x*`, the
  gradient at the optimum in both units, the active set, the checks at x₀ and
  x\*, the SLSQP counters and options, and the bounds used.
- `iterates.json`: the 36 iterates.
- `optimised_blade.png`: chord and twist, x₀ vs the optimum. This figure
  predates the shared figure style and is kept as the run produced it.
