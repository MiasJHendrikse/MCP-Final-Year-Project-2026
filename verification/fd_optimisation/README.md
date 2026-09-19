# FD-driven SLSQP, end to end (Phase 2, Stage A4)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. The optimum below is an
optimum *of this provisional box*; it is not a design. Bounds come from
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` and are not in `config/`.

**Re-run 2026-09-19 under the fixed generator rating.** The objective now
holds power above rated at `operating.rated_power_w` from
`config/rotor_design.yaml` — provisionally the baseline's own `P_aero(11 m/s;
x0) = 3822.19 W`, pending the nameplate (outstanding input B2) — instead of
at `P_aero(V_rated; d)`, which floated with the design and credited the
optimiser with a bigger generator (`docs/AEP_GAIN_AUDIT.md` §1.2, §3.2).
The floating-cap run of 2026-09-13 reached **+0.217 %** at an interior
point in 34 iterations; that number is history. What follows is the
fixed-rating run. It reproduces the audit's own fixed-rating re-optimisation
(`verification/aep_gain_audit/opt_none_fixed.json`: +0.121 %, 27
iterations, tip-twist floor active) through `src/` proper.

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
| `nit` / `nfev` / `njev` | 27 / 29 / 27 |
| objective evaluations in total | 1111 (29 + 27 × 20 for the FD gradients + 26 callback re-evaluations × 20 + post-checks) |
| wall time | 230.3 s |
| AEP at `x0` | 10.270157 MWh/yr |
| AEP at optimum | 10.282564 MWh/yr |
| **ΔAEP vs x0** | **+0.121 %** |
| active bounds (`|u|` or `|1-u|` < 1e-6) | **`twist_4` at its lower bound, −2°** (`u = 0.0`) |
| active envelope rows (row value < 1e-6 m) | none; tightest row is the floor at `r = 1.966 m` with 14.9 mm of slack |
| α range at optimum (post-check) | 0.64° … 6.52°, inside the cache band [-8°, 18°] |
| Re range at optimum | 60.1 k … 474 k, inside [40 k, 1 M] — the §4.5 cache extension is not triggered |
| `PolarDomainError` during the run | none; margin stayed at 0.05 |

Optimum (physical, `x_star` in `result.json`):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x0 | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x* | 0.285 | 0.173 | 0.107 | 0.078 | 0.046 m | 23.50 | 10.36 | 3.52 | 2.01 | **−2.00** ° |

The optimiser thinned the tip (chord_4 to 0.046 m, 1 mm above the grounded
0.045 m floor, `u = 0.002`, not active) and **twisted it down onto the
grounded −2° floor**, which is active: with the capped bins a constant, the
tip is unloaded as far as the bound allows. The root and mid-span moves are
smaller than in the floating-cap run (root chord +9 mm against +41 mm).
`optimised_blade.png` overlays the two distributions with their control
points at the Greville abscissae.

Gradient at the optimum (MWh/yr per unit `u`): 7.7e-3 on the active
`twist_4` (pointing into the bound, as a KKT point requires — AEP would
rise if the tip could twist below −2°) and at most 1.9e-4 on the nine free
variables, down from 0.56 at `x0`. SLSQP stopped on `ftol`: the last
accepted steps changed `fun` by less than 1e-8, i.e. less than 1e-7 MWh/yr.

## Iterates (`iterates.json`)

28 records, `k = 0 … 27`, each with `u_k`, `J_k` (MWh/yr), `fun_k`, `g_k`
(both units) and wall time. `g_k` is the gradient at the accepted iterate;
where SLSQP's last gradient call was elsewhere it is re-evaluated
(`callback_extra_evaluations = 26`). **The mid-run point for the adjoint's
Tier 3 check is `k = 13` (`nit / 2`).**

The first iterate overshoots, as before: SLSQP's initial step takes AEP
*down* to 10.119 MWh/yr at `k = 1` and the inexact line search accepts it.
It passes `x0` at `k = 2` (10.273). Everything after `k = 6` (10.282) is
refinement in the fourth decimal and beyond.

## The gain is 0.12 %, and why

The 2026-09-13 run's README flagged its 0.22 % as below the 2–6 % the brief
expected; `docs/AEP_GAIN_AUDIT.md` then showed why: under fixed λ = 6.5 at
every bin and an ideal cap, AEP is `Cp(6.5)` up to Reynolds, and the
polar-consistent Schmitz `x0` is the analytic maximiser of that. The audit
also found that 0.096 of the 0.217 % was the *cap* rising with the design,
not the blade improving. That part is now gone by construction, and
**+0.121 % is the honest current-model figure**: not a defect (neither gate
— > 15 % or negative — is tripped), and the global optimum of the problem as
posed (`verification/fd_optimisation_multistart/`, re-run the same day).
A larger honest number needs the machine model to contain a feature the real
machine has (rotor-speed ceiling, nameplate, above-rated mechanism — inputs
B1–B3, `docs/OUTSTANDING-INPUTS.md` §9); the audit's §3 measured what each
is worth.

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
