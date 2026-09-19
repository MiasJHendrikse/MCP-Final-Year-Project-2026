# FD-driven SLSQP, end to end (Phase 2, Stage A4)

**Re-run 2026-09-19 under the 300 rpm operating law and the project bounds in
`config/rotor_design.yaml`**: `chord_min_m = 0.045`, `chord_max_m = 0.30`,
`twist_min_deg = -2`, `twist_max_deg = 35`, and the local solidity cap
σ ≤ 0.5. The optimum below is an optimum of *that* box under *that* law; it
is not a design. The bounds moved out of
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` and into `config/` on
2026-09-19 (outstanding inputs B1 and O4, resolved); `PROVISIONAL_BOUNDS` is
retired. What changed and in which order every artefact was re-run is in
`verification/README.md`.

**Also re-run earlier the same day under the fixed generator rating.** The
objective holds power above rated at `operating.rated_power_w` from
`config/rotor_design.yaml` — provisionally the baseline's own `P_aero(11 m/s;
x0) = 3822.19 W`, pending the nameplate (outstanding input B2) — instead of
at `P_aero(V_rated; d)`, which floated with the design and credited the
optimiser with a bigger generator (`docs/AEP_GAIN_AUDIT.md` §1.2, §3.2).
Three runs are now in the git history: 2026-09-13 floating cap, +0.217 % in
34 iterations; 2026-09-19 fixed rating under the old law and the 0.45 m box,
+0.121 % in 27 iterations, tip-twist floor active (it reproduces the audit's
own `verification/aep_gain_audit/opt_none_fixed.json` through `src/` proper);
and this run. What follows is the current one.

## What was run

    python verification/fd_optimisation/run_fd_slsqp.py

    minimise   fun(u) = J(u) / |J(u0)|,   J = -AEP [MWh/yr]
    over       u in [0, 1]^10             (5 chord + 5 twist control points, scaled)
    subject to the polar-cache Reynolds envelope: 50 linear rows, margin 5 %
               the root/tip chord envelope and the local solidity cap
               (sigma_i <= 0.5), both from ScaledProblem.constraints()
    gradient   central finite differences, h* = 3e-6 (global, from
               verification/fd_step_size/sweep.json)
    solver     scipy.optimize.minimize, method="SLSQP", ftol=1e-8, maxiter=200

Starting point `x0` = `verification/baseline/x0.json` (fitted Schmitz). The
constraint set comes from `problem.constraints()`, which is the same call the
adjoint run makes, so no script can quietly run a different problem.

## Result

| quantity | value |
|---|---|
| termination | `Optimization terminated successfully` (exit mode 0) |
| `nit` / `nfev` / `njev` | 35 / 36 / 35 |
| objective evaluations in total | 1438 (36 + 35 × 20 for the FD gradients + 34 callback re-evaluations × 20 + post-checks) |
| wall time | 367.1 s |
| AEP at `x0` | 10.247707 MWh/yr |
| AEP at optimum | 10.262736 MWh/yr |
| **ΔAEP vs x0** | **+0.147 %** |
| active bounds (`|u|` or `|1-u|` < 1e-6) | **`chord_0` at its upper bound, 0.30 m** (`u = 1.0`) |
| active envelope rows (row value < 1e-6 m) | none; tightest row is the floor at `r = 1.966 m` with 17.7 mm of slack |
| active solidity rows | none; σ_max = 0.4089 against a cap of 0.5, tightest row `r = 0.334 m` with 0.091 of slack |
| α at optimum (post-check) | 0.58° … 28.39° over the 17 schedule bins; the 8 at or below the rating are 0.58° … 6.54°, inside the cache band [−8°, 18°]. The 10 above it are on the ceiling, hold a constant power and contribute nothing to the gradient |
| Re at optimum | 63.5 k … 338.9 k, inside [40 k, 1 M] — the §4.5 cache extension is not triggered |
| `PolarDomainError` during the run | none; margin stayed at 0.05 |

Optimum (physical, `x_star` in `result.json`):

| | c0 | c1 | c2 | c3 | c4 | θ0 | θ1 | θ2 | θ3 | θ4 |
|---|---|---|---|---|---|---|---|---|---|---|
| x0 | 0.276 | 0.182 | 0.097 | 0.079 | 0.067 m | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 ° |
| x* | **0.300** | 0.180 | 0.130 | 0.077 | 0.049 m | 24.35 | 10.42 | 5.60 | 2.10 | −1.85 ° |

The optimiser ran the root chord up onto the 0.30 m upper bound and moved
chord *inboard* (c2 0.097 → 0.130 m), thickening the blade near the root,
while leaving the tip twist just short of its −2° floor (unlike the previous
run, where the twist floor was the active one). With the ceiling dropping the
design point to λ = 5.712 the blade is no longer at Schmitz's own optimum, so
it has something to win and the answer is no longer interior.
`optimised_blade.png` overlays the two distributions with their control
points at the Greville abscissae.

Gradient at the optimum (MWh/yr per unit `u`): −7.0e-4 on the active
`chord_0` (negative, which is what an upper bound requires — AEP would rise
if the chord could pass 0.30 m) and at most 3.7e-4 on the nine free
variables, down from 0.21 at `x0`. SLSQP stopped on `ftol`: the last accepted
step changed `fun` by less than 1e-8, i.e. less than 1e-7 MWh/yr.

## Iterates (`iterates.json`)

36 records, `k = 0 … 35`, each with `u_k`, `J_k` (MWh/yr), `fun_k`, `g_k`
(both units) and wall time. `g_k` is the gradient at the accepted iterate;
where SLSQP's last gradient call was elsewhere it is re-evaluated
(`callback_extra_evaluations = 34`). **The mid-run point for the adjoint's
Tier 3 check is `k = 17` (`nit / 2`, rounded down).**

The first iterate no longer overshoots (the previous run's took AEP *down* to
10.119 MWh/yr and the inexact line search accepted it): `k = 1` is 10.2535,
already past `x0`. `k = 3` dips by 7.9e-4 and then the run climbs; everything
after `k = 6` (10.2623) is refinement in the fourth decimal and beyond,
reaching 10.262736 at `k = 35`. The last accepted step is numerically zero.

## The gain is 0.147 %, and why

`docs/AEP_GAIN_AUDIT.md` explained the old 0.12 %: under a fixed λ = 6.5 at
every bin and an ideal cap, AEP is `Cp(6.5)` up to Reynolds, and the
polar-consistent Schmitz `x0` is the analytic maximiser of that; the audit
also found that 0.096 of the original 0.217 % was the *cap* rising with the
design, not the blade improving, and that part is gone by construction. The
300 rpm ceiling removes the rest of the degeneracy — the design point is no
longer at the TSR Schmitz optimises for — which is why this run reaches
**+0.147 %** and pushes onto a bound.

Neither defect gate is tripped (> 15 % or negative), and the multi-start
study (`verification/fd_optimisation_multistart/`, re-run the same day) shows
this to be the global optimum of the problem as posed. A larger honest
number needs the machine model to contain a feature the real machine has
(nameplate, above-rated mechanism — inputs B1–B3, `docs/OUTSTANDING-INPUTS.md`
§9); the audit's §3 measured what each is worth.

## Reproduce

    python verification/fd_optimisation/run_fd_slsqp.py            # ~6 min
    python verification/fd_optimisation/run_fd_slsqp.py --step 1e-5 --maxiter 50

`--step` overrides `h*`; the default reads `h_star_global` from the
step-size study.

## Files

- `run_fd_slsqp.py` — the script (margin escalation 0.05 → 0.10 on a
  `PolarDomainError`, once, recorded; not triggered).
- `result.json` — everything in the table above, plus `u*`, `x*`, the
  gradient at the optimum in both units, the active set (bounds, envelope
  rows, solidity rows, `sigma_max`), the post-checks at `x0` and `x*`, the
  SLSQP counters and options, the bounds label.
- `iterates.json` — the 36 iterates.
- `optimised_blade.png` — chord and twist, `x0` vs optimum.
