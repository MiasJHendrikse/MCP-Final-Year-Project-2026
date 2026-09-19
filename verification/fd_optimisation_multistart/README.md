# Multi-start FD-driven SLSQP

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. Every optimum here is an
optimum of that provisional box.

## Question

Is A4's single-start optimum (`verification/fd_optimisation/`, +0.121 % AEP
from `x0`) the optimum of the problem, or the nearest local optimum to the
Schmitz blade?

**Re-run 2026-09-19 under the fixed generator rating** (see A4's README).
Same seed, same eight accepted starts (the acceptance test does not depend
on the cap), same solver settings; the 2026-09-13 floating-cap result —
nine interior optima at +0.217 % within 8.6e-4 of each other in `u` — is
in the git history.

## What was run

    python verification/fd_optimisation_multistart/run_multistart.py --sample
    python verification/fd_optimisation_multistart/run_multistart.py --start K   # K = 0..7, in parallel
    python verification/fd_optimisation_multistart/run_multistart.py --collect

Identical to A4 in every solver respect: `fun = J / |J(x0)|` (the same
normalisation, so `ftol = 1e-8` means the same thing), central FD at
`h* = 3e-6`, the 50-row polar-cache envelope constraint with margin 0.05,
`Bounds(0, 1)`, SLSQP with `ftol = 1e-8`, `maxiter = 200`. Only the starting
point differs.

**Starting points.** `u ~ U[0, 1]^10`, seed 20260913, accepted only where
the objective is evaluable: the envelope holds (else `CachedPolar` raises
before the first step), every station converges at all 17 operating points,
and `J` is finite. 8 accepted from 1,259 draws: 1,211 rejected for the
envelope, 40 for a non-converged station. The envelope ceiling is what makes
the box mostly infeasible — at 19.5 m/s the tip chord must stay below
0.148 m (`u_4 ≲ 0.24`) and the 60 %-span control point below ~0.25 m, so
"sampled across the provisional bounds" means sampled uniformly across the
*feasible* part of the box. All eight starts are far from `x0` (see
`multistart.png`, left panels: root chords 0.16–0.32 m, twist humps at
mid-span) and far below it in AEP (8.24–9.83 MWh/yr against 10.27; the fixed cap
raises a poor blade's AEP relative to the floating one, since it no longer
caps it at its own low `P_aero(11 m/s)`).

## Result

| start | AEP at start | `nit` | `nfev` | `njev` | AEP* [MWh/yr] | ΔAEP % vs x0 | `‖u* − u*_A4‖∞` | active set | wall [s]† |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 9.7364 | 41 | 45 | 41 | 10.2825655 | +0.121 | 1.1e-2 | `twist_4` lower | 598 |
| 1 | 9.6828 | 38 | 39 | 38 | 10.2825655 | +0.121 | 1.1e-2 | `twist_4` lower | 548 |
| 2 | 9.8261 | 42 | 44 | 42 | 10.2825655 | +0.121 | 1.1e-2 | `twist_4` lower | 610 |
| 3 | 9.2261 | 43 | 45 | 43 | 10.2825655 | +0.121 | 1.2e-2 | `twist_4` lower | 634 |
| 4 | 8.5675 | 42 | 44 | 42 | 10.2825655 | +0.121 | 1.1e-2 | `twist_4` lower | 622 |
| 5 | 8.2423 | 44 | 46 | 44 | 10.2825655 | +0.121 | 1.1e-2 | `twist_4` lower | 639 |
| 6 | 9.7194 | 40 | 42 | 40 | 10.2825655 | +0.121 | 1.1e-2 | `twist_4` lower | 574 |
| 7 | 9.4080 | 32 | 35 | 32 | 10.2825623 | +0.121 | 5.4e-3 | `twist_4` lower | 502 |
| **A4 (from x0)** | 10.2702 | 27 | 29 | 27 | **10.2825642** | **+0.121** | 0 | `twist_4` lower | 230 |

† the eight runs shared 12 cores concurrently, so their wall times are not
comparable with A4's; iteration counts are.

All eight terminated with `Optimization terminated successfully`, no
`PolarDomainError`, no active envelope row, and **every one on the same
active bound**: the tip twist control point at its grounded −2° floor.

**Best start vs A4:** start 4, AEP* 10.2825655 vs A4's 10.2825642 —
1.4e-6 MWh/yr (1.3e-5 %) *above* A4. Spread across all nine optima (eight
starts and A4): 3.2e-6 MWh/yr in AEP, 0.017 in `u` (∞-norm). That
`u`-spread is larger than the floating-cap study's 8.6e-4, and it is
entirely in the two root control points — `chord_0` (7.0 mm across the nine
optima, 0.283–0.290 m) and `twist_0` (0.44°, 23.36–23.80°) — with every
other variable agreeing to ≤ 5e-3 in `u`. Those are the two near-inert
directions the step-size study found (`|∂J/∂u| ≈ 1e-3` and 1e-4 MWh/yr per
`u` at `x0`): the seven starts that cluster at 10.2825655 are within
1.4e-3 of each other in `u`, while A4 and start 7 stopped 1.4e-6 and
3.2e-6 MWh/yr short along the root direction, where `ftol = 1e-8` on `fun`
is met before the last few millimetres of root chord are settled. That is
SLSQP's tolerance on a basin that is very flat in two of ten directions,
not a set of distinct optima: the AEP spread is 3e-7 relative.

**Conclusion.** Every start from anywhere in the feasible box converges to
the same point as A4 — the same active bound, the same blade to 7 mm at the
root and better than a millimetre elsewhere, the same AEP to 3e-6 MWh/yr.
Under this objective (AEP at fixed `λ = 6.5` with power held at the fixed
rating above 11 m/s), these provisional bounds and the polar-cache
envelope, the +0.121 % optimum is the global optimum of the problem as
posed, not a local one near Schmitz. The small gain over `x0` is therefore
not a starting-point artefact; taken with `verification/spline_fit_error/`
(the baseline's representation difference is −0.018 %, a seventh of the
gain and of the opposite sign) and `docs/AEP_GAIN_AUDIT.md`, it is a
property of the objective.

Two observations from the runs, for the record. Convergence from a poor
start costs 5–17 more iterations than from `x0` (32–44 vs 27): the first
six iterations do all the work (`multistart.png`, right panel) and the rest
is the same slow terminal phase A4 showed. And the tip-twist floor is the
one bound that binds, from every start: with the capped bins a constant,
unloading the tip is the direction the objective keeps asking for, and the
grounded −2° bound is what stops it. Where that bound sits is therefore
part of the answer, which is worth knowing before `chord_max_m` and the
root cut-out are settled (`docs/OUTSTANDING-INPUTS.md` §2, §5).

## Files

- `run_multistart.py` — the script (three modes above); reuses A4's
  `run()`, problem construction and `h*` unchanged.
- `starts.json` — the eight accepted starts (`u`, `x`, AEP, post-check) and
  every rejected draw with its reason.
- `start_K.json` — one per start: counters, optimum, distance to A4, active
  set, post-check, every iterate `(k, u_k, J_k, g_k, wall)`.
- `results.json` — the table, the best-vs-A4 comparison, the spreads.
- `multistart.png` — chord/twist of starts, optima, `x0` and A4; AEP vs
  iteration for every start.
