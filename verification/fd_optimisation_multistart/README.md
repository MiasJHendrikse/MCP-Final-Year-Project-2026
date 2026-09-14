# Multi-start FD-driven SLSQP

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. Every optimum here is an
optimum of that provisional box.

## Question

Is A4's single-start optimum (`verification/fd_optimisation/`, +0.217 % AEP
from `x0`) the optimum of the problem, or the nearest local optimum to the
Schmitz blade?

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
before the first step), every station converges at all 18 operating points,
and `J` is finite. 8 accepted from 1,259 draws: 1,211 rejected for the
envelope, 40 for a non-converged station. The envelope ceiling is what makes
the box mostly infeasible — at 19.5 m/s the tip chord must stay below
0.148 m (`u_4 ≲ 0.24`) and the 60 %-span control point below ~0.25 m, so
"sampled across the provisional bounds" means sampled uniformly across the
*feasible* part of the box. All eight starts are far from `x0` (see
`multistart.png`, left panels: root chords 0.16–0.32 m, twist humps at
mid-span) and far below it in AEP (7.35–9.56 MWh/yr against 10.27).

## Result

| start | AEP at start | `nit` | `nfev` | `njev` | AEP* [MWh/yr] | ΔAEP % vs x0 | `‖u* − u*_A4‖∞` | active set | wall [s]† |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 9.4369 | 41 | 45 | 41 | 10.2924816 | +0.217 | 2.0e-4 | none | 640 |
| 1 | 9.2615 | 37 | 40 | 37 | 10.2924816 | +0.217 | 3.3e-4 | none | 555 |
| 2 | 9.5608 | 40 | 43 | 40 | 10.2924816 | +0.217 | 3.1e-4 | none | 618 |
| 3 | 8.6106 | 43 | 45 | 43 | 10.2924816 | +0.217 | 5.3e-4 | none | 672 |
| 4 | 7.7758 | 39 | 41 | 39 | 10.2924815 | +0.217 | 7.6e-4 | none | 612 |
| 5 | 7.3484 | 45 | 49 | 45 | 10.2924816 | +0.217 | 1.8e-4 | none | 678 |
| 6 | 9.3993 | 38 | 40 | 38 | 10.2924816 | +0.217 | 3.3e-4 | none | 578 |
| 7 | 8.8915 | 40 | 44 | 40 | 10.2924816 | +0.217 | 9.9e-5 | none | 631 |
| **A4 (from x0)** | 10.2702 | 34 | 37 | 34 | **10.2924816** | **+0.217** | 0 | none | 361 |

† the eight runs shared 12 cores concurrently, so their wall times are not
comparable with A4's; iteration counts are.

All eight terminated with `Optimization terminated successfully`, no
`PolarDomainError`, no active bound, no active envelope row.

**Best start vs A4:** start 0, AEP* 10.292481588 vs A4's 10.292481589 —
1.8e-9 MWh/yr *below* A4. Spread across all nine optima (eight starts and
A4): 4.8e-8 MWh/yr in AEP, 8.6e-4 in `u` (∞-norm; 0.35 mm in root chord,
0.014° in root twist). That spread is SLSQP's `ftol` resolution on a
quadratic basin, not a set of distinct optima.

**Conclusion.** Every start from anywhere in the feasible box converges to
the same interior point as A4, to the solver's tolerance. Under this
objective (AEP at fixed `λ = 6.5` with simple power limiting), these
provisional bounds and the polar-cache envelope, the +0.217 % optimum is
the global optimum of the problem as posed, not a local one near Schmitz.
The small gain over `x0` is therefore not a starting-point artefact; taken
with `verification/spline_fit_error/` (the baseline's representation
difference is −0.031 %, a seventh of the gain and of the opposite sign), it
is a property of the objective. The 2–6 % expectation in
`PROJECT_DIRECTION_v2.md` §7.4 should be revisited against that.

Two observations from the runs, for the record. Convergence from a poor
start costs only 3–11 more iterations than from `x0` (37–45 vs 34): the
first six iterations do all the work (`multistart.png`, right panel) and the
remaining thirty are the same slow terminal phase A4 showed. And no run ever
touched a bound or the envelope, even starting near the ceiling: the
optimum is well interior to both.

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
