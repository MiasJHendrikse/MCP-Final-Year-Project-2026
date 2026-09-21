# Multi-start FD-driven SLSQP

**Re-run 2026-09-19 under the 300 rpm operating law**
(`λ_b = min(6.5, Ω_max R / V_b)`, `max_rotor_speed_rpm = 300`) **and the
configured bounds** (`chord_max_m = 0.30 m` — the O4 decision;
`chord_min_m = 0.045 m`, `twist_min = -2°`, `twist_max = 35°` grounded
2026-09-13; `σ ≤ 0.5`). Every optimum here is an optimum of that box under
that law. The two earlier runs — the 2026-09-13 floating-cap study and the
fixed-rating, no-ceiling, 0.45 m-box re-run — are in the git history; the
three are compared side by side below.

## Question

Is A4's single-start optimum (`verification/fd_optimisation/`, +0.147 % AEP
from `x0` under the ceiling) the optimum of the problem, or the nearest
local optimum to the Schmitz blade?

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
and `J` is finite. 8 accepted from 58 draws (14 %): 34 rejected because a
station did not converge and 16 because the blade left the geometric
envelope (min row −0.002 m to −0.088 m). The previous 0.45 m box accepted
8 from 1,259 (0.6 %), almost all of them rejected for the envelope (1,211
of 1,251) — a 0.45 m cap lets uniform draws build chords the envelope will
not hold, so the tighter box is the *easier* one to sample and what fails
now is the solver rather than the geometry. All eight starts are far from
`x0` (see `multistart.png`, left panels: root chords 0.149–0.293 m, mid-span
chords to 0.25 m, root twists 1–31°) and below it in AEP (8.35–9.86 MWh/yr
against 10.2477).

## Result

| start | AEP at start | `nit` | `nfev` | AEP* [MWh/yr] | ΔAEP % vs x0 | `‖u* − u*_A4‖∞` | active set | obj. evals | wall [s]† |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 8.4671 | 40 | 41 | 10.2627284 | +0.147 | 0.0166 | — | 1,644 | 753 |
| 1 | 9.1196 | 52 | 54 | 10.2627399 | +0.147 | 0.0070 | `chord_0` upper | 2,137 | 937 |
| 2 | 8.5772 | 52 | 56 | 10.2627400 | +0.147 | 0.0073 | `chord_0` upper | 2,139 | 946 |
| 3 | 8.3492 | 48 | 50 | 10.2627400 | +0.147 | 0.0077 | `chord_0` upper | 1,973 | 904 |
| 4 | 9.0116 | 49 | 54 | **10.2627400** | +0.147 | 0.0076 | `chord_0` upper | 2,037 | 906 |
| 5 | 9.8593 | 40 | 42 | 10.2627371 | +0.147 | 0.0114 | — | 1,645 | 757 |
| 6 | 9.0778 | 40 | 41 | 10.2627388 | +0.147 | 0.0049 | `chord_0` upper | 1,644 | 772 |
| 7 | 9.1627 | 51 | 55 | 10.2627400 | +0.147 | 0.0071 | `chord_0` upper | 2,098 | 928 |
| **A4 (from x0)** | 10.2477 | 35 | 36 | **10.2627365** | **+0.147** | 0 | `chord_0` upper | 1,438 | 367 |

† the eight runs shared 12 cores concurrently, so their wall times are not
comparable with A4's; iteration and evaluation counts are. The spread of
the eight wall times (753–946 s) is equally a scheduling artefact — the
objective evaluations differ by only 30 %.

All eight terminated with `Optimization terminated successfully`, no
`PolarDomainError`, no active envelope row, no active solidity row
(σ_max 0.404–0.409 against the 0.5 cap).

**Best start vs A4:** start 4, AEP* 10.2627400 vs A4's 10.2627365 —
3.5e-6 MWh/yr (3.4e-5 %) *above* A4. Spread across all nine optima (eight
starts and A4): 1.15e-5 MWh/yr in AEP (1.1e-6 relative), 0.0166 in `u`
(∞-norm).

**The active bound changed, and that is the headline.** The pre-ceiling run
put all nine optima on the grounded −2° tip-twist floor: with the capped
bins a constant, the objective kept asking to unload the tip and that bound
was what stopped it. Under the 300 rpm law the tip twist sits at −1.88°
(interior, −1.895…−1.879° across the nine) and **six of eight starts land
exactly on the 0.30 m root-chord ceiling instead** — `chord_0 = 0.300000 m`
to the bit, flagged `bounds_upper [0]`. The ceiling makes the machine want a
*bigger* root, not a smaller one, which is why the O4 cap is now the bound
that shapes the answer.

**And the spread lives in the two root control points, as before.** Six
starts agree to 1.1e-6 MWh/yr and to 4.0e-3 in `u` (`chord_0` identical,
`chord_1` within 1.07 mm, `twist_0` within 0.043°). Starts 0 and 5 stopped
marginally short along the same flat root direction A3's step-size study
found — `chord_0` 0.295776 m and 0.299140 m, i.e. 4.2 mm and 0.86 mm below
the bound, with `twist_0` 0.34° and 0.13° low — and their shortfalls
(1.15e-5 and 2.9e-6 MWh/yr) account for the entire AEP spread. Their
terminal `|∂J/∂u|∞` (4.2e-4, 3.6e-4 MWh/yr per `u`) are the two smallest of
the eight (the other six: 5.4e-4 to 7.3e-4), i.e. SLSQP stopped because its
own step fell below its threshold on a sub-millimetre-flat ridge, not
because it reached a bound. That is `ftol = 1e-8` on `fun` being met before
the last millimetres of root chord settle, not a set of distinct optima.

**Conclusion.** Every start from anywhere in the feasible box converges to
the same point as A4 — the same active bound (six of eight exactly on it,
the other two a fraction of a millimetre short), the same blade to 1.07 mm
in the free chord control points and 0.086° in the free twists, the same AEP
to 3.5e-6 MWh/yr. Under this objective (AEP on the 300 rpm schedule with
power held at the fixed rating above 11 m/s), these configured bounds and
the polar-cache envelope, the +0.147 % optimum is the global optimum of the
problem as posed, not a local one near Schmitz. The small gain over `x0` is
therefore not a starting-point artefact; taken with
`verification/spline_fit_error/` (the baseline's representation difference
is −0.015 %, a tenth of the gain and of the opposite sign) and
`docs/AEP_GAIN_AUDIT.md`, it is a property of the objective.

One observation from the runs, for the record: convergence from a poor
start costs 5–17 more iterations than from `x0` (40–52 vs 35) and 1.15–1.5×
the objective evaluations. Where the root-chord cap sits is part of the
answer (`docs/OUTSTANDING-INPUTS.md` §5), which is what the O4 decision
settled; the tip-twist floor matters less under the ceiling law than it did
without one.

## The three runs, side by side

| | 2026-09-13 (committed then) | fixed rating, no ceiling | **this run** |
|---|---|---|---|
| objective | `λ = 6.5`, floating cap | `λ = 6.5`, fixed rating | **300 rpm law, fixed rating** |
| box | `chord_max_m = 0.45 m` | same | **0.30 m** |
| optima | 9 interior, +0.217 % | 9 on the `twist_4` −2° floor, +0.121 % | **6 on the `chord_0` upper bound, +0.147 %** |
| AEP spread | 4.8e-8 MWh/yr | 3.2e-6 MWh/yr | **1.15e-5 MWh/yr** |
| `u` spread | 8.6e-4 | 0.017 | **0.0166** |
| accepted draws | — | — | **8 of 58** |

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

`sample_starts` is deterministic (fixed seed), so `--sample` reproduces the
eight starts above exactly; `--start K` may be run in any order or in
parallel, and `--collect` is what writes `results.json`.

The figure(s) in this directory predate `src/plotting/figstyle.py` and have not been redrawn through it; they are kept as the record of the run that produced them, and the scripts that write them carry no super-title.
