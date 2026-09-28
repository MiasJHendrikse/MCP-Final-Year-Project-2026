# Multi-start: is the energy optimum global?

The finite-difference-driven run in `verification/fd_optimisation/` finds
+0.147 % AEP starting from the Schmitz blade. The question here is whether
that is the optimum of the problem, or just the nearest local optimum to
Schmitz. To find out, the same optimisation is run from eight random starting
points.

**Answer: every start converges to the same optimum**, so +0.147 % is the
global optimum of the problem as posed.

The runs use the 300 rpm operating law (`λ_b = min(6.5, Ω_max R / V_b)`) and
the design bounds (chord 0.045–0.30 m, twist −2° to 35°, local solidity at
most 0.5).

## Method

    python verification/fd_optimisation_multistart/run_multistart.py --sample
    python verification/fd_optimisation_multistart/run_multistart.py --start K   # K = 0..7, can run in parallel
    python verification/fd_optimisation_multistart/run_multistart.py --collect

Everything is identical to the single-start run except the starting point:
`fun = J / |J(x0)|` (so `ftol = 1e-8` means the same thing), central
differences at `h* = 3e-6`, the 50-row polar-table envelope with margin 0.05,
`Bounds(0, 1)`, and SLSQP with `ftol = 1e-8`, `maxiter = 200`.

**Starting points.** Draws of `u ~ U[0, 1]^10` (seed 20260913), accepted only
where the objective can be evaluated: the envelope holds, every station
converges at all 17 operating points, and `J` is finite. Eight were accepted
from 58 draws (14 %). Of the rejected ones, 34 had a station that didn't
converge (explained below) and 16 put the blade outside the geometric
envelope (by 0.002 to 0.088 m). With the older 0.45 m chord bound, only 8 of
1,259 draws were accepted (0.6 %), almost all rejected by the envelope,
because a wider box lets random draws build chords the envelope can't hold.
All eight accepted starts are far from x₀ (see the left panels of
`multistart.png`: root chords 0.149–0.293 m, mid-span chords up to 0.25 m,
root twists 1–31°) and all have lower AEP (8.35–9.86 MWh/yr, against
10.2477).

## Result

| start | AEP at start | `nit` | `nfev` | AEP at optimum [MWh/yr] | change vs x₀ | `‖u* − u*_single‖∞` | active bounds | objective evaluations | time [s]† |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 8.4671 | 40 | 41 | 10.2627284 | +0.147 % | 0.0166 | – | 1,644 | 753 |
| 1 | 9.1196 | 52 | 54 | 10.2627399 | +0.147 % | 0.0070 | `chord_0` upper | 2,137 | 937 |
| 2 | 8.5772 | 52 | 56 | 10.2627400 | +0.147 % | 0.0073 | `chord_0` upper | 2,139 | 946 |
| 3 | 8.3492 | 48 | 50 | 10.2627400 | +0.147 % | 0.0077 | `chord_0` upper | 1,973 | 904 |
| 4 | 9.0116 | 49 | 54 | **10.2627400** | +0.147 % | 0.0076 | `chord_0` upper | 2,037 | 906 |
| 5 | 9.8593 | 40 | 42 | 10.2627371 | +0.147 % | 0.0114 | – | 1,645 | 757 |
| 6 | 9.0778 | 40 | 41 | 10.2627388 | +0.147 % | 0.0049 | `chord_0` upper | 1,644 | 772 |
| 7 | 9.1627 | 51 | 55 | 10.2627400 | +0.147 % | 0.0071 | `chord_0` upper | 2,098 | 928 |
| **single start (from x₀)** | 10.2477 | 35 | 36 | **10.2627365** | **+0.147 %** | 0 | `chord_0` upper | 1,438 | 367 |

† The eight runs shared 12 cores, so their wall times can't be compared with
the single-start run. The iteration and evaluation counts can.

All eight finished with `Optimization terminated successfully`, with no
polar-range errors, no active envelope row and no active solidity row
(largest solidity 0.404–0.409 against 0.5).

**Best start vs the single start.** Start 4 reaches 10.2627400 against
10.2627365, just 3.5e-6 MWh/yr (3.4e-5 %) higher. Across all nine optima the
spread is 1.15e-5 MWh/yr in AEP (1.1e-6 relative) and 0.0166 in `u`
(max norm).

**Six of the eight end exactly on the 0.30 m root-chord bound**
(`chord_0 = 0.300000 m`, flagged `bounds_upper [0]`). The tip twist ends
inside the box at about −1.88° (−1.895° to −1.879° across the nine). Under
the 300 rpm law the machine wants a *bigger* root, so the chord bound is what
shapes the answer.

**The spread is along the root.** Six starts agree to 1.1e-6 MWh/yr and to
4.0e-3 in `u` (`chord_0` identical, `chord_1` within 1.07 mm, `twist_0` within
0.043°). Starts 0 and 5 stopped just short along the same flat root direction
the step-size study found: `chord_0` at 0.295776 m and 0.299140 m (4.2 mm and
0.86 mm below the bound), with `twist_0` 0.34° and 0.13° low. Their shortfalls
(1.15e-5 and 2.9e-6 MWh/yr) account for all of the AEP spread. Their final
gradients (4.2e-4 and 3.6e-4 MWh/yr per `u`) are the two smallest of the eight
(the others range from 5.4e-4 to 7.3e-4). SLSQP stopped because its step fell
below its threshold on a very flat ridge, before the last millimetres of root
chord settled; these aren't separate optima.

**Conclusion.** Every start, from anywhere in the feasible box, converges to
the same point as the single-start run: the same active bound (six exactly on
it, two a fraction of a millimetre short), the same blade to within 1.07 mm
in the free chord control points and 0.086° in the free twists, and the same
AEP to within 3.5e-6 MWh/yr. So for this objective, these bounds and the
polar-table envelope, +0.147 % is the global optimum, not a local one near
Schmitz. The small gain over x₀ isn't caused by the starting point. Together
with `verification/spline_fit_error/` (the spline fit changes AEP by
−0.015 %, a tenth of the gain and of the opposite sign) and
`docs/DESIGN-BASIS.md` §2, it's a property of the objective.

Starting far from x₀ costs 5–17 more iterations than starting at x₀ (40–52
against 35) and 1.15–1.5 times as many objective evaluations.

## Compared with earlier configurations

| | floating rating | fixed rating, no ceiling | **current (300 rpm, fixed rating)** |
|---|---|---|---|
| objective | `λ = 6.5`, rating rising with the design | `λ = 6.5`, fixed rating | **300 rpm law, fixed rating** |
| chord bound | 0.45 m | 0.45 m | **0.30 m** |
| optima | 9 inside the box, +0.217 % | 9 on the −2° `twist_4` floor, +0.121 % | **6 on the `chord_0` upper bound, +0.147 %** |
| AEP spread | 4.8e-8 MWh/yr | 3.2e-6 MWh/yr | **1.15e-5 MWh/yr** |
| `u` spread | 8.6e-4 | 0.017 | **0.0166** |
| accepted draws | – | – | **8 of 58** |

## Why some random draws don't converge

`rejection_diagnosis.json`, from `diagnose_rejections.py`, re-solves the 34
draws rejected because a station didn't converge, at all seventeen operating
points, and records the solver's own failure message for every failed
station. There is exactly one failure mode:

- **All 3,120 failed station-points report "no sign change in the momentum
  region".** `bem.station.momentum_region_bracket` finds no windmill root.
  These stations are in the propeller-brake state, which the solver
  deliberately doesn't implement: such a station is reported rather than
  solved (see its docstring). There are no iteration-limit failures and no
  unreduced residuals.
- **Every failing station is twisted past its zero-induction inflow angle**,
  `atan(1 / lambda_r)`, by 1.1° to 25.1° (median 10.4°). The section meets the
  wind at a negative angle of attack before any induction, so it can't
  extract energy there. The failures are outboard of r/R = 0.37 and come from
  drawing the twist control points uniformly over the whole −2° to 35° range
  (rejected tip twists go up to 34°, while accepted ones are −1.3° to 13.1°).

So the solver is verified to converge in the windmill state
(`verification/phase_vi/`), and a blade twisted into the propeller-brake state
is outside that by construction. The random starts are drawn from the part of
the box where every station windmills, and the optimisation runs never visit
the other region (no run logged a failed evaluation).

    python verification/fd_optimisation_multistart/diagnose_rejections.py   # ~1 min

## Files

- `run_multistart.py`: the script (with the three modes above). It reuses the
  single-start run's `run()`, problem setup and `h*` unchanged.
- `starts.json`: the eight accepted starts (`u`, `x`, AEP, checks) and every
  rejected draw with its reason.
- `start_K.json`: one per start, with counters, the optimum, the distance to
  the single-start optimum, the active set, checks and every iterate.
- `results.json`: the table, the best-vs-single comparison and the spreads.
- `multistart.png`: chord and twist of the starts, the optima, x₀ and the
  single-start optimum, plus AEP against iteration for every start.
- `diagnose_rejections.py`, `rejection_diagnosis.json`: the rejected-draw
  analysis.

Sampling uses a fixed seed, so `--sample` reproduces the same eight starts.
`--start K` can run in any order or in parallel, and `--collect` writes
`results.json`.
