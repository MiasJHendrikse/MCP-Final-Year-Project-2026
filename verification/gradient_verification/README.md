# Gradient verification: the adjoint against finite differences

The adjoint gradients are checked in four tiers. Tiers 1 and 2 are unit
tests:

- `tests/test_adjoint_partials.py` checks every partial derivative against a
  complex step of the code's own residual, to `1e-13 max(1, |partial|)`.
- `tests/test_adjoint_transpose.py` checks the transpose identity to `1e-14`.

This folder holds the other two:

- **Tier 3** compares the assembled gradient with central finite differences,
  judged against the finite-difference noise floor.
- **Tier 4** works out where any remaining disagreement comes from.

**Summary.** At three points and for all ten design variables, the adjoint
agrees with finite differences to the precision the arithmetic allows. 29 of
the 30 comparisons pass the acceptance test. The one that fails (`chord_4` at
x₀) is caused by the acceptance scale for that variable being too small, not
by the adjoint. I kept the test failing on purpose so the case stays visible.

All runs use the 300 rpm operating law, the fixed generator rating (the
adjoint system has 17 × 25 = 425 states, with the above-rated bins as a
constant; `docs/adjoint_derivation.md` §7) and the design bounds in
`config/rotor_design.yaml`. The finite-difference references come from the
step-size study (`verification/fd_step_size/sweep.json`, called A3 below) and
the finite-difference-driven optimisation run (`verification/fd_optimisation/`,
A4). The three points are x₀, the run's mid-point iterate k = 17 (of 35), and
its optimum u\*.

## Reproducing it

    python verification/gradient_verification/run_tier3.py                 # ~1 min
    python verification/gradient_verification/run_tier3.py --mid-iterate 17
    python verification/gradient_verification/run_tier4.py                 # ~23 s, needs tier3.json

`run_tier3.py` measures ε_j itself at k = 17 and u\*, since the step-size
study only covers x₀. At x₀ it uses the committed values and also re-derives
them, reproducing them to 5e-7 relative. `--mid-iterate` must match half the
iteration count in `fd_optimisation/result.json` (currently 17 of 35).

## Tier 3: method

For each point and variable:

| quantity | source |
|---|---|
| `adjoint_j` | `ScaledProblem.jac_adjoint(u)`, converted to MWh/yr per unit `u` |
| `FD_j` | central difference at the variable's own step `h*_j` (from `sweep.json`) |
| `eps_j` | the noise floor `max(|g(h*) − g(h*/√10)|, |g(h*) − g(h*√10)|)`; at x₀ the committed value, elsewhere re-measured with the same three-step local sweep, since the floor depends on the point |

A comparison passes if `|adjoint_j − FD_j| ≤ 3 eps_j`, with a separate check
for gross errors, `|adjoint_j − FD_j| / |FD_j| ≤ 1e-3`. There is also a
Taylor-remainder test of the whole chain at each point. For three random
unit directions `v`, the remainder `r(ε) = |fun(u + εv) − fun(u) − ε gᵀv|` for
`ε = 1e-2, 1e-3, 1e-4` must fall by at least 30 times per decade. A correct
gradient gives about 100, a wrong one about 10.

## Tier 3: result

| point | `J` [MWh/yr] | worst variable | worst `|diff|/eps_j` | worst relative error | Taylor ratio (minimum over 3 directions × 2 decades) |
|---|---|---|---|---|---|
| x₀ | −10.247707 | `chord_4` | **16.30, fails** | 4.8e-8 | 98.9 |
| iterate k = 17 | −10.262676 | `chord_0` | **0.67** | 5.0e-6 | 98.5 |
| optimum u\* | −10.262736 | `chord_2` | **2.03** | 1.7e-5 | 98.7 |

29 of the 30 pairs are within `3 eps_j`, and 28 within `1 eps_j`. The two above
`1 eps_j` are `chord_2` at u\* (2.03) and `chord_4` at x₀ (16.30). Every
relative error is between 1.8e-10 and 1.7e-5, at least three decades below
the 1e-3 gross-error threshold.

The Taylor test passes at all three points in every direction, with the
remainder falling 98.5–100.9 times per decade. So the objective really is
smooth along these directions with the gradient the adjoint reports, and the
`chord_4` disagreement comes from the finite-difference *reference*, not the
adjoint.

`tier3_agreement.png` shows `|diff|/eps_j` and the relative error for every
variable at the three points, with the acceptance lines.

### The one failure

At x₀, `chord_4` (the tip chord control point, `h* = 3e-6`):

| quantity | value |
|---|---|
| adjoint | 0.1277659593 MWh/yr per `u` |
| finite difference at `h*` | 0.1277659601 MWh/yr per `u` |
| `|diff|` | 8.82e-10 |
| `eps_j` from the step-size study (the acceptance scale) | 5.41e-11 |
| `|diff|/eps_j` | **16.30** (limit 3) |
| round-off floor of `J` measured by Tier 4 at this `h*` | 4.00e-10 |
| `|diff|` / that floor | 2.20 |

The disagreement is 2.2 times the *measured* round-off floor, which is in line
with every other variable and far too small to affect the optimiser. The
problem is the scale it's divided by. The step-size study's `eps_j` is a
flatness estimate, not a measurement of round-off. For this variable it lands
on a degenerate pair of neighbouring steps and comes out 7.4 times too small.
It is the smallest of the ten `eps_j` at x₀ (5.4e-11, against a median of
6.5e-10), and `chord_4`'s own step-to-step jitter across the plateau is about
8.5e-10, roughly the size of the disagreement itself.

Under an earlier configuration this test passed, only because `chord_4`'s
`h*` happened to fall where `eps_j` was inflated. So the problem wasn't
introduced by the later re-run; the re-run exposed it. I left both the
estimator and the test's tolerance unchanged. The right fix is to measure
round-off directly instead of inferring it from flatness, and loosening the
tolerance would simply delete the evidence. As a result,
`tests/test_adjoint_gradient.py::test_adjoint_agrees_with_the_committed_fd_reference_at_x0`
fails, by design.

The full comparison at x₀ (MWh/yr per unit `u`; the other two points are in
`tier3.json`):

| | adjoint | FD(`h*_j`) | `|diff|` | `h*_j` | `eps_j` | `|diff|/eps_j` |
|---|---|---|---|---|---|---|
| chord_0 | −0.0114596782 | −0.0114596779 | 2.8e-10 | 3e-6 | 3.1e-10 | 0.90 |
| chord_1 | −0.1239682220 | −0.1239682219 | 5.2e-11 | 1e-5 | 9.7e-11 | 0.54 |
| chord_2 | −0.2073983263 | −0.2073983261 | 1.2e-10 | 1e-5 | 1.0e-9 | 0.11 |
| chord_3 | −0.0117105459 | −0.0117105453 | 5.6e-10 | 3e-6 | 6.5e-10 | 0.86 |
| chord_4 | 0.1277659593 | 0.1277659601 | 8.8e-10 | 3e-6 | 5.4e-11 | **16.30** |
| twist_0 | 0.0045530398 | 0.0045530398 | 2.3e-11 | 3e-5 | 7.5e-11 | 0.30 |
| twist_1 | 0.0091551900 | 0.0091551901 | 7.2e-11 | 3e-6 | 1.7e-10 | 0.44 |
| twist_2 | 0.0192045338 | 0.0192045337 | 8.4e-11 | 3e-6 | 1.1e-9 | 0.08 |
| twist_3 | −0.1076622868 | −0.1076622871 | 2.8e-10 | 3e-6 | 1.1e-9 | 0.27 |
| twist_4 | −0.0392018262 | −0.0392018262 | 7.0e-12 | 3e-6 | 1.8e-9 | 0.00 |

The re-measured `eps_j` at x₀ matches the committed values to 5e-7 relative
(the code is deterministic), including the small one.

**Cost at x₀.** `J` takes 0.33 s, a finite-difference gradient (20
evaluations) 6.5 s, and an adjoint gradient (one forward solve, the partials
and assembly) **0.27 s**, or 0.8 times one evaluation of `J`. The whole Tier 3
run took 241 objective evaluations and 3 adjoint evaluations. (These wall
times were measured while the multi-start runs were also running, so the
evaluation counts are the better comparison; unloaded, one evaluation takes
about 0.25 s.)

### The V-curves against the adjoint

`v_curve_vs_adjoint.png` re-plots the step-size study at x₀ (no re-run) with
`|g_j(h) − adjoint_j|` on the y-axis, using `run_sweep.py --reference`
(`adjoint_x0.json`). With an independent reference, the curves become the
textbook V the self-referenced plot couldn't show: `h²` truncation error from
1e-2 down to 1e-5, a minimum of 7.0e-12 to 8.8e-10 (MWh/yr per `u`) at
`h = 3e-6` to `3e-5`, and `1/h` round-off below that. The adjoint sits at the
bottom of every variable's V. That is Tier 3 in one picture: the finite
difference is limited by its own noise floor, and the adjoint is inside it.
That includes `chord_4`, whose V minimum is 8.8e-10 at `h = 1e-5`, the same
size as its Tier 3 disagreement.

## Tier 4: where the remaining disagreement comes from

    python verification/gradient_verification/run_tier4.py        # ~21 s

Tier 4 asks what the central-difference stencil `[u − h, u + h]` crosses,
using the adjoint system's own linearisation to get each station's
sensitivities:

- **Angle-of-attack knots** (the interpolant's C² breaks at every 0.5° node):
  `|dα_{b,i}/du_j| h` against the distance to the nearest node, with
  `dα/du_j = (dφ/du_j − N_θ[i,j]) span_j` and `dφ/du_j = −(∂R/∂d_j)/(∂R/∂φ)`.
- **Reynolds-number rows** (the same C² break in the Re direction, PCHIP in
  log Re): `|dRe_{b,i}/du_j| h = (Re/c_i) N_c[i,j] span_j h` against the
  nearest row in the polar table.
- **The Buhl blend at `a = 0.4`**: `|da_{b,i}/du_j| h` against `|a − 0.4|`,
  with `da/du_j` from the kernel's induction partials.
- **The round-off floor of `J`**, measured rather than assumed: `J` is
  sampled along `u + t e_j` for `t = −4e-12 … 4e-12`, a line is fitted, and the
  residual standard deviation is `δJ`. The central-difference floor is `δJ/h`.

### At the chosen steps, nothing is crossed and the rest is round-off

| point | worst variable | knot / Re-row / Buhl crossings at `h*_j`, any variable | nearest break in stencil half-widths (α / Re / Buhl) | `δJ` [MWh/yr] | round-off floor `δJ/h*` | observed `|diff|` | `|diff|` / floor |
|---|---|---|---|---|---|---|---|
| x₀ | `chord_4` (`h*` 3e-6) | **0 / 0 / 0** | 4.7 / 26 / 131 | 1.27e-15 | 4.00e-10 | 8.82e-10 | **2.20** |
| iterate k = 17 | `chord_0` (`h*` 3e-6) | **0 / 0 / 0** | 3.3 / 21 / 212 | 1.34e-15 | 4.25e-10 | 7.81e-10 | **1.84** |
| optimum u\* | `chord_2` (`h*` 1e-5) | **0 / 0 / 0** | 2.6 / 168 / 172 | 1.00e-15 | 1.00e-10 | 9.36e-11 | **0.94** |

At every point and for all ten variables, the stencil at `h*_j` moves each
station's α by at most 2.4e-3°, with the nearest knot at least 2.6 stencil
half-widths away. It moves Re by at most 2.88, with the nearest row at least
21 half-widths away, and moves `a` by at most 4.0e-5 against
`|a − 0.4| ≥ 9.1e-4` (at least 131 half-widths). **No stencil crosses a knot,
a row or the Buhl blend at any `h*_j`.** So agreement at `h*` is pure
round-off: the three worst disagreements are **2.20, 1.84 and 0.94 times** the
measured floor, obtained by fitting nine samples of `J`.

The 2.20 at x₀ is `chord_4`, the variable behind the failing test. Its
measured floor is `δJ/h* = 4.00e-10` and the disagreement is 8.82e-10, or 2.2
floors. The step-size study's `eps_j` for this variable is 5.41e-11, **7.4
times smaller than the true floor**, which is why the ratio comes out at 16.3
instead of about 2. In short, the adjoint and the finite difference agree as
closely as the arithmetic allows; it's the acceptance scale that is wrong for
this variable.

### The whole V-curve explained

`tier4_attribution.png` (x₀ and u\*), `tier4_attribution_iterate.png`
(iterate k = 17).

For the worst variable at each point, these plot `|FD(h) − adjoint|` over the
15 steps of the study, alongside the crossing counts and the floor. A smooth
`C h²` law is fitted on the crossing-free steps that sit well above the floor
(at least 30 times it, so the fit isn't fitted to noise), and subtracting it
gives the departure at every step:

| point | worst variable | slope on the truncation side (h ≥ 1e-4) | `C h²` fitted on | crossing-free down to | smallest step with a crossing | largest `|excess|` over `C h²` at a crossing step |
|---|---|---|---|---|---|---|
| x₀ | `chord_4` | 1.91 | `[3e-4, 1e-4, 3e-5]`, `C = 9.43` | `h ≤ 3e-4` | `1e-3` (4 knots) | 3.4e-4 at `h = 1e-2` (37 knots + 4 Re + 3 Buhl) |
| k = 17 | `chord_0` | 1.91 | `[1e-4]`, `C = 0.049` | `h ≤ 1e-4` | `3e-4` (1 knot) | 1.9e-6 at `h = 1e-2` (27 knots + 2 Re) |
| u\* | `chord_2` | 1.58 | **none fitted** (no step rises clear of its floor) | `h ≤ 3e-5` | `1e-4` (2 knots) | not defined, see below |

This shows three things:

- **The C² breaks are where theory says they should be.** Every stencil that
  crosses a knot or row has `h ≥ 1e-4`, two to four decades above any
  `h*_j`, and the crossing-free steps at `h ≤ 3e-5` follow pure `h²`
  truncation.
- **The departures from `C h²` come mostly from the higher-order tail of the
  smooth function, not the crossings.** At x₀ the excess falls from 3.4e-4 at
  `h = 1e-2` to 1.2e-5 at `3e-3` and 1.6e-8 at `1e-3`, a fall of 30–740 times
  per decade, faster than `h²`. At `h = 1e-3` the departure is inside the
  fit's own 0.3 % band (± 2.8e-8). The departures that are significant are at
  most 0.2 % of the error, about 1e-7 of the gradient there, and their sign
  isn't consistent (at x₀: −3.4e-4 at 1e-2, −1.2e-5 at 3e-3, +1.2e-9 at
  3.16e-4). A crossing term of the form `O(h)·Δf″` would be positive, of one
  sign, and wouldn't shrink faster than `h²`. Nothing here looks like that.
- **At u\*, the fit itself isn't possible.** No crossing-free step rises 30
  times above its floor (the cleanest, `h = 3e-5`, is at 1.5 times), so the
  `h²` regime can't be observed in double precision there. The JSON records
  `smooth_h2_fit_rule = "no law fitted: … the O(h^2) truncation regime is not
  observable at this point"` and stores `null` for the prediction and excess,
  rather than a made-up number.

So nothing is left unexplained at a Tier 3 point. At `h*_j` the disagreement
is round-off (0.94–2.20 times the measured floor), with **no** knot, Reynolds
row or Buhl crossings inside any stencil. Above 1e-4, the departures from
`h²` are a higher-order truncation tail. The Buhl blend is only crossed at
`h = 1e-2` (at x₀), nowhere near an `h*_j`.

## Files

- `run_tier3.py`: the Tier 3 script.
- `tier3.json`: for each point, the adjoint and finite-difference gradients
  in both units, the local three-step sweep, `h*_j`, `eps_j` (committed and
  re-measured), the errors, the pass flags, the Taylor tests with their
  directions and remainders, and the timings.
- `adjoint_x0.json`: the adjoint gradient at x₀, in the format
  `run_sweep.py --reference` reads.
- `v_curve_vs_adjoint.png`: the step-size study against the adjoint.
- `tier3_agreement.png`: the acceptance figure.
- `run_tier4.py`: the Tier 4 script (reads `tier3.json` and `sweep.json`).
- `tier4.json`: for each point, every variable's crossing counts and margins
  at `h*_j`, `δJ` and its samples, the round-off floor, and for the worst
  variable the full 15-step sweep with counts, the `C h²` law and the excess
  at every step. When no step clears 30 times its floor, the coefficient and
  everything derived from it are `null`.
- `tier4_attribution.png`, `tier4_attribution_iterate.png`: the explained
  V-curves.
