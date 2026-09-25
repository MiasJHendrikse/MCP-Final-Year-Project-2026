# FD step-size study at x0 (Phase 2, Stage A3)

**Under the configured bounds.** `chord_min_m = 0.045 m`, `chord_max_m = 0.30 m`,
`twist_min = -2°`, `twist_max = 35°`, all from `config/`. The scaled variable
`u = (d - lo) / span` and therefore every gradient below is stated in those
bounds. `chord_max_m = 0.30 m` was resolved by MJ on 2026-09-19 (B1/O4);
before that it was `0.45 m` in `tests/test_parameterisation.py::PROVISIONAL_BOUNDS`.

**Re-run 2026-09-19, second change of the day: the operating law and the
bounds.** `λ_b = min(6.5, Ω_max R / V_b)` with `Ω_max = 300 rpm` replaced the
fixed-`λ` law of the run below, and the box went from 0.45 m to 0.30 m of chord.
Every number in this file moved: `J(x0)` = -10.247707 MWh/yr and the gradient is
smaller and differently shaped than the fixed-`λ` result (outboard chord is no
longer the dominant lever, and `chord_3` changed sign). The **step-size
conclusions are unchanged** — the V-curve shape, the `10⁻⁸` relative accuracy at
`h*`, and the global `h* = 3.16e-06` all survive the change. The paragraphs
below the numbers are the fixed-rating run's, kept as history; the tables and
the reading of the gradient are the current run's.

**Re-run 2026-09-19 under the fixed generator rating.** The objective now
holds power above rated at `operating.rated_power_w` from
`config/rotor_design.yaml` (provisionally the baseline's own `P_aero(11 m/s;
x0) = 3822.19 W`, pending the nameplate -- outstanding input B2) instead of
at `P_aero(V_rated; d)`, which floated with the design
(`docs/AEP_GAIN_AUDIT.md` §3.2). `J(x0)` is unchanged to the bit; the
gradient is not: the capped bins no longer move with the blade, so every
component is smaller (roughly halved outboard) and two near-inert ones
(`chord_1`, `twist_0`) changed sign. The step-size conclusions are the same.

## What was run

    python verification/fd_step_size/run_sweep.py

At `u0 = bounds.to_scaled(x0)` (`verification/baseline/x0.json`), the central
difference of `fun(u) = J(u) / |J(u0)|`, `J = -AEP`, for all 10 scaled design
variables at 15 half-decade steps `h = logspace(-2, -9, 15)`. Every gradient
at every step is in `sweep.json`; nothing was discarded. No stencil left the
polar cache (`out_of_cache` is empty), so no step is recorded as
"out of cache".

| quantity | value |
|---|---|
| `J(x0)` | -10.247707 MWh/yr (AEP 10.247707) |
| one objective evaluation | 0.241 s (17 solves; the rated-speed solve is gone) |
| evaluations | 300 (15 steps × 2 × 10) |
| wall time | 79.6 s |

## Selection rule

`h*_j` is the step minimising the local flatness
`|g_j(h) - g_j(h/√10)| + |g_j(h) - g_j(h√10)|`. The global `h*` is the grid
step nearest the median of `log10 h*_j`. The FD **disagreement scale** is
`ε_j = max(|g_j(h*_j) - g_j(h*_j/√10)|, |g_j(h*_j) - g_j(h*_j√10)|)`; Tier 3
(adjoint vs FD) is measured against it, never against a looser number.

**`ε_j` is a flatness estimate, not a measured noise floor** — the name it had
until 2026-09-19 was wrong, and the 2026-09-19 re-run showed how wrong. It is
two differences between neighbouring grid steps, so it inherits the smoothness
assumption it is supposed to check: if both neighbours of `h*_j` sit on the same
interpolation kink, the flatness minimiser picks that step and `ε_j` collapses
instead of rising. `chord_4` is exactly that case here — `ε_4 = 5.4e-11` is the
smallest of the ten, while the median step-to-step jitter on `chord_4`'s own
plateau is `8.5e-10`, sixteen times larger. The adjoint agrees with the FD
gradient to `8.8e-10` there, which is O(1) against that jitter and 16.3× `ε_4`;
Tier 4 closes the argument by *measuring* the round-off floor of `J` at that
variable (`δJ/h*`, `δJ` from nine samples of `J` along `u + t e_4`,
`t = −4e-12 … 4e-12`): the floor is `4.00e-10`, i.e. **7.4 × larger than
`ε_4`**, and the adjoint–FD disagreement is `2.20 ×` that floor — the same
order as every other variable. So `chord_4`'s failure is a statement about
`ε_j`, not about the gradient. Tier 3 therefore reports one red variable at
`x0` (see
`gradient_verification/README.md` and the 2026-09-19 journal entry). The
estimator is unchanged from `master` — this was exposed by the re-run, not
caused by it — and MJ's decision on 2026-09-19 was to leave it and the test
alone and record the fragility rather than tune a threshold to hide it.

## Result

Gradient in MWh/yr per unit `u` (i.e. `fun` units × |J0|); `ε_j` likewise.

| variable | `h*_j` | `g_j(h*_j)` | `ε_j` | `ε_j / |g_j|` |
|---|---|---|---|---|
| chord_0 | 3e-06 | -0.011460 | 3.1e-10 | 2.7e-08 |
| chord_1 | 1e-05 | -0.123968 | 9.7e-11 | 7.8e-10 |
| chord_2 | 1e-05 |  -0.207398 | 1.0e-09 | 5.0e-09 |
| chord_3 | 3e-06 |  -0.011711 | 6.5e-10 | 5.6e-08 |
| chord_4 | 3e-06 |   0.127766 | 5.4e-11 | 4.2e-10 |
| twist_0 | 3e-05 |   0.004553 | 7.5e-11 | 1.7e-08 |
| twist_1 | 3e-06 |   0.009155 | 1.7e-10 | 1.8e-08 |
| twist_2 | 3e-06 |   0.019205 | 1.1e-09 | 5.8e-08 |
| twist_3 | 3e-06 |  -0.107662 | 1.1e-09 | 9.8e-09 |
| twist_4 | 3e-06 |  -0.039202 | 1.8e-09 | 4.6e-08 |

**Global `h* = 3.162e-06`.** This is the step the FD-driven SLSQP run (A4) uses.
Seven of the ten `h*_j` sit there too; `chord_1` and `chord_2` prefer `1e-05` and
`twist_0` prefers `3e-05`. All three of those are on the flat floor of their V,
which is the estimator's known weak spot (see *Selection rule* below).

Reading the gradient: under this law AEP **falls** with the inboard and mid-span
chord control points, most strongly `chord_2` (the 40–60 % span region), and
**rises** only with the tip control point `chord_4`; the outboard twist points
(`twist_3`, `twist_4`) push AEP down and are the only twist levers of any size.
Nine of the seventeen bins (**11.5–19.5 m/s**) are capped at the rating and
contribute nothing to the gradient; the 10.5 m/s bin sits on the ceiling but is
uncapped and does contribute. That is why the
largest lever is 0.207 MWh/yr per unit `u` against 0.558 before, and why
`chord_3` changed sign. For the record, the fixed-`λ` gradient at `x0` was
`chord_3 = +0.558`, `chord_4 = +0.386`, `twist_3 = -0.252` (0.558 MWh yr⁻¹ per
unit `u`, the fixed-`λ` form's largest lever).

## What the V-curve shows (`v_curve.png`)

`|g_j(h) - g_j(h*_j)|` vs `h`, log–log, per variable, with `h²` and `1/h`
guide lines through the geometric centre of the data.

The brief expected a *plateau* around `1e-4 – 1e-3` (J is C¹ not C²: the
Buhl blend at `a = 0.4`, and the polar interpolant's C² breaks at every
0.5° α-knot and every Reynolds row). **That is not what was measured.** The
curves are a clean V:

- from `h = 1e-2` down to about `1e-5` a decade, the deviation falls at the `h²`
  slope — the O(h²) truncation term of a central difference on a smooth
  function, several decades of it with no visible plateau;
- below that it turns and rises as `1/h` — round-off, from the objective's own
  precision (`brentq` at `xtol = 1e-14` on every station);
- the minimum sits at `1e-6 – 1e-5` (`3e-06` for seven of the ten variables,
  `1e-05` for two, `3e-05` for `twist_0`, whose gradient is the smallest),
  where `ε_j / |g_j|` is between `4e-10` and `6e-08`; at `h = 1e-8` the
  deviation is back up to `10⁻⁶ – 10⁻⁴`.

Why the C² breaks do not show at this point: their effect is O(h)·Δf″ *only
when the stencil straddles one*. At `x0` the stencil half-width in α is
`|dα/du_j| h`, of order a few degrees × `h`, and at `h ≤ 1e-3` it is well
below the 0.5° knot spacing, so almost no `(bin, station)` pair straddles a
knot for the steps that matter. The stations above `a = 0.4` in the
unlimited bins are likewise not within `h` of the blend. Tier 4 of the adjoint
plan is where this is counted explicitly; at `x0` the count is evidently
small. At `h = 1e-2` the deviation is `10⁻⁴ – 10⁻¹` relative — the largest
step is genuinely truncation-limited, not knot-limited. Under the 300 rpm law
the knots still do not show at `x0`, but they do at `u*`, where Tier 4 counts
a first crossing at `h = 1e-4` against `1e-3` here — the crossing step is a
property of the design point, not of the scheme.

Consequence for A4: `h* = 3e-6` gives a gradient with relative error
`≲ 1e-8` in every direction that matters, which is far below SLSQP's
`ftol = 1e-8` on an O(1) objective. Any `h` between `1e-6` and `1e-4` would
have served; the choice is not delicate.

## Re-plotting against the adjoint (later)

    python verification/fd_step_size/run_sweep.py --reference gradient.json

re-draws `v_curve.png` from `sweep.json` with `|g_j(h) - ref_j|` on the
y-axis, no re-run. `gradient.json` carries either `gradient_scaled` (fun
units per unit `u`) or `gradient_mwh_per_u` (converted with the sweep's
stored `J0`) and an optional `label`. `--replot` alone re-draws the committed
figure.

## Files

- `run_sweep.py` — the script.
- `sweep.json` — every gradient at every step; `h*_j`, `h*`, `ε_j` in both
  units; the flatness table; `x0`, `u0`, `J0`, timings, the bounds label.
- `v_curve.png` — the figure above.
