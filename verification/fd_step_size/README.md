# FD step-size study at x0 (Phase 2, Stage A3)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. The scaled variable
`u = (d - lo) / span` and therefore every gradient below is stated in those
bounds. Bounds come from `tests/test_parameterisation.py::PROVISIONAL_BOUNDS`
and are not in `config/`.

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
| `J(x0)` | -10.270157 MWh/yr (AEP 10.2702) |
| one objective evaluation | 0.207 s |
| evaluations | 300 (15 steps × 2 × 10) |
| wall time | 62.1 s |

## Selection rule

`h*_j` is the step minimising the local flatness
`|g_j(h) - g_j(h/√10)| + |g_j(h) - g_j(h√10)|`. The global `h*` is the grid
step nearest the median of `log10 h*_j`. The FD **noise floor** is
`ε_j = max(|g_j(h*_j) - g_j(h*_j/√10)|, |g_j(h*_j) - g_j(h*_j√10)|)`; Tier 3
(adjoint vs FD) is measured against it, never against a looser number.

## Result

Gradient in MWh/yr per unit `u` (i.e. `fun` units × |J0|); `ε_j` likewise.

| variable | `h*_j` | `g_j(h*_j)` | `ε_j` | `ε_j / |g_j|` |
|---|---|---|---|---|
| chord_0 | 1e-05 | -0.000145 | 2.3e-10 | 1.6e-06 |
| chord_1 | 3e-06 |  0.019356 | 9.9e-10 | 5.1e-08 |
| chord_2 | 3e-06 |  0.183099 | 2.2e-09 | 1.2e-08 |
| chord_3 | 1e-06 |  1.055235 | 5.8e-09 | 5.5e-09 |
| chord_4 | 3e-07 |  0.711802 | 1.4e-09 | 2.0e-09 |
| twist_0 | 3e-06 | -0.004424 | 1.5e-10 | 3.4e-08 |
| twist_1 | 1e-05 | -0.105533 | 3.9e-10 | 3.7e-09 |
| twist_2 | 1e-05 | -0.252586 | 1.1e-09 | 4.2e-09 |
| twist_3 | 1e-06 | -0.670088 | 1.2e-09 | 1.9e-09 |
| twist_4 | 3e-06 | -0.287430 | 1.7e-09 | 5.9e-09 |

**Global `h* = 3e-06`.** This is the step the FD-driven SLSQP run (A4) uses.

Reading the gradient: AEP rises with the outboard chord control points
(`chord_3`, `chord_4`, the 60–100 % span region) and falls with every twist
control point, most strongly `twist_3`. The root chord control point
(`chord_0`) is almost inert (-1.4e-4 MWh/yr per unit `u`), which is why its
relative noise floor is three orders worse than the rest: `ε` is absolute,
`g` is tiny.

## What the V-curve shows (`v_curve.png`)

`|g_j(h) - g_j(h*_j)|` vs `h`, log–log, per variable, with `h²` and `1/h`
guide lines through the geometric centre of the data.

The brief expected a *plateau* around `1e-4 – 1e-3` (J is C¹ not C²: the
Buhl blend at `a = 0.4`, and the polar interpolant's C² breaks at every
0.5° α-knot and every Reynolds row). **That is not what was measured.** The
curves are a clean V:

- from `h = 1e-2` down to about `1e-5` the deviation falls at the `h²`
  slope — the O(h²) truncation term of a central difference on a smooth
  function, four decades of it with no visible plateau;
- below about `1e-6` it rises as `1/h` — round-off, from the objective's own
  precision (`brentq` at `xtol = 1e-14` on every station);
- the minimum sits at `1e-6 – 1e-5` with relative error `~1e-9 – 1e-8`.

Why the C² breaks do not show at this point: their effect is O(h)·Δf″ *only
when the stencil straddles one*. At `x0` the stencil half-width in α is
`|dα/du_j| h`, of order a few degrees × `h`, and at `h ≤ 1e-3` it is well
below the 0.5° knot spacing, so almost no `(bin, station)` pair straddles a
knot for the steps that matter. The four stations above `a = 0.4` at the
rated point are likewise not within `h` of the blend. Tier 4 of the adjoint
plan is where this is counted explicitly; at `x0` the count is evidently
small. At `h = 1e-2` the deviation is `10⁻³ – 10⁻¹` relative — the largest
step is genuinely truncation-limited, not knot-limited.

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
  units; the flatness table; `x0`, `u0`, `J0`, timings, the provisional label.
- `v_curve.png` — the figure above.
