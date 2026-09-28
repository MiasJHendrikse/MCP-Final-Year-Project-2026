# Spline-fit error of the Schmitz reference blade

The reference blade x₀ isn't the analytic Schmitz blade itself. It is the
analytic blade fitted to the 5 + 5 control-point spline, because that's the
only kind of shape the optimiser can represent. This study measures how much
the fit changes the blade, and checks whether the optimiser's gain could just
be it undoing the fitting error.

**Answer: it can't.** Fitting costs 0.015 % of AEP, while the optimiser gains
0.147 %, nearly ten times more, and it gains that relative to the analytic
blade too.

The fit error depends on x₀ alone and uses no bounds. The AEP comparison uses
the finite-difference optimum found under the design bounds and the 300 rpm
operating law.

## Method

    python verification/spline_fit_error/run_fit_error.py          # ~5 s

x₀ (`verification/baseline/x0.json`) is rebuilt with
`design.build_schmitz_baseline()`, reproducing it exactly
(`x0_rebuild_max_abs_diff = 0`). The analytic Schmitz chord and twist it was
fitted from (design point Re = 200 k, α = 5.36°, Cl = 1.282, L/D = 98.5) are
then compared with the spline through x₀ (cubic, 5 + 5 clamped control points,
25 BEM strips).

## Fit error

At the 25 BEM stations, which are the only points the solver evaluates:

| | max abs | RMS | max relative |
|---|---|---|---|
| chord | **1.52 mm** (innermost station, r = 0.334 m) | **0.52 mm** | 0.58 % |
| twist | **0.265°** (innermost station) | **0.106°** | – |

On a 401-point grid over the whole span, hub to tip:

| | max abs | RMS |
|---|---|---|
| chord | 4.9 mm (at the hub, r = 0.30 m) | 0.65 mm |
| twist | 0.59° (at the hub) | 0.115° |

The largest errors on the fine grid are at the hub (r = 0.30 m), *inboard of
the first station* (0.334 m). The spline is fitted at the stations only, so it
extrapolates over the last 3 cm to the hub, where the Schmitz chord rises
steeply, but no BEM strip sits there. Between stations, the error is the
wiggle of a 5-point cubic against the Schmitz shape: about ±1 mm in chord and
±0.2° in twist, changing sign roughly every 0.4 m of span (bottom row of
`fit_error.png`). These match the RMS values recorded in `x0.json` when it was
built (0.51 mm, 0.106°).

## Does the fit error explain the optimiser's gain?

Three blades, evaluated with the same solver, resource and corrections:

| blade | AEP [MWh/yr] |
|---|---|
| analytic Schmitz (station values straight from the formulas) | 10.249287 |
| x₀, its spline fit | 10.247707 |
| finite-difference optimum | 10.262736 |

| difference | MWh/yr | % |
|---|---|---|
| fitting: x₀ − analytic | **−0.001580** | **−0.0154 %** |
| optimisation: optimum − x₀ | **+0.015029** | **+0.1467 %** |
| ratio | | **9.5×** |

**No, it doesn't.** Fitting costs the Schmitz blade 0.015 % of AEP, while the
optimiser gains 0.147 % over the fitted blade. Even if all the optimiser did
was undo the fitting error, it could only recover 0.015 %. The other 0.131 %
is an improvement over the *analytic* Schmitz blade, evaluated with the same
solver.

The shapes say the same thing. The optimiser moved the chord by up to 21.3 mm
(RMS 12.3 mm) and the twist by up to 2.02° (RMS 0.97°), against fitting errors
of 1.5 mm and 0.27°. So the optimum isn't the analytic blade recovered; it's a
different blade, with a thicker root on the 0.30 m bound, a thinner outboard
chord and the tip twisted down to −1.85° (see
`verification/fd_optimisation/README.md`).

This also shows the baseline is fair. The only difference between the
baseline and the optimum should be the optimisation itself. The fitting
difference is 0.015 %, a tenth of the effect being measured, and it works
*against* the baseline (the spline is slightly worse than the analytic curve),
so it can't have inflated the gain. The gain is small because of the
objective, not because of how the baseline was built (`docs/DESIGN-BASIS.md`
§2). `verification/fd_optimisation_multistart/` checks that it doesn't depend
on the starting point either.

Under earlier versions of the model the numbers differed slightly (a fitting
difference of −0.031 % when the rated power still rose with the design, and
−0.018 % against a +0.121 % gain once the rating was fixed), with the same
conclusion each time.

## Files

- `run_fit_error.py`: the script.
- `fit_error.json`: every number above, the signed errors at the stations and
  on the fine grid, the three chord and twist distributions, and the design
  point.
- `fit_error.png`: chord and twist (analytic, spline and optimum) and the
  signed fitting errors along the span. This figure predates the shared figure
  style and is kept as the run produced it.
