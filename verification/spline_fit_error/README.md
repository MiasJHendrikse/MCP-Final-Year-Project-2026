# Spline-fit error of the Schmitz baseline (PROJECT_DIRECTION_v2 §7.3 point 4)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision). The A4 gain quoted below is
under those bounds. This study itself uses no bound.

**Re-run 2026-09-19 under the fixed generator rating.** The fit error is a
property of `x0` and is unchanged. The AEP comparison is not: with the cap
floating, the analytic blade was credited with its *own* slightly higher
`P_aero(11 m/s)` as its rating, which was 0.013 % of the 0.031 %
representation difference recorded on 2026-09-13; under one rating for all
three blades the representation difference is −0.018 % and A4's gain is
+0.121 %. Same conclusion, smaller numbers.

## What was run

    python verification/spline_fit_error/run_fit_error.py          # ~5 s

`x0` (`verification/baseline/x0.json`) is rebuilt with
`design.build_schmitz_baseline()` — bit-identical (max |diff| = 0) — and the
analytic Schmitz chord and twist it was fitted from (design point
Re = 200 k, α = 5.36°, Cl = 1.282, L/D = 98.5) are compared with the spline
through `x0` (cubic, 5 + 5 clamped control points, 25 BEM strips).

## Fit error

At the 25 BEM stations (the only points the solver ever evaluates):

| | max abs | RMS | max relative |
|---|---|---|---|
| chord | **1.52 mm** (innermost station, r = 0.334 m) | **0.52 mm** | 0.58 % |
| twist | **0.265°** (innermost station) | **0.106°** | — |

On a 401-point grid over the full span, hub to tip:

| | max abs | RMS |
|---|---|---|
| chord | 4.9 mm (at the hub, r = 0.30 m) | 0.65 mm |
| twist | 0.59° (at the hub) | 0.115° |

The fine-grid extremes are at `r = r_hub = 0.30 m`, *inboard of the first
station* (0.334 m): the spline is fitted at the stations only and is
extrapolating over the last 3 cm to the hub, where the Schmitz chord is
rising steeply. No BEM strip sits there. Between the stations the error is
the oscillation of a 5-point cubic against `sin²(φ/2)·r`: ±1 mm in chord and
±0.2° in twist, both changing sign every ~0.4 m of span (`fit_error.png`,
bottom row). These agree with the RMS figures recorded in `x0.json` at
construction (0.51 mm, 0.106°).

## Does the fit error explain A4's gain?

Three blades, one solver, one resource, one set of corrections:

| blade | AEP [MWh/yr] |
|---|---|
| analytic Schmitz (station values straight from the formulas) | 10.271999 |
| `x0` = its spline projection | 10.270157 |
| A4 FD-SLSQP optimum | 10.282564 |

| difference | MWh/yr | % |
|---|---|---|
| representation: `x0` − analytic | **−0.001842** | **−0.018 %** |
| optimisation: optimum − `x0` | **+0.012407** | **+0.121 %** |
| ratio | | **6.7×** |

**No — the fit error does not explain the gain.** The projection onto the
spline costs the Schmitz blade 0.018 % of AEP; the optimiser gained 0.121 %
over the projected blade, nearly seven times more. Even if the optimiser had
done nothing but undo the projection error it could only have recovered
0.018 %; the remaining 0.103 % is above the *analytic* Schmitz blade
evaluated by the same solver. That is also visible in the shapes: the
optimiser moved the chord by up to 19 mm (RMS 6 mm) and the twist by up to
2.2° (RMS 0.62°), against fit errors of 1.5 mm and 0.27° — the optimum is
not the analytic blade found again, it is a different blade (thinner tip
twisted down onto the −2° bound; see
`verification/fd_optimisation/README.md`).

What the numbers *do* say about the gain being small (0.12 % against the
2–6 % expected in §7.4): the baseline is a fair one. §7.3's requirement was
that the only difference between baseline and optimum be the optimisation
itself; the representation difference is 0.02 %, a seventh of the effect
being measured, and it works *against* the baseline (the spline is slightly
worse than the analytic curve), so it cannot have inflated the gain. The
small gain is a property of the objective and starting point, not of the
baseline's construction. `verification/fd_optimisation_multistart/` tests
whether it is also a property of the starting point.

## Files

- `run_fit_error.py` — the script.
- `fit_error.json` — every number above, the per-station and fine-grid
  signed errors, the three chord/twist distributions, the design point.
- `fit_error.png` — chord and twist (analytic, spline, optimum) and the
  signed fit errors vs span.
