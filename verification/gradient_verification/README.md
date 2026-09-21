# Gradient verification: the discrete adjoint against finite differences (Phase 3, Tiers 3 and 4)

**Re-run 2026-09-19 under the 300 rpm operating law and the project bounds in
`config/rotor_design.yaml`** (`chord_max_m = 0.30 m`, resolved 2026-09-19;
`chord_min_m`, `twist_min_deg`, `twist_max_deg` grounded 2026-09-13). The
scaling `u = (d − lo)/span` and therefore every gradient here is stated under
*these* bounds; `tests/test_parameterisation.py::PROVISIONAL_BOUNDS` is
retired, so there is now one box and it lives in `config/`. The ordering of
the re-runs is in `verification/README.md`.

**Also re-run earlier the same day under the fixed generator rating.** The
objective holds power above rated at `operating.rated_power_w` from
`config/rotor_design.yaml` (provisionally the baseline's own `P_aero(11 m/s;
x0)`, pending input B2) instead of at `P_aero(V_rated; d)`; the adjoint
system lost its eighteenth operating point (the rated solve) and is 17 × 25 =
425 states, with the capped bins a constant in `J`
(`docs/adjoint_derivation.md` §7, revision note). Tiers 1 and 2 pass
unchanged on the 425-state system.

Tiers 3 and 4 below are measured against the current A3 sweep
(`verification/fd_step_size/sweep.json`), the mid-run iterate
`k = 17 = nit/2` of the current 35-iteration A4 run, and that run's optimum
`u*` (`verification/fd_optimisation/result.json`). The 2026-09-13
floating-cap numbers (30/30 within 3 ε, worst 2.78 ε, zero crossings at
`h*`, `k = 17` of a 34-iteration run) and the fixed-rating numbers
(`k = 13` of 27) are in the git history and in
`docs/adjoint_derivation.md` §9.

**One measured value fails its acceptance test, and it is reported as found:
`x0`, `chord_4`, `|diff|/eps_j = 16.30` against the limit of 3.** It does not
fail because the adjoint is wrong — Tier 4 measures the round-off floor of
`J` at that variable and finds the disagreement is 2.2 × that floor, while
A3's `eps_j` for the same variable is 7.4 × *below* it. The acceptance scale
is the fragile part; the details are under "One failure, and what it is"
below and in `verification/fd_step_size/README.md`. MJ's decision of
2026-09-19 was to leave both the estimator and the test unedited and record
this as a known fragility, so `pytest -q` stays red on
`tests/test_adjoint_gradient.py::test_adjoint_agrees_with_the_committed_fd_reference_at_x0`
until Phase 4 re-derives the floor.

Tiers 1 and 2 are tests, not artefacts: `tests/test_adjoint_partials.py`
(every partial against a complex step of the code's own residual, to
`1e-13 max(1, |partial|)`) and `tests/test_adjoint_transpose.py` (the
transpose identity to `1e-14`). This directory holds Tier 3 (the assembled
gradient against the FD noise floor) and Tier 4 (attribution of whatever
disagreement is left).

## Tier 3 — what was run

    python verification/gradient_verification/run_tier3.py        # ~1 min

At three points — `x0`, the FD-driven SLSQP run's mid-run iterate `k = 17`
(`verification/fd_optimisation/iterates.json`, `nit/2`), and that run's
optimum `u*` (`result.json`) — for all 10 scaled design variables:

| quantity | source |
|---|---|
| `adjoint_j` | `ScaledProblem.jac_adjoint(u)`, undone to MWh/yr per unit `u` |
| `FD_j` | central difference at `h*_j` (per variable, from `verification/fd_step_size/sweep.json`) |
| `eps_j` | the FD noise floor `max(|g(h*) − g(h*/√10)|, |g(h*) − g(h*√10)|)`: at `x0` the committed A3 value; at the other two points re-measured by the same three-step local sweep, since the floor is a property of the point |

Acceptance, derived not invented: `|adjoint_j − FD_j| ≤ 3 eps_j` and the
gross-error tripwire `|adjoint_j − FD_j| / |FD_j| ≤ 1e-3`. Plus the
whole-chain Taylor-remainder test at each point: three random unit
directions `v`, `r(ε) = |fun(u + εv) − fun(u) − ε gᵀv|` for
`ε = 1e-2, 1e-3, 1e-4`, ratio per decade asserted `≥ 30`.

## Tier 3 — result: two points pass, one variable at `x0` does not

| point | `J` [MWh/yr] | worst variable | worst `|diff|/eps_j` | worst rel. error | Taylor ratios (min over 3 draws × 2 decades) |
|---|---|---|---|---|---|
| `x0` | −10.247707 | `chord_4` | **16.30 — FAILS** | 4.8e-8 | 98.9 |
| iterate `k = 17` | −10.262676 | `chord_0` | **0.67** | 5.0e-6 | 98.5 |
| optimum `u*` | −10.262736 | `chord_2` | **2.03** | 1.7e-5 | 98.7 |

29 of the 30 (point, variable) pairs are within `3 eps_j`, and 28 within
`1 eps_j`. The two above `1 eps_j` are `chord_2` at `u*` (2.03) and
`chord_4` at `x0` (16.30). Every relative error is 1.8e-10 … 1.7e-5, three to
nine decades below the `1e-3` tripwire.

The Taylor-remainder test passes at all three points and all three draws:
the remainder falls by 98.5–100.9 × per decade, against an assertion of
≥ 30 × (a wrong gradient gives ≈ 10 ×). So `fun` really is smooth along
these directions with the gradient the adjoint reports; the `chord_4`
disagreement is a property of the FD *reference*, not of the adjoint.

## One failure, and what it is

At `x0`, `chord_4` (the tip chord control point, `h* = 3e-6`):

| quantity | value |
|---|---|
| adjoint | 0.1277659593 MWh/yr per `u` |
| FD at `h*` | 0.1277659601 MWh/yr per `u` |
| `|diff|` | 8.82e-10 |
| A3 `eps_j` (the acceptance scale) | 5.41e-11 |
| `|diff|/eps_j` | **16.30** (limit 3) |
| round-off floor of `J` measured by Tier 4 at this `h*` | 4.00e-10 |
| `|diff|` / that floor | 2.20 |

The disagreement is 2.2 × the *measured* round-off floor — the same order as
every other variable, and two decades below tolerances that would matter to
the optimiser. What is wrong is the scale it is divided by: A3's `eps_j` is a
*flatness* estimate (`max(|g(h*) − g(h*/√10)|, |g(h*) − g(h*√10)|)`), not a
measurement of round-off, and for this particular variable it lands on a
degenerate neighbour pair and comes out 7.4 × too small. It is the smallest
of the ten `eps_j` at `x0` (5.4e-11, against a median of 6.5e-10), and
`chord_4`'s own step-to-step jitter across the plateau is ≈ 8.5e-10, i.e.
about the `|diff|` itself.

Under the fixed-rating law the same test passed only because `chord_4`'s
`h*` happened to fall where `eps_j` was inflated; the failure is not
introduced by this re-run, it is exposed by it. Both the estimator and the
test are deliberately untouched (MJ, 2026-09-19): the estimator's fix is
"actually measure round-off instead of inferring it from flatness", which is
Phase 4 work, and editing the test's tolerance now would delete the evidence.
Recorded in `docs/journal/Session Notes/2026-09-19.md`.

At `x0`, per variable (MWh/yr per unit `u`; the full table for the other
two points is in `tier3.json`):

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
| twist_3 | −0.1076622868 | −0.1076622871 | 2.9e-10 | 3e-6 | 1.1e-9 | 0.27 |
| twist_4 | −0.0392018262 | −0.0392018262 | 7.0e-12 | 3e-6 | 1.8e-9 | 0.00 |

The re-measured `eps_j` at `x0` reproduces A3's committed values to 5e-7
relative (same evaluations, deterministic code) — including the small one, so
the number above is A3's, not a fresh draw.

**Wall time at `x0`:** `J` 0.33 s; FD gradient (20 evaluations) 6.5 s;
adjoint gradient (one forward solve + partials + assembly) **0.27 s** —
0.8 × `J`, against the brief's "≲ 2 × J". The whole Tier 3 run cost 241
objective evaluations and 3 adjoint evaluations. (Wall times in this file
were taken on the working desktop while the eight multi-start runs were
re-running; the evaluation counts are the comparable figures. The unloaded
per-evaluation cost is ≈ 0.25 s.)

## Tier 3 — the V-curves re-plotted against the adjoint

`v_curve_vs_adjoint.png` is A3's step-size sweep at `x0` (`sweep.json`,
no re-run) re-plotted with `|g_j(h) − adjoint_j|` on the y-axis, via
`run_sweep.py`'s `--reference` machinery (`adjoint_x0.json`); A3's own
figure is untouched. With an independent reference the curves become the
textbook V that A3's self-referenced plot could not show: `h²` truncation
from `1e-2` down to `1e-5`, a minimum of 7.0e-12 … 8.8e-10 (MWh/yr per `u`)
at `h = 3e-6 … 3e-5`, and `1/h` round-off below. The adjoint sits at the
bottom of every variable's V, which is the statement Tier 3 makes in one
picture: the FD reference is limited by its own noise floor, and the
adjoint is inside it — `chord_4` included, whose V minimum is 8.8e-10 at
`h = 1e-5`, i.e. the same magnitude as the Tier 3 disagreement it is
measured with.

`tier3_agreement.png` shows `|diff|/eps_j` and the relative error per
variable at the three points, with the acceptance lines.

## Tier 4 — degradation at the polar interpolation

    python verification/gradient_verification/run_tier4.py        # ~21 s

Tier 3 left 30 disagreements, 29 of them inside `3 eps_j` and two above
`1 eps_j`. Tier 4 explains all of them by asking what the central
stencil `[u − h, u + h]` crosses, using the adjoint system's own
linearisation for the per-station sensitivities:

- **α knots** (the interpolant's C² breaks at every 0.5° node):
  `|dα_{b,i}/du_j| h` against the distance to the nearest node, with
  `dα/du_j = (dφ/du_j − N_θ[i,j]) span_j` and `dφ/du_j = −(∂R/∂d_j)/(∂R/∂φ)`.
- **Reynolds rows** (the same C² break in the Re direction, PCHIP in log Re):
  `|dRe_{b,i}/du_j| h = (Re/c_i) N_c[i,j] span_j h` against the nearest cached row.
- **Buhl `a = 0.4`**, reported on its own line: `|da_{b,i}/du_j| h` against
  `|a − 0.4|`, with `da/du_j` from the kernel's induction partials.
- **The round-off floor of `J`**, measured, not assumed: `J` sampled along
  `u + t e_j`, `t = −4e-12 … 4e-12`, linear fit, residual std = `δJ`; the
  central-difference floor is `δJ/h`.

### Result at `h*_j`: no stencil crosses anything at any of the three points; the rest is round-off

| point | worst variable | knot / Re-row / Buhl crossings at `h*_j`, any variable | nearest break in stencil half-widths, any variable (α / Re / Buhl) | `δJ` [MWh/yr] | round-off floor `δJ/h*` | observed `|diff|` | `|diff|` / floor |
|---|---|---|---|---|---|---|---|
| `x0` | `chord_4` (`h*` 3e-6) | **0 / 0 / 0** | 4.7 / 26 / 131 | 1.27e-15 | 4.00e-10 | 8.82e-10 | **2.20** |
| iterate `k = 17` | `chord_0` (`h*` 3e-6) | **0 / 0 / 0** | 3.3 / 21 / 212 | 1.34e-15 | 4.25e-10 | 7.81e-10 | **1.84** |
| optimum `u*` | `chord_2` (`h*` 1e-5) | **0 / 0 / 0** | 2.6 / 168 / 172 | 1.00e-15 | 1.00e-10 | 9.36e-11 | **0.94** |

At every point and for all ten variables, the stencil at `h*_j` moves every
station's α by at most 2.4e-3° against knots a minimum of 2.6 stencil
half-widths away, moves Re by at most 2.88 against rows a minimum of 21 away,
and moves `a` by at most 4.0e-5 against `|a − 0.4| ≥ 9.1e-4` (a minimum of
131 half-widths). **No stencil crosses a knot, a row or the Buhl blend at any
`h*_j`** (the previously reported single `twist_4` knot crossing at `k = 13`
does not survive the re-run — the iterates moved and the stencil no longer
straddles a node). Agreement at `h*` is therefore pure round-off: the three
worst-variable disagreements are **2.20 ×, 1.84 × and 0.94 ×** the *measured*
`δJ/h*`, where `δJ` is not assumed but obtained by sampling `J` along
`u + t e_j` for `t = −4e-12 … 4e-12` and fitting a line (nine samples).

The 2.20 at `x0` is `chord_4`, and it is the variable behind the one red test.
Tier 4 measures its round-off floor at `δJ/h* = 4.00e-10` MWh/yr per `u`; the
adjoint–FD disagreement is 8.82e-10, i.e. 2.2 floors. A3's `eps_j` for the
same variable is 5.41e-11 — **7.4 × smaller than the true floor** — which is
why the ratio comes out at 16.3 instead of ≈ 2. This is the quantitative
statement of the failure: the adjoint and the FD reference agree to the
precision the arithmetic allows; the acceptance scale is wrong at this
variable. See "One failure, and what it is" above.

### The whole V attributed (`tier4_attribution.png`, `tier4_attribution_iterate.png`)

For the worst variable at each point, `|FD(h) − adjoint|` over A3's 15-step
grid next to the crossing counts and the floor. `tier4_attribution.png` carries
the reference and the optimum, and `tier4_attribution_iterate.png` the mid-run
iterate on its own (the report places it in Appendix D). Fitting the smooth law
`C h²` on the crossing-free steps that sit clear of the floor (at least
30 × the floor, so the fit is not itself fitted to noise) and subtracting it
gives the departure at every step:

| point | worst variable | truncation-side slope (h ≥ 1e-4) | `C h²` fitted on | crossing-free down to | smallest step with any crossing | largest `|excess|` over `C h²` at a crossing step |
|---|---|---|---|---|---|---|
| `x0` | `chord_4` | 1.91 | `[3e-4, 1e-4, 3e-5]`, `C = 9.43` | `h ≤ 3e-4` | `1e-3` (4 knots) | 3.4e-4 at `h = 1e-2` (37 knots + 4 Re + 3 Buhl) |
| `k = 17` | `chord_0` | 1.91 | `[1e-4]`, `C = 0.049` | `h ≤ 1e-4` | `3e-4` (1 knot) | 1.9e-6 at `h = 1e-2` (27 knots + 2 Re) |
| `u*` | `chord_2` | 1.58 | **none fitted** (no step rises clear of its floor) | `h ≤ 3e-5` | `1e-4` (2 knots) | not defined — see below |

Three things this shows. First, the C² breaks are where the theory says: the
steps whose stencils cross a knot or a row are all at `h ≥ 1e-4`, two to four
decades above any `h*_j`, and the crossing-free staircase at `h ≤ 3e-5` is
pure `h²` truncation. Second, the departures from `C h²` are *dominated by the
higher-order tail of the smooth function*, not by the crossings: at `x0` the
excess falls from 3.4e-4 at `h = 1e-2` to 1.2e-5 at `3e-3` to 1.6e-8 at
`1e-3`, i.e. by 30–740 × per decade, faster than the fitted `h²`; and at
`h = 1e-3` the departure is smaller than the fit's own 0.3 % extrapolation
band (`± 2.8e-8`), so the smallest crossing-populated step shows no departure
that is distinguishable from the fit. The departures that *are* significant
are ≤ 0.2 % of the error and ≈ 1e-7 of the gradient there, and their sign is
not stable: negative at all three points now, positive at `x0`'s worst
variable before the re-run, and different in sign between the two crossing
steps at `x0` (`−3.4e-4` at `1e-2`, `−1.2e-5` at `3e-3`, `+1.2e-9` at
`3.16e-4`). A crossing term of the form `O(h)·Δf″` would be positive, single
signed, and would not shrink faster than `h²`; nothing here has that
signature. Third, the `u*` panel is the honest failure of the *fit*, not of the
artefact: at that point no crossing-free step rises 30 × clear of its own
round-off floor (the cleanest, `h = 3e-5`, sits at 1.5 × its floor), so the
`C h²` regime is not observable in the double-precision `J` there. The artefact
records `smooth_h2_fit_rule = "no law fitted: … the O(h^2) truncation regime is
not observable at this point"` and emits `None` for the prediction and the
excess rather than a fabricated number (this used to propagate a NaN into
`points[2]` — fixed 2026-09-19, see `docs/journal/Session Notes/2026-09-19.md`).

Nothing is left unexplained at a Tier 3 point: at `h*_j` the disagreement is
round-off (0.94–2.20 × the measured floor, all three points), with **zero**
knot, Reynolds-row or Buhl crossings inside any `h*_j` stencil. Above `1e-4`
the departures from `h²` are a higher-order truncation tail that carries
counted crossings; the Buhl blend is crossed only at `h = 1e-2` (`x0`), never
near an `h*_j`.

## Reproduce

    python verification/gradient_verification/run_tier3.py                 # ~1 min
    python verification/gradient_verification/run_tier3.py --mid-iterate 17
    python verification/gradient_verification/run_tier4.py                 # ~23 s, needs tier3.json

`run_tier3.py` re-measures `eps_j` at `k = 17` and `u*` itself (there is no
A3 sweep at those points); at `x0` it takes `eps_j` from
`verification/fd_step_size/sweep.json` and additionally re-derives it, which
is what reproduces A3's values to 5e-7 relative. `--mid-iterate` must be kept
in step with `nit/2` of whatever `fd_optimisation/result.json` currently
holds (it is 17 of 35 since 2026-09-19).

## Files

- `run_tier3.py` — the script.
- `tier3.json` — everything above: per point, the adjoint and FD gradients
  in both units, the local three-step sweep, `h*_j`, `eps_j` (source and
  re-measured), absolute and relative errors, the pass flags, the Taylor
  draws with their directions and remainders, the timings.
- `adjoint_x0.json` — the adjoint gradient at `x0` in the format
  `run_sweep.py --reference` reads.
- `v_curve_vs_adjoint.png` — A3's sweep against the adjoint.
- `tier3_agreement.png` — the acceptance figure (disagreement over the
  acceptance scale | relative error, the three points named reference,
  mid-run iterate and optimum).
- `run_tier4.py` — the Tier 4 script (reads `tier3.json` and A3's `sweep.json`).
- `tier4.json` — per point: every variable's crossing counts and margins at
  `h*_j`, `δJ` and its samples, the round-off floor, the clean-step FD, and
  for the worst variable the full 15-step sweep with counts, the `C h²`
  law and the excess at every step. The fit is recorded with its rule and
  the steps it used (`smooth_h2_fit_rule`, `smooth_h2_fit_steps`); when no
  crossing-free step clears 30 × its round-off floor the coefficient and
  every derived prediction/excess are `null` rather than a fabricated or NaN
  value.
- `tier4_attribution.png` — the attributed V-curves for the reference and the
  optimum; `tier4_attribution_iterate.png` — the same for the mid-run iterate,
  alone.
