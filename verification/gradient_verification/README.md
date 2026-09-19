# Gradient verification: the discrete adjoint against finite differences (Phase 3, Tiers 3 and 4)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. The scaling
`u = (d - lo)/span` and therefore every gradient here is stated under
provisional bounds. Bounds come from
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` and are not in `config/`.

**Re-run 2026-09-19 under the fixed generator rating.** The objective now
holds power above rated at `operating.rated_power_w` from
`config/rotor_design.yaml` (provisionally the baseline's own `P_aero(11 m/s;
x0)`, pending input B2) instead of at `P_aero(V_rated; d)`; the adjoint
system lost its eighteenth operating point (the rated solve) and is 17 × 25
= 425 states, with the capped bins a constant in `J` (`docs/adjoint_derivation.md`
§7, revision note). Tiers 1 and 2 pass unchanged on the 425-state system.
Tiers 3 and 4 below are re-measured against the re-run A3 sweep and the
re-run A4 iterates (`k = 13 = nit/2` of a 27-iteration run, was 17 of 34).
The 2026-09-13 floating-cap numbers (30/30 within 3 ε, worst 2.78 ε, zero
crossings at `h*`) are in the git history and in
`docs/adjoint_derivation.md` §9.

Tiers 1 and 2 are tests, not artefacts: `tests/test_adjoint_partials.py`
(every partial against a complex step of the code's own residual, to
`1e-13 max(1, |partial|)`) and `tests/test_adjoint_transpose.py` (the
transpose identity to `1e-14`). This directory holds Tier 3 (the assembled
gradient against the FD noise floor) and Tier 4 (attribution of whatever
disagreement is left).

## Tier 3 — what was run

    python verification/gradient_verification/run_tier3.py        # ~55 s

At three points — `x0`, the FD-driven SLSQP run's mid-run iterate `k = 13`
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

## Tier 3 — result: all pass

| point | `J` [MWh/yr] | worst variable | worst `|diff|/eps_j` | worst rel. error | Taylor ratios (min over 3 draws × 2 decades) |
|---|---|---|---|---|---|
| `x0` | −10.270157 | `twist_3` | **2.04** | 1.4e-6 | 98.4 |
| iterate `k = 13` | −10.282512 | `twist_0` | **2.34** | 3.6e-6 | 98.7 |
| optimum `u*` | −10.282564 | `chord_4` | **1.56** | 9.1e-6 | 98.7 |

Every one of the 30 (point, variable) pairs is within `3 eps_j`; 26 of 30
are within `1 eps_j`. The four above `1 eps_j` — `twist_3` at `x0`,
`chord_0` and `twist_0` at `k = 13`, `chord_4` at `u*` — are all variables
whose `eps_j` is smallest in absolute terms (9.1e-11 … 7e-10 MWh/yr per
`u`), so the same round-off disagreement everywhere else sees is divided by
a smaller floor; Tier 4 below attributes them. The relative errors are
7e-10 … 9.1e-6, two to six decades below the tripwire; the largest, at
`u*`, is a gradient of 1.5e-5 MWh/yr per `u` (`twist_1`, nearly inert)
disagreeing by 1.4e-10.

At `x0`, per variable (MWh/yr per unit `u`; the full table for the other
two points is in `tier3.json`):

| | adjoint | FD(`h*_j`) | `|diff|` | `h*_j` | `eps_j` | `|diff|/eps_j` |
|---|---|---|---|---|---|---|
| chord_0 | −0.0010836794 | −0.0010836795 | 6.1e-11 | 1e-5 | 4.3e-10 | 0.14 |
| chord_1 | −0.0001335633 | −0.0001335631 | 1.9e-10 | 3e-6 | 4.0e-10 | 0.47 |
| chord_2 | 0.0741384500 | 0.0741384504 | 3.5e-10 | 3e-6 | 8.7e-10 | 0.40 |
| chord_3 | 0.5579270983 | 0.5579270973 | 1.0e-9 | 1e-6 | 1.4e-9 | 0.70 |
| chord_4 | 0.3864641320 | 0.3864641336 | 1.6e-9 | 3e-7 | 2.3e-9 | 0.69 |
| twist_0 | 0.0042873197 | 0.0042873194 | 2.7e-10 | 3e-6 | 4.1e-10 | 0.66 |
| twist_1 | −0.0146783030 | −0.0146783030 | 1.0e-11 | 1e-5 | 3.6e-10 | 0.03 |
| twist_2 | −0.0545534449 | −0.0545534448 | 5.6e-11 | 1e-5 | 9.9e-10 | 0.06 |
| twist_3 | −0.2516472640 | −0.2516472642 | 1.9e-10 | 1e-5 | 9.1e-11 | 2.04 |
| twist_4 | −0.0967123733 | −0.0967123736 | 2.5e-10 | 3e-6 | 1.7e-9 | 0.14 |

The Taylor remainders fall by 98.4–102.3× per decade at every point and
every draw (a wrong gradient gives ≈10×; the assertion is ≥ 30). The
re-measured `eps_j` at `x0` reproduces A3's committed values to 5e-7
relative (same evaluations, deterministic code).

**Wall time at `x0`:** `J` 0.203 s; FD gradient (20 evaluations) 4.04 s;
adjoint gradient (one forward solve + partials + assembly) **0.233 s** —
1.15 × `J`, against the brief's "≲ 2 × J". The whole Tier 3 run cost 241
objective evaluations and 3 adjoint evaluations, 50.6 s.

## Tier 3 — the V-curves re-plotted against the adjoint

`v_curve_vs_adjoint.png` is A3's step-size sweep at `x0` (`sweep.json`,
no re-run) re-plotted with `|g_j(h) − adjoint_j|` on the y-axis, via
`run_sweep.py`'s `--reference` machinery (`adjoint_x0.json`); A3's own
figure is untouched. With an independent reference the curves become the
textbook V that A3's self-referenced plot could not show: `h²` truncation
from `1e-2` down to `1e-5`, a minimum of 1e-11 … 1e-10 (fun units per `u`)
at `h = 1e-5 … 3e-6`, and `1/h` round-off below. The adjoint sits at the
bottom of every variable's V, which is the statement Tier 3 makes in one
picture: the FD reference is limited by its own noise floor, and the
adjoint is inside it.

`tier3_agreement.png` shows `|diff|/eps_j` and the relative error per
variable at the three points, with the acceptance lines.

## Tier 4 — degradation at the polar interpolation

    python verification/gradient_verification/run_tier4.py        # ~21 s

Tier 3 left 30 disagreements, all inside `3 eps_j`, four of them above
`1 eps_j`. Tier 4 explains every one of them by asking what the central
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

### Result at `h*_j`: one stencil crosses one knot, to no measurable effect; the rest is round-off

| point | worst variable | knot / Re-row / Buhl crossings at `h*_j`, any variable | nearest break in stencil half-widths, any variable (α / Re / Buhl) | `δJ` [MWh/yr] | round-off floor `δJ/h*` | observed `|diff|` | `|diff|` / floor |
|---|---|---|---|---|---|---|---|
| `x0` | `twist_3` (`h*` 1e-5) | 0 / 0 / 0 | 3.4 / 53 / 1067 | 8.4e-16 | 8.4e-11 | 1.9e-10 | **2.2** |
| iterate `k = 13` | `twist_0` (`h*` 3e-6) | **1** / 0 / 0 | **0.90** / 66 / 2018 | 1.3e-15 | 4.2e-10 | 9.7e-11 | **0.23** |
| optimum `u*` | `chord_4` (`h*` 3e-7) | 0 / 0 / 0 | 2.1 / 57 / 1839 | 1.4e-15 | 4.4e-9 | 1.1e-9 | **0.25** |

At every point and for all ten variables, the stencil at `h*_j` moves every
station's α by at most 1.6e-4° against knots 7.9e-5 … 3.4e-4° away, moves
Re by at most 5.5 against rows 24 … 59 away, and moves `a` by at most
9.8e-6 against `|a − 0.4| ≥ 2.7e-3`. **One stencil straddles a knot:**
`twist_4` at `k = 13`, `h* = 3e-6`, where one station's α sits 7.9e-5°
from a 0.5° node and the stencil half-width is 8.9e-5° (margin 0.90
half-widths). Its measured degradation is **6.6e-13 MWh/yr per `u`, 0.0004
`eps_j`** — the O(h)·Δf″ term of a single knot crossing at `h = 3e-6` is
three decades below the round-off floor, and this pair is the *best*
agreeing of the ten at that point. Buhl crossings contribute zero
everywhere; Reynolds-row crossings zero everywhere. Elsewhere agreement is
round-off: the three worst-variable disagreements are 2.2 ×, 0.23 × and
0.25 × the measured `δJ/h*`. The 2.2 at `x0` is `twist_3`, whose A3
flatness rule picked `h* = 1e-5` and whose `eps_j` (9.1e-11) is the
smallest of the thirty — a disagreement of 1.9e-10 that is 0.03–0.7 `eps`
for every other variable at `x0` is 2.04 `eps` for it. At the largest
crossing-free step on the grid (`h = 1e-4`) its FD is 2.3e-9 from the
adjoint, on the smooth `C h²` law (`C = 0.23`), and at `h = 3e-6` it is
9.5e-11. There is no shrinking of `h` to do.

### The whole V attributed (`tier4_attribution.png`)

For the worst variable at each point, `|FD(h) − adjoint|` over A3's 15-step
grid next to the crossing counts and the floor. Fitting the smooth law
`C h²` on the crossing-free steps clear of the floor and subtracting it
gives the polar / Buhl contribution at every crossing-populated step, with
its sign:

| point | truncation-side slope (h ≥ 1e-4) | crossing-free down to | largest excess over `C h²` at a crossing step | steps with a Buhl crossing |
|---|---|---|---|---|
| `x0`, `twist_3` | 2.69 | `h ≤ 1e-4` (2 knots at 3e-4, 3 at 1e-3) | +2.7e-4 at `h = 1e-2` (60 knots) | none |
| `k = 13`, `twist_0` | 2.05 | `h ≤ 1e-4` (2 knots at 3e-4) | +1.2e-6 at `h = 1e-2` (61 knots + 1 Buhl) | `h = 1e-2` only |
| `u*`, `chord_4` | 2.30 | `h ≤ 3e-5` (1 knot at 1e-4) | −3.2e-3 at `h = 1e-2` (68 knots + 23 Re + 1 Buhl); −1.0e-5 at `h = 3e-4` (3 knots + 2 Re) | `h = 1e-2` only |

Three things this shows. First, the C² breaks are where the theory says:
the *only* steps whose error departs from `C h²` are steps whose stencil
crosses a knot, a row or the Buhl blend, and every such step is at
`h ≥ 1e-4` — two to four decades above any `h*_j` — with the single
exception at `k = 13` recorded above, whose effect is unmeasurable. Second,
the polar contribution is small even where it exists: at `x0` and `k = 13`
the crossing steps sit *above* the `h²` line (positive excess, slopes 2.69
and 2.05 — the O(h)·Δf″ terms from a few knot crossings adding to smooth
truncation); at `u*` the 23 Reynolds-row crossings and the Buhl crossing at
`h = 1e-2` pull the error *below* the `h²` line (excess −3.2e-3, the
crossing terms cancelling truncation), and the 3e-4 step shows the
non-monotone dip visible in the right-hand panel — the signature the brief
predicted for FD at a C² break, and one the adjoint does not have. Third,
the Buhl blend is crossed only at `h = 1e-2`, never inside an `h*_j`
stencil, so its line in the attribution reads zero at every Tier 3 point.

Nothing is left unexplained: at `h*_j` the disagreement is round-off
(0.23–2.2 × the measured floor) plus one counted knot crossing worth
0.0004 `eps`; above `1e-4` the departures from `h²` are counted knot / row
/ Buhl crossings; in between, smooth truncation.

## Reproduce

    python verification/gradient_verification/run_tier3.py                 # ~55 s
    python verification/gradient_verification/run_tier3.py --mid-iterate 13
    python verification/gradient_verification/run_tier4.py                 # ~21 s, needs tier3.json

## Files

- `run_tier3.py` — the script.
- `tier3.json` — everything above: per point, the adjoint and FD gradients
  in both units, the local three-step sweep, `h*_j`, `eps_j` (source and
  re-measured), absolute and relative errors, the pass flags, the Taylor
  draws with their directions and remainders, the timings.
- `adjoint_x0.json` — the adjoint gradient at `x0` in the format
  `run_sweep.py --reference` reads.
- `v_curve_vs_adjoint.png` — A3's sweep against the adjoint.
- `tier3_agreement.png` — the acceptance figure.
- `run_tier4.py` — the Tier 4 script (reads `tier3.json` and A3's `sweep.json`).
- `tier4.json` — per point: every variable's crossing counts and margins at
  `h*_j`, `δJ` and its samples, the round-off floor, the clean-step FD, and
  for the worst variable the full 15-step sweep with counts, the `C h²`
  law and the excess at every step.
- `tier4_attribution.png` — the attributed V-curves.
