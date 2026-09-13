# Gradient verification: the discrete adjoint against finite differences (Phase 3, Tiers 3 and 4)

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional** (pending
the hub-radius / root-attachment decision); `chord_min_m = 0.045 m`,
`twist_min = -2°`, `twist_max = 35°` are grounded. The scaling
`u = (d - lo)/span` and therefore every gradient here is stated under
provisional bounds. Bounds come from
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` and are not in `config/`.

Tiers 1 and 2 are tests, not artefacts: `tests/test_adjoint_partials.py`
(every partial against a complex step of the code's own residual, to
`1e-13 max(1, |partial|)`) and `tests/test_adjoint_transpose.py` (the
transpose identity to `1e-14`). This directory holds Tier 3 (the assembled
gradient against the FD noise floor) and Tier 4 (attribution of whatever
disagreement is left).

## Tier 3 — what was run

    python verification/gradient_verification/run_tier3.py        # ~55 s

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

## Tier 3 — result: all pass

| point | `J` [MWh/yr] | worst variable | worst `|diff|/eps_j` | worst rel. error | Taylor ratios (min over 3 draws × 2 decades) |
|---|---|---|---|---|---|
| `x0` | −10.270157 | `chord_4` | **2.78** | 1.5e-7 | 98.4 |
| iterate `k = 17` | −10.292371 | `chord_3` | **1.22** | 1.3e-5 | 98.6 |
| optimum `u*` | −10.292482 | `chord_4` | **1.14** | 7.6e-5 | 98.6 |

Every one of the 30 (point, variable) pairs is within `3 eps_j`; 27 of 30
are within `1 eps_j`. The three above `1 eps_j` are all chord control
points with the smallest `h*_j` (`chord_4` at `3e-7`, `chord_3` at `1e-6`),
where the FD estimate is deepest in the round-off side of its V-curve —
Tier 4 below attributes them. The relative errors are 1e-9 … 7.6e-5,
three to six decades below the tripwire; the largest, at `u*`, is a
gradient of 6e-5 MWh/yr per `u` disagreeing by 5e-9, i.e. the same absolute
noise as everywhere else divided by a small number.

At `x0`, per variable (MWh/yr per unit `u`; the full table for the other
two points is in `tier3.json`):

| | adjoint | FD(`h*_j`) | `|diff|` | `h*_j` | `eps_j` | `|diff|/eps_j` |
|---|---|---|---|---|---|---|
| chord_0 | −0.0001447580 | −0.0001447580 | 2.1e-11 | 1e-5 | 2.3e-10 | 0.09 |
| chord_1 | 0.0193557676 | 0.0193557680 | 3.2e-10 | 3e-6 | 9.9e-10 | 0.33 |
| chord_2 | 0.1830988925 | 0.1830988932 | 6.2e-10 | 3e-6 | 2.2e-9 | 0.29 |
| chord_3 | 1.0552353984 | 1.0552353972 | 1.1e-9 | 1e-6 | 5.8e-9 | 0.20 |
| chord_4 | 0.7118018781 | 0.7118018742 | 3.9e-9 | 3e-7 | 1.4e-9 | 2.78 |
| twist_0 | −0.0044235196 | −0.0044235197 | 9.6e-11 | 3e-6 | 1.5e-10 | 0.64 |
| twist_1 | −0.1055328836 | −0.1055328836 | 6.9e-11 | 1e-5 | 3.9e-10 | 0.18 |
| twist_2 | −0.2525859665 | −0.2525859665 | 6.4e-11 | 1e-5 | 1.1e-9 | 0.06 |
| twist_3 | −0.6700875210 | −0.6700875205 | 4.9e-10 | 1e-6 | 1.2e-9 | 0.39 |
| twist_4 | −0.2874304744 | −0.2874304747 | 3.1e-10 | 3e-6 | 1.7e-9 | 0.19 |

The Taylor remainders fall by 98.4–100.8× per decade at every point and
every draw (a wrong gradient gives ≈10×; the assertion is ≥ 30). The
re-measured `eps_j` at `x0` reproduces A3's committed values to 1e-8
relative (same evaluations, deterministic code).

**Wall time at `x0`:** `J` 0.213 s; FD gradient (20 evaluations) 4.27 s;
adjoint gradient (one forward solve + partials + assembly) **0.244 s** —
1.15 × `J`, against the brief's "≲ 2 × J". The whole Tier 3 run cost 241
objective evaluations and 3 adjoint evaluations.

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

*(section written by `run_tier4.py`'s run; see below)*

## Reproduce

    python verification/gradient_verification/run_tier3.py                 # ~55 s
    python verification/gradient_verification/run_tier3.py --mid-iterate 17

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
