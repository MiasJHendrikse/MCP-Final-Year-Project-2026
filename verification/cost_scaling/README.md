# `verification/cost_scaling/` — gradient cost against the number of design variables

Phase 4, Step 3 (plan §8.6). Produced 2026-09-19 by `run_scaling.py` from the
working tree of the commit that carries this directory (`scaling.json::src_commit`
is that commit's parent, `ceb7dd3`, because the script and its output were
committed together). Wall time 11.6 min, almost all of it the FD gradient at
`n = 160` (320 objective evaluations, five times).

**Claim it backs.** Under the operating law and the configured box, one adjoint
gradient costs **1.09–1.12 objective evaluations, independent of `n`**
(`p = −0.003` on `t = a·nᵖ` over `n = 10 … 160`), while central FD costs
exactly `2n` evaluations (`p = 1.000`). The Phase 4 constraint's gradient — the
same chain over the nine-point load set — costs **0.55 ×** the objective
adjoint (9/17 = 0.53 points). The FD-to-adjoint ratio is 18 at the production
`n = 10` and 297 at `n = 160`.

## What was run

For `n_chord = n_twist = k ∈ {5, 10, 20, 40, 80}` (`n = 2k`), 25 BEM strips
throughout, a fresh `BladeParameterisation`, `DesignBounds.from_config(k, k)`,
`ScaledProblem`, `BEMSystem` and `RootMomentSystem` per `k`, all evaluated at
the least-squares projection of the analytic Schmitz blade onto that block
count (`build_schmitz_baseline(parameterisation)` — the same construction as
`x0`; `k = 5` *is* `x0`). Timed, best of 5 (`time.perf_counter`):

| measured | how |
|---|---|
| `J` | `ScaledProblem.J(u)`: 17 forward solves, one per bin |
| central FD | `ScaledProblem.jac_fd(u, h*)` at the committed global `h* = 3.162277660168379e-06` in `u`; `2n` evaluations of `fun` |
| tangent (direct), total | one `solve` + one `partials`, then `n` products `BEMSystem.tangent_from_parts(parts, phi, d, e_j)` — the helper was added for this study so the partials are computed once, not per direction; `BEMSystem.tangent` now wraps it and computes what it always computed |
| tangent, the `n` products alone | the same `n` products on partials already in hand — the only part of the tangent that grows with `n` |
| adjoint | `BEMSystem.gradient(d)`: one solve, partials over 17 × 25 stations, `ψ`, one assembly |
| moment adjoint | `RootMomentSystem.gradient(d)` with `m_ref_nm = 177.3755406092969` passed explicitly (`ScaledProblem.load_system()` only builds for the 5 + 5 parameterisation, so the constant is passed and cross-checked once against the `k = 5` problem's `m_ref`) |
| forward solves per gradient | counted by wrapping `solve` (adjoint, tangent, moment adjoint) or from `n_fun_evals` (FD); not timed |

Every timed gradient was also checked to be the same gradient: tangent vs
adjoint agree to `≤ 4.4e-17` and FD vs adjoint to `≤ 9.0e-11` (max abs
difference relative to `max |dfun/du|`) at every `k`.

## The numbers (`scaling.json`)

| k | n | AEP at the point [MWh/yr] | vs x0 | `J` [s] | FD [s] | tangent total [s] | tangent products [s] | adjoint [s] | moment adjoint [s] | FD / adjoint | adjoint / `J` | solves FD / tan / adj / mom |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 10 | 10.247707 | 0 | 0.206 | 4.25 | 0.242 | 0.00035 | 0.231 | 0.128 | 18.4 | 1.12 | 20 / 1 / 1 / 1 |
| 10 | 20 | 10.249287 | +0.0154 % | 0.208 | 8.47 | 0.231 | 0.00069 | 0.229 | 0.126 | 36.9 | 1.10 | 40 / 1 / 1 / 1 |
| 20 | 40 | 10.249287 | +0.0154 % | 0.209 | 17.16 | 0.233 | 0.00140 | 0.230 | 0.130 | 74.7 | 1.10 | 80 / 1 / 1 / 1 |
| 40 | 80 | 10.249287 | +0.0154 % | 0.212 | 34.09 | 0.235 | 0.00283 | 0.231 | 0.127 | 147.6 | 1.09 | 160 / 1 / 1 / 1 |
| 80 | 160 | 10.249287 | +0.0154 % | 0.209 | 67.80 | 0.236 | 0.00589 | 0.228 | 0.128 | 297.4 | 1.09 | 320 / 1 / 1 / 1 |

Fits `t = a·nᵖ` on the five points:

| method | `a` [s] | `p` | expectation |
|---|---|---|---|
| central FD | 0.425 | **1.000** | 1 (`2n` evaluations of a constant-cost `J`) |
| tangent, total | 0.239 | −0.004 | flat: one solve + partials dominates; affine fit `c = 0.2357 s`, `b = −1.8e-6 s` per variable — a slope indistinguishable from zero |
| tangent, the `n` products alone | 3.3e-5 | **1.018** | 1, with a slope four orders below FD's (`3.3e-5` vs `0.425` s per variable) |
| adjoint | 0.233 | **−0.003** | 0: the cost is the partials over the stations, not the variables; the `N_cᵀ` assembly is 25 × n and invisible |
| moment adjoint | 0.127 | 0.001 | 0, at 9/17 of the objective adjoint |
| `J` | 0.204 | 0.006 | 0 |

`scaling.png` plots the five points per method on log–log axes with the fits.

## The same blade throughout

The AEP column shows it: `k = 5` is `x0` to the bit (10.247707391810723); from
`k = 10` the 25-station Schmitz targets are reproduced to `2.7e-5 m` rms and the
AEP is 10.249287 (+0.0154 % vs `x0` — the fit error the 5 + 5 projection carries,
`verification/spline_fit_error/`), identical to nine digits from `k = 20` on.
For `k > 25` the fit is underdetermined and `lstsq` returns the minimum-norm
control points: the station chords and twists are reproduced to round-off
(`1e-16` rms), and 2 (`k = 40`) and 32 (`k = 80`) of the control points fall
outside the configured box. Neither matters here — the BEM sees the stations,
the timing sees the operation count — and neither `k > 5` point is a design.
`KS` at the point is 1.00044 at `k = 5` and 1.00012 from `k = 10` on, the same
blade's moment to 0.03 %.

## The honest note (derivation §8, restated with these numbers)

With a diagonal `∂R/∂x` and one scalar state per station, the adjoint does
**not** save a linear solve — `ψ` is 425 divisions. Its advantage over central
FD is the **`2n` factor and accuracy**: at the production `n = 10` it is 18 ×
faster per gradient and free of the FD noise floor and of the `O(h)·Δf″`
contamination at the polar interpolant's C² breaks; at `n = 160` it is 297 ×
faster. The tangent (direct) mode is *also* flat in `n` here, because one
`partials` (Python loops over 17 × 25 stations, ~0.23 s) dominates the `n`
cheap numpy products (`3.3e-5 s` each): with this residual structure forward
mode is as cheap as reverse mode up to `n ≈ 7000`, and the adjoint's case
against it is the accuracy argument and the independence of the two code
paths (Tier 2), not the cost. What the study does establish is that the
project's gradient cost is a constant ~1.1 objective evaluations however many
control points the blade carries, so the representation study's choice of
5 + 5 was made on representation grounds alone, never on gradient cost.

## Tests

`tests/test_cost.py` carries the structural half of this study: the adjoint
makes exactly one forward solve at `k ∈ {5, 10, 20, 40}`
(`test_adjoint_forward_solve_count_is_independent_of_n`) and `jac_fd` costs
`2n` objective evaluations at `k ∈ {5, 10}` through the real `ScaledProblem`
(`test_fd_forward_solve_count_is_2n`). Wall times are machine-dependent and
are not asserted.

## Reproduce

    python verification/cost_scaling/run_scaling.py            # ~12 min, rewrites scaling.json and scaling.png
    python verification/cost_scaling/run_scaling.py --replot   # redraw scaling.png from the committed scaling.json

Machine: Windows 11, AMD64 Family 25 Model 33 (Ryzen 5000 class), Python 3.12.10,
NumPy 2.2.5. The ratios and exponents are what to compare across machines; the
absolute seconds are not.

## What would make it wrong

A change to `BEMSystem.partials` that made its cost depend on `n` (it must not:
the loops are over `(b, i)`); a change to `central_difference` that stopped
costing `2n`; or a `tangent_from_parts` that recomputed the partials per
direction. The two tests catch the first and second; the third shows up as
the tangent total climbing with `n`.
