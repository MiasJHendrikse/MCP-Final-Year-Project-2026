# Gradient cost against the number of design variables

**The claim.** One adjoint gradient costs **1.09–1.12 objective evaluations,
however many design variables there are** (fitted exponent `p = −0.003` for
`t = a·nᵖ` over `n = 10 … 160`), while central finite differences cost exactly
`2n` evaluations (`p = 1.000`). The root-moment constraint's adjoint, the same
calculation over the nine load cases, costs **0.55 times** the energy adjoint
(9/17 = 0.53 of the operating points). Finite differences are 18 times slower
than the adjoint at the `n = 10` used in the project, and 297 times slower at
`n = 160`.

All timings use the 300 rpm operating law and the configured bounds. The full
run takes 11.6 minutes, almost all of it the finite-difference gradient at
`n = 160` (320 objective evaluations, repeated five times).

## Method

For `n_chord = n_twist = k ∈ {5, 10, 20, 40, 80}` (`n = 2k`), with 25 BEM
strips throughout, the script builds a fresh `BladeParameterisation`,
`DesignBounds.from_config(k, k)`, `ScaledProblem`, `BEMSystem` and
`RootMomentSystem` for each `k`. Everything is evaluated at the least-squares
fit of the analytic Schmitz blade to that number of control points
(`build_schmitz_baseline(parameterisation)`, the same construction as x₀; at
`k = 5` it *is* x₀). Each timing is the best of five (`time.perf_counter`):

| measured | how |
|---|---|
| `J` | `ScaledProblem.J(u)`: 17 forward solves, one per bin |
| central finite differences | `ScaledProblem.jac_fd(u, h*)` at the committed `h* = 3.162277660168379e-06`; `2n` evaluations of `fun` |
| tangent (direct mode), total | one `solve` and one `partials`, then `n` products with `BEMSystem.tangent_from_parts(parts, phi, d, e_j)`. The helper computes the partials once rather than per direction; `BEMSystem.tangent` wraps it and returns the same result as before |
| tangent, the `n` products alone | the same products on partials already computed, which is the only part of the tangent that grows with `n` |
| adjoint | `BEMSystem.gradient(d)`: one solve, partials over 17 × 25 stations, `ψ`, one assembly |
| moment adjoint | `RootMomentSystem.gradient(d)` with `m_ref_nm = 177.3755406092969` passed explicitly (`ScaledProblem.load_system()` only builds for 5 + 5, so the constant is passed in and checked once against the `k = 5` value) |
| forward solves per gradient | counted by wrapping `solve` (adjoint, tangent, moment adjoint) or from `n_fun_evals` (finite differences); not timed |

Each timed gradient is also checked to be the right gradient: the tangent and
adjoint agree to within 4.4e-17, and finite differences and adjoint to within
9.0e-11 (the largest absolute difference relative to `max |dfun/du|`), at
every `k`.

## Results

`scaling.json`:

| k | n | AEP at the point [MWh/yr] | vs x₀ | `J` [s] | FD [s] | tangent total [s] | tangent products [s] | adjoint [s] | moment adjoint [s] | FD / adjoint | adjoint / `J` | solves FD / tan / adj / mom |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 10 | 10.247707 | 0 | 0.206 | 4.25 | 0.242 | 0.00035 | 0.231 | 0.128 | 18.4 | 1.12 | 20 / 1 / 1 / 1 |
| 10 | 20 | 10.249287 | +0.0154 % | 0.208 | 8.47 | 0.231 | 0.00069 | 0.229 | 0.126 | 36.9 | 1.10 | 40 / 1 / 1 / 1 |
| 20 | 40 | 10.249287 | +0.0154 % | 0.209 | 17.16 | 0.233 | 0.00140 | 0.230 | 0.130 | 74.7 | 1.10 | 80 / 1 / 1 / 1 |
| 40 | 80 | 10.249287 | +0.0154 % | 0.212 | 34.09 | 0.235 | 0.00283 | 0.231 | 0.127 | 147.6 | 1.09 | 160 / 1 / 1 / 1 |
| 80 | 160 | 10.249287 | +0.0154 % | 0.209 | 67.80 | 0.236 | 0.00589 | 0.228 | 0.128 | 297.4 | 1.09 | 320 / 1 / 1 / 1 |

Power-law fits `t = a·nᵖ` over the five points:

| method | `a` [s] | `p` | expected |
|---|---|---|---|
| central finite differences | 0.425 | **1.000** | 1 (`2n` evaluations of a constant-cost `J`) |
| tangent, total | 0.239 | −0.004 | flat, because one solve plus the partials dominates; an affine fit gives `c = 0.2357 s` and `b = −1.8e-6 s` per variable, a slope indistinguishable from zero |
| tangent, the `n` products alone | 3.3e-5 | **1.018** | 1, but with a slope four orders of magnitude below finite differences (3.3e-5 against 0.425 s per variable) |
| adjoint | 0.233 | **−0.003** | 0: the cost is the partials over the stations, not the variables; the `N_cᵀ` assembly is 25 × n and negligible |
| moment adjoint | 0.127 | 0.001 | 0, at 9/17 of the energy adjoint |
| `J` | 0.204 | 0.006 | 0 |

`scaling.png` plots the five points for each method on log-log axes, with the
fits.

## It's the same blade throughout

The AEP column confirms it. `k = 5` is x₀ exactly (10.247707391810723). From
`k = 10` on, the 25 station chords and twists of the Schmitz blade are
reproduced to 2.7e-5 m RMS, and the AEP is 10.249287 (+0.0154 % vs x₀, which is
the fitting error of the 5 + 5 spline; see `verification/spline_fit_error/`),
identical to nine digits from `k = 20` on. For `k > 25` the fit is
underdetermined and `lstsq` returns the minimum-norm control points. The
station values are still reproduced to round-off (1e-16 RMS), but 2 (at
`k = 40`) and 32 (at `k = 80`) of the control points fall outside the
configured box. That doesn't matter here, because the BEM only sees the
stations and the timing only sees the operation count, and none of the `k > 5`
points is meant as a design. KS at the point is 1.00044 at `k = 5` and 1.00012
from `k = 10` on, the same blade's moment to within 0.03 %.

## What the adjoint actually saves

With a diagonal `∂R/∂x` and one unknown per station, the adjoint doesn't save
a linear solve: `ψ` is just 425 divisions. Its advantage over central finite
differences is the `2n` factor and accuracy. At `n = 10` it is 18 times faster
per gradient, and free of the finite-difference noise floor and of the
contamination at the polar interpolant's C² breaks; at `n = 160` it is 297
times faster.

The tangent (direct) mode is *also* flat in `n` here, because one `partials`
call (Python loops over 17 × 25 stations, about 0.23 s) dominates the `n` cheap
numpy products (3.3e-5 s each). With this residual structure, forward mode is
as cheap as reverse mode up to about `n ≈ 7000`. So the case for the adjoint
over forward mode rests on accuracy and on having two independent code paths to
check against each other (Tier 2), not on cost.

What the study does establish is that the gradient costs a constant 1.1 or so
objective evaluations however many control points the blade has. The choice
of 5 + 5 control points was therefore made purely on how well they represent
the blade, never on gradient cost.

## Tests

`tests/test_cost.py` covers the structural part of this study. The adjoint
makes exactly one forward solve at `k ∈ {5, 10, 20, 40}`
(`test_adjoint_forward_solve_count_is_independent_of_n`), and `jac_fd` costs
`2n` objective evaluations at `k ∈ {5, 10}` through the real `ScaledProblem`
(`test_fd_forward_solve_count_is_2n`). Wall times depend on the machine, so
they aren't asserted.

## Reproducing it

    python verification/cost_scaling/run_scaling.py            # ~12 min, rewrites scaling.json and scaling.png
    python verification/cost_scaling/run_scaling.py --replot   # redraws scaling.png from the committed scaling.json

Measured on Windows 11 with an AMD Ryzen 5000-series CPU, Python 3.12.10 and
NumPy 2.2.5. The ratios and exponents are what to compare across machines; the
absolute times aren't.

## What would make it wrong

- A change to `BEMSystem.partials` that made its cost depend on `n`. It
  shouldn't, since its loops are over bins and stations.
- A change to `central_difference` so that it no longer costs `2n`.
- A `tangent_from_parts` that recomputed the partials for each direction.

The two tests catch the first and second; the third would show up as the
tangent total rising with `n`.
