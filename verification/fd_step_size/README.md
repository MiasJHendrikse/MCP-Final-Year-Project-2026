# Choosing the finite-difference step size

Finite-difference gradients need a step size: too large and truncation error
dominates, too small and round-off does. This study sweeps the step for all
ten design variables at the Schmitz blade x₀ and picks the best one. It also
records the noise scale `ε_j` that the adjoint gradients are later judged
against.

**Result: `h* = 3.16e-6`** (in the scaled variable `u`) gives gradients with a
relative error of about 1e-8, far below anything that matters to the
optimiser.

The design variables are scaled as `u = (d − lo) / span` using the bounds in
`config/rotor_design.yaml` (chord 0.045–0.30 m, twist −2° to 35°), so every
gradient here is in those units. The objective uses the 300 rpm operating law
and the fixed generator rating.

## Method

    python verification/fd_step_size/run_sweep.py

At `u0 = bounds.to_scaled(x0)`, the script takes the central difference of
`fun(u) = J(u) / |J(u0)|` (`J = −AEP`) for all ten variables at 15 steps
spaced half a decade apart, `h = logspace(−2, −9, 15)`. Every gradient at every
step is saved in `sweep.json`, and nothing was discarded. No stencil left the
polar table.

| quantity | value |
|---|---|
| `J(x0)` | −10.247707 MWh/yr (AEP 10.247707) |
| one objective evaluation | 0.241 s (17 rotor solves) |
| evaluations | 300 (15 steps × 2 × 10) |
| wall time | 79.6 s |

## How the step is chosen

`h*_j` is the step that minimises the local flatness
`|g_j(h) − g_j(h/√10)| + |g_j(h) − g_j(h√10)|`, and the global `h*` is the
grid step nearest the median of `log10 h*_j`. The noise scale is

    ε_j = max(|g_j(h*_j) − g_j(h*_j/√10)|, |g_j(h*_j) − g_j(h*_j√10)|)

and the adjoint-versus-finite-difference comparison (Tier 3) is judged against
it, never against anything looser.

**`ε_j` is a flatness estimate, not a measured noise floor**, and that
difference turns out to matter. Because it is built from differences between
neighbouring steps, it relies on the same smoothness it is meant to check. If
both neighbours of `h*_j` happen to sit on the same small kink, the flatness
test selects that step and `ε_j` comes out too small. `chord_4` is exactly that
case: `ε_4 = 5.4e-11` is the smallest of the ten, while the typical
step-to-step jitter on `chord_4`'s own plateau is 8.5e-10, sixteen times larger.

The adjoint agrees with the finite difference to 8.8e-10 for that variable,
which is in line with the jitter but 16.3 times `ε_4`. Tier 4 settles it by
measuring the round-off floor of `J` directly: it samples `J` nine times along
`u + t e_4` for `t = −4e-12 … 4e-12`, fits a line, and uses the residual `δJ`.
The floor `δJ/h*` is 4.00e-10, **7.4 times larger than `ε_4`**, and the
adjoint–finite-difference disagreement is 2.20 times that floor, the same as
every other variable. So `chord_4`'s failure says something about `ε_j`, not
about the gradient. Tier 3 therefore reports one failing variable at x₀ (see
`gradient_verification/README.md`). I left the estimator and the test
unchanged, so the weakness stays visible instead of being hidden by a looser
threshold.

## Result

Gradients in MWh/yr per unit `u` (`fun` units times `|J0|`); `ε_j` likewise.

| variable | `h*_j` | `g_j(h*_j)` | `ε_j` | `ε_j / |g_j|` |
|---|---|---|---|---|
| chord_0 | 3e-06 | −0.011460 | 3.1e-10 | 2.7e-08 |
| chord_1 | 1e-05 | −0.123968 | 9.7e-11 | 7.8e-10 |
| chord_2 | 1e-05 | −0.207398 | 1.0e-09 | 5.0e-09 |
| chord_3 | 3e-06 | −0.011711 | 6.5e-10 | 5.6e-08 |
| chord_4 | 3e-06 | 0.127766 | 5.4e-11 | 4.2e-10 |
| twist_0 | 3e-05 | 0.004553 | 7.5e-11 | 1.7e-08 |
| twist_1 | 3e-06 | 0.009155 | 1.7e-10 | 1.8e-08 |
| twist_2 | 3e-06 | 0.019205 | 1.1e-09 | 5.8e-08 |
| twist_3 | 3e-06 | −0.107662 | 1.1e-09 | 9.8e-09 |
| twist_4 | 3e-06 | −0.039202 | 1.8e-09 | 4.6e-08 |

**Global `h* = 3.162e-06`**, which the finite-difference-driven optimisation
uses. Seven of the ten `h*_j` are at that step too. `chord_1` and `chord_2`
prefer 1e-05 and `twist_0` prefers 3e-05, all on the flat bottom of their
V-curves, where the flatness estimate is weakest.

**What the gradient says.** AEP falls as the inboard and mid-span chord
control points grow, most strongly `chord_2` (40–60 % of the span), and rises
only with the tip chord `chord_4`. The outboard twist points (`twist_3`,
`twist_4`) reduce AEP and are the only twist variables with much effect. Nine
of the seventeen bins (11.5–19.5 m/s) are at rated power and contribute
nothing to the gradient. The 10.5 m/s bin is on the rotor-speed ceiling but
below rated, so it does contribute. That's why the largest gradient is 0.207
MWh/yr per unit `u`, compared with 0.558 under the earlier fixed-λ law
(where the largest were `chord_3 = +0.558`, `chord_4 = +0.386` and
`twist_3 = −0.252`), and why `chord_3` has changed sign. The step-size
conclusions were the same under both laws.

## The V-curve

`v_curve.png` plots `|g_j(h) − g_j(h*_j)|` against `h` on log-log axes for each
variable, with `h²` and `1/h` guide lines.

Because the objective is C¹ but not C² (the Buhl blend at `a = 0.4`, and the
polar interpolant's C² breaks at every 0.5° angle-of-attack knot and every
Reynolds row), I expected a plateau around 1e-4 to 1e-3. **That isn't what
the sweep shows.** The curves are a clean V:

- From `h = 1e-2` down to about 1e-5, the deviation falls with slope `h²`, the
  truncation error of a central difference on a smooth function, over several
  decades with no plateau.
- Below that it rises as `1/h`, which is round-off from the objective's own
  precision (`brentq` with `xtol = 1e-14` at every station).
- The minimum is between 1e-6 and 1e-5 (3e-06 for seven variables, 1e-05 for
  two, and 3e-05 for `twist_0`, whose gradient is the smallest), where
  `ε_j / |g_j|` is between 4e-10 and 6e-08. By `h = 1e-8` the deviation is
  back up to 1e-6 to 1e-4.

The C² breaks don't show here because they only matter when the stencil
*straddles* one. At x₀ the stencil spans `|dα/du_j| h` in angle of attack,
roughly a few degrees times `h`. For `h ≤ 1e-3` that's far smaller than the
0.5° knot spacing, so almost no bin-and-station pair straddles a knot at the
steps that matter, and the stations above `a = 0.4` aren't within `h` of the
blend either. At `h = 1e-2` the deviation is 1e-4 to 1e-1 relative, which is
genuine truncation error, not knot crossings. Tier 4 counts crossings
explicitly: at x₀ the first one appears at `h = 1e-3`, and at the optimum u\*
at `h = 1e-4`, so where they appear depends on the design point.

For the optimisation, this means `h* = 3e-6` gives gradients with relative
error of about 1e-8 or better in every direction that matters, far below
SLSQP's `ftol = 1e-8` on an objective of order one. Any step between 1e-6 and
1e-4 would have worked, so the choice isn't delicate.

## Re-plotting against another gradient

    python verification/fd_step_size/run_sweep.py --reference gradient.json

redraws `v_curve.png` from `sweep.json` with `|g_j(h) − ref_j|` on the y-axis,
without re-running anything. `gradient.json` holds either `gradient_scaled`
(`fun` units per unit `u`) or `gradient_mwh_per_u` (converted using the stored
`J0`), and an optional `label`. `--replot` on its own redraws the committed
figure. `verification/gradient_verification/v_curve_vs_adjoint.png` is this
plot with the adjoint as the reference.

## Files

- `run_sweep.py`: the script.
- `sweep.json`: every gradient at every step; `h*_j`, `h*` and `ε_j` in both
  units; the flatness table; x₀, `u0`, `J0`, timings and the bounds used.
- `v_curve.png`: the figure.
