# Solver convergence on the NREL Phase VI rotor

This checks that the BEM solver converges reliably, including on blade shapes
well away from the one it was built for.

## Residual histories

`generate_residual_histories.py` → `residual_histories.json`,
`residual_convergence.png`.

The requirement was for the solver to converge across the whole operating
envelope, for the baseline **and** perturbed geometries, with the residual
histories committed. These are those histories: not a summary, the actual sequence of `(φ, R)` the
root-finder evaluated at every station, so the reduction can be read off
rather than taken on trust.

### Coverage

13 cases × 19 stations = **247 stations**, all converged.

| | |
|---|---|
| baseline operating points | v = 5, 7, 10, 15, 20 m/s at 71.63 rpm |
| tip-speed-ratio sweep | λ = 2.0, 4.0, 6.0, 7.5 at v = 7 m/s |
| perturbed geometry | chord ×0.7, chord ×1.3, twist −6°, twist +6° |

The perturbed cases are the point. A fixed-geometry sweep only ever visits the
blade as built; an optimiser line search visits blades that are 30 % off it,
and those are where a bracketing scheme fails if it is going to.

### Result

| | |
|---|---|
| stations failed | **0 of 247** |
| worst relative residual `\|R\|/R₀` | **2.94e-15** |
| maximum root-finder iterations | **13** |
| maximum residual evaluations | **14** |

The convergence criterion is `|R| ≤ 1e-9 · R₀` — relative to the initial
residual norm, never absolute. What is achieved is
six orders tighter than that, and a returned pole would sit around `|R| ~ 1e15`,
so the criterion has roughly twenty-four orders of magnitude of separation
either side. That is what makes it a check rather than a tuned threshold.

`residual_convergence.png` plots best-`|R|`-so-far against evaluation count for
all 247 stations on one axis. The shape is the claim: every station drops
superlinearly to machine precision inside ~14 evaluations, with no stalls, no
plateaus, and no station taking a materially different path from the rest.

### Context: what this replaced

Originally the bracket was located by sampling the residual at **2000**
trial φ per station and taking the sign change nearest `atan(1/λ_r)` — a
heuristic that could select a pole (both cross the residual's sign) and that
raised `ValueError` when it found nothing, aborting an entire sweep over one
station. The bracket is now derived from the momentum region's boundary and
checked at run time; see `bem/station.momentum_region_bracket`.

## Regenerating

```
python verification/phase_vi/generate_residual_histories.py
```

Takes a few seconds. Regenerate when the solver changes; the summary block at
the top of the JSON is what to compare.
