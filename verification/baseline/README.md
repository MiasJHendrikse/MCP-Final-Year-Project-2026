# Schmitz baseline — `x₀` and the reference numbers

Plan step 1.7. The plan describes the evaluated performance of this blade as
**"the reference numbers for the entire results chapter"**, so this directory
is a fixed artefact: `x0.json` is what Phases 2–5 start from, and
`baseline_reference.json` is what every later result is compared against.

```
python verification/baseline/generate_baseline.py
```
→ `x0.json`, `baseline_reference.json`, `baseline_geometry.png`,
`baseline_operating_line.png`

(`baseline.png`, the earlier three-panel figure, is no longer written by this
script; the two files above replace it, and it is kept only because the
dashboard still points at it.)

Regeneration is deliberate, not routine — the same rule as the golden files. If
these numbers move, something in the solver or the polar layer moved, and the
commit message has to say what and why.

**Re-run 2026-09-19 under the 300 rpm operating law and the configured
bounds.** `x0.json`'s design vector is **byte-identical** to the committed
one (the bounds do not move a point that was already inside them); the only
change in that file is its `feasibility` block, which flipped from
`checked: false` to `checked: true` with no violations now that
`config/rotor_design.yaml` carries the box (below). The *reference numbers*
moved for the reason recorded in the commit that introduced the ceiling: the
rotor now runs λ = 6.5 only up to V_c = 9.67 m/s and holds 300 rpm above it,
so the 11 m/s design point is no longer at the design TSR (see "History"
below for the pre-ceiling figures, which are kept deliberately).

## Construction

1. Analytic **Schmitz** distribution at the SG6043 maximum-L/D point.
2. Projected onto the **same spline parameterisation the optimiser uses** (5+5
   control points, the count the representation study selected), by linear
   least squares. `x₀` is the *fitted* vector, not the analytic blade — the
   optimiser can only start somewhere it can represent.
3. Feasibility checked against the design-variable bounds.
4. Evaluated with the same solver and settings as everything else.

| design point | |
|---|---|
| Reynolds | 200,000 |
| α_design | 5.36° |
| C_L | 1.2825 |
| max L/D | 98.5 |

## Fitting error

| | RMS |
|---|---|
| chord | **0.515 mm** (0.389 % of the 132 mm mean chord) |
| twist | **0.106°** |

Reported rather than absorbed, per the plan. Well below any plausible
manufacturing tolerance, so `x₀` is a faithful representation of the Schmitz
blade rather than an approximation of it.

## Reference numbers

| | |
|---|---|
| Cp at design point (V = 11 m/s, λ = 5.711986642890533, 300 rpm) | **0.4516** |
| Ct at design point | **0.7029** |
| aerodynamic power at design point | **3656.9 W** |
| peak Cp on the fixed-λ curve | **0.4720 at λ = 6.5** |
| root bending moment | **177.38 N·m per blade** (was 206.7) |
| peak thrust | **517.51 N at 11 m/s** (was 597.2) |
| AEP | **10.2477 MWh/yr** (was absent) |

The earlier edition of this file said AEP was "outstanding", feasibility was
unchecked and the above-rated line was unswept. All three are now answered
(the sites' Weibull inputs arrived; the bounds are in `config/`; the 300 rpm
ceiling *is* the above-rated operating law) and the numbers above are the
answers. The paragraphs in "What is outstanding" below are kept, struck
through in substance, so that the change is visible rather than silently
overwritten.

**The peak Cp landing exactly at the design λ is the strongest available check
on the construction.** The blade is designed for λ = 6.5; if the Schmitz
formulae, the twist convention, or the solver wiring were wrong, the peak would
sit somewhere else. `tests/test_baseline.py` asserts it. Note that the peak is
on the *fixed-λ* curve, which is what the construction check needs; the
schedule the machine actually runs caps λ at 300 rpm from V = 9.67 m/s, so the
design point at 11 m/s sits at λ = 5.71 and Cp = 0.4516, below the peak — the
blade is deliberately no longer at its design TSR at rated wind speed.

Cp = 0.472 against the Betz limit of 0.593 is what a real blade with drag
achieves — the ~20 % shortfall is profile drag plus tip loss, not an error.

Root bending moment is **per blade**, not per rotor: the root attachment
carries one blade's load, and a rotor-summed figure would overstate it by a
factor of B.

### History

| | 2026-09-13 (no rotor-speed ceiling, λ = 6.5 everywhere) | 2026-09-19 (300 rpm) |
|---|---|---|
| Cp at 11 m/s | 0.4720 | 0.4516 (λ = 5.712) |
| root bending moment | 206.7 N·m | 177.38 N·m |
| peak thrust | 597.2 N | 517.51 N |
| AEP | not computable | 10.2477 MWh/yr |

The drop is arithmetic, not regression: capping the rotor speed lowers λ above
9.67 m/s, which lowers Cp and Ct at the same wind speed, hence both loads and
the AEP the baseline extracts. `tests/test_baseline.py` pins the current
numbers; the pre-ceiling ones live in the git history and in the 2026-09-19
journal entry.

## What is outstanding

**Nothing in this artefact is blocked any more.** The section is kept as a
record of what the earlier edition said and what answered it.

- *"AEP is absent, not estimated."* It was, pending the Weibull parameters
  (`docs/OUTSTANDING-INPUTS.md` §1). MJ supplied them; `aep_mwh_per_year` is
  now **10.2477 MWh/yr**, inside the artefact's own `[8, 12]` sanity band
  (`aep_in_sanity_band: true`).
- *"Feasibility was not checked."* It is now: `x0.json`'s `feasibility` block
  reads `checked: true`, `clipped: false`, `violations: []`, because
  `config/rotor_design.yaml` carries the box (`chord_min_m = 0.045`,
  `chord_max_m = 0.30`, `twist_min_deg = −2`, `twist_max_deg = 35`). The
  bounds question is closed as of 2026-09-19. Context that made it a real
  question at the time: the blade runs **263 mm chord at the root to 69 mm at
  the tip**, twist **23.07° → 0.57°**, and a 263 mm root chord on a 2.0 m
  blade is wide — characteristic of Schmitz. It turns out to sit inside the
  0.30 m cap with 36 mm to spare. Had it not, `DesignBounds.clip_physical`
  would have clipped it and `feasibility.violations` would say so rather than
  leaving the caller to reconstruct the violation.
- *"Above-rated operating line not swept."* It is specified now: the 300 rpm
  ceiling in `config/rotor_design.yaml` is the operating law, and the sweep
  runs 3.0 → 11.0 m/s in 0.5 m/s steps (17 points, all converged) with
  λ = 6.5 up to V_c = 9.67 m/s and 300 rpm above it. Three knots of the sweep
  sit on the ceiling (10.0, 10.5, 11.0 m/s at λ = 6.28, 5.98, 5.71) — visible
  in `baseline_operating_line.png` as the break in the λ curve.

**`baseline_reference.json`'s `outstanding.above_rated_operating_line` string
is still the old text** ("… is not yet specified"). It is a hard-coded
sentence in `generate_baseline.py`, it is now inaccurate, and it was left
alone deliberately on 2026-09-19: rewriting it would edit a published
artefact's provenance field in a re-run whose subject is the artefacts'
numbers, and the correction belongs in the same commit as whatever fixes the
generator's own idea of "outstanding". Recorded here and in the journal rather
than silently dropped.

## Note on the interpolant fix

The first version of this baseline had a station fail to converge at
v = 3 m/s, Re = 70,995 — which turned out to be a genuine value discontinuity
in the polar interpolant at alpha knots, for any Reynolds number between cache
rows. That is fixed (see `tests/golden/README.md`, 2026-09-10). The numbers
above are post-fix. Every operating point and every λ in the sweep now
converges, which `test_every_operating_point_converges` asserts.
