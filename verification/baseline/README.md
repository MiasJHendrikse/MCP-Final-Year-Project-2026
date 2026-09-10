# Schmitz baseline — `x₀` and the reference numbers

Plan step 1.7. The plan describes the evaluated performance of this blade as
**"the reference numbers for the entire results chapter"**, so this directory
is a fixed artefact: `x0.json` is what Phases 2–5 start from, and
`baseline_reference.json` is what every later result is compared against.

```
python verification/baseline/generate_baseline.py
```
→ `x0.json`, `baseline_reference.json`, `baseline.png`

Regeneration is deliberate, not routine — the same rule as the golden files. If
these numbers move, something in the solver or the polar layer moved, and the
commit message has to say what and why.

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
| Cp at design point (λ = 6.5, V = 11 m/s) | **0.4720** |
| Ct at design point | **0.8111** |
| peak Cp | **0.4720 at λ = 6.5** |
| root bending moment | **206.7 N·m per blade** |
| peak thrust | **597.2 N at 11 m/s** |
| AEP | **outstanding — see below** |

**The peak Cp landing exactly at the design λ is the strongest available check
on the construction.** The blade is designed for λ = 6.5; if the Schmitz
formulae, the twist convention, or the solver wiring were wrong, the peak would
sit somewhere else. `tests/test_baseline.py` asserts it.

Cp = 0.472 against the Betz limit of 0.593 is what a real blade with drag
achieves — the ~20 % shortfall is profile drag plus tip loss, not an error.

Root bending moment is **per blade**, not per rotor: the root attachment
carries one blade's load, and a rotor-summed figure would overstate it by a
factor of B.

## What is outstanding

**AEP is absent, not estimated.** It needs the Weibull parameters, which are
`TODO` (see [`docs/OUTSTANDING-INPUTS.md`](../../docs/OUTSTANDING-INPUTS.md)
§1). `baseline_reference.json` carries an `outstanding` block naming it as
blocked. An estimated AEP sitting in the reference numbers would be
indistinguishable from a real one to anyone reading this later, which is
exactly the failure ground rule 3 exists to prevent.

**Feasibility was not checked**, and the artefact says so — `feasibility.checked`
is `false` with the reason attached. The bounds are `TODO`
(`OUTSTANDING-INPUTS.md` §2). A baseline that was never checked must not read
as one that passed.

Context for whoever sets those bounds: this blade runs **263 mm chord at the
root to 69 mm at the tip**, twist **23.07° → 0.57°**. A 263 mm root chord on a
2.0 m blade is wide — characteristic of Schmitz — so whether the baseline is
feasible is a genuinely open question, not a formality. If it is not, step 1.7
requires it to be clipped and the violation *recorded*, which
`DesignBounds.clip_physical` returns rather than leaving the caller to
reconstruct.

**Above-rated operating line not swept.** Power limiting above 11 m/s is part
of the AEP model (step 1.5) and is not yet specified, so the operating line
stops at rated rather than extrapolating a strategy that does not exist.

## Note on the interpolant fix

The first version of this baseline had a station fail to converge at
v = 3 m/s, Re = 70,995 — which turned out to be a genuine value discontinuity
in the polar interpolant at alpha knots, for any Reynolds number between cache
rows. That is fixed (see `tests/golden/README.md`, 2026-09-10). The numbers
above are post-fix. Every operating point and every λ in the sweep now
converges, which `test_every_operating_point_converges` asserts.
