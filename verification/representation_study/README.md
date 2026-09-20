# Blade representation study

Plan section 4.3 and step 1.6's final bullet: the control-point count is to be
*justified* by a documented comparison across 6, 8, 10, 12 and 16, not
asserted.

**Conclusion: 10 total — 5 chord + 5 twist.** That is the plan's intended
configuration; this study is what makes it a justified choice rather than a
stated one. The reasoning is below, including why the count that fits best is
not the one selected.

```
python verification/representation_study/run_study.py
```
→ `representation_study.json`, `representation_study.png`

**Not re-run in the 2026-09-19 sweep of the optimisation artefacts.** The
study is a geometric comparison — how well each control-point count holds the
analytic Schmitz chord and twist — evaluated at the design point in
`verification/baseline/`, and that blade has not changed. The bounds and the
operating law do not enter it. Its committed output is from 2026-09-10.

## Method

The reference is the **analytic Schmitz blade** for this rotor (R = 2.0 m,
B = 3, λ = 6.5) at the SG6043 maximum-L/D point — L/D = 98.5 at α = 5.36°,
Cl = 1.2825, at Re = 200,000. That reference is the honest target: it is the
shape the optimiser is started from, so what matters is whether the
parameterisation can hold *it*, not whether it can hold a generic taper chosen
to be easy to fit.

The reference blade runs from 263 mm chord at the root to 69 mm at the tip
(mean 132 mm), with twist 23.07° → 0.57°.

Because the parameterisation is linear in its control points, fitting is a
single linear least-squares solve — exact, no starting guess, no local minima.
The error reported is a property of the **basis**, not of an optimiser that
fitted it.

## Results

| total | chord/twist | degree | chord RMS | (% of mean) | twist RMS | cond(N) | sawtooth reproduced |
|---|---|---|---|---|---|---|---|
| 6 | 3 + 3 | 2 | 4.112 mm | 3.111 % | 0.8626° | 3.17 | 9.8 % |
| 8 | 4 + 4 | 3 | 0.642 mm | 0.486 % | 0.2926° | 5.98 | 9.8 % |
| **10** | **5 + 5** | **3** | **0.515 mm** | **0.389 %** | **0.1061°** | **5.59** | **14.5 %** |
| 12 | 6 + 6 | 3 | 0.427 mm | 0.323 % | 0.0301° | 5.43 | 15.6 % |
| 16 | 8 + 8 | 3 | 0.122 mm | 0.093 % | 0.0041° | 5.60 | 20.7 % |

"Sawtooth reproduced" is the fraction of a pure strip-to-strip oscillation the
basis can represent — 0 % means the null space of §4.2 is invisible to the
optimiser, 100 % is what per-station design variables would give.

## Why 10

**6 is excluded outright.** Three control points per block cannot carry a
cubic, so the degree drops to 2, and it shows: 3.1 % chord error, 0.86° of
twist error, and a visible miss at both the root and the tip in the left-hand
panel of the figure. It cannot represent the shape it would be asked to start
from.

**Fitting accuracy stops meaning anything after 8.** At 10 control points the
chord RMS error is **0.515 mm on a 132 mm mean chord**. No composite blade of
this size is manufactured to a tolerance anywhere near that, so the additional
accuracy from 12 (0.427 mm) or 16 (0.122 mm) buys nothing physical. The
error curve in the middle panel is steep from 6 to 8 and flat from 8 to 12,
which is the signature of a basis that has already captured the shape.

**Conditioning does not discriminate here**, which is worth saying plainly
because plan 4.3 lists it first. cond(N) sits between 5.4 and 6.0 across every
count from 8 upward — the basis functions are well separated at all of these
resolutions on 25 strips. The low value at 6 is an artefact of the reduced
degree, not a virtue. So the criterion the plan expected to be decisive is not,
and the decision rests on the other three.

**16 fits best and is still not chosen.** This is the one genuine trade in the
table, and it goes the other way for two reasons: it reproduces **20.7 %** of a
strip-to-strip oscillation against 14.5 % at 10, i.e. it hands the optimiser
appreciably more of the exact sawtooth freedom §4.2 says the parameterisation
exists to remove; and it costs 60 % more chain-rule work through the
parameterisation in every adjoint evaluation, for a fitting improvement that is
already far below manufacturing tolerance.

**10 sits at the knee.** Twist error at 10 (0.106°) is an order of magnitude
better than at 8 (0.293°) for one extra control point per block, while chord
error is essentially unchanged — twist is the more demanding of the two
distributions here, and 10 is where it becomes negligible.

## What this study does not settle

- **Demonstrable optimisation improvement**, the fifth criterion in plan 4.3.
  That needs the optimiser, which is Phase 5. If a converged optimum at 10
  control points turns out to be materially worse than at 12, this conclusion
  should be revisited — the study is cheap to re-run.
- **Ease of imposing bounds**, the third criterion. The bounds are still `TODO`
  (see `docs/OUTSTANDING-INPUTS.md`), and a manufacturability envelope that
  constrains root chord tightly might favour more control points inboard, or an
  uneven split between the chord and twist blocks rather than the even one
  assumed throughout here.
- **An uneven chord/twist split.** Every row splits the total evenly. Given
  that twist error dominates at low counts, `4 chord + 6 twist` might beat
  `5 + 5` at the same cost. Not explored, and worth a look before Phase 5.
