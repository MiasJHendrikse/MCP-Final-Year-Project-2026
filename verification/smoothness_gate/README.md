# Is the objective smooth enough to differentiate?

This is the check that had to pass before any gradient work: is `J` smooth enough to
differentiate? For each design variable in turn, hold the others at `x0` and
sweep that one; plot `J` against it and the **first difference**, which is
where the defects actually show.

```
python verification/smoothness_gate/run_gate.py [--points N]
```
→ `smoothness_gate.json`, `smoothness_gate_objective.png`,
`smoothness_gate_difference.png` (~12.6 min at 300 points)

**Verdict: PASSED, with localised C² defects.** `J` is **C¹ but not C²**.
That is sufficient for the adjoint and for finite-difference verification,
which is what the gradient work needs. The defects are single-step changes of
slope, never jumps in value, and they trace to two structural features of the
model (below), not to an implementation error.

## What was run

10 design variables × 300 points = **3000 objective evaluations**, 12.6
minutes. Each `J` is 17 wind-speed bins, so ~51,000 rotor solves.

The objective is `objective.annual_energy_mwh`: the Weibull-weighted AEP in
MWh/yr, under the committed operating law (fixed tip-speed-ratio tracking, the
300 rpm rotor-speed ceiling, the configured rating). This is the functional
the gradients differentiate.

## The sweep range

The stated range is the **design box** of `config/rotor_design.yaml`
(`parameterisation.bounds`, read with `DesignBounds.from_config`): chord
0.045–0.30 m, twist −2 to 35 degrees. Four of the ten variables cannot be
swept across the whole box: below some chord a station Reynolds number falls
under the polar cache's 40 000 floor and `J` is undefined (`PolarDomainError`),
and above some twist the BEM solver no longer converges. Each variable is
therefore swept across the largest interval containing its `x0` value on which
`J` is defined **and** every station converges, found by bisection; the box,
the swept interval and the fraction are recorded per sweep in the JSON.

| variable | box | swept | fraction |
|---|---|---|---|
| chord CP 0 | 45.0 – 300.0 mm | 140.2 – 300.0 mm | 63 % |
| chord CP 1 | 45.0 – 300.0 mm | 45.0 – 300.0 mm | 100 % |
| chord CP 2 | 45.0 – 300.0 mm | 45.0 – 300.0 mm | 100 % |
| chord CP 3 | 45.0 – 300.0 mm | 45.0 – 300.0 mm | 100 % |
| chord CP 4 | 45.0 – 300.0 mm | 45.0 – 300.0 mm | 100 % |
| twist CP 0 | −2.0 – 35.0° | −2.0 – 35.0° | 100 % |
| twist CP 1 | −2.0 – 35.0° | −2.0 – 35.0° | 100 % |
| twist CP 2 | −2.0 – 35.0° | −2.0 – 28.5° | 83 % |
| twist CP 3 | −2.0 – 35.0° | −2.0 – 19.4° | 58 % |
| twist CP 4 | −2.0 – 35.0° | −2.0 – 12.2° | 38 % |

Two consequences worth stating. The cache floor is what makes the low-chord
corner of the box unreachable, so the polar cache is part of the design
problem's boundary rather than only a numerical table. And the high-twist
corner is unreachable because the solver will not converge there, which is a
statement about the model's validity range, not about the box.

## Results

| variable | converged | undefined | distinct | flat steps | max Δ/med | max Δ²/med |
|---|---|---|---|---|---|---|
| chord CP 0 | ✓ | 0 | 1.000 | 0 | 2.98 | 2.16 |
| chord CP 1 | ✓ | 0 | 1.000 | 0 | 4.16 | 35.37 |
| chord CP 2 | ✓ | 0 | 1.000 | 0 | 2.37 | 198.83 |
| chord CP 3 | ✓ | 0 | 1.000 | 0 | 2.19 | 436.95 |
| chord CP 4 | ✓ | 0 | 1.000 | 0 | 1.43 | 12.31 |
| twist CP 0 | ✓ | 0 | 1.000 | 0 | 8.79 | 397.84 |
| twist CP 1 | ✓ | 0 | 1.000 | 0 | 4.15 | 17.86 |
| twist CP 2 | ✓ | 0 | 1.000 | 0 | 3.06 | 12.54 |
| twist CP 3 | ✓ | 0 | 1.000 | 0 | 3.00 | 15.65 |
| twist CP 4 | ✓ | 0 | 1.000 | 0 | 5.33 | 3.58 |

Checking for the four kinds of defect:

- **Staircasing — absent.** `distinct = 1.000` means every one of the 300 `J`
  values is unique on every sweep, and there is not a single zero step in any
  first difference. A piecewise-constant interpolant anywhere in the chain
  would show here immediately.
- **High-frequency noise — absent.** The first differences are smooth curves;
  `max Δ/median` between 1.4 and 8.8 is ordinary curvature.
- **Discontinuous jumps — absent.** All 3000 evaluations converged. A jump
  would give a `max Δ/median` in the hundreds.
- **Kinks — present, on several variables, and explained below.**

## Where the kinks are, and what kind they are

`max Δ²/median` is the detector. The largest localised ratios are:

| variable | at | ratio | what the first difference does |
|---|---|---|---|
| chord CP 3 | 53.5 mm | 437 | continuous, slope changes |
| twist CP 0 | 0.35 – 0.48° | 398 and 396 | continuous, slope changes sign |
| chord CP 2 | 62.1 mm | 199 | continuous, slope changes |
| chord CP 1 | 118.3 mm | 35 | continuous, slope changes |
| twist CP 1 | 27.2° | 18 | continuous, slope changes |
| twist CP 3 | 9.9° | 16 | continuous, slope changes |
| twist CP 2 | 16.6° | 13 | continuous, slope changes |
| chord CP 4 | 46.7 mm | 12 | continuous, slope changes |
| twist CP 4 | 11.3° | 4 | continuous, slope changes |

In every case the first difference is **continuous** across the step — there
is no jump in the derivative — and only its slope changes. That is a **C²
discontinuity, not a C¹ one**, which is why the adjoint is unaffected.

The ratio is a detector, not a magnitude, and it is inflated where `J` is
nearly flat: twist CP 0's whole sweep spans 0.099 MWh/yr and its derivative
changes sign, so its median second difference is small and the ratio is large.
The table should be read with that in mind.

## Why they are there, and why they are not bugs

Two structural features of the model are C¹ only:

- **Buhl's correction**, constructed to match momentum theory in **value and
  slope** at `a = 0.4` — C⁰ and C¹, and nothing further.
  `tests/test_corrections.py` asserts both of those exactly. Curvature is not
  matched and was never claimed to be.
- **The polar interpolant**, which has a C² break at every 0.5° angle-of-attack
  knot and every Reynolds row (recorded in `run_sweep.py`'s docstring: the
  interpolant is C¹ with these knots).

Both are properties of the cited methods, not of the implementation. A
C²-continuous blend would require a different empirical curve (quintic rather
than quadratic) and would no longer be Buhl's, so "fixing" it would mean
departing from the cited method to remove something that does not obstruct the
work.

**Expected effect on the gradient work:**

- **Adjoint: none.** The adjoint needs `J` to be C¹. It is. The derivative is
  well-defined and continuous everywhere in the swept range.
- **Finite-difference verification: none at usable step sizes.** FD converges
  for a C¹ function. A step straddling a break picks up a second-order error,
  the same order as FD truncation error anyway.
- **Quasi-Newton (SLSQP's BFGS update): mild.** The Hessian approximation sees
  a curvature jump as a station crosses a break. Expect slightly degraded
  superlinear convergence near such a crossing, not incorrect results.

## What this check has already caught

A station failed to converge at Re = 70,995 while the Schmitz baseline was
being built — a genuine value discontinuity in the polar interpolant at every
α knot, for any Reynolds number between cache rows. That is precisely a
*staircasing/jump* defect of the kind this gate exists to find, and it is
fixed (see `tests/golden/README.md`): the clean `distinct = 1.000` column above
is the post-fix state.

## The objective the gate sweeps

The gate sweeps the committed objective, `objective.annual_energy_mwh`. An
earlier pass of it swept `objective.energy_surrogate`, the same per-bin powers
summed with **unit weights**, while the wind resource was unresolved — the
argument being that the Weibull weights are a fixed convex combination over
the bins, independent of `d`, so both forms exercise the identical
`d`-dependent chain and a positive fixed weighting can neither create nor
remove a discontinuity.

That argument was tested rather than merely asserted, sweep by sweep: every
sweep agreed (3000/3000 converged, `distinct = 1.000` throughout, zero flat
first-difference steps, the difference ratios tracking within a few percent
with the same variables elevated). It was a real check, not a formality: had
the verdict changed, the weighting would have been coupled to `d` somewhere it
should not be, and that is a defect worth finding before any gradient work.

In `run_gate.py` the weights are computed **once, outside the sweep loop**,
which makes their independence from `d` structural rather than a claim in a
docstring. The same single `bin_powers` call feeds both the objective value and
the convergence flag, so the rotor is not solved twice.

## The figures

`smoothness_gate_objective.png` shows `J` along each chord sweep and along each
twist sweep, each curve normalised to its own range so that ten curves of
different magnitude can be read on one axis. `smoothness_gate_difference.png`
shows the first difference of the same sweeps, in MWh/yr per unit physical
variable, which is where the C² defects appear. Both are written by
`src/plotting/figstyle.py`.
