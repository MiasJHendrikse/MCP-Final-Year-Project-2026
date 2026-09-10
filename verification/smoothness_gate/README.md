# Objective smoothness gate

Plan step 1.8 — the readiness test for all gradient work, and the step the
Phase 1 brief calls the single most valuable in Phase 1.

```
python verification/smoothness_gate/run_gate.py [--points N]
```
→ `smoothness_gate.json`, `smoothness_gate.png` (~11 min at 300 points)

**Verdict: PASSED, with one documented non-smoothness.** `J` is **C¹ but not
C²**. That is sufficient for the adjoint and for finite-difference
verification, which is what Phase 2 and Phase 3 need. The C² defect is inherent
to the Buhl correction rather than an implementation error, and its expected
effect is recorded below.

## What was run

10 design variables × 300 points = **3000 objective evaluations**, 10.6 minutes.
Each `J` is 17 wind-speed bins, so ~51,000 rotor solves.

That this is feasible at all is the Task 5 cost gate: the audit costed this
sweep at **~5.5 days per pass** at 9.0 s per operating point.

The objective swept is `objective.energy_surrogate` — per-bin power with **unit
weights** rather than Weibull weights, because the wind resource is still
`TODO`. The Weibull weights are a fixed convex combination over the bins and do
not depend on `d`, so the surrogate exercises the identical d-dependent chain
and cannot create or remove a discontinuity. Sweep ranges are provisional
(chord ±40 % of its `x₀` value, twist ±8°) pending the manufacturability
bounds. See [`docs/OUTSTANDING-INPUTS.md`](../../docs/OUTSTANDING-INPUTS.md).

## Results

| variable | converged | distinct | flat steps | max Δ/med | max Δ²/med |
|---|---|---|---|---|---|
| chord CP 0 | ✓ | 1.000 | 0 | 3.47 | 2.95 |
| chord CP 1 | ✓ | 1.000 | 0 | 2.38 | 1.39 |
| chord CP 2 | ✓ | 1.000 | 0 | 2.01 | 1.08 |
| chord CP 3 | ✓ | 1.000 | 0 | 1.73 | 1.19 |
| chord CP 4 | ✓ | 1.000 | 0 | 1.58 | **5.42** |
| twist CP 0 | ✓ | 1.000 | 0 | 2.19 | 1.28 |
| twist CP 1 | ✓ | 1.000 | 0 | 1.79 | 1.07 |
| twist CP 2 | ✓ | 1.000 | 0 | 1.68 | 1.11 |
| twist CP 3 | ✓ | 1.000 | 0 | 1.71 | 1.27 |
| twist CP 4 | ✓ | 1.000 | 0 | 2.49 | **2.91** |

Against the four defects the plan names:

- **Staircasing — absent.** `distinct = 1.000` means every one of 300 `J`
  values is unique on every sweep, and there is not a single zero step in any
  first difference. A piecewise-constant interpolant anywhere in the chain
  would show here immediately. (This is the check that would have failed
  before the Task 3 interpolant fix.)
- **High-frequency noise — absent.** The first differences are smooth monotone
  curves; `max Δ/median` between 1.6 and 3.5 is ordinary curvature, and every
  sweep is a clean unimodal `J`.
- **Discontinuous jumps — absent.** All 3000 evaluations converged. A jump
  would give a `max Δ/median` in the hundreds.
- **Kinks — present, in two variables, and explained below.**

## The kink: stations crossing the Buhl blend

`chord CP 4` and `twist CP 4` — both **tip** control points — show visible
corners in their first-difference plots and elevated `max Δ²/median`.

Traced to the cause. As the tip control point moves, stations cross the
turbulent-wake threshold `a = 0.4` one at a time and switch onto Buhl's
empirical branch:

| chord CP 4 | stations with a > 0.4 |
|---|---|
| 40.2 – 45.8 mm | 0 |
| **47.6 mm** | 1 ← crossing |
| **66.1 mm** | 2 ← crossing |
| **80.9 mm** | 3 ← crossing |

Refining to 81 points across the first crossing settles what kind of
non-smoothness it is:

```
d = 46.625 mm   dJ/dd = +2658.4
d = 46.675 mm   dJ/dd = +2571.1     steps of -87
d = 46.725 mm   dJ/dd = +2483.8
d = 46.775 mm   dJ/dd = +2396.3
d = 46.825 mm   dJ/dd = +2308.8
d = 46.875 mm   dJ/dd = +2244.1     steps of -33
d = 46.925 mm   dJ/dd = +2211.3
d = 46.975 mm   dJ/dd = +2180.3
```

The derivative is **continuous** — no jump, relative discontinuity 0.046 over
one step, which is ordinary variation. What changes is the *slope* of the
derivative, from −87 to −33 per step. That is a **C² discontinuity, not a C¹
one**.

**Why it is there, and why it is not a bug.** Buhl's correction is constructed
to match momentum theory in **value and slope** at `a = 0.4` — C⁰ and C¹, and
nothing further. `tests/test_corrections.py` asserts both of those exactly.
Curvature is not matched and was never claimed to be. A C²-continuous blend
would require a different empirical curve (quintic rather than quadratic) and
would no longer be Buhl's, so "fixing" it would mean departing from the cited
method to remove something that does not obstruct the work.

**Expected effect on gradient work**, which is what the exit criterion asks to
be recorded:

- **Adjoint (Phase 3): none.** The adjoint needs `J` to be C¹. It is. The
  derivative is well-defined and continuous everywhere in the swept range,
  including across the crossings.
- **Finite-difference verification (Phase 2): none at usable step sizes.** FD
  converges for a C¹ function. A step straddling a crossing picks up a second-
  order error, which is the same order as FD truncation error anyway.
- **Quasi-Newton (SLSQP's BFGS update): mild.** The Hessian approximation sees
  a curvature jump when a station crosses. Expect slightly degraded
  superlinear convergence near such a crossing, not incorrect results. Worth
  remembering if Phase 5 reports slow convergence in a region where stations
  are entering the turbulent-wake state.

## Why the tip control points and not the others

The tip is where the Prandtl loss factor `F` is smallest, so the momentum-side
induction is highest and stations are closest to the `a = 0.4` threshold to
begin with. The inboard control points move stations that are nowhere near it,
which is why `max Δ²/median` sits at 1.07–1.39 for the middle of the blade and
rises only at the two tip variables.

## What this gate has already caught

Before it was even built. While constructing the Schmitz baseline (step 1.7), a
station failed to converge at Re = 70,995 — which turned out to be a genuine
value discontinuity in the polar interpolant at every α knot, for any Reynolds
number between cache rows. That is precisely a *staircasing/jump* defect of the
kind this gate exists to find, and it is fixed (see
`tests/golden/README.md`, 2026-09-10). The clean `distinct = 1.000` column
above is the post-fix state.

## Re-running against the real objective

Once the Weibull parameters land, swap `energy_surrogate` for
`annual_energy_mwh` and re-run. Roughly 13 minutes. The conclusions above
should be unchanged — the weights do not depend on `d` — and if they are not,
that is itself a finding worth chasing.
