# Golden-file regression snapshot

Task 0 of `docs/repo-audit-instructions-REVISED.md`. These files record the
BEM solver's behaviour **exactly as it stood before any remediation work
began** — commit `77653fa`, branch `task-0-golden-regression`, captured
2026-08-23.

They are the regression net under the whole work order. Tasks 1, 4 and 5 each
carry "Phase VI cross-check figures unchanged" as an acceptance criterion, and
the only previously recorded figures were two summary percentages. A refactor
can shift spanwise `a`, `a'` and `φ` materially while barely moving a mean, so
"unchanged" needs to mean *checked to a tolerance*.

## Contents

| File | What it holds |
|---|---|
| `phase_vi_cp_lambda.json` | `Cp(λ)` and `Ct(λ)` for the NREL Phase VI rotor, λ = 1.5 … 7.5 in steps of 0.5, at `v_inf` = 7.0 m/s, ρ = 1.225 kg/m³. 13 points spanning the cross-validated λ range. |
| `phase_vi_spanwise.json` | Per-station `r`, `φ`, `a`, `a'`, `α`, `Cl`, `Cd`, `Cn`, `Ct`, `Cq`, `F`, `Re`, `W` and the local loads `dT/dr`, `dQ/dr`, at four operating points on the Sequence S fixed-RPM line (71.63 rpm). |
| `cross_tool_summary.json` | Our solver's side of all three cross-tool comparisons, plus the per-station Reynolds bucket assignment. |

### The four spanwise operating points

Chosen to span the envelope, per Task 0's requirement for at least one
attached-flow and one near-stall condition:

| `v_inf` | λ | α range | Regime |
|---|---|---|---|
| 5.0 m/s | 7.545 | 1.9 – 4.3° | Fully attached |
| 7.0 m/s | 5.389 | 3.4 – 10.2° | Attached, stall approaching inboard |
| 10.0 m/s | 3.772 | 5.9 – 19.7° | Stall onset / near-stall |
| 15.0 m/s | 2.515 | 10.6 – 31.9° | Deep stall |

The 15 m/s point is past the S809 cache's +18° upper bound, so
`bem.airfoil.S809Polar`'s α-clamp was live there. That was deliberate: **Task 4
deleted that clamp**, and this point is where the deletion showed up — Cp moved
−25.2 %, with 16 of 19 stations previously clamped. See the 2026-09-10 entry in
the change log below.

## Tolerance

`rtol = 1e-10`, with a per-field `atol = 1e-12 × max|golden value in that
field|`. The reasoning is written out in full in the docstring of
`tests/test_golden_regression.py`; in short, the floor is set by `brentq`'s
own `xtol = rtol = 1e-12` on the root location φ, and the ceiling by the fact
that the smallest physically meaningful change (0.01 % of `Cp`) is 1e-4 —
eight orders of magnitude above the tolerance.

**Measured sensitivity.** Perturbing every station's chord by a relative
1e-9 moves `Cp` by 9.2e-10 and is caught; 1e-10 moves it by 9.2e-11 and is
not. The response is linear with gain ≈ 0.92 across 1e-12 … 1e-8, so the
effective detection threshold sits at the stated tolerance with no dead zone.
Verified end-to-end: a 1e-6 relative change to one station's chord fails four
of the five tests with messages naming the station and quantity.

## Important: the CCBlade and pyBEMT reference values are stale

`cross_tool_summary.json` records mean |Cp| deviations of **13.32 %**
(CCBlade) and **15.79 %** (pyBEMT) against the stored external values — not
the **0.51 %** and **1.83 %** published in
`docs/validation/bem-cross-validation.md`.

**This is not a solver regression.** It is a stale-reference problem, found
while capturing this snapshot:

- `docs/validation/ccblade_case/results.json` and `pybemt_case/results.json`
  were committed in `05a6f45` (2026-07-25).
- The S809 cache was rebuilt at `Ncrit = 5` in `4593b11` and extended to 1.3M
  in `5b5574f` (2026-07-28). Neither is an ancestor of `05a6f45`.
- At that time `compare_pybemt.RE_BUCKETS` stopped at 500k, so **every station
  clamped to the 500k ceiling** — the exact defect `build_polar_cache.py`'s own
  comment records. The bucket list now runs to 1.3M and the stations resolve to
  600k–900k.

So the external columns in those files were produced by CCBlade and pyBEMT fed
from polar tables the repo no longer contains. Both sides of that comparison
are stale; only our side has been updated. Re-establishing the published
figures requires **re-running the external tools against the current cache**,
which needs their separate virtualenvs and is outside Task 0's scope.

The two percentages are pinned here as **tracking values**, so a refactor
cannot move them unnoticed. They are not validation results and must not be
quoted as such.

The third comparison is unaffected: `data/qblade/bem_solver_phase_vi_result_Re*.csv`
postdate the cache rebuild and reproduce **bitwise** from today's tree.
`test_qblade_pinned_matches_committed_csv` asserts that exactly.

### Consequence for the work order

Tasks 1, 4 and 5 list "the Phase VI cross-checks reproduce their existing
figures — CCBlade 0.51 %, pyBEMT 1.83 %" as an acceptance criterion. That
criterion is **not satisfiable as written** from the current tree, by anything
those tasks do. Until the external tools are re-run, the golden files here are
the operative check for "the Phase VI cross-check figures are unchanged".

## Change log

Every regeneration, what moved, and why. A task that changes these files
appears here.

### 2026-08-29 — Task 2 (gap-free cache, ±180° extension)

**What moved:** 12 of 41 `phase_vi_cp_lambda` values (λ = 3.5 … 6.0) and 81 of
1160 `phase_vi_spanwise` values, all at relative 1e-5 … 3e-4.

**Why:** Task 2 closed the eight interior holes in the S809 cache. Before, a
query landing on a hole read a linear bridge between the converged neighbours
either side; now it reads the documented fill. The spanwise movement is
entirely at the **v = 7 m/s** point, stations **12–14** (α = 7.17–7.66°,
Re = 938k–945k) — the three whose lookup straddles α = 7.5° between the 900k
and 1.0M curves, which is the (1.0M, 7.5°) hole. Station 11 (α = 8.18°) and
station 15 (α = 6.64°) do not straddle it and did not move. The difference at
the filled node is ~0.001 in Cl — the same figure `validate_polars` check 1 had
already measured there.

**What did not move:** `cross_tool_summary.json` in full, including both
tracking percentages and the QBlade pinned per-station CSVs, and
`test_snapshot_is_bitwise_reproducible`.

**Not yet moved by Task 2:** the 15 m/s deep-stall point. The cache now extends
to ±180°, but `S809Polar`'s α clamp is deliberately pinned to the
XFOIL-converged band, so the clamp path is still live there. Task 4 deletes it;
that is when this point should change.

### 2026-09-10 — Task 4 (polar adapter replaced; clamping and substitution removed)

**What moved:** 26 of 41 `phase_vi_cp_lambda` values, 992 of 1160
`phase_vi_spanwise` values, and 420 values in `cross_tool_summary` — all of the
latter inside `qblade_pinned`.

**What did not move, and why that matters.** The `sweep` block —
our side of the CCBlade and pyBEMT comparisons — is **bitwise unchanged**, all
28 values, and both tracking percentages are identical to every digit
(13.324537038059685 % and 15.787123761875879 %). So is
`reynolds_buckets_per_station`. Those paths run on `TableS809Polar`, a
deliberately discretised table over the already-±180° Viterna export, which
Task 4 does not touch. That is the control: it shows the movement below comes
from the adapter being replaced, not from anything shifting underneath the
whole snapshot.

**Why, in two separable parts.** Both causes are present, and they are
distinguishable by where they appear:

*1. The interpolant (bilinear → cubic-in-α, PCHIP-in-log(Re)).* Small,
everywhere, and the only cause wherever the whole span is attached:

| point | Cp before | Cp after | |
|---|---|---|---|
| `cp_lambda`, λ = 4.0 … 7.5 | — | — | **+0.007 % … +0.099 %** |
| spanwise v = 5.0 m/s (α ≤ 4.3°) | 0.406172 | 0.406212 | +0.010 % |
| spanwise v = 7.0 m/s (α ≤ 10.2°) | 0.367031 | 0.367103 | +0.020 % |

The largest single interpolant-only move in the whole snapshot is **+4.1 %**,
on `qblade_pinned[Re=100k].Cd` at station 9 — Cd is convex in α, and linear
interpolation between 0.5° knots systematically over-reads it.

*2. The α clamp deleted.* Large, and confined to points with stations past the
+18° XFOIL-converged band. This is the change `tests/golden_reference.py` and
this file both said in writing should happen here:

| point | stations past 18° | Cp before | Cp after | |
|---|---|---|---|---|
| spanwise v = 10.0 m/s | 7 / 19 | 0.258726 | 0.256597 | −0.823 % |
| spanwise v = 15.0 m/s | 16 / 19 | 0.150104 | 0.112335 | **−25.16 %** |
| `cp_lambda` λ = 3.5 | — | 0.220746 | 0.212980 | −3.518 % |
| `cp_lambda` λ = 2.5 | — | 0.128481 | 0.093871 | −26.94 % |
| `cp_lambda` λ = 1.5 | — | 0.060918 | 0.023875 | **−60.81 %** |

The mechanism is clearest in drag. At v = 15 m/s the six innermost stations sit
at α = 25–32°; the old adapter returned `Cd(18°) ≈ 0.079` for every one of
them, and the ±180° cache gives 0.34–0.43. That is a **+340 % to +449 %**
correction on the largest `Cd` values in the snapshot, and it is what drives
Cp down. Clamped lift does not fall off past stall and clamped drag does not
rise, so every previously-clamped point was reporting too much power.

**None of these numbers are a regression.** The old values were computed from
lift and drag the airfoil was not producing at those angles. The new ones are
read from the Viterna extrapolation that Task 2 built, verified and recorded —
the same extrapolation that was exported to the `.plr` QBlade has been running
against all along, so our solver and QBlade now share a post-stall model as
well as a measured band.

**Also re-anchored in the same commit,** for the same cause and with the same
before/after recorded at the site:
`validation/validate_powercurve.PHASE_VI_SEQUENCE_S_REFERENCE`, and
`data/qblade/bem_solver_phase_vi_result_Re{100000,500000}.csv`.

**Coverage limits this exposed.** Removing the Re clamp made two sweeps stop
where the S809 cache actually stops, rather than running on a clamped Reynolds
number: the Cp-λ sweep in `validate_powercurve` check 4 (λ ≤ 7.5; λ = 8.0 puts
the tip station at Re = 1,384,108 against a 1,300,000 ceiling) and its
fixed-speed sweep in check 2 (≤ 20 m/s; at 25 m/s the *root* station, which has
the largest chord and so runs out of cache first at fixed RPM, reaches
1,312,857). Neither is a Task 4 defect — both are the cache's real coverage,
previously hidden.

**What did not move:** `test_snapshot_is_bitwise_reproducible` still passes,
and the four fixed operating conditions (wind speeds, λ values, ρ, ν, RPM) are
unchanged, so this is a like-for-like comparison.

### 2026-09-10 — Task 6 (Ning γ-form): the golden JSON files were **not** regenerated

Recorded because the absence of a change is itself the result.

Task 6 rewrote `corrections.corrected_axial_induction` in Ning's γ₁γ₂γ₃ form.
That is the same equation reparameterised, not a different model, so the only
movement it can produce is floating-point reassociation. Measured across the
whole snapshot:

| file | values moved | largest move |
|---|---|---|
| `phase_vi_cp_lambda.json` | 4 of 41 (all `Ct`) | 2.3e-15 relative |
| `phase_vi_spanwise.json` | 63 of 1160 | **2.4e-15 relative** |
| `cross_tool_summary.json` | 19 of 868 (all `qblade_pinned`) | 1.9e-15 relative |
| `cross_tool_summary.sweep` | **0** | — |

The largest move anywhere is 2.4e-15, five orders of magnitude below the
`rtol = 1e-10` these files are compared at. All three JSON regression tests
**pass unchanged**, so per this file's own rule — regeneration happens only
when the test fails, and never as a routine step — the golden files were left
exactly as Task 4 wrote them.

**What did have to be regenerated:**
`data/qblade/bem_solver_phase_vi_result_Re{100000,500000}.csv`. Those are
compared **bitwise** by `test_qblade_pinned_matches_committed_csv`, not to a
tolerance, so a 4-ulp change breaks them by design. One row moved in each file
(station 18 at Re = 100k: `phi_deg` 5.010656781762976 → 5.010656781762972).

The distinction is worth keeping straight when reading the next task's diff: a
tolerance test staying green and a bitwise test going red are consistent
outcomes for the same change, and here they bracket the size of it.

### 2026-09-10 — Task 5 (Ning bracket and residual): golden JSONs again **not** regenerated

Same outcome as Task 6, for a different reason, and with one instructive
difference.

Task 5 replaced the residual with Ning's form and the 2000-point bracket scan
with a region-classified bracket. Both changes alter *how* the root is located,
not *which* root it is, so the movement is convergence noise around the same
solution rather than a different answer:

| file | values moved | largest move |
|---|---|---|
| `phase_vi_cp_lambda.json` | 25 of 41 | 1.1e-13 relative |
| `phase_vi_spanwise.json` | 760 of 1160 | **8.0e-12 relative** |
| `cross_tool_summary.json` | 347 of 868 | 1.2e-12 relative |

The largest move anywhere is 8.0e-12, still comfortably under `rtol = 1e-10`,
so all three JSON tests pass and the files were left alone. Only the bitwise
`data/qblade/bem_solver_phase_vi_result_Re*.csv` needed regenerating.

**The instructive difference:** `cross_tool_summary.sweep` **did** move this
time — 14 values, and both tracking percentages shifted in their last three
digits (13.324537038059685 → 13.324537038062505). In Tasks 4 and 6 that block
was bitwise frozen and served as the control. It is not frozen here, and that
is correct: `TableS809Polar` pins the *polar* those runs see, not the *solver*.
Task 4 and Task 6 changed the polar layer, so the pinned-table paths were
insulated from them; Task 5 changed the solver itself, which nothing insulates.
A reader comparing the three change-log entries should expect exactly this, and
its absence would have been the thing to investigate.

**What did not move:** `reynolds_buckets_per_station`, the fixed operating
conditions, and `test_snapshot_is_bitwise_reproducible`.

## Regenerating

```
python tests/generate_golden.py
```

Regeneration is **not** a routine step. Per the work order's ground rule 6:

> The golden regression from Task 0 is the net under all of this. It runs after
> every task. Any change to it is deliberate, separately committed, and
> justified in the commit message — never absorbed silently into a code change.

So: commit the code change, run the test, and only if it fails for a reason you
can state does the regeneration happen — as a **second commit** whose message
names the task and the reason.

Real behaviour changes are expected from Tasks 2, 3, 4, 5 and 6 (the polar
layer and the residual form). Task 1 is a pure rewiring of ρ and ν through a
config layer and should be numerically inert — if it moves these numbers, the
rewiring is wrong.

## Running

```
python -m pip install pytest      # or: pip install -e ".[test]"
pytest tests/test_golden_regression.py
```

from the repo root. `pyproject.toml` puts `src/` and `tests/` on the path, so
no install step and no `cd src` is needed. Takes about 3½ minutes on unmodified code: 17 Phase VI
operating points at the ~9 s each the audit measured, plus the fast pinned-table
sweeps. Task 5's cost target (under 0.1 s per operating point) brings this down
to seconds.

### 2026-09-10 — interpolant knot-continuity fix (found by plan step 1.7)

The first regeneration in this log that is a **correctness fix rather than a
refactor**, so unlike the Task 5 and Task 6 entries the golden files did have
to move.

**What moved:** 26 of 41 `phase_vi_cp_lambda`, 992 of 1160 `phase_vi_spanwise`,
362 of 868 `cross_tool_summary` (including 14 in `sweep`), and both QBlade CSVs.

**Largest moves:** `Cd` up to **2.6e-4** relative, `Cp` up to **4.2e-5**. Above
the `rtol = 1e-10` these files are compared at, hence a real regeneration.

**Why.** The polar interpolant had a value discontinuity at every alpha knot,
for any Reynolds number *between* cache rows — Cl jumping up to 5.4e-3 on
SG6043. Blending cubic-spline *coefficients* across log(Re) with PCHIP does not
preserve continuity at the alpha knots, because continuity there is a linear
constraint on the coefficient channels and PCHIP is nonlinear in its data. The
old golden numbers were computed on a polar with those jumps in it; the new
ones are not.

It surfaced as a BEM station on the Schmitz baseline failing to converge at
Re = 70,995 with a genuinely discontinuous residual — which is Task 5's
reported-non-convergence doing exactly the job it was built for. Under the
Phase 0 solver this would have been a `ValueError` aborting the sweep, or worse,
a quietly wrong number.

**Direction of the change:** these values are *more* correct, not merely
different. On a cache row at a cache node the surface now reproduces the cached
number bitwise, so the underlying data is unchanged — only the interpolation
between data points moved, and it moved onto a surface that is actually C¹.
