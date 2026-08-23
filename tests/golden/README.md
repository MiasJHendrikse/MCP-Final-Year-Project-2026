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
`bem.airfoil.S809Polar`'s α-clamp is live there. That is deliberate: **Task 4
deletes that clamp**, and this point is where the deletion will show up.
It should change visibly, not quietly.

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
