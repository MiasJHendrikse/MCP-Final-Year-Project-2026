# Golden-file regression snapshot

These files record the BEM solver's behaviour **exactly as it was before a
round of restructuring and correctness work** on the solver and polar layer
(captured at commit `a121b37`). They're the safety net under that work, and
under everything since.

The only earlier record of the Phase VI cross-checks was two summary
percentages. A refactor can change the spanwise `a`, `a'` and `φ` materially
while barely moving a mean, so "unchanged" has to mean *checked to a
tolerance*, station by station.

## Contents

| File | What it holds |
|---|---|
| `phase_vi_cp_lambda.json` | `Cp(λ)` and `Ct(λ)` for the NREL Phase VI rotor, λ = 1.5 … 7.5 in steps of 0.5, at `v_inf` = 7.0 m/s and ρ = 1.225 kg/m³: 13 points spanning the cross-validated range |
| `phase_vi_spanwise.json` | per-station `r`, `φ`, `a`, `a'`, `α`, `Cl`, `Cd`, `Cn`, `Ct`, `Cq`, `F`, `Re`, `W` and the local loads `dT/dr`, `dQ/dr`, at four operating points on the Sequence S fixed-speed line (71.63 rpm) |
| `cross_tool_summary.json` | this solver's side of all three cross-tool comparisons, plus the Reynolds-number bucket each station is assigned to |

### The four spanwise operating points

Chosen to span the envelope, with at least one fully attached and one
near-stall condition:

| `v_inf` | λ | α range | regime |
|---|---|---|---|
| 5.0 m/s | 7.545 | 1.9 – 4.3° | fully attached |
| 7.0 m/s | 5.389 | 3.4 – 10.2° | attached, stall approaching inboard |
| 10.0 m/s | 3.772 | 5.9 – 19.7° | stall onset |
| 15.0 m/s | 2.515 | 10.6 – 31.9° | deep stall |

The 15 m/s point is beyond the S809 table's +18° upper limit, so when the
snapshot was taken, the old polar adapter was clamping α there. That was
deliberate: removing the clamp was one of the planned changes, and this point
is where it showed up (Cp moved −25.2 %, with 16 of 19 stations previously
clamped; see the change log below).

## Tolerance

`rtol = 1e-10`, with a per-field `atol = 1e-12 × max|golden value in that
field|`. The full reasoning is in the docstring of
`tests/test_golden_regression.py`. In short, the floor is set by `brentq`'s own
`xtol = rtol = 1e-12` on the root φ, and the ceiling by the fact that the
smallest physically meaningful change (0.01 % of `Cp`) is 1e-4, eight orders of
magnitude above the tolerance.

**Measured sensitivity.** Perturbing every station's chord by a relative 1e-9
moves `Cp` by 9.2e-10 and is caught; 1e-10 moves it by 9.2e-11 and isn't. The
response is linear, with gain ≈ 0.92 from 1e-12 to 1e-8, so the detection
threshold sits at the stated tolerance with no dead zone. Checked end to end: a
relative change of 1e-6 to one station's chord fails four of the five tests,
with messages naming the station and quantity.

## Important: the stored CCBlade and pyBEMT values are out of date

`cross_tool_summary.json` records mean |Cp| differences of **13.32 %**
(CCBlade) and **15.79 %** (pyBEMT) against the stored external results, not
the **0.51 %** and **1.83 %** in `docs/validation/bem-cross-validation.md`.

**This isn't a solver regression.** It's an out-of-date reference, found while
capturing this snapshot:

- `docs/validation/ccblade_case/results.json` and `pybemt_case/results.json`
  were committed in `3f19e93`.
- The S809 table was later rebuilt at `Ncrit = 5` (`3827f04`) and extended to
  Re = 1.3M (`2f2f017`).
- At that time `compare_pybemt.RE_BUCKETS` stopped at 500k, so **every station
  was clamped to the 500k table**, which is the problem noted in
  `build_polar_cache.py`. The buckets now run to 1.3M, and the stations use
  600k–900k.

So the external results in those files came from CCBlade and pyBEMT running on
polar tables the repository no longer contains, and only this solver's side has
been updated since. The published 0.51 % and 1.83 % are genuine for the
comparison as it was run, but reproducing them against the current polar table
would mean **re-running the external tools**, which need their own virtual
environments.

The two percentages above are pinned here as **tracking values**, so a change
can't move them unnoticed. They aren't validation results and shouldn't be
quoted as such.

The QBlade comparison isn't affected:
`data/qblade/bem_solver_phase_vi_result_Re*.csv` were produced after the table
rebuild and reproduce **bitwise** from the current code, which
`test_qblade_pinned_matches_committed_csv` checks exactly. Until the external
tools are re-run, the golden files here are the working check that the Phase VI
cross-check results haven't changed.

## Change log

Every regeneration, what moved, and why.

### Gap-free polar table, extended to ±180°

**What moved:** 12 of 41 `phase_vi_cp_lambda` values (λ = 3.5 … 6.0) and 81 of
1160 `phase_vi_spanwise` values, all by a relative 1e-5 to 3e-4.

**Why:** the eight interior holes in the S809 table were closed. Before, a
lookup landing on a hole read a linear bridge between the converged
neighbours; now it reads the documented fill. The spanwise change is entirely
at the **7 m/s** point, stations **12–14** (α = 7.17–7.66°, Re = 938k–945k),
the three whose lookup straddles α = 7.5° between the 900k and 1.0M curves,
which is where the (1.0M, 7.5°) hole was. Stations 11 (α = 8.18°) and 15
(α = 6.64°) don't straddle it and didn't move. The difference at the filled
point is about 0.001 in Cl, the same figure `validate_polars` check 1 had
already measured there.

**What didn't move:** `cross_tool_summary.json` in full (including both
tracking percentages and the QBlade pinned values), and
`test_snapshot_is_bitwise_reproducible`. The 15 m/s deep-stall point also
didn't move yet, because the old adapter still clamped α to the converged
range.

### Polar adapter replaced (clamping and substitution removed)

**What moved:** 26 of 41 `phase_vi_cp_lambda` values, 992 of 1160
`phase_vi_spanwise` values, and 420 values in `cross_tool_summary`, all within
`qblade_pinned`.

**What didn't move, and why that matters.** The `sweep` block, this solver's
side of the CCBlade and pyBEMT comparisons, is **bitwise unchanged** (all 28
values), and both tracking percentages are identical to every digit
(13.324537038059685 % and 15.787123761875879 %). So is
`reynolds_buckets_per_station`. Those paths run on `TableS809Polar`, a
discretised table over the ±180° Viterna export that this change doesn't
touch. That's the control: it shows the movement below comes from replacing
the adapter, not from something shifting under the whole snapshot.

**There are two separate causes**, and they show up in different places.

*1. The interpolant changed from bilinear to cubic in α and PCHIP in log(Re).*
This is small, everywhere, and the only cause wherever the whole span is
attached:

| point | Cp before | Cp after | |
|---|---|---|---|
| `cp_lambda`, λ = 4.0 … 7.5 | – | – | **+0.007 % … +0.099 %** |
| spanwise 5.0 m/s (α ≤ 4.3°) | 0.406172 | 0.406212 | +0.010 % |
| spanwise 7.0 m/s (α ≤ 10.2°) | 0.367031 | 0.367103 | +0.020 % |

The largest interpolant-only change in the snapshot is **+4.1 %**, in
`qblade_pinned[Re=100k].Cd` at station 9: Cd is convex in α, and linear
interpolation between 0.5° points systematically overestimates it.

*2. The α clamp was removed.* This is large, and only affects points with
stations beyond the +18° converged range:

| point | stations past 18° | Cp before | Cp after | |
|---|---|---|---|---|
| spanwise 10.0 m/s | 7 / 19 | 0.258726 | 0.256597 | −0.823 % |
| spanwise 15.0 m/s | 16 / 19 | 0.150104 | 0.112335 | **−25.16 %** |
| `cp_lambda` λ = 3.5 | – | 0.220746 | 0.212980 | −3.518 % |
| `cp_lambda` λ = 2.5 | – | 0.128481 | 0.093871 | −26.94 % |
| `cp_lambda` λ = 1.5 | – | 0.060918 | 0.023875 | **−60.81 %** |

Drag shows the mechanism most clearly. At 15 m/s the six innermost stations
sit at α = 25–32°. The old adapter returned `Cd(18°) ≈ 0.079` for all of them,
while the ±180° table gives 0.34–0.43, a **+340 % to +449 %** correction on the
largest drag values in the snapshot, and that's what brings Cp down. Clamped
lift doesn't fall off past stall and clamped drag doesn't rise, so every point
that used to be clamped was reporting too much power.

**None of this is a regression.** The old values were computed from lift and
drag the airfoil doesn't produce at those angles. The new ones use the Viterna
extrapolation built into the table, the same one QBlade has been using through
the exported `.plr` all along, so this solver and QBlade now share a post-stall
model as well as the measured range.

**Also updated in the same commit**, for the same reason and with the
before/after recorded where they live:
`validation/validate_powercurve.PHASE_VI_SEQUENCE_S_REFERENCE` and
`data/qblade/bem_solver_phase_vi_result_Re{100000,500000}.csv`.

**Coverage limits this exposed.** Removing the Reynolds clamp made two sweeps
stop where the S809 table genuinely ends, instead of running on a clamped
Reynolds number: the Cp–λ sweep in `validate_powercurve` check 4 (λ ≤ 7.5;
λ = 8.0 puts the tip station at Re = 1,384,108 against a 1,300,000 limit) and
its fixed-speed sweep in check 2 (≤ 20 m/s; at 25 m/s the *root* station, with
the largest chord, reaches 1,312,857). Neither is a defect; both are the
table's real coverage, previously hidden.

**What didn't move:** `test_snapshot_is_bitwise_reproducible` still passes, and
the fixed operating conditions (wind speeds, λ, ρ, ν, rpm) are unchanged, so
this is a like-for-like comparison.

### High-induction correction rewritten in Ning's γ-form: not regenerated

Recorded because the absence of a change is itself the result.

`corrections.corrected_axial_induction` was rewritten in Ning's γ₁γ₂γ₃ form.
That's the same equation reparameterised, not a different model, so the only
possible change is floating-point reordering. Measured over the whole snapshot:

| file | values moved | largest change |
|---|---|---|
| `phase_vi_cp_lambda.json` | 4 of 41 (all `Ct`) | 2.3e-15 relative |
| `phase_vi_spanwise.json` | 63 of 1160 | **2.4e-15 relative** |
| `cross_tool_summary.json` | 19 of 868 (all `qblade_pinned`) | 1.9e-15 relative |
| `cross_tool_summary.sweep` | **0** | – |

The largest change anywhere is 2.4e-15, five orders of magnitude below the
`rtol = 1e-10` these files are compared at. All three JSON regression tests
**pass unchanged**, so under this file's own rule (regenerate only when a test
fails, never as a routine step) the golden files were left as they were.

**What did have to be regenerated:**
`data/qblade/bem_solver_phase_vi_result_Re{100000,500000}.csv`. Those are
compared **bitwise** by `test_qblade_pinned_matches_committed_csv`, so a 4-ulp
change breaks them by design. One row changed in each file (station 18 at
Re = 100k: `phi_deg` 5.010656781762976 → 5.010656781762972). A tolerance test
staying green and a bitwise test going red are consistent outcomes for the same
change, and together they bracket its size.

### Ning's residual and a new root bracket: not regenerated

The residual was replaced with Ning's form, and the 2000-point bracket scan
with a bracket derived from the momentum region. Both change *how* the root is
found, not *which* root it is, so the movement is convergence noise around the
same solution:

| file | values moved | largest change |
|---|---|---|
| `phase_vi_cp_lambda.json` | 25 of 41 | 1.1e-13 relative |
| `phase_vi_spanwise.json` | 760 of 1160 | **8.0e-12 relative** |
| `cross_tool_summary.json` | 347 of 868 | 1.2e-12 relative |

The largest change is 8.0e-12, still well under `rtol = 1e-10`, so all three
JSON tests pass and the files were left alone. Only the bitwise QBlade CSVs
needed regenerating.

**One difference worth noting:** this time `cross_tool_summary.sweep` **did**
move, by 14 values, and both tracking percentages changed in their last three
digits (13.324537038059685 → 13.324537038062505). In the two previous entries
that block was frozen and served as the control. It isn't here, and that's
correct: `TableS809Polar` fixes the *polar* those runs see, not the *solver*.
The previous changes were to the polar layer, so the fixed-table paths were
insulated from them; this one changed the solver itself, which nothing
insulates. If this block *hadn't* moved, that would have been worth
investigating.

**What didn't move:** `reynolds_buckets_per_station`, the fixed operating
conditions, and `test_snapshot_is_bitwise_reproducible`.

### Interpolant continuity fix: regenerated

This is the first regeneration that is a **correctness fix rather than a
refactor**, so the golden files did have to change.

**What moved:** 26 of 41 `phase_vi_cp_lambda`, 992 of 1160
`phase_vi_spanwise`, 362 of 868 `cross_tool_summary` (including 14 in `sweep`),
and both QBlade CSVs.

**Largest changes:** `Cd` by up to **2.6e-4** relative, `Cp` by up to
**4.2e-5**, above the `rtol = 1e-10`, hence a real regeneration.

**Why.** The polar interpolant had a jump in value at every angle-of-attack
knot, for any Reynolds number *between* table rows, with Cl jumping by up to
5.4e-3 on SG6043. Blending cubic-spline *coefficients* across log(Re) with PCHIP
doesn't preserve continuity at the α knots, because continuity there is a
linear constraint on the coefficients and PCHIP is nonlinear in its data. The
old golden numbers were computed on a polar with those jumps in it; the new
ones aren't.

It surfaced when a BEM station on the Schmitz baseline failed to converge at
Re = 70,995 with a genuinely discontinuous residual, which is the solver's
reported-non-convergence doing exactly what it was built for. Under the
original solver this would have been a `ValueError` aborting the sweep, or
worse, a quietly wrong number.

**These values are more correct, not just different.** At a table node on a
table row, the surface now reproduces the stored value bitwise, so the
underlying data is unchanged; only the interpolation between data points moved,
onto a surface that is genuinely C¹.

## Regenerating

```
python tests/generate_golden.py
```

Regeneration is **not** a routine step. Any change to these files should be
deliberate, committed separately, and explained in the commit message, never
absorbed quietly into a code change. So: commit the code change, run the test,
and only if it fails for a reason you can state, regenerate in a **second
commit** whose message gives that reason.

## Running

```
python -m pip install pytest      # or: pip install -e ".[test]"
pytest tests/test_golden_regression.py
```

Run from the repository root. `pyproject.toml` puts `src/` and `tests/` on the
path, so no install step or `cd src` is needed.
