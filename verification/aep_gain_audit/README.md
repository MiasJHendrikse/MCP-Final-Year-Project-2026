# AEP-gain audit — measurements behind `docs/AEP_GAIN_AUDIT.md`

**History: under the retired provisional bounds (`chord_max_m = 0.45 m`).**
Nothing here is a project result. It is the evidence for the audit of why the
optimiser gained +0.217 % under the bounds and operating law in force on
2026-09-13, and where an honest larger number would come from. Read the audit,
not this file, for the argument.

**One exception: `cross_evaluate_xstar.json` was re-run on 2026-09-19** under
the 300 rpm operating law and the configured bounds, because it evaluates
blades rather than optimising them and its input `x*` now *is* the
ceiling-constrained optimum. Everything else here still reflects the retired
box.

`reoptimise.py` and `diagnostics.py` deliberately still import
`PROVISIONAL_BOUNDS` from `tests/test_parameterisation.py`, which was retired
on 2026-09-19 (the box now lives in `config/rotor_design.yaml`); the two
scripts and every `opt_*.json` are kept as written, as the record behind the
audit, and are **not** expected to reproduce under today's configuration
without being ported. The re-run artefacts that supersede them live in
`verification/fd_optimisation/`, `verification/adjoint_optimisation/` and
`verification/fd_optimisation_multistart/`.

## What is here

| file | what |
|---|---|
| `diagnostics.py` → `diagnostics.json` | Evaluation-only sweeps through the repo's own objective chain: bin-by-bin gain attribution, fixed-rating re-evaluation, Cp–λ, design-TSR and operating-TSR sweeps, design-Re, rotor-speed-ceiling losses, the literature's baselines (linear, Betz, datasheet-α Schmitz), root clipping, spanwise loading, and a cross-evaluation of the re-optimised blades under every strategy |
| `reoptimise.py` → `opt_<ceiling>_<rating>.json` | A4's SLSQP (central FD, `h = 3e-6`, `ftol = 1e-8`, the envelope constraint, provisional bounds, from `x0`) on a modified objective: generator rating fixed at `x0`'s 3822.2 W or floating; optional tip-speed ceiling (Region 2½) |
| `opt_none_fixed.json` | fixed rating, no ceiling: **+0.121 %** |
| `opt_60_fixed.json` | 60 m/s tip-speed ceiling (286 rpm), fixed rating: **+0.34 %** |
| `opt_55_fixed.json` | 55 m/s (263 rpm), fixed rating: **+2.33 %** |
| `opt_50_fixed.json` | 50 m/s (239 rpm), fixed rating: **+7.84 %** |
| `opt_55_float.json` | 55 m/s, rating floating with the design: +7.40 % |
| `cross_evaluate_xstar.py` → `cross_evaluate_xstar.json` | **re-run 2026-09-19.** No re-optimisation: two fixed blades (`x0`, `x*`) plus the audit's ceiling-optimised context blades, evaluated bin by bin under each operating strategy. `x*` is now the 300 rpm-law optimum, so the 300 rpm case is its own law, not a counterfactual |

Gains are over the Schmitz `x0` evaluated *under the same strategy*.

## What the re-run says

| ceiling | `x*` over `x0` | `opt_60` | `opt_55` | `opt_50` | first bin below λ = 6.5 | `V_c` |
|---|---|---|---|---|---|---|
| none | +0.058 % | −0.069 % | −0.430 % | −0.985 % | — | — |
| **300 rpm (config)** | **+0.147 %** | +0.104 % | −0.154 % | −0.610 % | 10.5 m/s | 9.67 m/s |
| 60 m/s | +0.270 % | +0.337 % | +0.205 % | −0.170 % | 9.5 m/s | 9.23 m/s |
| 55 m/s | +0.898 % | +1.883 % | +2.327 % | +2.150 % | 8.5 m/s | 8.46 m/s |
| 50 m/s | +2.470 % | +4.727 % | +7.573 % | +7.844 % | 8.5 m/s | 7.69 m/s |

The `opt_*` columns are the audit's provisional-bound blades, re-evaluated
under each strategy. Two readings:

- **The ceiling is what costs AEP, and the blade gives back most of it.** Under
  the config law the ceiling removes 0.0225 MWh/yr from `x0` (10.2702 →
  10.2477) and the re-optimisation recovers 0.0150 of it (+0.147 %). The
  interesting part is the `none` row: with the ceiling *removed*, the same
  blade's advantage over `x0` falls to +0.058 %, and the blade optimised for
  no ceiling at all (`opt_none_fixed`) gains +0.121 %. So `x*` is tuned to the
  constrained law in force, not uniformly better than `x0`; the ceiling is
  also what the cap pushes against — under it the optimiser wants *more* root
  chord, not less, which is why `chord_0` sits on the configured 0.30 m upper
  bound.
- **Ceiling-specific optimisation is worth more than transfer.** `x*`
  under-performs the blade optimised for a 55 m/s or 50 m/s ceiling by a
  factor of 2.6 and 3.2 respectively (0.898 vs 2.327 %, 2.470 vs 7.844 %).
  A tighter ceiling changes the optimal shape, not just the operating point —
  so if the machine eventually tightens below 300 rpm, `x*` must be re-run
  under the new law rather than re-evaluated. (Those context blades are
  infeasible under the configured 0.30 m box: `opt_55`'s root chord is
  0.441 m. They are context, not candidates.)

The bin listing shows where the gain lives: from 10.5 m/s up — the bins the
ceiling bends below λ = 6.5 — `x*` beats `x0` by +0.86 % rising to +12.27 % at
19.5 m/s, while the first two bins lose 1.17 % and 0.32 %. That is consistent
with the extra root chord `x*` carries: it pays a little in the low-wind bins
and collects a lot in the high-wind ones.

## Reproduce

    python verification/aep_gain_audit/reoptimise.py none fixed verification/aep_gain_audit/opt_none_fixed.json
    python verification/aep_gain_audit/reoptimise.py 55 fixed verification/aep_gain_audit/opt_55_fixed.json
    python verification/aep_gain_audit/diagnostics.py
    python verification/aep_gain_audit/cross_evaluate_xstar.py

Wall time: 3–7 min per `reoptimise.py` case, ~4 min for `diagnostics.py`,
~1 min for `cross_evaluate_xstar.py`. The last is the only one expected to
reproduce as committed; it reads `x*` from
`verification/adjoint_optimisation/result.json`.
