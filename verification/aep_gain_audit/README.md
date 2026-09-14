# AEP-gain audit — measurements behind `docs/AEP_GAIN_AUDIT.md`

**Under provisional bounds: `chord_max_m = 0.45 m` is provisional.** Nothing
here is a project result; it is the evidence for the audit of why the
optimiser gains +0.217 % and where an honest larger number would come from.
Read the audit, not this file, for the argument.

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

Gains are over the Schmitz `x0` evaluated *under the same strategy*.

## Reproduce

    python verification/aep_gain_audit/reoptimise.py none fixed verification/aep_gain_audit/opt_none_fixed.json
    python verification/aep_gain_audit/reoptimise.py 55 fixed verification/aep_gain_audit/opt_55_fixed.json
    python verification/aep_gain_audit/diagnostics.py

Wall time: 3–7 min per `reoptimise.py` case, ~4 min for `diagnostics.py`.
