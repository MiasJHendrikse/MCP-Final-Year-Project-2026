# Why maximising energy gained so little: the measurements

**This is a record of an earlier investigation, not a current result.** When
the first energy optimisation found only +0.217 % over the Schmitz blade, I
investigated why. The conclusion is summarised in `docs/DESIGN-BASIS.md` §2.
This folder holds the measurements behind it, made with the bounds and
operating law in use at the time (including a 0.45 m chord bound that has
since changed).

`reoptimise.py` and `diagnostics.py` still import a set of bounds
(`PROVISIONAL_BOUNDS`) that has since been removed, so they won't run on the
current code without changes. I kept them and their `opt_*.json` outputs as
they were. The current versions of these results are in
`verification/fd_optimisation/`, `verification/adjoint_optimisation/` and
`verification/fd_optimisation_multistart/`.

**The one exception is `cross_evaluate_xstar.py`**, which is kept up to date
under the current 300 rpm law and bounds. It only evaluates blades and doesn't
optimise them, and its input x\* is now the current energy optimum.

## What's here

| file | contents |
|---|---|
| `diagnostics.py` → `diagnostics.json` | evaluation-only sweeps through the project's own objective: the gain bin by bin, re-evaluation with a fixed rating, Cp–λ, design and operating tip-speed-ratio sweeps, design Reynolds number, losses from rotor-speed ceilings, the literature's usual baselines (linear, Betz, Schmitz at a datasheet angle of attack), root clipping, spanwise loading, and every re-optimised blade under every operating strategy |
| `reoptimise.py` → `opt_<ceiling>_<rating>.json` | the finite-difference SLSQP run (central differences, `h = 3e-6`, `ftol = 1e-8`, the envelope constraint and the bounds of the time, from x₀) on a modified objective: generator rating fixed at x₀'s 3822.2 W or rising with the design, with an optional tip-speed ceiling |
| `opt_none_fixed.json` | fixed rating, no ceiling: **+0.121 %** |
| `opt_60_fixed.json` | 60 m/s tip-speed ceiling (286 rpm), fixed rating: **+0.34 %** |
| `opt_55_fixed.json` | 55 m/s (263 rpm), fixed rating: **+2.33 %** |
| `opt_50_fixed.json` | 50 m/s (239 rpm), fixed rating: **+7.84 %** |
| `opt_55_float.json` | 55 m/s, rating rising with the design: +7.40 % |
| `cross_evaluate_xstar.py` → `cross_evaluate_xstar.json` | kept up to date. No re-optimisation: x₀ and x\*, plus the ceiling-optimised blades above for context, evaluated bin by bin under each operating strategy. x\* is now the optimum under the 300 rpm law, so the 300 rpm row is its own law, not a hypothetical |

Gains are measured against the Schmitz blade x₀ evaluated *under the same
strategy*.

## What the up-to-date cross-evaluation shows

| ceiling | x\* over x₀ | `opt_60` | `opt_55` | `opt_50` | first bin below λ = 6.5 | `V_c` |
|---|---|---|---|---|---|---|
| none | +0.058 % | −0.069 % | −0.430 % | −0.985 % | – | – |
| **300 rpm (current)** | **+0.147 %** | +0.104 % | −0.154 % | −0.610 % | 10.5 m/s | 9.67 m/s |
| 60 m/s | +0.270 % | +0.337 % | +0.205 % | −0.170 % | 9.5 m/s | 9.23 m/s |
| 55 m/s | +0.898 % | +1.883 % | +2.327 % | +2.150 % | 8.5 m/s | 8.46 m/s |
| 50 m/s | +2.470 % | +4.727 % | +7.573 % | +7.844 % | 8.5 m/s | 7.69 m/s |

The `opt_*` columns are the older blades, re-evaluated under each strategy.
Two things stand out:

- **The ceiling costs energy, and the optimised blade wins most of it back.**
  Under the current law, the ceiling takes 0.0225 MWh/yr from x₀ (10.2702 →
  10.2477), and re-optimising recovers 0.0150 of it (+0.147 %). With the
  ceiling *removed*, the same blade's advantage over x₀ falls to +0.058 %,
  while the blade optimised for no ceiling (`opt_none_fixed`) gains +0.121 %.
  So x\* is tuned to the current law rather than better than x₀ everywhere.
  Under the ceiling the optimiser also wants *more* root chord, which is why
  `chord_0` sits on the 0.30 m upper bound.
- **Optimising for a specific ceiling beats transferring a blade.** x\* falls
  short of the blades optimised for a 55 m/s or 50 m/s ceiling by factors of
  2.6 and 3.2 (0.898 against 2.327 %, and 2.470 against 7.844 %). A tighter
  ceiling changes the best shape, not just the operating point, so if the
  machine's ceiling were ever lower than 300 rpm, x\* would need re-optimising
  rather than just re-evaluating. (Those comparison blades don't fit the
  current 0.30 m bound: `opt_55`'s root chord is 0.441 m. They're context, not
  candidates.)

The bin-by-bin listing shows where the gain comes from. From 10.5 m/s up,
the bins where the ceiling pushes λ below 6.5, x\* beats x₀ by +0.86 %, rising
to +12.27 % at 19.5 m/s, while the two lowest bins lose 1.17 % and 0.32 %.
That fits the extra root chord x\* carries: it gives up a little in light winds
and gains a lot in strong ones.

## Reproducing it

    python verification/aep_gain_audit/reoptimise.py none fixed verification/aep_gain_audit/opt_none_fixed.json
    python verification/aep_gain_audit/reoptimise.py 55 fixed verification/aep_gain_audit/opt_55_fixed.json
    python verification/aep_gain_audit/diagnostics.py
    python verification/aep_gain_audit/cross_evaluate_xstar.py

Each `reoptimise.py` case takes 3–7 minutes, `diagnostics.py` about 4 minutes
and `cross_evaluate_xstar.py` about 1 minute. Only the last one is expected to
reproduce the committed output; it reads x\* from
`verification/adjoint_optimisation/result.json`.
