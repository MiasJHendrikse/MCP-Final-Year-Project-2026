# Polar sensitivity of the minimum-material result

**New 2026-09-26**, review roadmap item 1 (the Devil's Advocate CRITICAL
point). The mass optimum `x_m` matches AEP(x0) to ~1e-12 on the production
SG6043 cache (XFOIL, n_crit = 9), whose lift differs from the UIUC
measurements by 0.11-0.15 RMS at the design Reynolds numbers. Does the
saving survive that uncertainty?

## What was run

    python verification/polar_sensitivity/run_polar_sensitivity.py    # ~5 min

Six perturbed polar surfaces, each built from the production grid with the
change applied over the XFOIL band (-8 .. 18 deg) and a fresh C1
interpolant, and routed through every polar lookup (forward solver, adjoint
system, problem), so forward values and gradients come from the same surface:

- **uniform lift offset `+/- s(Re)`**, `s` the measured lift RMS difference
  against UIUC per Reynolds number (0.144 at 100 k .. 0.024 at 500 k),
  interpolated in log Re and held at its end values outside that range;
- **tilted lift offset `+/- s(Re) w(Re)`**, `w` falling linearly in log Re
  from +1 at 62 000 to -1 at 339 000 (the blades' Reynolds range): the shape
  of error that does not cancel between blades, sized by the measured error
  (an adversarial test, not a model of the true error);
- **drag scaled by 0.8 and 1.2**, an assumed +/- 20 %.

Under each surface: (a) the three fixed blades are re-evaluated, and (b) the
`delta = 0` mass problem is **re-optimised** from `x0`, every reference
re-measured on that surface. The production case is re-optimised too, as a
check of the path: it returns -3.379 %, the committed result.

## Result (`polar_sensitivity.json`)

| surface | AEP(x_m) vs AEP(x0), fixed blades | re-optimised shell saving |
|---|---|---|
| production (n_crit = 9) | -0.0000 % | **-3.379 %** |
| lift + s(Re) | +0.188 % | -6.353 % |
| lift - s(Re) | **-0.262 %** | **-2.938 %** |
| lift + s w: low Re up, high Re down | -0.032 % | -3.092 % |
| lift - s w: low Re down, high Re up | +0.042 % | -4.114 % |
| drag x 0.8 | -0.005 % | -3.330 % |
| drag x 1.2 | +0.005 % | -3.443 % |

- **The fixed blade's energy parity does not survive the polar error.** A
  lift error of the measured size moves `x_m` against `x0` by -0.26 % to
  +0.19 %: it is the *uniform* offset, not the tilted one, that does not
  cancel, because it moves the lift at which Schmitz's planform is optimal
  and the two planforms respond differently. At the reference-energy
  exchange rate, -0.26 % of energy is worth about 2.9 points of material.
- **A saving survives on every surface.** Re-optimised, the minimum-material
  blade saves between 2.94 % and 6.35 % of shell material at equal energy.
  The sign of the result is robust to the polar uncertainty; its size, and
  the exact blade, are not (-3.38 % is one value in a band of about 3-6 %).
- Drag of +/- 20 % barely matters (+/- 0.005 %; saving 3.33-3.44 %).

## Not done: XFOIL at n_crit = 7 and 11

The first plan was to rebuild the cache at the n_crit study's neighbouring
settings, 7 and 11. Both caches were built on 2026-09-26, but XFOIL's
gap-closure retries (`xfoil.close_polar_gaps`) failed to converge at
Re = 100 000; by MJ's instruction the run was stopped and the change
reverted (registry entries, configs and data removed; nothing committed).
The lift perturbations above bound the lift error directly instead.
