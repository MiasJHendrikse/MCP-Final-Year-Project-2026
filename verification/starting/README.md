# Starting torque of the committed blades (relative)

**New 2026-09-26**, review roadmap (minor item on starting). The design takes
the 3 m/s cut-in as a machine requirement; this artefact asks whether the
minimum-material blade `x_m` starts as readily as the Schmitz reference `x0`.

## What was run

    python verification/starting/run_starting_torque.py     # ~2 s

At rest every section meets the wind at `alpha = 90 deg - theta(r)` with no
induction, so the parked rotor torque is

    Q_0 = B int (1/2) rho V^2 c(r) C_l(90 deg - theta(r)) r dr

(the drag term vanishes at `phi = 90 deg`), evaluated at the cut-in speed on
the 25 BEM stations.

## Result (`starting_torque.json`)

| | parked torque at 3 m/s | vs `x0` | per unit shell mass, vs `x0` | share from the inboard half |
|---|---|---|---|---|
| `x0` | 0.377 N m | -- | -- | 0.80 |
| `x_c` | 0.445 N m | +18.0 % | +11.2 % | 0.78 |
| **`x_m`** | **0.351 N m** | **-6.9 %** | **-3.7 %** | 0.80 |

`x_m` starts with about 7 % less aerodynamic torque than the reference (its
inboard chord is 3-6 % smaller, and 80 % of the parked torque comes from the
inboard half); per unit of shell mass, a rough stand-in for the rotor's
inertia, the deficit is about 4 %. The minimum-material blade therefore
starts somewhat less readily than Schmitz, which is the trade the small-turbine
literature warns about.

## What it cannot say

- A ratio between blades, not a prediction: every angle here (67-91 deg) is
  in the cache's Viterna extension, and the station Reynolds numbers at 3 m/s
  (about 10 000 to 45 000) lie below the cache's 40 000 floor, so the
  coefficients are read on the 40 000 row. Viterna barely depends on Reynolds
  number at these angles, and both simplifications apply to all three blades
  alike.
- No generator cogging or friction torque and no rotor-inertia model is
  specified, so whether any of the three blades actually starts at 3 m/s is
  not answered; shell mass is only a proxy for inertia (inertia weights mass
  by `r^2`).
