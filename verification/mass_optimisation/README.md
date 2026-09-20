# Phase 5: the minimum-material blade (AEP-constrained mass optimisation)

**New 2026-09-20** (plan `docs/PLAN-mass-objective-2026-09-20.md`, Step 4).
Run under the 300 rpm operating law
(`lambda(V) = min(6.5, Omega_max R / V)`, `V_c = 9.67 m/s`), the fixed rating
`3822.189755449124 W` (B2 provisional), the project bounds
(`chord_max_m = 0.30`, `chord_min_m = 0.045`, twist `-2 .. 35 deg`, local
solidity cap 0.5) and the polar-cache envelope at margin 0.05. Every reference
is the committed Schmitz `x0` (`verification/baseline/x0.json`): AEP(x0),
`KS0`/`M_ref`, the root chord `c00`, `delta_ref`/`D0`, and the shell area
`0.43705266291483963 m^2` the objective is normalised by.

## What was run

    python verification/mass_optimisation/run_mass_slsqp.py --delta 0     # ~20 s, the production optimum
    python verification/mass_optimisation/run_multistart.py               # ~4 min, 10 starts at delta = 0
    python verification/mass_optimisation/run_sweep.py                     # ~4 min, the delta sweep and the ablation
    python verification/mass_optimisation/run_mass_checks.py               # ~1 min, Tiers 1-3, Taylor, the reference table
    python verification/mass_optimisation/run_cross_evaluation.py          # ~1 min, x0 / x_c / x_m under other resources and ceilings

All of them read `verification/baseline/x0.json`,
`verification/load_constraint/result_eps0.json` (`x_c`, the Phase 4 energy
optimum), `verification/fd_optimisation_multistart/results.json` (the
`0.0166` agreement spread in `u`) and `verification/fd_step_size/sweep.json`
(`h*`); `run_mass_checks.py` and `run_cross_evaluation.py` read
`result_delta0.json`, so the order above is the order. `mass_optimisation`
runs after `load_constraint` in `verification/README.md`. The shared code is
`_common.py` (problem construction, the recorder, the KKT estimate, the blade
record, the figures) -- lifted from `load_constraint/run_constrained_slsqp.py`,
not copied a fourth time.

## The problem

    minimise    m(u) = k_P int c dr / (k_P int c dr)(x0)        shell material, geometric, exact gradient
    subject to  AEP(u) / AEP(x0) - (1 - delta)   >= 0            energy floor, objective adjoint
                KS0 - KS(u)                      >= 0            Phase 4 moment cap (eps = 0), moment adjoint
                KS0 / c00^2 - KS(u) / c0(u)^2    >= 0            root-stress proxy, moment adjoint
                D0 - D(u)                        >= 0            tip-deflection proxy, deflection adjoint
                c_i - c_{i+1} >= 0 (4), theta_i - theta_{i+1} >= 0 (4)   monotone control points, constant Jacobian
                envelope (50 rows), solidity (25 rows), u in [0, 1]^10

`ScaledProblem.constraints_for_mass_problem(delta)`; SLSQP, `ftol = 1e-10`,
`maxiter = 400`. One nine-point load solve per accepted iterate serves the
moment, stress and deflection rows (`deflection_solves_total = 0` in every
run: the deflection row always found the moment row's state). The state rows
are guarded: a trial point a row cannot evaluate returns `-1e3` and is logged
in `evaluation_failures` -- **no run logged one**.

The rows are what `docs/adjoint_derivation.md` sections 10 and 11 derive;
the material model is `src/objective/mass.py` (`k_P = 2.048504996756389`,
`k_A = 0.06850241963050001`, exact from the SG6043 coordinates); the solid
proxy `k_A int c^2 dr` is reported beside the shell figure, never optimised.

## The result (`result_delta0.json`, `blade_delta0.png`)

**`x_m`: -4.238 % shell material (-8.532 % on the solid proxy) at exactly
AEP(x0)** -- the floor's slack is `-1.2e-12`. SLSQP exit 0 in 22 iterations,
17 s; 43 objective solves, 23 objective adjoints, 30 load solves.

| | `x0` (Schmitz) | `x_c` (energy optimum, Phase 4) | **`x_m` (mass optimum)** |
|---|---|---|---|
| AEP vs `x0` | 0 | +0.1465 % | **-0.0000 %** (slack -1.2e-12) |
| shell material `k_P int c dr` vs `x0` | 0 | +6.171 % | **-4.238 %** |
| solid proxy `k_A int c^2 dr` vs `x0` | 0 | +14.196 % | **-8.532 %** |
| `KS / KS0` (moment cap) | 1.0000 | 1.0000 (active) | 0.9562 (inactive, slack 0.044) |
| stress proxy ratio | 1.0000 | 0.8451 | 1.0000 (active) |
| deflection proxy ratio | 1.0000 | 0.7618 | 1.0000 (active) |
| rated-point root moment [N m] | 177.38 | 177.41 | 169.55 (-4.41 %) |
| rated-point rotor thrust [N] | 517.5 | 525.5 | 500.9 (-3.2 %) |
| cut-in `Re_min` | 73 170 | 62 149 | 58 380 (floor 40 000) |
| chord control points [mm] | 276 / 182 / 97 / 79 / 67 | 300 / 182 / 126 / 78 / 48 | **270 / 152 / 124 / 72 / 45** |
| twist control points [deg] | 24.5 / 9.7 / 4.0 / 1.0 / 0.6 | 24.3 / 10.6 / 5.2 / 2.2 / -1.9 | **23.4 / 8.5 / 5.0 / 1.2 / -1.8** |

(`reference_blades.json`, `reference_blades.png`; the `x_c` row is the
"+6.2 % material for +0.15 % energy" motivation of the re-pitch, measured on
the same code.)

**Where the material comes off.** The tip: `chord_4` sits on its 45 mm lower
bound (multiplier 0.073) and the outer third of the blade is thinner than
Schmitz; the root stays close to Schmitz (270 vs 276 mm) because the stress
row pins `KS / c0^2` and the deflection row pins `int M (R - r) / c^3 dr`,
which the `1 / c^3` panel of `blade_delta0.png` shows is dominated by the
outer stations. The twist moves with it (about -1 deg inboard, -2.4 deg at
the tip control point) to hold the energy. The monotone rows are slack
(minimum 0.027) at `delta = 0`.

**Active set and KKT** (`kkt`): `aep_floor`, `stress`, `deflection`,
`lower:chord_4`; all four multipliers positive; residual after projection
`2.3e-6` on a gradient of norm 0.93. The AEP row's multiplier is the
**exchange rate at the margin: 14.5 % material per 1 % energy** (both the
objective and the row are fractions of `x0`'s, so the multiplier needs no
unit conversion). The moment cap is inactive: at the mass optimum the load
is 4.4 % below Schmitz's without being asked, because a thinner blade at the
same energy carries less thrust.

**Post-checks.** `alpha` on the uncapped points within `[1.49, 8.83] deg`
(band `[-8, 18]`); `reynolds_min = 58 380` against the 40 000 floor, with the
tightest envelope row (`floor r = 1.966`) 13.5 mm slack; every station
converged; no `PolarDomainError`, no margin escalation.

**Sanity band** (stated in the plan before the run): shell saving 2-10 %, solid
5-20 %, AEP floor active with a non-negative multiplier, moment cap inactive,
no failed trials, `alpha` within, `Re_min` above the floor. All met;
`-4.24 %` is below the scratch's deflection-free `-7.5 %`, as it must be (the
ablation below reproduces `-7.59 %` without the two proxies).

## Multi-start (`multistart_delta0.json`, `multistart_delta0.png`)

Ten starts -- `x0`, `x_c`, and the eight committed random starts of
`fd_optimisation_multistart/starts.json` (drawn for the energy problem; most
are infeasible for this one and SLSQP finds the feasible set first). **All
ten exit 0 at the same point**: spread in `u` `1.7e-5` (criterion `0.0166`),
spread in the objective `3.9e-11`, worst constraint violation `-7.7e-12`, no
failed trial point in any run. 216 s.

## The energy-floor sweep (`pareto.json`, `pareto.png`)

Cold from `x0` and warm from the previous floor at each `delta`; cold and
warm agree to `<= 2.4e-6` in `u` everywhere (criterion `0.0166`).

| `delta` | energy given up | shell saved | solid saved | exchange rate (KKT) | secant to the next point | active rows beyond AEP / stress / deflection |
|---|---|---|---|---|---|---|
| 0 | 0 | **4.238 %** | 8.532 % | 14.54 | 8.72 | `lower:chord_4` |
| 0.25 % | 0.2500 % | 6.417 % | 12.080 % | 5.94 | 4.94 | `lower:chord_4` |
| 0.5 % | 0.5000 % | 7.651 % | 14.014 % | 4.22 | 2.88 | `lower:chord_4` |
| 1 % | 1.0000 % | 9.090 % | 16.209 % | 2.14 | 1.49 | `lower:chord_4`, monotone `chord_1-chord_2`, `twist_3-twist_4` |
| 2 % | 2.0000 % | 10.576 % | 17.925 % | 1.09 | -- | `lower:chord_3`, `lower:chord_4`, monotone `chord_3-chord_4`, `twist_3-twist_4` |

The front is concave, the saving is monotone in `delta`, and **every secant
lies between the KKT exchange rates at its two ends** (8.72 between 14.54
and 5.94, and so on) -- the adjoint multiplier and the front's own slope
check each other, which is what the tangents in `pareto.png` show. The
stress and deflection rows stay active along the whole front; from `delta =
1 %` the monotone rows start to bind (the optimiser would otherwise put a
chord bump at the second control point), and at `2 %` two tip control
points sit on the 45 mm bound. The first quarter-percent of energy buys the
most material (2.2 points for 0.25 %); after 1 % the rate is near 1:1.

Read against decision 2 of the plan (the reference stays at `lambda = 6.5`):
a `lambda = 6.0` Schmitz would raise the floor by about 0.16 %, which at the
first secant (8.7 points per %) costs roughly 1.4 of the 4.2 points of shell
material. Stated, not re-run.

## Ablation at `delta = 0` (`ablation.json`, `pareto.png` middle panel)

| rows | shell saved | solid saved | stress ratio | deflection ratio | `KS / KS0` |
|---|---|---|---|---|---|
| full set | **4.238 %** | 8.532 % | 1.000 | 1.000 | 0.956 |
| without the deflection row | 4.824 % | 10.165 % | 1.000 | 1.064 | 0.968 |
| without the stress row | 5.654 % | 13.609 % | 1.806 | 1.000 | 0.955 |
| without both proxies | 7.586 % | 18.650 % | 2.137 | 1.129 | 0.972 |
| without the monotone rows | 4.238 % | 8.532 % | 1.000 | 1.000 | 0.956 |
| without the moment cap | 4.238 % | 8.532 % | 1.000 | 1.000 | 0.956 |

Removing a row never shrinks the saving (to the solver tolerance; the
`no monotone` and `no moment cap` optima are the full-set optimum to
`2e-9` points -- those rows are slack there). **The stress row is doing the
work**: without it the optimiser halves the root chord's section modulus
(stress ratio 1.81) for another 1.4 points of shell; without the deflection
row it accepts 6.4 % more tip deflection for 0.6 points. Without both, the
scratch's `-7.6 %` returns with a root that carries 2.1x Schmitz's stress
proxy -- the number the re-pitch proposal warned was not a blade. The
saving that survives both proxies, `-4.2 %`, is the defensible one.

## Checks (`checks.json`)

At `x0` and at `x_m`:

| check | `x0` | `x_m` |
|---|---|---|
| shell / solid gradient vs central FD (relative) | 1.7e-11 / 4.0e-12 | 1.7e-11 / 2.0e-12 |
| manufacturing rows exact to the bit; min slack | yes; 0.0072 | yes; 0.0267 |
| forward twins: deflection / moment (relative) | 2.2e-16 / 2.2e-16 | 2.2e-16 / 2.2e-16 |
| Tier 1 `dD/dphi`, `dD/dd`, `dKS/dx`, `dKS/dd` (mixed) | 1.1e-15, 6.7e-16, 8.3e-16, 5.8e-16 | 2.1e-15, 1.1e-15, 7.8e-16, 6.1e-16 |
| Tier 2 tangent vs adjoint; assembly identity | 3.9e-16; 0 | 3.8e-16; 0 |
| Tier 3 stress row, worst `abs(adj - fd) / eps_j` | 0.366 (`chord_2`) | 0.372 (`chord_1`) |
| Tier 3 deflection row, worst `abs(adj - fd) / eps_j` | 0.812 (`twist_1`) | 0.908 (`twist_0`) |
| Taylor remainder ratios per decade (six directions) | 97-100 | 96-100 |

Acceptance `3 eps_j`; every variable of both rows passes at both points, so
the round-off floor was not measured (the Phase 4 rule: measured only on
failure). The `eps_j` fragility of `gradient_verification/` (the red
`chord_4` at `x0`) did not reappear on these rows at `x_m`.

**Re-evaluation of `x_m` from its JSON** through the forward path only
(`objective()`, `root_moment`, `tip_deflection_at`, `shell_area`): AEP, the
rated root moment, the rated tip deflection and the shell area agree with
the recorded numbers to `0.0` relative (plan, Verification 5).

## Cross-evaluation (`cross_evaluation.json`)

No optimisation: `x0`, `x_c`, `x_m` evaluated under other conditions. The
energy parity of `x_m` holds under one resource and one operating law; this
is how much of it survives elsewhere.

| condition | `x_c` AEP vs `x0` | **`x_m` AEP vs `x0`** | `x_m` rated moment vs `x0` |
|---|---|---|---|
| central resource (`k = 1.709, c = 7.274`), 300 rpm | +0.147 % | **-0.000 %** | -4.41 % |
| `k = 1.5, c = 7.049` | +0.135 % | +0.000 % | (unchanged) |
| `k = 1.5, c = 7.633` | +0.132 % | -0.007 % | |
| `k = 2.0, c = 7.049` | +0.163 % | +0.015 % | |
| `k = 2.0, c = 7.633` | +0.162 % | -0.003 % | |
| no ceiling | +0.063 % | +0.098 % | -3.41 % |
| ceiling 60 m/s | +0.262 % | -0.183 % | -4.85 % |
| ceiling 55 m/s | +0.835 % | -0.590 % | -5.41 % |
| ceiling 50 m/s | +2.342 % | -2.583 % | -5.80 % |

**Resource: robust** -- over the four corners `x_m` stays within
`[-0.007, +0.015] %` of `x0`: the per-bin energy differences between the
two blades are small and of both signs, so re-weighting the bins moves the
total very little. **Ceiling: not robust** -- a lower rotor-speed ceiling pushes the
rated point to lower `lambda`, where the thin tip of `x_m` loses more than
Schmitz's, and by 50 m/s it gives up 2.6 % of energy (while `x_c`, with
its fat tip, gains 2.3 %). The mass optimum is an optimum *for the 300 rpm
machine*; a different ceiling wants a different blade, and the honest
statement of the headline is "-4.2 % material at equal energy under the
configured law". The resource band used here (`c` from the log-law
cross-check, `k +/- 0.25`) is an assumption of this artefact -- plan 1.3's
sensitivity band is still not chosen and this does not choose it.

## What would make it wrong

- **The material proxy.** `k_P int c dr` is a constant-thickness shell of one
  section with no spar, no root insert, no web, no ply drop -- mass is
  proportional to planform area. A real laminate schedule would weight the
  root more (and a solid or spar-dominated blade would follow `c^2` or
  `c^3 t`); the solid proxy is reported beside it as the other end of the
  range (-8.5 % vs -4.2 %). No kg figure is quoted: `structure.*` in
  `config/rotor_design.yaml` is `TODO` and only multiplies this number.
- **The loads are operating-range only (B3).** `L` is the nine points with
  `P <= P_rated` plus the rated point; no parked, gust, yawed or fault case.
  A survival load case would set the root, not the rated point, and the
  stress row would be a different row.
- **The three structural rows are relative proxies**, all on the same thin-
  shell assumption (`Z ~ c^2 t`, `I ~ c^3 t`, mass `~ c`), pinned to
  Schmitz's values. They say the blade is no worse than Schmitz on each
  measure; they do not say Schmitz was adequate. No fatigue, no buckling
  (plan, cut list).
- **The Reynolds floor.** `Re_min = 58 380` against the 40 000 cache floor,
  13.5 mm of envelope slack at the tip station; the tip control point sits
  on the 45 mm chord bound. Both are numbers the polar cache and the bounds
  set, not physics, and a thinner cache or a lower `chord_min_m` would let
  the optimiser go further.
- **The monotone rows** are on the control points (variation diminishing),
  which is sufficient for a monotone spline but not necessary; a mould
  constraint might allow shapes these forbid, or forbid shapes these allow
  (the `delta = 1 %` optimum has two equal chord control points).
- **The parameterisation**: 5 + 5 cubic B-spline control points, selected in
  `representation_study/` for the energy objective. A finer chord basis
  might find more material to remove at the same energy; a coarser one
  less.
- **One resource, one law** -- see the cross-evaluation.
