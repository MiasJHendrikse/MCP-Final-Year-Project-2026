# Phase 5: the minimum-material blade (AEP-constrained mass optimisation)

**New 2026-09-20** (plan `docs/PLAN-mass-objective-2026-09-20.md`, Step 4);
**re-run the same afternoon with the 60 mm buildable-tip floor**
(`manufacturing.min_chord_m`, five linear rows `c_i >= 0.060` of the mass
problem; MJ's decision, plan amendment "Decision 6"; the morning run
without it is in the journal and in this directory's git history at
`5227ba6`, not here). Run under the 300 rpm operating law
(`lambda(V) = min(6.5, Omega_max R / V)`, `V_c = 9.67 m/s`), the fixed
rating `3822.189755449124 W` (B2 provisional), the project bounds
(`chord_max_m = 0.30`, `chord_min_m = 0.045` -- the box is unchanged, see
below -- twist `-2 .. 35 deg`, local solidity cap 0.5) and the polar-cache
envelope at margin 0.05. Every reference is the committed Schmitz `x0`
(`verification/baseline/x0.json`): AEP(x0), `KS0`/`M_ref`, the root chord
`c00`, `delta_ref`/`D0`, and the shell area `0.43705266291483963 m^2` the
objective is normalised by. `x0`'s tip is 67 mm, above the floor, so
nothing about the reference moved.

## Why the floor was added (2026-09-20, afternoon)

The morning optimum sat on the 45 mm chord bound at the tip (`-4.24 %`
shell, 6:1 taper). The 2026-09-13 basis of that bound ("SG6043 10 % t/c ->
~4.5 mm laminate minimum") was the thickness of the *whole* section: a 45 mm
chord is 4.5 mm thick in total, with no room for two skins, a bond line and
a web, and the optimiser had found the hole in the constraint set.
`misc/blade_shape_experiment/` (scratch, gitignored) tried a 60 / 80 /
100 mm floor, a taper-ratio row `c_root <= k c_tip` and convexity rows on
the chord control points:

- **60 mm** (a 6 mm section: two 2 mm skins + bond + web) is feasible at
  `delta = 0`, saves 3.4 %, and lies below `x0`'s tip, so the comparison is
  unchanged. **Taken.**
- **80 mm** is *infeasible* at `delta = 0` (the best energy a blade with an
  80 mm tip can reach under the load cap is 0.002 % short of Schmitz's) and
  saves 0.8 % at `delta = 0.5 %`; because it sits above `x0`'s 67 mm tip, a
  fair 80 mm study would need Schmitz re-fitted under the same floor --
  decision 2's reference chain reopened. Not taken.
- The taper-ratio and convexity rows were **slack** at every optimum with
  the floor in place: the floor alone shapes the tip. Not added.

**Why a row and not the box.** The box bound `chord_min_m` sets the scaling
of the design variables, and every committed Phase 1-4 gradient and
optimisation artefact -- `h*_j`, `eps_j`, the red test's FD reference,
`x*`, `x_c` -- was measured in that scaling. Raising the box to 0.060 was
tried first: the red test's committed FD reference then disagreed with the
adjoint on every chord variable (the scaling had moved, not the gradient),
and the honest fix would have been to re-run Phases 2-4, which the plan
forbids. So the floor is five rows `c_i - 0.060 >= 0` beside the monotone
rows (`ScaledProblem.min_chord_rows`, `manufacturing.min_chord_m` in
config with the basis above), the box stays 0.045, and nothing outside this
directory moved. `x_c`, the Phase 4 energy optimum with its 47.7 mm tip,
violates the row by 12 mm and is kept as the energy reference, evaluated as
it is (and used as an infeasible start, which SLSQP handles).

## What was run

    python verification/mass_optimisation/run_mass_slsqp.py --delta 0     # ~20 s, the production optimum
    python verification/mass_optimisation/run_multistart.py               # ~3 min, 10 starts at delta = 0
    python verification/mass_optimisation/run_sweep.py                     # ~4 min, the delta sweep and the ablation
    python verification/mass_optimisation/run_mass_checks.py               # ~1 min, Tiers 1-3, Taylor, the reference table, the renders
    python verification/mass_optimisation/run_cross_evaluation.py          # ~1 min, x0 / x_c / x_m under other resources and ceilings

Every figure of this directory is redrawn from the committed JSON, without
solving anything, by

    python verification/mass_optimisation/replot_figures.py

(and each runner has its own `--replot`: `run_mass_checks.py`,
`run_multistart.py`, `run_sweep.py`).

All of them read `verification/baseline/x0.json`,
`verification/load_constraint/result_eps0.json` (`x_c`),
`verification/fd_optimisation_multistart/results.json` (the `0.0166`
agreement spread in `u`) and `verification/fd_step_size/sweep.json` (`h*`);
`run_mass_checks.py` and `run_cross_evaluation.py` read `result_delta0.json`,
so the order above is the order. `mass_optimisation` runs after
`load_constraint` in `verification/README.md`. The shared code is
`_common.py` (problem construction, the recorder, the KKT estimate, the
blade record, the figures) -- lifted from
`load_constraint/run_constrained_slsqp.py`, not copied a fourth time.

## The problem

    minimise    m(u) = k_P int c dr / (k_P int c dr)(x0)        shell material, geometric, exact gradient
    subject to  AEP(u) / AEP(x0) - (1 - delta)   >= 0            energy floor, objective adjoint
                KS0 - KS(u)                      >= 0            Phase 4 moment cap (eps = 0), moment adjoint
                KS0 / c00^2 - KS(u) / c0(u)^2    >= 0            root-stress proxy, moment adjoint
                D0 - D(u)                        >= 0            tip-deflection proxy, deflection adjoint
                c_i - c_{i+1} >= 0 (4), theta_i - theta_{i+1} >= 0 (4)   monotone control points, constant Jacobian
                c_i - 0.060 >= 0 (5)                                     the buildable-tip floor, constant Jacobian
                envelope (50 rows), solidity (25 rows), u in [0, 1]^10

`ScaledProblem.constraints_for_mass_problem(delta)`; SLSQP, `ftol = 1e-10`,
`maxiter = 400`. One nine-point load solve per accepted iterate serves the
moment, stress and deflection rows (`deflection_solves_total = 0` in every
run). The state rows are guarded: a trial point a row cannot evaluate
returns `-1e3` and is logged in `evaluation_failures` -- **no run logged
one**. The rows are what `docs/adjoint_derivation.md` sections 10 and 11
derive; the material model is `src/objective/mass.py`
(`k_P = 2.048504996756389`, `k_A = 0.06850241963050001`, exact from the
SG6043 coordinates); the solid proxy `k_A int c^2 dr` is reported beside the
shell figure, never optimised.

## The result (`result_delta0.json`, `blade_delta0_geometry.png`, `blade_delta0_loads.png`, `blades_rendered_iso.png`, `blades_rendered_plan.png`)

**Two headline numbers, from one sweep** (MJ, 2026-09-20, plan "Decision 7"):

- **at Schmitz's energy exactly (`delta = 0`): `x_m` carries -3.379 % shell
  material (-7.278 % on the solid proxy)** -- floor slack `-3.7e-13`;
- **within 0.5 % of Schmitz's energy (`delta = 0.5 %`): -5.271 % shell
  (-9.666 % solid)**.

SLSQP exit 0 in 18 iterations, 14 s; 31 objective solves, 19 objective
adjoints, 22 load solves.

| | `x0` (Schmitz) | `x_c` (energy optimum, Phase 4) | **`x_m` (mass optimum, delta = 0)** |
|---|---|---|---|
| AEP vs `x0` | 0 | +0.1465 % | **-0.0000 %** (slack -3.7e-13) |
| shell material `k_P int c dr` vs `x0` | 0 | +6.171 % | **-3.379 %** |
| solid proxy `k_A int c^2 dr` vs `x0` | 0 | +14.196 % | **-7.278 %** |
| `KS / KS0` (moment cap) | 1.0000 | 1.0000 (active) | 0.9685 (inactive, slack 0.032) |
| stress proxy ratio | 1.0000 | 0.8451 | 1.0000 (active) |
| deflection proxy ratio | 1.0000 | 0.7618 | 1.0000 (active) |
| rated-point root moment [N m] | 177.38 | 177.41 | 171.74 (-3.18 %) |
| rated-point rotor thrust [N] | 517.5 | 525.5 | 505.5 (-2.3 %) |
| cut-in `Re_min` | 73 170 | 62 149 | 71 127 (floor 40 000) |
| chord control points [mm] | 276 / 182 / 97 / 79 / 67 | 300 / 182 / 126 / 78 / 48 | **271 / 150 / 135 / 60 / 60** |
| twist control points [deg] | 24.5 / 9.7 / 4.0 / 1.0 / 0.6 | 24.3 / 10.6 / 5.2 / 2.2 / -1.9 | **23.7 / 8.1 / 6.2 / -0.0 / -0.1** |
| taper ratio root / tip | 4.1 | 6.3 | **4.5** |

(`reference_blades.json`, `reference_blades.png`; the `x_c` row is the
"+6.2 % material for +0.15 % energy" motivation of the re-pitch, measured on
the same code. `blades_rendered_iso.png` and `blades_rendered_plan.png` draw
the three blades in isometric and plan view, and `blades_rendered_edge.png`
the edge-on view; the root cylinder and
the tip rounding are drawn only: neither is modelled, and the blade the
numbers describe starts at the 0.30 m cut-out.)

**Where the material comes off.** The outer third: the two tip control
points sit on their 60 mm rows (multipliers 0.34 and 0.97) and their
monotone row is tight (equal chords), so the blade is a straight 60 mm from
r ~ 1.5 m out; the root stays close to Schmitz (271 vs 276 mm) because the
stress row pins `KS / c0^2` and the deflection row pins
`int M (R - r) / c^3 dr`. The mid-span control point rises (135 vs 97 mm)
to make back the energy the tip gives up, which is what the `1 / c^3` panel
of `blade_delta0_loads.png` says the deflection row wants. The tip twist goes to
about 0 deg (Schmitz 0.6 deg).

**Active set and KKT** (`kkt`, non-negative least squares -- the active set
is degenerate, `e_3`, `e_4` and their monotone row are dependent, so plain
least squares can return a spurious sign): `aep_floor`, `stress`,
`deflection`, `min chord chord_3`, `min chord chord_4`,
`monotone chord_3-chord_4` (multiplier 0); residual after projection
`1.3e-6` on a gradient of norm 0.93. The AEP row's multiplier is the
**exchange rate at the margin:
11.0 % material per 1 % energy** (both the objective and the row are
fractions of `x0`'s, so the multiplier needs no unit conversion). The
moment cap is inactive: the load is 3.2 % below Schmitz's without being
asked.

**Post-checks.** `alpha` on the uncapped points within `[1.41, 8.80] deg`
(band `[-8, 18]`); `reynolds_min = 71 127` against the 40 000 floor, with the
tightest envelope row (`floor r = 1.966`) 25.5 mm slack; every station
converged; no `PolarDomainError`, no margin escalation.

**Sanity band** (stated in the plan before the first run): shell saving
2-10 %, solid 5-20 %, AEP floor active with a non-negative multiplier,
moment cap inactive, no failed trials, `alpha` within, `Re_min` above the
floor. All met (`sanity_notes` empty); `-3.38 %` is below the scratch's
deflection-free `-7.5 %`, as it must be.

## Multi-start (`multistart_delta0.json`, `multistart_delta0.png`)

Ten starts -- `x0`, `x_c` (below the floor, an infeasible start), and the
eight committed random starts of `fd_optimisation_multistart/starts.json`
(drawn for the energy problem; most are infeasible for this one and SLSQP
finds the feasible set first). **All ten exit 0 at the same point**: spread
in `u` `1.1e-6` (criterion `0.0166`), spread in the objective `3.4e-11`, no
failed trial point in any run. 185 s.

## The energy-floor sweep (`pareto.json`, `pareto_front.png`)

Cold from `x0` and warm from the previous floor at each `delta`; cold and
warm agree to `<= 5.7e-6` in `u` everywhere (criterion `0.0166`).

| `delta` | energy given up | shell saved | solid saved | exchange rate (KKT) | secant to the next point | tip chords [mm] | rows beyond AEP / stress / deflection / the two min-chord rows |
|---|---|---|---|---|---|---|---|
| 0 | 0 | **3.379 %** | 7.278 % | 11.02 | 4.90 | 60 / 60 | monotone `chord_3-chord_4` |
| 0.25 % | 0.2500 % | 4.605 % | 8.557 % | 3.11 | 2.66 | 60 / 60 | + monotone `twist_3-twist_4` |
| **0.5 %** | 0.5000 % | **5.271 %** | **9.666 %** | 2.31 | 1.83 | 60 / 60 | + monotone `twist_2-twist_3` |
| 1 % | 1.0000 % | 6.187 % | 11.490 % | 1.56 | 1.34 | 60 / 60 | same |
| 2 % | 2.0000 % | 7.530 % | 14.294 % | 1.22 | -- | 60 / 60 | same |

The front is concave, the saving is monotone in `delta`, and **every secant
lies between the KKT exchange rates at its two ends** (4.90 between 11.02
and 3.11, and so on) -- the adjoint multiplier and the front's own slope
check each other, which is what the tangents in `pareto_front.png` show. The
stress and deflection rows and the two min-chord rows stay active along the
whole front; from `delta = 0.25 %` the monotone twist rows bind too (the
outer twist control points equalise). The first quarter-percent of energy
buys the most material (1.2 points for 0.25 %); after 1 % the rate is
near 1:1. With the 45 mm floor the front ran 4.2 -> 10.6 %; the 60 mm floor
costs 0.9 points at `delta = 0` and 2.4 at `delta = 0.5 %` -- the price of
a buildable tip.

Read against decision 2 of the plan (the reference stays at `lambda = 6.5`):
a `lambda = 6.0` Schmitz would raise the floor by about 0.16 %, which at the
first secant (4.9 points per %) costs roughly 0.8 of the 3.4 points of shell
material. Stated, not re-run.

## Ablation at `delta = 0` (`ablation.json`, `pareto_ablation.png`)

| rows | shell saved | solid saved | stress ratio | deflection ratio | `KS / KS0` |
|---|---|---|---|---|---|
| full set | **3.379 %** | 7.278 % | 1.000 | 1.000 | 0.968 |
| without the deflection row | 4.303 % | 9.750 % | 1.000 | 1.078 | 0.982 |
| without the stress row | 4.527 % | 11.997 % | 1.702 | 1.000 | 0.971 |
| without both proxies | 6.532 % | 17.086 % | 1.981 | 1.131 | 0.989 |
| without the manufacturing block (monotone + min chord) | 4.238 % | 8.532 % | 1.000 | 1.000 | 0.956 |
| without the moment cap | 3.379 % | 7.278 % | 1.000 | 1.000 | 0.968 |

Removing a row never shrinks the saving. The `no moment cap` optimum is the
full-set optimum to `1e-9` points (the cap is slack). Without the
manufacturing block the morning's 45 mm optimum returns exactly (`-4.238 %`,
tip on the box bound): **the min-chord floor costs 0.86 points**, and the
monotone rows on their own cost nothing at `delta = 0` (the morning
ablation showed `no monotone == full` with the 45 mm box).
**The stress row is doing the work**: without it the optimiser thins the
root (stress ratio 1.70) for another 1.1 points of shell; without the
deflection row it accepts 7.8 % more tip deflection for 0.9 points. Without
both, `-6.5 %` returns with a root that carries 2.0x Schmitz's stress proxy
-- the number the re-pitch proposal warned was not a blade. The saving that
survives both proxies and the buildable tip, `-3.4 %`, is the defensible one.

## Checks (`checks.json`)

At `x0` and at `x_m`:

| check | `x0` | `x_m` |
|---|---|---|
| shell / solid gradient vs central FD (relative) | 1.7e-11 / 4.0e-12 | 2.6e-11 / 2.8e-12 |
| manufacturing rows exact to the bit; min slack | yes; 0.0070 (`x0`'s tip clears the floor by 7 mm) | yes; -2.4e-17 (the tight rows) |
| forward twins: deflection / moment (relative) | 2.2e-16 / 2.2e-16 | 2.2e-16 / 2.2e-16 |
| Tier 1 `dD/dphi`, `dD/dd`, `dKS/dx`, `dKS/dd` (mixed) | 1.1e-15, 6.7e-16, 8.3e-16, 5.8e-16 | 7.2e-16, 5.7e-16, 7.2e-16, 3.9e-16 |
| Tier 2 tangent vs adjoint; assembly identity | 3.9e-16; 0 | 3.4e-16; 0 |
| Tier 3 stress row, worst `abs(adj - fd) / eps_j` | 0.366 (`chord_2`) | 1.014 (`twist_0`) |
| Tier 3 deflection row, worst `abs(adj - fd) / eps_j` | 0.812 (`twist_1`) | 0.291 (`twist_2`) |
| Taylor remainder ratios per decade (six directions) | 97-100 | 97-100 |

Acceptance `3 eps_j`; every variable of both rows passes at both points, so
the round-off floor was not measured (the Phase 4 rule: measured only on
failure). The `x0` numbers are the morning's exactly, because the box -- and
so the scaling `h*` and `eps_j` live in -- did not move; the `eps_j`
fragility of `gradient_verification/` (the red `chord_4` at `x0`) did not
reappear on these rows at `x_m`.

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
| central resource (`k = 1.709, c = 7.274`), 300 rpm | +0.147 % | **-0.000 %** | -3.18 % |
| `k = 1.5, c = 7.049` | +0.135 % | -0.000 % | (unchanged) |
| `k = 1.5, c = 7.633` | +0.132 % | -0.006 % | |
| `k = 2.0, c = 7.049` | +0.163 % | +0.011 % | |
| `k = 2.0, c = 7.633` | +0.162 % | -0.002 % | |
| no ceiling | +0.063 % | +0.084 % | -2.49 % |
| ceiling 60 m/s | +0.262 % | -0.100 % | -3.48 % |
| ceiling 55 m/s | +0.835 % | -0.444 % | -3.89 % |
| ceiling 50 m/s | +2.342 % | -1.840 % | -4.27 % |

**Resource: robust** -- over the four corners `x_m` stays within
`[-0.006, +0.011] %` of `x0`: the per-bin energy differences between the
two blades are small and of both signs, so re-weighting the bins moves the
total very little. **Ceiling: not robust** -- a lower rotor-speed ceiling
pushes the rated point to lower `lambda`, where the thinner tip of `x_m`
loses more than Schmitz's, and by 50 m/s it gives up 1.8 % of energy (while
`x_c`, with its fat mid-span, gains 2.3 %); the 60 mm floor made this
milder than the 45 mm blade's -2.6 %. The mass optimum is an optimum *for
the 300 rpm machine*; the honest statement of the headline is "-3.4 %
material at equal energy under the configured law". The resource band used
here (`c` from the log-law cross-check, `k +/- 0.25`) is an assumption of
this artefact -- plan 1.3's sensitivity band is still not chosen and this
does not choose it.

## What would make it wrong

- **The material proxy.** `k_P int c dr` is a constant-thickness shell of one
  section with no spar, no root insert, no web, no ply drop -- mass is
  proportional to planform area. A real laminate schedule would weight the
  root more (and a solid or spar-dominated blade would follow `c^2` or
  `c^3 t`); the solid proxy is reported beside it as the other end of the
  range (-7.3 % vs -3.4 %). No kg figure is quoted: `structure.*` in
  `config/rotor_design.yaml` is `TODO` and only multiplies this number.
- **The tip floor is a judgement.** 60 mm = "a 6 mm section can carry two
  2 mm skins, a bond and a web" -- one sentence, no laminate behind it. The
  whole outer third of the optimum sits on it, so it is the number that
  sets the headline: without it (the 45 mm box) -4.2 %, at 80 mm ~-0.8 %
  at `delta = 0.5 %` against an unfloored Schmitz (`misc/blade_shape_experiment`).
  It is a row, not the box, so the Phase 1-4 scaling and artefacts stand.
- **The loads are operating-range only (B3).** `L` is the nine points with
  `P <= P_rated` plus the rated point; no parked, gust, yawed or fault case.
  A survival load case would set the root, not the rated point, and the
  stress row would be a different row.
- **The three structural rows are relative proxies**, all on the same thin-
  shell assumption (`Z ~ c^2 t`, `I ~ c^3 t`, mass `~ c`), pinned to
  Schmitz's values. They say the blade is no worse than Schmitz on each
  measure; they do not say Schmitz was adequate. No fatigue, no buckling
  (plan, cut list).
- **The Reynolds floor.** `Re_min = 71 127` against the 40 000 cache floor,
  25.5 mm of envelope slack at the tip station -- not close, with the
  60 mm chord.
- **The monotone rows** are on the control points (variation diminishing),
  which is sufficient for a monotone spline but not necessary; a mould
  constraint might allow shapes these forbid, or forbid shapes these allow
  (the front's optima have equal tip chords and, from 0.25 %, equal outer
  twists).
- **The parameterisation**: 5 + 5 cubic B-spline control points, selected in
  `representation_study/` for the energy objective. A finer chord basis
  might find more material to remove at the same energy; a coarser one
  less.
- **The root cylinder in the renders** is drawn, not modelled: no mass, no
  load, no drag. The blade the numbers describe starts at r = 0.30 m.
- **One resource, one law** -- see the cross-evaluation.

## How much energy the shape can move, and the multipliers in common units (`energy_split.json`, added 2026-09-26)

Review roadmap item 6 and the minor item on multiplier scaling.
`run_energy_split.py` splits each blade's AEP into the bins capped at the
fixed rating (a constant whatever the shape) and the rest:

| | AEP [MWh/yr] | capped (fixed by the rating) | design-dependent | capped bins |
|---|---|---|---|---|
| `x0` | 10.2477 | 4.2902 (41.9 %) | 5.9575 (58.1 %) | 9 |
| `x_c` | 10.2627 | 4.2902 | 5.9725 | 9 |
| `x_m` | 10.2477 | 4.2902 | 5.9575 | 9 |

On the design-dependent basis the energy optimum's gain is **+0.252 %**
(+0.1465 % of total AEP), and the marginal exchange rate at the reference
energy is **6.41 % shell material per 1 % of design-dependent energy**
(11.02 % per 1 % of total AEP). All three blades cap the same nine bins, so
the split is the same for each.

The rating (3822 W) is `x0`'s power at 11 m/s and lambda = 6.5, a point the
300 rpm ceiling forbids. Under the operating law `x0` first reaches it at
**11.22 m/s** (lambda = 5.60), `x_c` at 11.13 m/s and `x_m` at 11.27 m/s
(`blades.*.rating_first_reached_at_ms`); 11 m/s is the load set's rated
point, not where rated power is first reached.

The KKT multipliers of `x_m`, each multiplied by its row's reference scale so
that they share a unit (percent of shell material per percent of the row's
reference value; per millimetre for the floor rows):

| row | raw multiplier | row scale | normalised |
|---|---|---|---|
| AEP floor | 11.02 | 1 (fraction of AEP(x0)) | **11.0 %/%** |
| stress proxy | 0.00331 | `KS0 / c00^2` = 13.15 m^-2 | **0.044 %/%** |
| deflection proxy | 0.198 | `D0` = 1.0004 | **0.198 %/%** |
| floor, chord 3 | 0.336 | per m of chord | **0.034 %/mm** |
| floor, chord 4 | 0.969 | per m of chord | **0.097 %/mm** |
| monotone chord 3-4 | 0 | -- | 0 (dependent row) |

These are local rates. The ablation measures the total: removing the stress
row saves 1.15 points only because the stress proxy then rises by 70 %, at a
marginal price of 0.044 % material per 1 % stress.

    python verification/mass_optimisation/run_energy_split.py     # ~10 s
