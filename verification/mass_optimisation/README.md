# The minimum-material blade

This is the project's main result: the blade with the least shell material
that still delivers the Schmitz reference blade's annual energy, without
increasing the root moment, root stress or tip deflection, and with a
buildable tip.

**The result.** At exactly the Schmitz blade's energy, the optimised blade
x_m uses **3.379 % less shell material** (7.278 % less on the solid-section
measure). If 0.5 % of the energy can be given up, the saving rises to
**5.271 %** (9.666 % solid).

All runs use the 300 rpm operating law (`λ(V) = min(6.5, Ω_max R / V)`,
ceiling reached at 9.67 m/s), the fixed rating of 3822.189755449124 W, the
design bounds (chord 0.045–0.30 m, twist −2° to 35°, local solidity at most
0.5) and the polar-table envelope with a 0.05 margin. Every reference value
comes from the Schmitz blade x₀ (`verification/baseline/x0.json`): its AEP,
its root moment (`KS0`, `M_ref`), its root chord `c00`, its tip deflection
(`delta_ref`, `D0`), and its shell area of 0.43705266291483963 m², which
normalises the objective.

## Reproducing it

    python verification/mass_optimisation/run_mass_slsqp.py --delta 0     # ~20 s, the main optimum
    python verification/mass_optimisation/run_multistart.py               # ~3 min, 10 starts at delta = 0
    python verification/mass_optimisation/run_sweep.py                    # ~4 min, the energy sweep and the ablation
    python verification/mass_optimisation/run_mass_checks.py              # ~1 min, gradient checks, reference table, renders
    python verification/mass_optimisation/run_cross_evaluation.py         # ~1 min, the three blades under other conditions
    python verification/mass_optimisation/run_energy_split.py             # ~10 s, the energy split and scaled multipliers

These read `verification/baseline/x0.json`,
`verification/load_constraint/result_eps0.json` (x_c),
`verification/fd_optimisation_multistart/results.json` (the 0.0166 agreement
tolerance in u) and `verification/fd_step_size/sweep.json` (h\*).
`run_mass_checks.py` and `run_cross_evaluation.py` also read
`result_delta0.json`, so run them in the order above. The shared code
(problem setup, the iteration recorder, the KKT estimate, the figures) is in
`_common.py`.

All the figures can be redrawn from the committed JSON, without solving
anything, with `python verification/mass_optimisation/replot_figures.py`.

## The problem

    minimise    m(u) = k_P ∫ c dr / (k_P ∫ c dr)(x0)          shell material (exact gradient)
    subject to  AEP(u) / AEP(x0) − (1 − δ)   ≥ 0             energy floor (energy adjoint)
                KS0 − KS(u)                  ≥ 0             root-moment cap (moment adjoint)
                KS0 / c00² − KS(u) / c0(u)²  ≥ 0             root-stress proxy (moment adjoint)
                D0 − D(u)                    ≥ 0             tip-deflection proxy (deflection adjoint)
                c_i − c_{i+1} ≥ 0 (4), θ_i − θ_{i+1} ≥ 0 (4) chord and twist decrease toward the tip
                c_i − 0.060 ≥ 0 (5)                          60 mm minimum buildable chord
                envelope (50 rows), solidity (25 rows), u ∈ [0, 1]¹⁰

The problem is built by `ScaledProblem.constraints_for_mass_problem(delta)`
and solved with SLSQP (`ftol = 1e-10`, `maxiter = 400`). One nine-point load
solve per accepted iterate serves the moment, stress and deflection
constraints. If a trial point can't be evaluated, the constraint returns
−1000 and the failure is logged in `evaluation_failures`; no run logged one.

The constraint derivations are in `docs/adjoint_derivation.md` §10 and §11.
The material model is in `src/objective/mass.py`, with
`k_P = 2.048504996756389` and `k_A = 0.06850241963050001` calculated exactly
from the SG6043 coordinates. The solid measure `k_A ∫ c² dr` is reported
alongside the shell figure but never optimised.

### Why there is a 60 mm minimum chord

Without it, the optimum sat on the 45 mm lower chord bound at the tip
(−4.24 % shell, a 6:1 taper). But 45 mm is the thickness of the *whole*
SG6043 section at that chord (about 4.5 mm), leaving no room for two skins,
a bond line and a web, so the optimiser had found a gap in the constraints.
In an exploratory study I tried a 60, 80 and 100 mm floor, a taper-ratio
limit and convexity constraints on the chord:

- **60 mm** (a 6 mm thick section: two 2 mm skins, a bond line and a web) is
  feasible at δ = 0, saves 3.4 %, and is below x₀'s 67 mm tip, so the
  comparison with Schmitz is unaffected. I used this one.
- **80 mm** is infeasible at δ = 0 (the best such blade falls 0.002 % short
  of the Schmitz energy) and saves 0.8 % at δ = 0.5 %. Because it is above
  x₀'s tip, a fair comparison would also need the Schmitz blade refitted
  with the same floor. Not used.
- The taper-ratio and convexity constraints were never active once the 60 mm
  floor was in place, so I didn't add them.

The floor is written as five constraints rather than by raising the box's
lower bound, because the box also sets the scaling of the design variables.
Every earlier gradient and optimisation result (the step sizes, the noise
floor, x\*, x_c) was measured in that scaling. When I tried raising the box
to 0.060, the committed finite-difference reference no longer matched the
adjoint on the chord variables, only because the scaling had changed. Keeping
the floor as constraints (`ScaledProblem.min_chord_rows`,
`manufacturing.min_chord_m` in the config) leaves everything else valid. The
energy optimum x_c has a 47.7 mm tip, so it violates the floor by 12 mm. It
is kept as the energy reference and also used as an infeasible starting
point, which SLSQP handles.

## The result

`result_delta0.json`, `blade_delta0_geometry.png`, `blade_delta0_loads.png`,
`blades_rendered_iso.png`, `blades_rendered_plan.png`

SLSQP exits normally (code 0) after 18 iterations and 14 s, using 31 objective
solves, 19 energy adjoints and 22 load solves.

| | x₀ (Schmitz) | x_c (energy optimum) | **x_m (minimum material, δ = 0)** |
|---|---|---|---|
| AEP vs x₀ | 0 | +0.1465 % | **−0.0000 %** (slack −3.7e-13) |
| shell material `k_P ∫ c dr` vs x₀ | 0 | +6.171 % | **−3.379 %** |
| solid measure `k_A ∫ c² dr` vs x₀ | 0 | +14.196 % | **−7.278 %** |
| `KS / KS0` (moment cap) | 1.0000 | 1.0000 (active) | 0.9685 (inactive, slack 0.032) |
| stress ratio | 1.0000 | 0.8451 | 1.0000 (active) |
| deflection ratio | 1.0000 | 0.7618 | 1.0000 (active) |
| rated-point root moment [N m] | 177.38 | 177.41 | 171.74 (−3.18 %) |
| rated-point rotor thrust [N] | 517.5 | 525.5 | 505.5 (−2.3 %) |
| lowest Reynolds number (cut-in) | 73 170 | 62 149 | 71 127 (floor 40 000) |
| chord control points [mm] | 276 / 182 / 97 / 79 / 67 | 300 / 182 / 126 / 78 / 48 | **271 / 150 / 135 / 60 / 60** |
| twist control points [°] | 24.5 / 9.7 / 4.0 / 1.0 / 0.6 | 24.3 / 10.6 / 5.2 / 2.2 / −1.9 | **23.7 / 8.1 / 6.2 / −0.0 / −0.1** |
| taper ratio (root / tip) | 4.1 | 6.3 | **4.5** |

(`reference_blades.json`, `reference_blades.png`.) The x_c column shows why
the objective changed: the energy optimum gains 0.15 % energy but needs 6.2 %
more material. The rendered figures show the three blades in isometric, plan
and edge-on views (`blades_rendered_edge.png`). The root cylinder and the
rounded tip are only drawn; the modelled blade starts at the 0.30 m cut-out.

**Where the material comes from.** The outer third. The two tip control
points sit on the 60 mm floor (multipliers 0.34 and 0.97) and are equal, so
the blade is a constant 60 mm chord from about r = 1.5 m outward. The root
stays close to Schmitz (271 vs 276 mm), because the stress constraint fixes
`KS / c0²` and the deflection constraint fixes `∫ M (R − r) / c³ dr`. The
mid-span control point grows (135 vs 97 mm) to recover the energy lost at the
tip, which is what the `1/c³` panel of `blade_delta0_loads.png` shows the
deflection constraint favours. The tip twist goes to about 0° (Schmitz has
0.6°).

**Active constraints and KKT multipliers.** The active set is the AEP floor,
stress, deflection, the minimum chord at control points 3 and 4, and the
chord monotonicity between 3 and 4 (with a zero multiplier). The multipliers
come from non-negative least squares (`kkt` in the JSON), because the last
three constraints are linearly dependent and plain least squares can give a
spurious sign. The residual is 1.3e-6 on a gradient of norm 0.93. The AEP
multiplier is the **marginal exchange rate: 11.0 % material per 1 % of
energy**. Both the objective and the constraint are fractions of x₀'s
values, so no unit conversion is needed. The moment cap isn't active: the
root moment ends up 3.2 % below Schmitz's without being asked to.

**Other checks.** The angle of attack at the operating points below rated
stays within 1.41°–8.80° (the polar's valid range is −8° to 18°). The lowest
Reynolds number is 71 127 against the 40 000 floor, and the tightest envelope
constraint has 25.5 mm of slack. Every station converged, with no polar-range
errors.

Before the first run I set out the range the result should fall in: shell
saving 2–10 %, solid 5–20 %, the AEP floor active with a non-negative
multiplier, the moment cap inactive, no failed trial points, and the angle of
attack and Reynolds number in range. All were met. The −3.38 % is also below
the −7.5 % the exploratory study found without a deflection constraint, as
it should be.

## Multi-start

`multistart_delta0.json`, `multistart_delta0.png`

Ten starting points: x₀, x_c (which violates the tip floor), and the eight
random starts committed in `fd_optimisation_multistart/starts.json`. The
random starts were drawn for the energy problem, and most are infeasible
here, so SLSQP first has to find the feasible region. **All ten converge to
the same point**: the spread is 1.1e-6 in u (tolerance 0.0166) and 3.4e-11 in
the objective, with no failed trial points. The runs take 185 s in total.

## The energy–material trade-off

`pareto.json`, `pareto_front.png`

Each energy floor δ is solved twice, from x₀ and warm-started from the
previous floor; the two agree to within 5.7e-6 in u everywhere (tolerance
0.0166).

| δ | energy given up | shell saved | solid saved | marginal rate (KKT) | secant to next point | tip chords [mm] | other active constraints |
|---|---|---|---|---|---|---|---|
| 0 | 0 | **3.379 %** | 7.278 % | 11.02 | 4.90 | 60 / 60 | chord monotonicity 3–4 |
| 0.25 % | 0.2500 % | 4.605 % | 8.557 % | 3.11 | 2.66 | 60 / 60 | plus twist monotonicity 3–4 |
| **0.5 %** | 0.5000 % | **5.271 %** | **9.666 %** | 2.31 | 1.83 | 60 / 60 | plus twist monotonicity 2–3 |
| 1 % | 1.0000 % | 6.187 % | 11.490 % | 1.56 | 1.34 | 60 / 60 | same |
| 2 % | 2.0000 % | 7.530 % | 14.294 % | 1.22 | – | 60 / 60 | same |

(The AEP, stress, deflection and both minimum-chord constraints are active
throughout.)

The curve is concave and the saving increases steadily with δ. **Every
secant lies between the KKT rates at its two ends** (4.90 lies between 11.02
and 3.11, and so on), so the adjoint multipliers and the curve's own slope
confirm each other; the tangents drawn in `pareto_front.png` show this. The
first quarter per cent of energy buys the most material (1.2 points), and
beyond 1 % the rate is close to 1:1. With the 45 mm floor instead, the
savings ran from 4.2 to 10.6 %, so the buildable tip costs 0.9 points at
δ = 0 and 2.4 points at δ = 0.5 %.

The reference blade uses λ = 6.5. A Schmitz blade designed for λ = 6.0 would
raise the energy floor by about 0.16 %, which at the first secant (4.9 points
per %) would cost roughly 0.8 of the 3.4 points. I estimated this but didn't
re-run it.

## Which constraints matter (ablation at δ = 0)

`ablation.json`, `pareto_ablation.png`

| constraints | shell saved | solid saved | stress ratio | deflection ratio | `KS / KS0` |
|---|---|---|---|---|---|
| all | **3.379 %** | 7.278 % | 1.000 | 1.000 | 0.968 |
| without deflection | 4.303 % | 9.750 % | 1.000 | 1.078 | 0.982 |
| without stress | 4.527 % | 11.997 % | 1.702 | 1.000 | 0.971 |
| without stress and deflection | 6.532 % | 17.086 % | 1.981 | 1.131 | 0.989 |
| without the manufacturing constraints (monotonicity and minimum chord) | 4.238 % | 8.532 % | 1.000 | 1.000 | 0.956 |
| without the moment cap | 3.379 % | 7.278 % | 1.000 | 1.000 | 0.968 |

Removing a constraint never reduces the saving. Without the moment cap the
optimum is the same to within 1e-9 points, since the cap isn't active.
Without the manufacturing constraints, the 45 mm-tip optimum returns exactly
(−4.238 %), so **the minimum chord costs 0.86 points**. The monotonicity
constraints on their own cost nothing at δ = 0.

**The stress constraint does the most work.** Without it the optimiser thins
the root (stress ratio 1.70) for another 1.1 points of material. Without the
deflection constraint it accepts 7.8 % more tip deflection for 0.9 points.
Without both, the saving reaches 6.5 %, but the root then carries twice the
Schmitz blade's stress, which wouldn't be a usable blade. The saving that
survives both structural constraints and the buildable tip, 3.4 %, is the
one I'd defend.

## Gradient and consistency checks

`checks.json`, at x₀ and at x_m:

| check | x₀ | x_m |
|---|---|---|
| shell / solid gradient vs central FD (relative) | 1.7e-11 / 4.0e-12 | 2.6e-11 / 2.8e-12 |
| manufacturing constraints exact; smallest slack | yes; 0.0070 (x₀'s tip clears the floor by 7 mm) | yes; −2.4e-17 (the active ones) |
| forward re-evaluation: deflection / moment (relative) | 2.2e-16 / 2.2e-16 | 2.2e-16 / 2.2e-16 |
| Tier 1: `dD/dphi`, `dD/dd`, `dKS/dx`, `dKS/dd` | 1.1e-15, 6.7e-16, 8.3e-16, 5.8e-16 | 7.2e-16, 5.7e-16, 7.2e-16, 3.9e-16 |
| Tier 2: tangent vs adjoint; assembly identity | 3.9e-16; 0 | 3.4e-16; 0 |
| Tier 3: stress constraint, worst `|adj − fd| / ε_j` | 0.366 (`chord_2`) | 1.014 (`twist_0`) |
| Tier 3: deflection constraint, worst `|adj − fd| / ε_j` | 0.812 (`twist_1`) | 0.291 (`twist_2`) |
| Taylor remainder ratio per decade (six directions) | 97–100 | 97–100 |

The acceptance threshold is 3 ε_j, and every variable passes for both
constraints at both points, so the round-off floor didn't need measuring
(it is only measured when a check fails). The single noise-floor failure in
`gradient_verification/` (`chord_4` at x₀) doesn't reappear here.

Re-evaluating x_m from its JSON through the forward path alone (the
objective, root moment, tip deflection and shell area) reproduces the
recorded AEP, rated root moment, rated tip deflection and shell area exactly.

## Robustness

`cross_evaluation.json`. There is no re-optimisation here: the three blades
are simply evaluated under other conditions, to see how much of x_m's energy
parity survives.

| condition | x_c AEP vs x₀ | **x_m AEP vs x₀** | x_m rated moment vs x₀ |
|---|---|---|---|
| design resource (`k = 1.709, c = 7.274`), 300 rpm | +0.147 % | **−0.000 %** | −3.18 % |
| `k = 1.5, c = 7.049` | +0.135 % | −0.000 % | (unchanged) |
| `k = 1.5, c = 7.633` | +0.132 % | −0.006 % | |
| `k = 2.0, c = 7.049` | +0.163 % | +0.011 % | |
| `k = 2.0, c = 7.633` | +0.162 % | −0.002 % | |
| no ceiling | +0.063 % | +0.084 % | −2.49 % |
| ceiling 60 m/s | +0.262 % | −0.100 % | −3.48 % |
| ceiling 55 m/s | +0.835 % | −0.444 % | −3.89 % |
| ceiling 50 m/s | +2.342 % | −1.840 % | −4.27 % |

**It is robust to the wind resource.** Over the four corners, x_m stays
within −0.006 to +0.011 % of x₀. The bin-by-bin energy differences between
the two blades are small and of both signs, so re-weighting the bins barely
moves the total.

**It is not robust to the rotor-speed ceiling.** A lower ceiling pushes the
rated point to a lower tip-speed ratio, where x_m's thinner tip loses more
than Schmitz's does. At 50 m/s it gives up 1.8 % of energy, while x_c, with
its wider mid-span, gains 2.3 %. (The 60 mm floor makes this milder than the
−2.6 % of the 45 mm-tip blade.) So x_m is an optimum *for the 300 rpm
machine*, and the fair statement of the result is "−3.4 % material at equal
energy under the configured operating law". The resource range used here
(`c` from the log-law cross-check, `k` ± 0.25) is itself an assumption of
this check.

## How much of the energy the blade shape can change

`energy_split.json` (`run_energy_split.py`) splits each blade's AEP into the
bins above rated, where power is fixed by the rating whatever the shape, and
the rest:

| | AEP [MWh/yr] | above rated (fixed) | shape-dependent | bins above rated |
|---|---|---|---|---|
| x₀ | 10.2477 | 4.2902 (41.9 %) | 5.9575 (58.1 %) | 9 |
| x_c | 10.2627 | 4.2902 | 5.9725 | 9 |
| x_m | 10.2477 | 4.2902 | 5.9575 | 9 |

Measured against the shape-dependent part only, the energy optimum's gain is
**+0.252 %** (+0.1465 % of total AEP), and the marginal exchange rate at the
reference energy is **6.41 % shell material per 1 % of shape-dependent
energy** (11.02 % per 1 % of total AEP). All three blades reach rated power
in the same nine bins, so the split is the same for each.

The rating (3822 W) is x₀'s power at 11 m/s and λ = 6.5, but the 300 rpm
ceiling doesn't allow that point. Under the operating law, x₀ first reaches
rated power at **11.22 m/s** (λ = 5.60), x_c at 11.13 m/s and x_m at
11.27 m/s (`blades.*.rating_first_reached_at_ms`). So 11 m/s is the rated
point of the load set, not where rated power is first reached.

The KKT multipliers of x_m, each scaled by its constraint's reference value
so they share a unit (per cent of shell material per per cent of the
constraint's reference value, or per millimetre for the floor):

| constraint | raw multiplier | scale | scaled |
|---|---|---|---|
| AEP floor | 11.02 | 1 (fraction of AEP(x0)) | **11.0 %/%** |
| stress | 0.00331 | `KS0 / c00²` = 13.15 m⁻² | **0.044 %/%** |
| deflection | 0.198 | `D0` = 1.0004 | **0.198 %/%** |
| floor, chord 3 | 0.336 | per m of chord | **0.034 %/mm** |
| floor, chord 4 | 0.969 | per m of chord | **0.097 %/mm** |
| chord monotonicity 3–4 | 0 | – | 0 (dependent) |

These are local rates. The ablation measures the total effect: removing the
stress constraint saves 1.15 points only because the stress then rises by
70 %, at a marginal price of 0.044 % material per 1 % of stress.

## What would make it wrong

- **The material measure.** `k_P ∫ c dr` treats the blade as a single
  constant-thickness shell of one section, with no spar, root insert, web or
  ply drops, so mass is proportional to planform area. A real laminate
  schedule would put more weight at the root, and a solid or spar-dominated
  blade would scale more like `c²` or `c³ t`. The solid measure is reported
  as the other end of that range (−7.3 % vs −3.4 %). The kilogram figures for
  one assumed laminate are in `verification/absolute_material/`.
- **The tip floor is a judgement.** "A 6 mm section can hold two 2 mm skins,
  a bond line and a web" is a single argument, not a laminate design. The
  whole outer third of the optimum sits on it, so it largely sets the result:
  −4.2 % without it (the 45 mm bound), and about −0.8 % at δ = 0.5 % with an
  80 mm floor. Because it is a constraint and not the box, the earlier
  gradient results stay valid.
- **Loads cover the operating range only.** The load set is the nine points
  at or below rated power; there are no parked, gust, yawed or fault cases.
  A survival load case would size the root instead of the rated point, and
  the stress constraint would look different.
- **The structural constraints are relative.** All three rest on the same
  thin-shell assumption (`Z ~ c² t`, `I ~ c³ t`, mass `~ c`) and are pinned
  to the Schmitz values. They say the blade is no worse than Schmitz on each
  measure, not that Schmitz was adequate. Fatigue and buckling aren't
  modelled.
- **The Reynolds floor.** The lowest Reynolds number is 71 127 against the
  40 000 limit of the polar table, with 25.5 mm of envelope slack at the tip,
  so this isn't close.
- **The monotonicity constraints** act on the control points. That is
  sufficient for a monotone spline but not necessary, so a real mould
  constraint might allow or forbid different shapes. The optima have equal
  tip chords and, from δ = 0.25 %, equal outer twists.
- **The parameterisation.** 5 + 5 cubic B-spline control points, chosen in
  `representation_study/` for the energy objective. A finer chord basis
  might find more material to remove at the same energy.
- **The root cylinder in the renders** is only drawn: it has no mass, load or
  drag. The modelled blade starts at r = 0.30 m.
- **One resource and one operating law.** See the robustness section above.
