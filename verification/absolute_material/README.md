# The blades in kilograms, megapascals and millimetres

The design problem in `verification/mass_optimisation/` works entirely in
percentages of the Schmitz blade, so it doesn't need any material data. This
folder applies one assumed laminate (`config/rotor_design.yaml::structure`,
sources in `docs/OUTSTANDING-INPUTS.md` §11) to the three blades and reports
their mass, root stress and tip deflection in engineering units.

Everything here is evaluated from `baseline/x0.json`,
`load_constraint/result_eps0.json` (x_c) and
`mass_optimisation/result_delta0.json` (x_m). The exception is the two
optimisation runs in section 4, which test whether an absolute stress
constraint would change the optimum.

## Reproducing it

    python verification/absolute_material/run_absolute_material.py     # 71 s
    python verification/absolute_material/run_spanwise_stress.py       # ~5 s

This reads the three blades, `mass_optimisation/ablation.json` (for the
comparison case), `fd_step_size/sweep.json` (h\*) and the config, and writes
`result.json` and `absolute_material.png`. It reuses
`mass_optimisation/_common.py` and the Tier 3 and Taylor helpers from
`mass_optimisation/run_mass_checks.py`, so it runs after `mass_optimisation`.
Nothing else reads its output.

## The inputs

All recorded in `result.json::structural_inputs`:

| field | value | basis |
|---|---|---|
| `laminate_density_kg_m3` | 1920 | E-LT-5500/EP-3 unidirectional E-glass/epoxy, 54 % fibre volume, VARTM; Griffith & Ashwill (2011), SAND2011-3779, p. 38 |
| `shell_thickness_m` | 0.002 | a design choice: the 60 mm tip is a 6 mm section, "two 2 mm skins, a bond line and a web" |
| `youngs_modulus_pa` | 41.8e9 | the same laminate's longitudinal modulus, Table 19 |
| `allowable_stress_pa` | 702e6 | the same laminate's compressive strength (95/95), Table 19; the lower of tension (972 MPa) and compression |
| `safety_factor` | 3.5725 | Germanischer Lloyd loads × materials factor for a wet hand-laid laminate that isn't post-cured: `1.35 × (1.35 × 1.35 × 1.1 × 1.2 × 1.1)`, Tables 9 and 12 of the same report |
| design allowable | **196.5 MPa** | `702 / 3.5725` |
| `tip_clearance_m` | not set | depends on the tower and hub geometry; there's no defensible number for it |
| `k_P`, `k_A`, `k_I`, `k_Z` | 2.0485, 0.0685, 0.003293, 0.05164 | calculated exactly from `data/airfoils/sg6043.dat` for a constant-thickness skin (`src/objective/mass.py`) |

These describe one skin of one laminate at constant thickness around the
whole section, with spanwise fibres, no spar cap, no web mass and no root
insert. They are representative values for that construction, not
measurements of a built blade.

## 1. The three blades

`result.json::blades`, `absolute_material.png`

| | x₀ (Schmitz) | x_c (energy optimum) | **x_m (minimum material)** |
|---|---|---|---|
| shell mass per blade | 1.678 kg | 1.782 kg (+6.17 %) | **1.622 kg (−3.38 %)** |
| shell mass per rotor (3 blades) | 5.035 kg | 5.346 kg | **4.865 kg** |
| solid-section mass per blade (reported only) | 4.290 kg | 4.899 kg | **3.977 kg (−7.28 %)** |
| thin-shell root stress `KS M_ref / (k_Z c0² t)` | 22.6 MPa | 19.1 MPa | **22.6 MPa** |
| as a share of the 196.5 MPa allowable | 11.5 % | 9.7 % | **11.5 %** |
| tip deflection `delta_ref D / (E k_I t)` | 111.8 mm (5.6 % of R) | 85.2 mm | **111.8 mm** |
| against a tip clearance | – | – | – (none specified) |

The saving is 1.678 − 1.622 = 0.057 kg per blade, or 0.17 kg per rotor. That
is the same −3.38 % as the main result, because the mass is the shell area
times a constant. x_m's stress and deflection match x₀'s exactly because both
relative constraints are active at x_m, which is what "no worse than the
reference" looks like in engineering units.

x₀'s deflection at the rated point alone (without the KS aggregate) is
111.79 mm (`references.rated_tip_deflection_x0_mm`); the KS value over the
load set is 111.84 mm.

## 2. Gradient checks on the absolute constraints

`result.json::verification`. The absolute stress constraint is
`1 − σ(u) / 196.5 MPa ≥ 0`. The absolute deflection constraint,
`1 − δ_tip(u) / clearance ≥ 0`, is checked with a **150 mm clearance
injected purely as a test value**. Both are checked at x₀ and x_m against the
3 ε_j threshold at the committed h\*:

| constraint | point | worst `|adj − fd| / ε_j` | Taylor ratios per decade |
|---|---|---|---|
| stress (absolute) | x₀ | 0.728 (`chord_1`) | 114.6, 102.4 |
| stress (absolute) | x_m | 1.314 (`twist_1`) | 100.9, 102.1 |
| deflection (absolute) | x₀ | 0.839 (`twist_1`) | 144.8, 94.2 |
| deflection (absolute) | x_m | 0.249 (`chord_0`) | 105.0, 100.3 |

All pass, so the round-off floor wasn't measured (it's only measured when a
check fails). Tiers 1 and 2 are inherited from the relative constraints: each
absolute constraint is the same cached adjoint result times a constant, which
`tests/test_mass_problem.py` checks to 1e-12.

## 3. Adding the absolute stress constraint

`runs.added_warm_from_x_m`. Here the absolute stress constraint is added to
the seven committed ones, starting from x_m. SLSQP stops after one iteration
at the same point (no change in u, shell −3.3788 %, same active set). The
new constraint's slack at x_m is 0.885, since the root carries only 11.5 % of
the allowable. **Adding it doesn't change the optimum.**

## 4. Replacing the relative stress constraint with the absolute one

`runs.swapped_*`. Here the relative stress constraint is replaced by the
absolute one. The problem is solved from x₀ (20 iterations) and from x_m
(14), and the two agree to 6.7e-8 in u (tolerance 0.0166).

| | x_m | **absolute stress instead of relative** | "no stress" ablation |
|---|---|---|---|
| shell material vs x₀ | −3.379 % | **−4.527 %** | −4.527 % |
| shell mass per blade | 1.622 kg | **1.602 kg** | – |
| root stress (thin shell) | 22.6 MPa (ratio 1.000) | **38.4 MPa (ratio 1.702)** | ratio 1.702 |
| tip deflection | 111.8 mm | **111.8 mm** (active) | ratio 1.000 |
| chord control points [mm] | 271 / 150 / 135 / 60 / 60 | **208 / 179 / 127 / 61 / 60** | – |
| active constraints | AEP, stress, deflection, minimum chord ×2, monotonicity | AEP, deflection, minimum chord at `chord_4`, twist monotonicity | – |
| distance from x_m in u (max norm) | 0 | **0.247** | – |

The result **is exactly the "no stress" optimum** from
`mass_optimisation/ablation.json` (distance 0.0 in u, the same −4.527 %).
Even after the root thins from 271 to 208 mm, the operating stress is only
38 MPa against a 196.5 MPa allowable, so the absolute constraint never gets
within a factor of five of binding and constrains nothing. The optimiser
takes an extra 1.15 points of material out of the root, moving 0.247 in u
(15 times the agreement tolerance), and ends at a root carrying 1.7 times the
reference's stress per unit load.

**Why the relative constraints stay the design constraints.** The load set
covers normal operation only. The loads that size a small blade's root are
the IEC 61400-2 cases, such as parked in the 50-year gust and fatigue, which
this model doesn't compute. So a 196.5 MPa allowable against a 22.6 MPa
operating stress isn't evidence that the root has a factor of nine in hand.
It shows that normal operation isn't the sizing case. The relative
constraint holds whatever the sizing load is: the new root carries no more
stress per unit load than the reference. The absolute constraint is kept as a
verified check, available by name and inactive at every blade, rather than as
a replacement. The moment cap `KS ≤ KS0` is kept for a different reason: it
limits the load on the hub, shaft and tower, and isn't a strength check at
all.

## 5. Stress along the span

`spanwise_stress.json`. The stress constraint only looks at the root section.
`run_spanwise_stress.py` evaluates the same thin-shell stress,
`M(r) / (k_Z t c(r)²)`, at every BEM station at the rated point, where `M(r)`
is the moment about each station from the load outboard of it (as in Table 9.1
of the report):

| | root (about `r_hub`) | peak along the span | where | peak / root |
|---|---|---|---|---|
| x₀ | 22.6 MPa | **39.8 MPa** | r = 0.95 m (r/R 0.48) | 1.76 |
| x_c | 19.1 MPa | 33.1 MPa | r = 0.88 m (r/R 0.45) | 1.73 |
| x_m | 22.6 MPa | **40.2 MPa** | r = 0.81 m (r/R 0.41) | 1.78 |

With a constant-thickness skin, the root isn't the most stressed section of
any of the three blades. The peak is near mid-span, at about 1.75 times the
root value. x_m's peak is 1.0 % above x₀'s and 0.14 m further inboard; x_c's
is 16.9 % below. The peak is still only 20 % of the 196.5 MPa allowable, so
the conclusion above (normal operation isn't the sizing case) still holds.
But it does mean the stress constraint holds the root at the reference value,
not the blade's highest stress. A spanwise stress constraint, with one
adjoint right-hand side per station, would constrain the peak.

## What would make this wrong

- **A different laminate or thickness** rescales every kilogram, megapascal
  and millimetre here (mass and stress scale with `1/t`, deflection with
  `1/(E t)`), but none of the percentages in `mass_optimisation/`. The values
  used are recorded in `result.json::structural_inputs` and pinned in
  `tests/test_config.py`.
- **A spar cap.** The thin-shell constants put all the bending material in
  the skin. A blade with a spar cap would be stiffer and stronger than these
  figures say, by an amount only a laminate design can give. Until then,
  every absolute figure carries the label "thin-shell".
- **The load case.** Section 4 shows that at the operating loads the absolute
  stress constraint can't size a root. Quoting the 196.5 MPa margin as a
  design margin would be a mistake.
- **The bending axis.** `k_I` is taken about the chord-parallel axis through
  the contour's centroid. The 0.086° rotation of the principal axes is
  neglected, as is the section's twist relative to the rotor plane (as in the
  relative measures).
- **Any change to x₀, x_c or x_m** means this folder has to be re-run.
