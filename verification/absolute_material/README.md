# The absolute material figures: kilograms, megapascals, millimetres

**New 2026-09-20 (evening).** The committed mass problem
(`verification/mass_optimisation/`, the 60 mm floor, `delta = 0`) is
unchanged; this artefact puts the laminate recorded that evening in
`config/rotor_design.yaml::structure` (provenance:
`docs/MATERIALS-STRUCTURAL-INPUTS.md`) under the three committed blades and
reports what it turns the proxies into. Everything here is *evaluated* from
`baseline/x0.json`, `load_constraint/result_eps0.json` (`x_c`) and
`mass_optimisation/result_delta0.json` (`x_m`) except the two optimisation
runs of section 4, which exist to measure whether the absolute stress row
moves the optimum.

## What was run

    python verification/absolute_material/run_absolute_material.py     # 71 s

Reads the three blades above, `mass_optimisation/ablation.json` (the
comparison case), `fd_step_size/sweep.json` (`h*`) and the config. Writes
`result.json` and `absolute_material.png`. Shares `mass_optimisation/_common.py`
and the Tier 3 / Taylor helpers of `mass_optimisation/run_mass_checks.py`.
Runs after `mass_optimisation` in `verification/README.md`'s order.

## The inputs (all in `result.json::structural_inputs`)

| field | value | basis |
|---|---|---|
| `laminate_density_kg_m3` | 1920 | E-LT-5500/EP-3 unidirectional E-glass/epoxy, 54 % fibre volume, VARTM -- Griffith & Ashwill (2011) SAND2011-3779 p. 38 |
| `shell_thickness_m` | 0.002 | design decision: the 60 mm tip is a 6 mm section, "two 2 mm skins + bond + web" |
| `youngs_modulus_pa` | 41.8e9 | same laminate, `E_L`, Table 19 |
| `allowable_stress_pa` | 702e6 | same laminate, ultimate compressive strength `UCS_L` (95/95), Table 19; the smaller of tension (972) and compression |
| `safety_factor` | 3.5725 | GL combined loads x materials factor for a wet hand-laid, non-post-cured laminate: `1.35 x (1.35 x 1.35 x 1.1 x 1.2 x 1.1)`, Tables 9 and 12 of the same report |
| design allowable | **196.5 MPa** | `702 / 3.5725` |
| `tip_clearance_m` | **TODO** | machine geometry; no defensible number |
| `k_P`, `k_A`, `k_I`, `k_Z` | 2.0485, 0.0685, 0.003293, 0.05164 | exact from `data/airfoils/sg6043.dat` for a constant-thickness skin (`src/objective/mass.py`) |

The construction these describe: one skin of one laminate at constant
thickness round the whole section, fibres spanwise, no spar cap, no web
mass, no root insert. They are representative laminate values for that
construction, not measurements of a built blade.

## 1. The three blades (`result.json::blades`, `absolute_material.png`)

| | `x0` (Schmitz) | `x_c` (energy optimum) | **`x_m` (mass optimum)** |
|---|---|---|---|
| shell mass, per blade | 1.678 kg | 1.782 kg (+6.17 %) | **1.622 kg (-3.38 %)** |
| shell mass, per rotor (3 blades) | 5.035 kg | 5.346 kg | **4.865 kg** |
| solid-section mass, per blade (reported, never optimised) | 4.290 kg | 4.899 kg | **3.977 kg (-7.28 %)** |
| thin-shell root stress `KS M_ref / (k_Z c0^2 t)` | 22.6 MPa | 19.1 MPa | **22.6 MPa** |
| against the 196.5 MPa design allowable | 11.5 % | 9.7 % | **11.5 %** |
| tip deflection `delta_ref D / (E k_I t)` | 111.8 mm (5.6 % R) | 85.2 mm | **111.8 mm** |
| against a tip clearance | -- | -- | -- (none specified) |

The saving in kilograms is `1.678 - 1.622 = 0.057 kg` per blade, `0.17 kg`
per rotor: the same -3.38 % as the committed headline, because the mass is
the shell area times a constant. The stress and deflection of `x_m` equal
`x0`'s to the digit because both relative rows are active at `x_m` -- that
is what "no worse than the reference" means in engineering units.

The rated-point deflection of `x0` on its own (no KS aggregate) is
111.79 mm (`references.rated_tip_deflection_x0_mm`); the KS over the load
set is 111.84 mm.

## 2. The absolute rows' Tier 3 and Taylor (`result.json::verification`)

`stress_absolute` (`1 - sigma(u) / 196.5 MPa >= 0`) and, under an
**injected 150 mm clearance that is a test value and nothing else**,
`deflection_absolute` (`1 - delta_tip(u) / clearance >= 0`), at `x0` and
`x_m`, acceptance `3 eps_j` at the committed `h*`:

| row | point | worst `abs(adj - fd) / eps_j` | Taylor ratios per decade |
|---|---|---|---|
| stress_absolute | `x0` | 0.728 (`chord_1`) | 114.6, 102.4 |
| stress_absolute | `x_m` | 1.314 (`twist_1`) | 100.9, 102.1 |
| deflection_absolute | `x0` | 0.839 (`twist_1`) | 144.8, 94.2 |
| deflection_absolute | `x_m` | 0.249 (`chord_0`) | 105.0, 100.3 |

All pass; the round-off floor was not measured (the Phase 4 rule: only on
failure). Tiers 1 and 2 are the relative rows' -- each absolute row is the
same cached adjoint result times a constant, which
`tests/test_mass_problem.py` asserts to `1e-12`.

## 3. Adding the absolute stress row (`runs.added_warm_from_x_m`)

Rows: the committed seven plus `stress_absolute`, warm from `x_m`. SLSQP
exits at iteration 1 at the same point (`u_inf` shift `0.0`, shell
`-3.3788 %`, the same active set). The absolute row's slack at `x_m` is
0.885 -- the root carries 11.5 % of its allowable. **The committed optimum
is unchanged by the absolute row.**

## 4. Replacing the relative stress row by the absolute one (`runs.swapped_*`)

Rows: the committed seven with `stress` **replaced** by `stress_absolute`,
cold from `x0` (20 iterations) and warm from `x_m` (14); the two agree to
`6.7e-8` in `u` (criterion `0.0166`).

| | committed `x_m` | **stress -> stress_absolute** | `no stress` ablation |
|---|---|---|---|
| shell material vs `x0` | -3.379 % | **-4.527 %** | -4.527 % |
| shell mass per blade | 1.622 kg | **1.602 kg** | -- |
| root stress (thin shell) | 22.6 MPa (ratio 1.000) | **38.4 MPa (ratio 1.702)** | ratio 1.702 |
| tip deflection | 111.8 mm | **111.8 mm** (row active) | ratio 1.000 |
| chord control points [mm] | 271 / 150 / 135 / 60 / 60 | **208 / 179 / 127 / 61 / 60** | -- |
| active rows | AEP, stress, deflection, min chord x2, monotone | AEP, deflection, min chord `chord_4`, monotone twist | -- |
| `u_inf` from `x_m` | 0 | **0.247** | -- |

The swapped optimum **is the "no stress" ablation optimum of
`mass_optimisation/ablation.json`** (`u_inf` distance `0.0`, the same
-4.527 %): with the allowable at 196.5 MPa and the operating stress at 38 MPa
even after the root has thinned from 271 to 208 mm, the absolute row never
comes within a factor of five of binding, so it constrains nothing and the
optimiser takes the extra 1.15 points of material out of the root. That is
the measured answer to "does the optimum move when the relative row is
replaced by the absolute one": **it moves by 0.247 in `u` (15 x the
agreement criterion) and saves 1.15 points more, to exactly the design the
plan's ablation had already labelled "not a blade"** -- a root carrying
1.7 x the reference's stress per unit load.

**Why the relative rows stay the design rows.** The load set `L` is the
operating set; the loads that size a small blade's root are IEC 61400-2's
parked-at-`V_e50`, gust and fatigue cases, which this model does not
compute. A 196.5 MPa allowable against a 22.6 MPa operating stress is
therefore not evidence that the root has a factor of nine in hand -- it is
evidence that the operating case is not the sizing case. The relative row
is the statement that holds whatever the sizing load is: the optimum's root
carries no more stress per unit of it than the reference's. The absolute
row is kept as a *check* -- available by name, verified, slack at every
committed blade, and the number the report can quote -- not as a
replacement. The moment cap `KS <= KS0` is kept for the same reason in the
other direction: it caps the hub, shaft and tower load, and is not a
strength check at all.

## What would make this wrong

- **A different laminate or thickness** rescales every kilogram, MPa and
  mm here (mass and stress `~ 1/t`, deflection `~ 1/(E t)`) but not one
  percentage in `mass_optimisation/`; `result.json::structural_inputs`
  records what was used, and the fields are pinned in
  `tests/test_config.py`.
- **A spar cap.** The thin-shell `k_I`, `k_Z` put all the bending material
  on the skin; a blade with a cap is stiffer and stronger than these numbers
  say, by an amount only a laminate schedule can give. Until one exists the
  label "thin-shell" stays on every absolute figure.
- **The load case.** Section 4 is the demonstration: at the operating
  loads the absolute stress row cannot size a root. Quoting the 196.5 MPa
  margin as a design margin would be the error this artefact exists to
  forestall.
- **The section axis.** `k_I` is about the chord-parallel axis through the
  contour's centroid; the principal-axis rotation neglected is 0.086 deg
  and the twist of the section relative to the rotor plane is neglected
  exactly as in the relative proxies.
- **A change to `x0`, `x_c` or `x_m`** (any upstream re-run) re-runs this.

## Re-run order

After `mass_optimisation` (it reads `result_delta0.json` and
`ablation.json`); nothing reads this artefact's output.
