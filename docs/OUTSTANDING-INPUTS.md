# Model assumptions and the inputs that would replace them

Every value the design depends on is either measured, derived, or an
engineering assumption with a stated basis. This document goes through the
assumptions one at a time. For each, it gives the value used, the external
input that would replace it, where that input would go, and what would need
re-running. Where an input is missing, I built and tested the mechanism
around it, and the code raises an error instead of falling back to a default.

What the project designs, and what it found, is in `docs/DESIGN-BASIS.md`.
The section numbers here are kept stable because other documents refer to
them.

---

## 1. Wind resource: Weibull `k` and `c`

**Used.** A GASP extraction at the site at 50 m, for all directions:
`A = 8.8 m/s`, `k = 1.87`, `V̄ = 7.8 m/s`. This is extrapolated to the 20 m
hub height with the Justus & Mikhail (1976) correlation (shear exponent
0.2079), giving `k = 1.7092`, `c = 7.2738 m/s` and `V̄ = 6.4876 m/s`. The
extrapolation, with a log-law cross-check, is in `verification/wind_resource/`,
and the values are in `config/site.yaml` (`wind_resource:`).

**Assumptions.** The extraction is a single point rather than an area
average, so the `location.gwa_area` field is left unresolved. The longitude
is East (16°33′28.2″ E), which matches the site imagery. Nothing downstream
reads either field.

**What would replace it.** A measurement campaign at the site, or an
area-averaged atlas value at hub height. `config/site.yaml` would change and
everything from `verification/baseline/` onward would be re-run, in the order
given in `verification/README.md`. How much the final blades depend on the
resource is already checked by evaluation in
`verification/mass_optimisation/cross_evaluation.json`: over four Weibull
corners, x_m stays within −0.006 to +0.011 % of the reference energy.

---

## 2. Design-variable bounds and the minimum buildable chord

**Used** (`config/rotor_design.yaml`, `parameterisation.bounds` and
`manufacturing:`):

| field | value | basis |
|---|---|---|
| `chord_min_m` | 0.045 m | the lower bound of the box. It also sets the scaling of the design variables that every gradient check was done in, so it stays fixed |
| `manufacturing.min_chord_m` | 0.060 m | the minimum buildable chord, added as five constraints `c_i ≥ 0.060`. A 60 mm SG6043 section is about 6 mm thick, enough for two 2 mm skins, a bond line and a web. It is below the Schmitz tip chord (67 mm), so the reference blade is unaffected |
| `chord_max_m` | 0.30 m | limits the local solidity at the root cut-out (r = 0.3 m) to 0.48, about the edge of where BEM is valid. 0.45 m would give 0.72, with the blades nearly touching. c/R = 0.15 is at the top of the commercial small-turbine range. The Schmitz root chord is 276 mm, so the reference fits |
| `twist_min_deg` | −2° | the Schmitz tip twist (0.57°) with some margin |
| `twist_max_deg` | 35° | the Schmitz root twist (23.07°) with some margin |
| `constraints.max_local_solidity` | 0.5 | `σ_i = B c_i / (2π r_i) ≤ 0.5` at every station. With the 0.30 m bound it can't become active (the first station would need 350 mm) |

**What would replace it.** A hub and root-attachment design (a different
maximum chord) or a laminate design (a different minimum section). Changing
the box changes the variable scaling, so every gradient and optimisation
result would need re-running. Changing the minimum chord only affects
`verification/mass_optimisation/`. The moment-capped energy optimum x_c has a
47.7 mm tip, below the floor, so it is kept as a comparison blade rather than
a buildable design.

---

## 3. Buhl (2005), NREL/TP-500-36834: source check

**Used.** The Glauert/Buhl high-induction correction in Ning's γ-form, blended
at `a = 0.4`. `tests/test_corrections.py` checks that it joins momentum
theory at `a = 0.4` with matching value and slope (C⁰ and C¹), and that
`C_T(1) = 2`. Those three conditions fix all three coefficients uniquely.

**Assumption.** The constants were transcribed correctly from the original
report. The tests would catch a transcription mistake, but not constants
copied faithfully from a secondary source.

**What would replace it.** The primary report, with the equation number and
the stated `a_c` noted next to the provenance comment in
`src/bem/corrections.py`. Nothing would need re-running.

---

## 4. Ning (2014), *Wind Energy* 17(9), 1327–1345: source check

**Used.** The single-residual formulation and the momentum-region
classification in `src/bem/station.py`. These were derived from the structure
of the residual and verified numerically over 247 stations
(`verification/phase_vi/`). The propeller-brake region, `φ ∈ (−π/4, 0)`, is
not implemented. A station that would need it is reported rather than solved
by an untested branch.

**What would replace it.** The primary article, to check
`momentum_region_bracket` against Ning's own region definitions. Nothing
would need re-running.

---

## 5. Root cut-out

**Used.** `geometry.root_fraction: 0.15` in `config/rotor_design.yaml`, a
conventional cut-out for a rotor this size and a modelling choice. It set the
Reynolds range of the SG6043 polar table and the inboard end of the designed
span. The part inside `r_hub = 0.30 m` carries no chord, load or material in
any result; the root cylinder in the rendered blades is only drawn.

**What would replace it.** A real hub and root-attachment geometry, after
which the polar table's range and the reference blade would be re-checked.

---

## 6. NREL Phase VI measured performance

**Used.** None. Every Phase VI number in the repository is a prediction
compared with other predictions (mean |C_P| difference: 0.51 % against
CCBlade, 1.83 % against pyBEMT; `docs/validation/bem-cross-validation.md`).

**What would replace it.** A measured Cp–λ dataset from a primary source.
`src/validation/plot_bem_comparison.py` can plot measurements as a fourth
series. The design wouldn't need re-running.

---

## 7. The AEP sanity band

**Used.** `config/rotor_design.yaml` `sanity:` 8–12 MWh/yr of aerodynamic
shaft energy. This comes from the resource uncertainty (±10 % on `c` and `k`
gives 8.21–12.01; the log-law roughness range in the height extrapolation
gives 9.73–11.10), rounded outward. It corresponds to a capacity factor of
0.24–0.36 on the aerodynamic rating. The reference blade's 10.2477 MWh/yr is
inside it. A test,
`test_the_sanity_band_is_still_narrow_enough_to_catch_a_bug`, checks that the
band still catches unlimited power (14.6), sea-level density (13.0),
unnormalised bin weights (12.9) and any factor-of-two error.

**Assumption.** The band, and every AEP figure, is aerodynamic shaft energy.
No drivetrain efficiency is applied anywhere. Converting to electrical energy
would need an efficiency and a config field (for example
`10.25 × 0.90 = 9.2 MWh/yr`, still inside the band).

---

## 8. Justus & Mikhail (1976): source check

**Used.** `src/objective/height_extrapolation.py` implements the correlation
with the constants `0.37` and `0.0881` and the 10 m anchor height.
`tests/test_height_extrapolation.py` checks the identity at `z = z_ref`,
monotonicity in `z`, the sign of both height relations, and agreement with an
independent log-law calculation over the plausible roughness range.

**What would replace it.** The primary paper, to confirm the constants and
the form of the `k(z)` relation. This one matters more than §3 or §4, because
the two numbers it produces feed into every AEP figure.

---

## 9. The machine: rotor-speed ceiling, rating and behaviour above rated

Three facts about the machine shape the energy calculation. Each is a
modelling assumption with the basis given here and in the comments in
`config/rotor_design.yaml` (`operating:`). A real generator datasheet would
replace the first two.

### Maximum rotor speed: 300 rpm

`operating.max_rotor_speed_rpm: 300`, which is a 62.8 m/s tip speed at
R = 2 m. The ceiling takes over at 9.67 m/s, and the operating law is
`λ(V) = min(6.5, Ω_max R / V)` (`src/objective/power.py::tsr_schedule`).

The basis:
- Without a ceiling the rotor would reach 341 rpm (71.4 m/s tip speed) at
  rated wind, above the roughly 60–65 m/s limit usually applied to keep
  residential turbines quiet.
- The closest commercial machine in size and rating (Skystream 3.7, 2.4 kW,
  3.72 m rotor) runs up to about 330 rpm, or 64 m/s.
- Off-the-shelf direct-drive 3 kW permanent-magnet generators are typically
  rated at 250–300 rpm.

Lower ceilings (263 rpm / 55 m/s and 239 rpm / 50 m/s) would mean designing
for an under-sized generator, so I didn't use them. The effect of the ceiling
on the final blades is checked by evaluation
(`mass_optimisation/cross_evaluation.json`): x_m loses 0.10 / 0.44 / 1.84 %
of energy relative to the reference at 60 / 55 / 50 m/s.

**What would replace it.** The rated rpm from a generator datasheet, after
which every optimisation would be re-run. See also §10.

### Generator rating: 3 822.19 W

`operating.rated_power_w` is the reference blade's aerodynamic power at
11 m/s and λ = 6.5, at full precision. A test,
`tests/test_baseline.py::test_configured_rating_is_the_baselines_aerodynamic_rated_power`,
pins that basis. Above rated, `power_per_bin` holds the power at this value,
and the adjoint treats those bins as constants
(`adjoint.system.BEMSystem.J_capped`). Under the 300 rpm ceiling the
reference runs at λ = 5.71 at 11 m/s, so its power at rated wind is slightly
below the rating. That's expected, because a generator's rating doesn't
change with the operating schedule.

**What would replace it.** A nameplate rating. The number would be replaced,
the test above retired, and every AEP result re-run.

### Above rated: an ideal power hold

The model holds `P = P_rated` above 11 m/s without modelling how (stall,
torque control or furling). For energy this is reasonable: an above-rated bin
adds a constant to the energy and nothing to its gradient. It is **not** a
model of the loads above rated. At wind speeds above 11 m/s at 300 rpm, the
BEM solution gives more aerodynamic power than the rating, with nothing to
shed the difference. At 20 m/s the inboard stations sit at α ≈ 30°, in the
Viterna extrapolation. For that reason every load constraint is applied at
and below rated only. That's nine operating points, and the rated point
(11 m/s, 300 rpm, λ = 5.71) carries a KS weight of 0.957. The cut-out moment
is reported as a separate check only (`verification/load_constraint/README.md`).

**What would replace it.** The real machine's power-limiting mechanism. That
would add a trim equation to the adjoint and allow loads above rated to be
constrained. It is outside the scope of this project.

---

## 10. Kestrel e400 rotor speed: a second reference for §9

The Kestrel e400 (3 kW, 4.0 m rotor, made in Gqeberha) is the closest
regional machine in size and rating. Its datasheet rotor speed, and its tip
speed or noise rating if given, would be cited next to the Skystream 3.7
figure. If it doesn't support about 300 rpm, the value would be revisited;
results would only need re-running if the config value changed.

---

## 11. The laminate and the tip clearance (`structure:`)

The optimisation itself doesn't use any of these values. The objective is
the shell-material measure `k_P ∫ c dr` per blade, and every design
constraint is relative to the reference blade. So the main result is a
percentage (−3.38 % shell material at the reference energy, −7.28 % on the
solid measure, −5.27 % within 0.5 % of the reference energy), whatever the
laminate. The values below turn those percentages into kilograms,
megapascals and millimetres.

### The assumed construction

One skin of a single E-glass/epoxy laminate, of constant thickness `t` all the
way around the section, with the fibres running spanwise. There is no spar
cap, web, root insert or gelcoat. The web and bond line fit inside the 60 mm
tip's thickness ("two 2 mm skins, a bond line and a web") but aren't given
any mass or stiffness here.

Three section constants turn a chord into a mass, a section modulus and a
second moment of area. For a thin skin they can be calculated exactly from
the airfoil coordinates (`src/objective/mass.py::section_coefficients`,
`data/airfoils/sg6043.dat`):

| constant | definition | SG6043 | used for |
|---|---|---|---|
| `k_P` | perimeter per unit chord | 2.048505 | mass `m = ρ t k_P ∫c dr` |
| `k_A` | area per unit chord² | 0.068502 | the solid measure (reported only) |
| `k_I` | `∮ (y − ȳ)² ds` per unit chord³ | 0.003293 | `I = k_I c³ t`, deflection |
| `k_Z` | `k_I / max|y − ȳ|` per unit chord² | 0.051643 | `Z = k_Z c² t`, root stress |

Here `ȳ` is the centroid of the contour weighted by length (0.0389 c above
the chord line), and the extreme fibre is on the upper surface, 0.0638 c
above it. Bending is taken about a chord-parallel axis through `ȳ`
(flapwise). I neglected the 0.086° rotation of the principal axes of the
cambered section, and the section's twist relative to the rotor plane, as
the relative measures also do. The values are pinned in `tests/test_mass.py`
and checked against a contour subdivided 400 times.

**These are representative values for that construction, not measurements
of a built blade, and none of the absolute figures is a certification.**

### The values

| field | value | basis and source |
|---|---|---|
| `laminate_density_kg_m3` | 1920 kg/m³ | E-LT-5500/EP-3 unidirectional E-glass/epoxy, 54 % fibre volume, made by VARTM. Griffith & Ashwill (2011), *The Sandia 100-meter All-glass Baseline Wind Turbine Blade: SNL100-00*, SAND2011-3779, p. 38 |
| `shell_thickness_m` | 0.002 m | a design choice tied to the 60 mm minimum chord ("two 2 mm skins, a bond line and a web"). Constant along the span and around the section; about four plies of a heavy unidirectional fabric (ply thickness around 0.5 mm; Bir 2005, NREL/TP-500-38929, Fig. 5) |
| `youngs_modulus_pa` | 41.8 GPa | the same laminate's longitudinal modulus, Table 19 of the same report |
| `allowable_stress_pa` | 702 MPa | the same laminate's longitudinal compressive strength (a 95/95 value), Table 19. The tensile strength is 972 MPa; compression governs at a flapwise-loaded root, so the lower value is used |
| `safety_factor` | 3.5725 | Germanischer Lloyd (2010) partial safety factors for ultimate strength, as tabulated in Griffith & Ashwill (2011), Tables 9 and 12 (pp. 25–26): loads `γ_f = 1.35` (normal load case) times materials `γ_M = 1.35 × 1.35 (ageing) × 1.1 (temperature) × 1.2 (wet hand lay-up) × 1.1 (not post-cured) = 2.6463`. I chose the hand lay-up and non-post-cured factors because that matches the assumed construction; the source's own infused, post-cured laminate would give 2.977. The design allowable is `702 / 3.5725 = 196.5 MPa` |
| `tip_clearance_m` | not set | depends on the machine: tower diameter, hub overhang, cone and tilt. Germanischer Lloyd gives the allowed remaining clearance as a percentage of the unloaded clearance, and that belongs to the machine, so I couldn't justify a number from the literature. The deflection is reported in mm instead |

All four laminate values come from one row of one table, so that the density,
stiffness and strength belong to the same material. The same table's
triaxial laminate (`[±45]₂[0]₂`, `E_L = 27.7 GPa`, 1850 kg/m³) would be a more
usual skin, but no strength is tabulated for it, and mixing rows would amount
to inventing a laminate. This choice errs in two known directions. A skin
with off-axis plies would deflect `41.8 / 27.7 = 1.51` times more at the same
thickness. A wet hand lay-up at a lower fibre fraction would be lighter and
less stiff per unit thickness than the infused laminate in the table. I
haven't quantified either.

IEC 61400-2 (2013) is the certification standard for a machine like this. It
is cited for what it prescribes (its own load cases and partial factors), not
for any number used here.

### What the values give

    mass            m        = ρ t k_P ∫c dr                    per blade, kg
    root stress     σ(u)     = KS(u) M_ref / (k_Z c₀(u)² t)     Pa
    tip deflection  δ_tip(u) = δ_ref D(u) / (E k_I t)           m

`KS`, `M_ref`, `D` and `δ_ref` are the same quantities the relative
constraints use. The absolute constraints are written as
`1 − σ/(σ_allow/SF) ≥ 0` and `1 − δ_tip/tip_clearance ≥ 0`, divided through
so that SLSQP sees values of order one. Their gradients are the relative
constraints' adjoint gradients times a constant (`src/gradients/problem.py`).
The results are in `verification/absolute_material/`:

| | x₀ (Schmitz) | x_c (energy optimum) | x_m (minimum material) |
|---|---|---|---|
| shell mass per blade | 1.678 kg | 1.782 kg | **1.622 kg** (−3.38 %; −0.057 kg per blade, −0.17 kg per rotor) |
| shell mass per rotor | 5.035 kg | 5.346 kg | 4.865 kg |
| solid-section mass per blade (reported only) | 4.290 kg | 4.899 kg | 3.977 kg |
| thin-shell root stress | 22.6 MPa | 19.1 MPa | 22.6 MPa, 11.5 % of 196.5 MPa |
| tip deflection (rated point) | 111.8 mm (5.6 % of R) | 85.2 mm | 111.8 mm |

### Why the relative constraints stay the design constraints

The absolute stress constraint is built and verified (Tier 3 at x₀ and x_m).
Added to the design problem, it leaves x_m unchanged (slack 0.885). Used
instead of the relative stress constraint, however, it moves the optimum to
−4.53 % shell material, with a 208 mm root carrying 1.70 times the
reference's stress per unit load (still only 38 MPa). That happens because at
the operating loads the 196.5 MPa allowable is never within a factor of five
of binding.

The loads that size a small blade's root are the ones IEC 61400-2
prescribes, such as the parked case in a 50-year gust and fatigue, and this
model doesn't compute them. A nine-fold margin at the operating point shows
that the operating point isn't the sizing case; it isn't evidence of a real
margin. The relative constraint holds whatever the sizing load turns out to
be, so it stays as the design constraint, and the absolute figures are
reported as thin-shell checks. The moment cap `KS ≤ KS₀` is kept alongside
as a limit on the load passed to the hub, shaft and tower, not as a strength
check.

### The tip clearance

The absolute deflection constraint is built and passes Tier 3 with a test
clearance. `ScaledProblem.absolute_rows_available()` reports why it isn't
included. Once a clearance is known, setting the field re-runs
`verification/absolute_material/` and adds the constraint to the checks. It
would only change the optimum if the clearance were less than the thin-shell
blade's 111.8 mm deflection at the rated point, and in that case it's the
laminate, not the planform, that should change.

A real laminate schedule (a spar cap, or thickness varying along the span)
would replace the constants `t`, `k_I` and `k_Z` with functions of radius.
Until then, every absolute figure is labelled "thin-shell".

### Where each value lives

- `config/rotor_design.yaml::structure`: the values, with their basis in a comment.
- `src/config/schema.py::DesignRotorConfig`: the fields; the design allowable is derived.
- `src/objective/mass.py`: `k_I`, `k_Z`, `mass_kg`.
- `src/gradients/problem.py`: `stress_absolute_pa`, `deflection_absolute_m`, the two absolute constraints, `absolute_rows_available`.
- `tests/test_config.py`, `tests/test_mass.py`, `tests/test_mass_problem.py`: pinned values, Tier 3 at x₀ and x_m.
- `verification/absolute_material/`: the figures above and the two optimisation runs.
