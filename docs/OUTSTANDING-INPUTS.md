# Model assumptions and the inputs that would replace them

Every value the design rests on is either measured, derived, or an engineering
assumption with a stated basis. This file lists, section by section, the
assumptions that a specific external input would replace, what that input
is, where it goes, and what re-runs when it lands. No value here was invented
to fill a gap: where an input is absent, the mechanism around it is built and
tested, and the code raises rather than defaults.

The design purpose, the model and the settled results are in
`docs/DESIGN-BASIS.md`; this file is the companion list of replaceable
inputs. Section numbers are stable because other documents cite them.

---

## 1. Wind resource — Weibull `k` and `c`

**As used.** A GASP point extraction at the site at 50 m (omni-directional):
`A = 8.8 m/s`, `k = 1.87`, `V̄ = 7.8 m/s`; extrapolated to the 20 m hub height
by the Justus & Mikhail (1976) correlation with shear exponent 0.2079, giving
`k = 1.7092`, `c = 7.2738 m/s`, `V̄ = 6.4876 m/s`. Recorded with a log-law
cross-check in `verification/wind_resource/` and carried in
`config/site.yaml` (`wind_resource:`), with `latitude_deg` and
`longitude_deg`.

**Assumptions.** The extraction is a single point, not an area selection, so
`location.gwa_area` is unresolved and stays so; the source atlas is GASP. The
longitude is recorded as East (16°33′28.2″ E), consistent with the pinned site
imagery. Nothing downstream reads the area field or the longitude.

**What would replace it.** A site measurement campaign or an area-averaged
atlas extraction at hub height. Then `config/site.yaml` changes and the chain
from `verification/baseline/` onward re-runs in the order given in
`verification/README.md`. The sensitivity of the settled blades to the
resource is bounded by evaluation in
`verification/mass_optimisation/cross_evaluation.json` (four Weibull corners:
x_m within −0.006 … +0.011 % of the reference energy).

---

## 2. Design-variable bounds and the buildable-tip floor

**As used** (`config/rotor_design.yaml`, `parameterisation.bounds` and
`manufacturing:`):

| field | value | basis |
|---|---|---|
| `chord_min_m` | 0.045 m | the box bound; it sets the scaling of the design variables in which every gradient artefact was measured, and is therefore fixed |
| `manufacturing.min_chord_m` | 0.060 m | the buildable-tip floor, five linear rows `c_i ≥ 0.060` of the design problem: a 60 mm SG6043 section is about 6 mm thick, room for two 2 mm skins, a bond line and a web; below the Schmitz tip (67 mm), so the reference is unchanged |
| `chord_max_m` | 0.30 m | local solidity at the root cut-out (r = 0.3 m): 0.30 m gives σ = 0.48, the conventional edge of BEM validity (0.45 m would give 0.72, blades nearly touching); c/R = 0.15 is the top of the small-turbine commercial range; the Schmitz root is 276 mm, so the reference is feasible and unclipped |
| `twist_min_deg` | −2° | Schmitz tip twist 0.57° with margin |
| `twist_max_deg` | 35° | Schmitz root twist 23.07° with margin |
| `constraints.max_local_solidity` | 0.5 | `σ_i = B c_i / (2π r_i) ≤ 0.5` at every station; with the 0.30 m box it cannot bind (350 mm would be needed at the first station) and reports inactive |

**What would replace it.** A hub radius and root-attachment concept (a
different upper chord) or a laminate concept (a different minimum section).
Changing the box re-scales the design variables and re-runs every gradient
and optimisation artefact; changing the floor re-runs
`verification/mass_optimisation/` only. The energy-optimum companion blade
x_c has a 47.7 mm tip and is below the floor; it is kept as an evaluated
reference, not as a feasible design.

---

## 3. Buhl (2005), NREL/TP-500-36834 — provenance

**As used.** The Glauert/Buhl high-induction relation in Ning's γ-form, with
the blend at `a = 0.4`; `tests/test_corrections.py` asserts C⁰ and C¹
continuity against momentum theory at `a = 0.4` and the F-independent anchor
`C_T(1) = 2`, the three conditions that pin all three coefficients uniquely.

**Assumption.** The constants are transcribed correctly from the original
report. The tests would catch a transcription error but not constants
faithfully transcribed from a secondary source.

**What would replace it.** The primary report, filed under `docs/references/`,
with the equation number and the stated `a_c` recorded beside the PROVENANCE
block in `src/bem/corrections.py`. No re-run.

---

## 4. Ning (2014), *Wind Energy* 17(9), 1327–1345 — provenance

**As used.** The single-residual formulation and the momentum-region
classification in `src/bem/station.py`, derived from the residual's
structure and verified numerically over 247 stations
(`verification/phase_vi/`). The propeller-brake region on `φ ∈ (−π/4, 0)` is
deliberately not implemented: a station that would need it is reported, not
solved by an untested branch.

**What would replace it.** The primary article, filed under
`docs/references/`, with `momentum_region_bracket` reconciled against Ning's
region definitions. No re-run.

---

## 5. Root cut-out fraction

**As used.** `geometry.root_fraction: 0.15` in `config/rotor_design.yaml` — a
conventional root cut-out for a rotor of this size, and a modelling choice. It
sized the SG6043 cache's Reynolds bounds and sets the inboard end of the
parameterised span; the segment inside `r_hub = 0.30 m` carries no chord, no
load and no material in any result. The root cylinder in the rendered blades
is drawn only.

**What would replace it.** A hub and root-attachment geometry. Then the
SG6043 cache bounds and the reference blade are re-checked.

---

## 6. NREL Phase VI measured performance data

**As used.** None. Every Phase VI number in this repository is a prediction
cross-checked against other predictions (CCBlade 0.51 %, pyBEMT 1.83 % mean
|C_P| deviation; `docs/validation/bem-cross-validation.md`).

**What would replace it.** A numeric measured Cp–λ dataset from a primary
source. `src/validation/plot_bem_comparison.py` plots measurements as a
fourth series should one be obtained. No re-run of the design.

---

## 7. The AEP sanity band

**As used.** `config/rotor_design.yaml` `sanity:` 8–12 MWh/yr of aerodynamic
shaft energy, derived from the recorded resource uncertainty (±10 % on `c`
and `k` spans 8.21–12.01; the log-law roughness band on the height
extrapolation spans 9.73–11.10) and rounded outward; capacity factor
0.24–0.36 on the aerodynamic rating. The reference blade's 10.2477 MWh/yr
lies inside it. `test_the_sanity_band_is_still_narrow_enough_to_catch_a_bug`
asserts that the band still fails on unlimited power (14.6), sea-level
density (13.0), unnormalised bin masses (12.9) and any factor of two.

**Assumption.** The band and every AEP figure are aerodynamic shaft energy;
no drivetrain efficiency is applied anywhere in the chain. Stating energy
electrically would introduce an `η` and a config field
(`10.25 × 0.90 = 9.2 MWh/yr`, still inside the band).

---

## 8. Justus & Mikhail (1976) — provenance

**As used.** `src/objective/height_extrapolation.py` implements the
correlation with constants `0.37` and `0.0881` and the 10 m anchor height;
`tests/test_height_extrapolation.py` checks the identity at `z = z_ref`,
monotonicity in `z`, the sign of both height relations, and agreement with an
independent log-law calculation across the plausible roughness range.

**What would replace it.** The primary paper, filed under `docs/references/`,
confirming the constants and the form of the `k(z)` relation. This item
carries more weight than §3 or §4: the two numbers it sets multiply through
every AEP figure.

---

## 9. The machine — rotor-speed ceiling, rating, above-rated limiting

Three machine facts set the structure of the energy functional. Each is an
assumption of the model with a stated basis (`config/rotor_design.yaml`
`operating:` and its comments); a specified generator replaces the first two.

### Maximum rotor speed — 300 rpm

`operating.max_rotor_speed_rpm: 300` (62.8 m/s tip speed at R = 2 m; the
ceiling bites at V_c = 9.67 m/s). The operating law is
`λ(V) = min(6.5, Ω_max R / V)` (`src/objective/power.py::tsr_schedule`).
Basis: without a ceiling the rotor would run 341 rpm = 71.4 m/s tip speed at
rated, above the ~60–65 m/s noise-conscious upper limit generally applied to
residential-scale small turbines; the nearest commercial analogue in size and
rating (Skystream 3.7, 2.4 kW, 3.72 m rotor) runs to ~330 rpm ≈ 64 m/s;
off-the-shelf direct-drive 3 kW permanent-magnet generators cluster at
250–300 rpm rated. Lower ceilings (263 rpm / 55 m/s, 239 rpm / 50 m/s) would
amount to designing for an under-sped generator and were not adopted. The
sensitivity of the settled blades to the ceiling is by evaluation
(`mass_optimisation/cross_evaluation.json`: x_m loses 0.10 / 0.44 / 1.84 % of
energy relative to the reference at 60 / 55 / 50 m/s).

**What would replace it.** A generator datasheet rated rpm. Then the value is
replaced and every optimisation artefact re-runs. See also §10.

### Generator rating — 3 822.19 W

`operating.rated_power_w` is the reference blade's aerodynamic power at
11 m/s and λ = 6.5, at full precision
(`tests/test_baseline.py::test_configured_rating_is_the_baselines_aerodynamic_rated_power`
pins that basis). Above rated, `power_per_bin` holds power at this value and
the adjoint carries the capped bins as constants
(`adjoint.system.BEMSystem.J_capped`). Under the 300 rpm ceiling the
reference runs λ = 5.71 at 11 m/s, so its power on the schedule is a little
below the rating; the rating is a nameplate and does not move with the
schedule.

**What would replace it.** A nameplate rating. Then the number is replaced,
that test is retired, and every AEP artefact re-runs.

### Above-rated limiting — ideal power hold

The model holds `P = P_rated` above 11 m/s with no mechanism (stall, torque
control or furling) modelled. A capped bin contributes a constant to the
energy and nothing to its gradient, which is defensible for the energy axis.
It is **not** a model of the loads above rated: at (V > 11 m/s, Ω_max) the
BEM state produces P_aero > P_rated with nothing shedding the difference, and
at 20 m/s / 300 rpm the inboard stations sit at α ≈ 30° on the Viterna
extrapolation. Consequently every load constraint is applied at and below
rated (nine operating points; the rated point at 11 m/s, 300 rpm, λ = 5.71
carries KS weight 0.957), and the cut-out moment is reported as a labelled
post-check only (`verification/load_constraint/README.md`).

**What would replace it.** The actual machine's limiting mechanism. It would
add a trim equation to the adjoint structure and allow loads above rated to
be constrained; it is outside the project's stated scope.

---

## 10. Kestrel e400 rated rotor speed — citation to strengthen §9

The Kestrel e400 (3 kW, 4.0 m rotor diameter, manufactured in Gqeberha) is
the nearest regional comparable in size and rating. Its datasheet rated rotor
speed and, if given, tip speed or noise rating would be cited beside the
Skystream 3.7 figure in the config comment and the report. If it does not
support ~300 rpm, the basis is recorded and the value revisited; artefacts
re-run only if the config value changes.

---

## 11. The laminate and the tip clearance — `structure:`

**As used.** The committed optimisation uses none of them: the objective is
the shell material proxy `k_P ∫ c dr` per blade and every production
constraint row is relative to the reference blade, so the headline is a
percentage (−3.38 % shell material at the reference energy, −7.28 % on the
solid proxy; −5.27 % within 0.5 % of it) and stays one whatever the laminate.

**Resolved: the laminate** (`laminate_density_kg_m3 = 1920`,
`shell_thickness_m = 0.002`, `youngs_modulus_pa = 41.8e9`,
`allowable_stress_pa = 702e6`, `safety_factor = 3.5725`). One laminate, one
stated construction — a 2 mm skin of the E-LT-5500/EP-3 unidirectional
E-glass/epoxy laminate of Griffith & Ashwill (2011) SAND2011-3779 Table 19,
tied to the 60 mm buildable tip ("two 2 mm skins + bond + web"), under the
GL combined safety factor for a wet hand-laid, non-post-cured laminate.
Provenance, value by value: `docs/MATERIALS-STRUCTURAL-INPUTS.md`. These
are representative laminate values for that construction, not measurements
of a built blade. What they give (`verification/absolute_material/`):

| | x₀ | x_c | x_m |
|---|---|---|---|
| shell mass per blade | 1.678 kg | 1.782 kg | **1.622 kg** (−0.057 kg; −0.17 kg per rotor) |
| thin-shell root stress `KS M_ref / (k_Z c₀² t)` | 22.6 MPa | 19.1 MPa | 22.6 MPa — 11.5 % of the 196.5 MPa design allowable |
| tip deflection `δ_ref D / (E k_I t)` | 111.8 mm | 85.2 mm | 111.8 mm |

The absolute stress row (`1 − σ/(σ_allow/SF) ≥ 0`) is built, verified
(Tier 3 at x₀ and x_m) and **additional**: added to the committed set it
leaves x_m unchanged (slack 0.885); put *in place of* the relative stress
row it returns the "no stress row" ablation optimum (−4.53 %, a root at 1.70×
the reference's stress per unit load, still only 38 MPa), because at the
operating loads of the load set the allowable is never within a factor of
five of binding. The sizing loads for a small blade's root are the IEC
61400-2 parked, gust and fatigue cases this model does not compute, so the
relative rows stay the design rows and the absolute figures are checks with
the label "thin-shell": no spar cap, no web, no root insert. The moment cap
`KS ≤ KS₀` is a load cap on the hub, shaft and tower, kept beside the stress
rows, which are strength checks.

**Outstanding: `tip_clearance_m`.** Machine geometry — tower diameter, hub
overhang, cone, tilt — of the same class as §9's nameplate and rotor speed.
GL states the allowed *remaining* clearance as a percentage of the unloaded
one, and the unloaded one is the machine's; no number for it can be
defended from the literature, so it stays `Unresolved`. The absolute
deflection row (`1 − δ_tip/tip_clearance ≥ 0`) is built and passes Tier 3
under an injected test clearance; `ScaledProblem.absolute_rows_available()`
reports why it is not assembled. When the clearance lands: set the field,
`verification/absolute_material/` re-runs, and the row joins the checks. It
changes the optimum only if the clearance is below the thin-shell blade's
111.8 mm at the rated point — in which case the laminate, not the planform,
is what would be revisited.

A laminate schedule (spar cap, variable thickness) replaces `t`, `k_I` and
`k_Z` by a schedule; until one exists the thin-shell label stays.
