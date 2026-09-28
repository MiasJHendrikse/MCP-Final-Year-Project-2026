# Design basis

This document sets out what the project designs, the assumptions behind it,
the model it uses and what it found. Every number is read from a file in the
repository, and the file is named beside it, so any figure can be traced back
to the run that produced it.

## 1. Purpose and research questions

The aim is to find the lightest buildable blade for a small wind turbine with
a 2.0 m radius, three blades and fixed pitch, at a low-wind, high-altitude
site in the Khomas Hochland, Namibia. The blade has to deliver the same annual
energy as the classical Schmitz design for the same rotor, without increasing
the root bending moment, the root stress or the tip deflection. It also has to
respect a minimum buildable chord and keep the chord and twist decreasing
from root to tip.

Annual energy is therefore a constraint rather than the objective. Section 2
explains why.

The two research questions are:

1. **Design.** What is the lightest buildable blade that matches the Schmitz
   blade's annual energy at this site, with no increase in root moment, root
   stress or tip deflection? And how much material can be saved for each
   per cent of energy given up?
2. **Method.** Are discrete-adjoint gradients of a BEM model accurate enough
   to drive this constrained design? How do they compare with finite
   differences for accuracy, cost and robustness?

## 2. Why the objective is material and not energy

The project started by maximising annual energy production (AEP). The
optimiser converged cleanly, and two independent gradient methods (finite
differences and the adjoint) found the same optimum from several starting
points. But the optimum was only about 0.15 % better than the Schmitz blade,
where the literature led me to expect 2–6 %.

The reason turned out to lie in the problem, not the optimiser:

- **The machine runs at a fixed tip-speed ratio.** Below the rotor-speed
  ceiling, every wind speed is the same operating point in λ. Only the
  Reynolds number differs between wind-speed bins, so AEP is essentially one
  number, C_P at λ = 6.5, weighted over the bins.
- **The Schmitz blade is the analytic maximiser of exactly that number.** It
  is built from the same SG6043 polar at the same λ. Only second-order effects
  are left to gain: Reynolds dependence, tip loss and the spline fit.
- **The 2–6 % figure compares against weaker baselines.** Under the same
  objective, a linear-taper, linear-twist blade scores 4.0–4.7 % below
  Schmitz. A Schmitz blade built at a datasheet angle of attack instead of the
  real polar's best point scores 2.1 % below. The literature's gains are
  measured against blades like those.
- **Part of the first number was not a real gain.** An early version let the
  rated power rise with the design. With the generator rating fixed, as it
  must be for a real machine, the gain fell from 0.217 % to 0.121 %.

I also tested lower rotor-speed ceilings. They do make blade shape matter
more (a 55 m/s tip-speed limit gives a 2.3 % gain), but they amount to
designing for an under-sized generator, and part of the gain comes from the
Schmitz blade being scored deep in stall. Picking a ceiling to get a bigger
number would not have been honest.

The useful finding is that energy is flat near the optimum. The energy
optimum is only 0.15 % better than Schmitz, yet it needs 6 % more shell
material. That flatness is design freedom, and this project spends it on
material: hold the energy at the Schmitz value and make the blade lighter.

The measurements behind this section are in `verification/aep_gain_audit/`,
`verification/aep_optimisation_experiment/` and
`verification/load_constraint/`.

## 3. Site and wind resource

| quantity | value | basis | file |
|---|---|---|---|
| Site | 22°25′40.7″ S, 16°33′28.2″ E; a ridge west of Windhoek at about 1 800 m | pinned satellite location; the elevation is the assumed ridge height | `config/site.yaml` `location:` |
| Hub height | 20 m | a typical free-standing tower for this size of rotor | `config/site.yaml` |
| Resource at 50 m | Weibull k = 1.87, A = 8.8 m/s, V̄ = 7.81 m/s | GASP extraction at the site (a single point, all directions) | `verification/wind_resource/wind_resource_20m.json` → `reference_level` |
| Resource at 20 m | k = 1.7092, c = 7.2738 m/s, V̄ = 6.4876 m/s | Justus & Mikhail height extrapolation, shear exponent 0.2079, cross-checked with the log law | `config/site.yaml` `wind_resource:`; `wind_resource_20m.json` → `extrapolation`, `cross_check_log_law` |
| Temperature | 293.15 K (annual mean) | design value for the site, not the ISA temperature | `config/rotor_design.yaml` `atmosphere:` |
| Pressure | 81 489 Pa | ISA barometric relation at 1 800 m | same |
| Air density | ρ = 0.9684 kg/m³ | ρ = p/(R_s T), about 21 % below sea level | same |
| Kinematic viscosity | ν = 1.869 × 10⁻⁵ m²/s | μ = 1.81 × 10⁻⁵ Pa s at 20 °C, divided by ρ | same |

The design isn't re-optimised for different resource assumptions. Instead,
the final blades are re-evaluated at four Weibull corners (section 7.5).

## 4. Machine

| quantity | value | basis | file |
|---|---|---|---|
| Rotor radius | R = 2.0 m, fixed | keeps the blade's Reynolds numbers inside the measured SG6043 range; not a design variable | `config/rotor_design.yaml` `geometry:` |
| Blades | B = 3 | | same |
| Root cut-out | 0.15 R = 0.30 m | conventional for this size; a modelling choice | same |
| Design tip-speed ratio | λ_d = 6.5 | | `operating:` |
| Rated / cut-in / cut-out wind speed | 11 / 3 / 20 m/s | rated at about 1.7 times the mean wind speed (11 / 6.49), low for a low-wind site | same |
| Rated power | P_rated = 3 822.19 W, held constant above rated | the reference blade's aerodynamic power at 11 m/s and λ = 6.5; a generator's rating doesn't change with the blade | `operating.rated_power_w` |
| Rotor-speed ceiling | Ω_max = 300 rpm (62.8 m/s tip speed; the ceiling takes over at 9.67 m/s) | a tip-speed limit of about 60–65 m/s is usual for quiet residential turbines; comparable machines run at 250–330 rpm (e.g. Skystream 3.7 at about 330 rpm, 64 m/s) | `operating.max_rotor_speed_rpm` |
| Operating law | λ(V) = min(6.5, Ω_max R / V); P = min(P_aero, P_rated) | variable-speed tracking at a fixed tip-speed ratio, with a rotor-speed ceiling and an ideal power hold above rated | `src/objective/power.py::tsr_schedule` |
| Airfoil | SG6043 along the whole span | good low-Reynolds performance at this scale; S809 is used only for validating the solver | `data/airfoils/README.md`, `data/polars/sg6043/README.md` |

These machine values are modelling assumptions, each with the basis given.
Real generator data would replace the ceiling and the rating, after which
every optimisation result would be re-run in the order given in
`verification/README.md`. Above rated, the model assumes an ideal power hold
with no limiting mechanism, so loads are only constrained at and below rated
(section 6).

## 5. Model

- **Solver.** Blade-element momentum in Ning's single-residual form, with one
  unknown (the inflow angle φ) per station and operating point. It includes
  Prandtl tip and hub loss and the Glauert/Buhl high-induction correction,
  blended smoothly (C¹) at a = 0.4. The root is found by bracketing, which
  always converges. The blade is split into 25 strips (`src/bem/`;
  `docs/adjoint_derivation.md` §1–§5; `verification/phase_vi/`).
- **Polars.** An XFOIL polar table for SG6043 covering Re 40 k – 1 M at
  N_crit = 9, extended past stall with Viterna, and read through a C¹
  (α, Re) interpolant with analytic derivatives. XFOIL is never called during
  optimisation (`data/polars/sg6043/`, `src/polars/`,
  `verification/polar_interpolant/`, `results/sg6043_uiuc_validation/`,
  `results/ncrit_sensitivity/`).
- **Energy.** AEP is summed over 17 wind-speed bins of 1 m/s, with Weibull
  bin weights taken from the cumulative distribution and the operating law
  above. Bins above rated are constants. The energy function is C¹, with one
  curvature break where stations cross the Buhl threshold
  (`src/objective/`; `verification/smoothness_gate/README.md`).
- **Loads.** The normal load per blade is ½ρw²cC_n. From it come the root
  flapwise moment, the spanwise moment and an Euler–Bernoulli unit-load tip
  deflection with section stiffness proportional to c³. These are combined
  over the nine operating points at or below rated with a KS aggregate
  (ρ_KS = 100) (`src/objective/loads.py`; `docs/adjoint_derivation.md`
  §10–§11).
- **Material.** The objective is the shell area, k_P ∫c dr, where
  k_P = 2.048505 is the SG6043 perimeter per unit chord. A solid-section
  measure, k_A ∫c² dr with k_A = 0.068502, is reported alongside but never
  optimised (`src/objective/mass.py`).
- **Laminate.** The kilogram and stress figures assume one construction: a
  2 mm skin of unidirectional E-glass/epoxy (E-LT-5500/EP-3: ρ = 1920 kg/m³,
  E_L = 41.8 GPa, compressive strength 702 MPa; Griffith & Ashwill 2011,
  SAND2011-3779, Table 19). The skin has constant thickness around the whole
  section, with no spar cap, web or root insert. The design allowable is
  702 / 3.5725 = 196.5 MPa, using the Germanischer Lloyd safety factors for a
  wet hand-laid laminate that hasn't been post-cured. The thin-shell section
  constants k_I = 0.003293 (I = k_I c³t) and k_Z = 0.051643 (Z = k_Z c²t) are
  computed exactly from the SG6043 contour. These are representative values
  for that construction, not measurements of a built blade
  (`docs/OUTSTANDING-INPUTS.md` §11, `config/rotor_design.yaml::structure`).
- **Parameterisation.** Clamped cubic B-splines with 5 chord and 5 twist
  control points (10 design variables). The variables are scaled to [0, 1] by
  the box: chord 0.045–0.30 m, twist −2° to 35°. The number of control points
  was chosen in `verification/representation_study/`; at 5 + 5 the spline
  matches the analytic Schmitz blade to 0.515 mm RMS in chord and 0.106° in
  twist.
- **Reference blade x₀.** The Schmitz distribution with wake rotation at the
  SG6043 best-L/D point (Re 200 k, α = 5.36°, C_l = 1.282, L/D = 98.5) and
  λ = 6.5, fitted to the spline. The fit is within 1.52 mm (0.52 mm RMS) in
  chord and 0.265° (0.106° RMS) in twist, which costs 0.0154 % of AEP
  (`verification/baseline/x0.json`, `verification/spline_fit_error/`).
- **Gradients.** A discrete adjoint on the diagonal ∂R/∂φ system, solved for
  three right-hand sides (energy, root moment, tip deflection) that share one
  forward solve. The material gradient is analytic. The gradients are
  checked in four tiers (`src/adjoint/`, `src/gradients/problem.py`;
  section 7.4).

## 6. The design problem

    minimise    m(u) = k_P ∫ c dr / (k_P ∫ c dr)(x0)          shell material (exact gradient)
    subject to  AEP(u) / AEP(x0) − (1 − δ)   ≥ 0             energy floor (δ = 0 for the main result)
                KS0 − KS(u)                  ≥ 0             root-moment cap
                KS0 / c00² − KS(u) / c0(u)²  ≥ 0             root-stress proxy (Z ∝ c² t)
                D0 − D(u)                    ≥ 0             tip-deflection proxy (I ∝ c³ t)
                c_i − c_{i+1} ≥ 0 (4), θ_i − θ_{i+1} ≥ 0 (4) chord and twist decrease toward the tip
                c_i − 0.060 ≥ 0 (5)                          60 mm minimum buildable chord
                polar envelope (50 rows), local solidity ≤ 0.5 (25 rows), u ∈ [0, 1]¹⁰

It is solved with SLSQP (ftol 1e-10, maxiter 400). The problem is assembled
in `src/gradients/problem.py::constraints_for_mass_problem`, and
`verification/mass_optimisation/README.md` describes it in full.

A few notes on the constraints:

- **The three structural constraints are relative to x₀.** They share one
  thin-shell assumption: mass ∝ c, section modulus ∝ c²t, stiffness ∝ c³t.
  Each one says the new blade is no worse than the Schmitz blade, so no
  material data enters them. The laminate in section 5 turns the same
  quantities into kilograms, megapascals and millimetres
  (`verification/absolute_material/`), and adds an absolute stress
  constraint, σ ≤ 196.5 MPa. That constraint never comes within a factor of
  five of binding at the operating loads. If it replaced the relative stress
  constraint, stress would effectively stop being constrained: the optimum
  moves to −4.53 % material with a root at 1.70 times the reference stress
  per unit load. The loads that actually size a small blade's root (parked,
  gust and fatigue cases) are outside this model. So the relative
  constraints do the design work and the absolute figures serve as checks.
  An absolute deflection constraint would need a tip clearance, which
  depends on the tower and hub geometry and isn't known.
- **The 60 mm minimum chord is a manufacturing judgement.** A 60 mm SG6043
  section is about 6 mm thick, enough for two 2 mm skins, a bond line and a
  web. The 45 mm lower bound of the box is the thickness of the whole section
  and couldn't actually be built. The floor sits below the reference tip
  chord (67 mm), so x₀ and everything normalised by it are unchanged. It is
  written as a constraint rather than a new box bound because the box also
  sets the scaling of the design variables, and changing it would invalidate
  the gradient checks.
- **Loads are constrained at and below rated only.** There are nine operating
  points. The rated point (11 m/s, 300 rpm, λ = 5.71) carries a KS weight of
  0.957. The cut-out condition is only reported, because the model holds
  rated power above 11 m/s with no real limiting mechanism
  (`verification/load_constraint/README.md`).
- **The monotonicity constraints act on the control points.** A B-spline is
  monotone whenever its control points are (the variation-diminishing
  property), so this is sufficient but a little stricter than necessary.

## 7. Results

### 7.1 The reference blade x₀

AEP 10.2477 MWh/yr. At the rated point (11 m/s, 300 rpm, λ = 5.712):
C_P = 0.4516, C_T = 0.7029, root flapwise moment 177.38 N m per blade, rotor
thrust 517.5 N. Chord control points 275.8 / 181.6 / 97.1 / 78.6 / 67.0 mm;
twist 24.53 / 9.73 / 4.04 / 1.00 / 0.59°; shell area 0.43705 m² per blade
(`verification/baseline/baseline_reference.json`, `x0.json`;
`verification/mass_optimisation/reference_blades.json`).

### 7.2 Energy is nearly insensitive to blade shape

| blade | AEP vs x₀ | shell vs x₀ | solid vs x₀ | file |
|---|---|---|---|---|
| x\*, unconstrained energy optimum | +0.147 % | | | `verification/adjoint_optimisation/result.json` (the adjoint and FD runs find the same point; eight random starts converge to it within 0.0166 in u, `fd_optimisation_multistart/`) |
| x_c, energy optimum with the root moment capped | +0.1465 % | +6.171 % | +14.196 % | `verification/load_constraint/result_eps0.json`; `mass_optimisation/reference_blades.json` |

Tightening the cap on the root moment by 2 / 5 / 10 % costs 0.0206 / 0.1282 /
0.5785 % of AEP. The shadow prices range from 0.0017 to 0.129 % of AEP per
1 % of moment (`load_constraint/pareto.json`, `result_eps*.json`).

### 7.3 The minimum-material blade x_m

| | δ = 0 (at the reference energy) | δ = 0.5 % | file |
|---|---|---|---|
| shell material vs x₀ | **−3.379 %** | **−5.271 %** | `result_delta0.json`, `pareto.json` (`steps[2]`) |
| solid measure vs x₀ | −7.278 % | −9.666 % | same |
| AEP vs x₀ | −1.4 × 10⁻¹⁰ % (floor active) | −0.500 % | same |
| KS / KS₀ (moment cap) | 0.9685, not active | | same |
| stress and deflection ratios | 1.000, 1.000 (both active) | | same |
| shell mass per blade, with the laminate above | **1.622 kg** (x₀ 1.678 kg; −0.057 kg per blade, −0.17 kg per rotor) | | `verification/absolute_material/result.json` |
| thin-shell root stress; tip deflection | 22.6 MPa (11.5 % of 196.5 MPa); 111.8 mm | | same |
| rated root moment; thrust | 171.74 N m (−3.18 %); 505.5 N (−2.3 %) | | `reference_blades.json` |
| chord control points [mm] | 271.4 / 149.7 / 135.2 / 60.0 / 60.0 | | `result_delta0.json` → `x_m` |
| twist control points [°] | 23.65 / 8.09 / 6.24 / −0.03 / −0.12 | | same |
| active constraints | AEP floor, stress, deflection, minimum chord at points 3 and 4, chord monotonicity between 3 and 4 | | `active_set`, `kkt` |
| material saved per 1 % of energy, at the margin | 11.02 % | 2.31 % | `kkt` multiplier of the AEP constraint |
| optimiser | SLSQP exit 0, 18 iterations, 31 objective solves, 19 energy adjoints, 22 load solves, no failed trial steps | | `result_delta0.json` |
| multi-start | ten starts (x₀, x_c and eight random) agree to 1.1 × 10⁻⁶ in u | | `multistart_delta0.json` |

The energy–material trade-off (`pareto.json`): for δ = 0 / 0.25 / 0.5 / 1 /
2 %, the shell saving is 3.379 / 4.605 / 5.271 / 6.187 / 7.530 % and the solid
saving 7.278 / 8.557 / 9.666 / 11.490 / 14.294 %. The marginal rates from the
KKT multipliers are 11.02 / 3.11 / 2.31 / 1.56 / 1.22, and every secant
between neighbouring points lies between the rates at its ends. Cold and warm
starts agree to within 5.7 × 10⁻⁶ in u.

Removing constraints one at a time at δ = 0 (`ablation.json`) shows what each
one costs. With all of them the saving is 3.379 %. Without the deflection
constraint it is 4.303 % (deflection ratio 1.078); without the stress
constraint, 4.527 % (stress ratio 1.702); without both, 6.532 % (ratios 1.981
and 1.131); without the manufacturing constraints, 4.238 % (the tip drops to
the 45 mm bound). Removing the moment cap changes nothing, because it isn't
active.

### 7.4 Gradient verification

- **Tier 1** (analytic partial derivatives against the complex step): within
  1 × 10⁻¹³ · max(1, |∂|) on the energy chain, and 1.1 × 10⁻¹⁵ on the moment
  and deflection chains at x₀ and x_m (`tests/test_adjoint_partials.py`,
  `mass_optimisation/checks.json`).
- **Tier 2** (transpose identity): within 1 × 10⁻¹⁴; 3.9 × 10⁻¹⁶ and
  3.4 × 10⁻¹⁶ for the load systems.
- **Tier 3** (assembled gradient against central finite differences at
  h\* = 3.16 × 10⁻⁶): every component is within 3 ε_j at the mid-run iterate
  and at x\*. At x₀, one component (chord_4) is 16.30 times the acceptance
  scale. Tier 4 traces its residual (8.8 × 10⁻¹⁰) to the round-off floor of
  the objective itself (measured at 4.0 × 10⁻¹⁰, so 2.2 floors), not to the
  adjoint. The fragile quantity is the acceptance estimate ε_j. I left the
  corresponding test failing so the case stays visible, rather than loosening
  the tolerance (`verification/gradient_verification/README.md`,
  `tier3.json`, `tier4.json`). The stress and deflection constraints pass at
  x₀ and x_m (worst case 1.014 ε_j).
- **Cost.** One adjoint gradient costs 1.09–1.12 objective evaluations,
  independent of the number of design variables (fitted exponent −0.003 over
  n = 10 to 160). Central finite differences cost 2n evaluations, so the
  adjoint is 18 times faster at n = 10 and 297 times faster at n = 160
  (`verification/cost_scaling/README.md`).

### 7.5 Robustness of x_m

These are evaluations only, with no re-optimisation
(`mass_optimisation/cross_evaluation.json`). At four resource corners
(k ∈ {1.5, 2.0}, c ∈ {7.049, 7.633} m/s), x_m stays within −0.006 to
+0.011 % of x₀'s energy. Under lower rotor-speed ceilings it loses energy
relative to x₀: −0.100 % at 60 m/s, −0.444 % at 55 m/s and −1.840 % at
50 m/s (and +0.084 % with no ceiling). The main result therefore holds for
the 300 rpm operating law. The resource range used here is an assumption of
the check.

### 7.6 Solver verification

The residual converges to |R|/R₀ ≈ 10⁻¹⁵ at 247 stations across the NREL
Phase VI operating envelope (`verification/phase_vi/`). On the Phase VI
geometry, the mean |C_P| difference is 0.51 % against CCBlade and 1.83 %
against pyBEMT (`docs/validation/bem-cross-validation.md`). Those two runs
used an earlier version of the S809 polar table and haven't been repeated
since; the QBlade comparison is current. These are comparisons between BEM
codes; measured Phase VI performance data isn't part of the verification.

## 8. Limitations

- The material result is a saving in shell area, in per cent. The kilogram
  figures apply that saving to one assumed laminate and a 2 mm skin; they are
  not a weighed blade.
- The structural constraints are relative measures based on a thin shell.
  There is no model of a spar, root insert, buckling or fatigue. The absolute
  stress and deflection are thin-shell checks at the operating loads, not a
  certification against the IEC 61400-2 load cases, and no tip clearance is
  specified.
- Loads are steady and aerodynamic, at and below rated. There are no gust,
  parked, yawed or fault cases.
- There is one wind resource and one operating law. Sensitivity to both is
  checked by evaluation only (section 7.5).
- The polars come from XFOIL at N_crit = 9, checked against UIUC wind-tunnel
  data at Re 100–500 k. There is no CFD.
- The solver is cross-validated against other BEM codes, not against measured
  rotor performance.
- The rotor radius is fixed. The root cylinder only appears in the renders.
  Nothing was manufactured or field-tested.
- The number of control points (5 + 5) was chosen using the energy function.

## 9. Where each number comes from

| quantity | file | key or section |
|---|---|---|
| site, atmosphere, resource | `config/site.yaml` | `location`, `atmosphere`, `wind_resource` |
| resource extrapolation | `verification/wind_resource/wind_resource_20m.json` | `extrapolation`, `cross_check_log_law` |
| rotor, operating law, bounds, constraints | `config/rotor_design.yaml` | `geometry`, `operating`, `parameterisation`, `manufacturing`, `structure` |
| reference blade | `verification/baseline/baseline_reference.json`, `x0.json` | `design_point`, `aep_mwh_per_year`, `root_bending_moment_nm`, `peak_thrust_n` |
| spline fit error | `verification/spline_fit_error/fit_error.json`, README | table |
| number of control points | `verification/representation_study/README.md` | table |
| smoothness | `verification/smoothness_gate/README.md` | verdict |
| finite-difference step size | `verification/fd_step_size/sweep.json`, README | `h*_j`, `eps_j` |
| Tiers 3 and 4 | `verification/gradient_verification/tier3.json`, `tier4.json`, README | the one failure; Tier 4 |
| cost scaling | `verification/cost_scaling/scaling.json`, README | table |
| energy optimum x\* | `verification/adjoint_optimisation/result.json`; `fd_optimisation_multistart/results.json` | `x_star`, `u_star`; spread |
| moment-capped energy optimum x_c, trade-off | `verification/load_constraint/result_eps0.json`, `pareto.json`, `checks.json` | `x_c`, `kkt` |
| x_m, trade-off curve, ablation, checks, robustness | `verification/mass_optimisation/result_delta0.json`, `pareto.json`, `ablation.json`, `checks.json`, `cross_evaluation.json`, `reference_blades.json` | as named |
| kilograms, MPa and mm | `verification/absolute_material/result.json`, README; `docs/OUTSTANDING-INPUTS.md` §11 | `blades[*].absolute`, `verification`, `runs`, `comparison` |
| adjoint derivation | `docs/adjoint_derivation.md` | §1–§8 energy, §10 moment, §11 deflection |
| cross-validation | `docs/validation/bem-cross-validation.md` | results |
