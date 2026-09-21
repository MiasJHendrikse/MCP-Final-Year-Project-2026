# Design basis — the minimum-material blade at the Schmitz reference energy

This is the one document that states what the project designs, under which
assumptions, with which model, and what it found. Every number below is read
from a committed artefact; the path is given beside it. The README points
here; the report is written from here and from the artefacts this file names.

## 1. Purpose and research question

**Design purpose.** Find the manufacturable blade of least material for a
2.0 m-radius, three-bladed, fixed-pitch small wind turbine at a low-wind,
high-altitude site in the Khomas Hochland, Namibia, subject to the blade
delivering the annual energy of the classical Schmitz design for the same
rotor under the machine's operating law, at no greater root flapwise moment,
root-stress proxy or tip-deflection proxy than the Schmitz blade, with a
buildable minimum chord and monotone planform and twist.

Annual energy is a constraint of the design problem, not its objective. The
reason is a property of the design basis: under a fixed tip-speed-ratio
operating law with a rotor-speed ceiling, the annual energy of a
polar-consistent Schmitz blade is within 0.15 % of the best that any blade of
the same rotor delivers (§6.2), so the design freedom that a verified
gradient method can spend lies in material, loads and buildability, not in
energy.

**Research questions.**

1. *Design.* What is the minimum-material manufacturable blade that delivers
   the Schmitz reference's annual energy at this site under the operating
   law, at no greater root flapwise moment, root-stress proxy or
   tip-deflection proxy — and what is the exchange rate between annual energy
   and material along the constrained front?
2. *Method.* Are discrete-adjoint gradients of a blade-element momentum
   energy functional, a root-moment functional and a tip-deflection
   functional accurate enough to drive that constrained design, and how do
   they compare with finite differences in accuracy, cost and robustness?

## 2. Site and resource (assumptions of the model, with their basis)

| quantity | value | basis | artefact |
|---|---|---|---|
| Site | 22°25′40.7″ S, 16°33′28.2″ E; highland ridge west of Windhoek; ≈ 1 800 m AMSL | pinned satellite location; the elevation is the assumed ridge height | `config/site.yaml` `location:` |
| Hub height | 20 m | typical free-standing tower for this rotor class | `config/site.yaml` |
| Resource at 50 m | Weibull k = 1.87, A = 8.8 m/s, V̄ = 7.81 m/s | GASP point extraction at the site (a point, not an area; omni-directional) | `verification/wind_resource/wind_resource_20m.json` → `reference_level` |
| Resource at 20 m | k = 1.7092, c = 7.2738 m/s, V̄ = 6.4876 m/s | Justus & Mikhail height extrapolation, shear exponent 0.2079; log-law cross-check in the same file | `config/site.yaml` `wind_resource:`; `wind_resource_20m.json` → `extrapolation`, `cross_check_log_law` |
| Temperature | 293.15 K (annual mean) | design value for the site, not the ISA temperature | `config/rotor_design.yaml` `atmosphere:` |
| Pressure | 81 489 Pa | ISA barometric relation at 1 800 m | same |
| Air density | ρ = 0.9684 kg/m³ | ρ = p/(R_s T); ≈ 21 % below sea level | same |
| Kinematic viscosity | ν = 1.869 × 10⁻⁵ m²/s | μ = 1.81 × 10⁻⁵ Pa s at 20 °C over ρ | same |

The resource uncertainty carried by the design is handled by evaluation, not
by re-optimisation: the settled blades are re-evaluated at four Weibull
corners (§6.5).

## 3. Machine (assumptions of the model, with their basis)

| quantity | value | basis | artefact |
|---|---|---|---|
| Rotor radius | R = 2.0 m, fixed | set by polar-validation quality (station Reynolds numbers inside the measured SG6043 range); not a design variable | `config/rotor_design.yaml` `geometry:` |
| Blades | B = 3 | | same |
| Root cut-out | 0.15 R = 0.30 m | conventional for this size; a modelling choice | same |
| Design tip-speed ratio | λ_d = 6.5 | | `operating:` |
| Rated / cut-in / cut-out wind speed | 11 / 3 / 20 m/s | ≈ 1.7 V̄ (11 / 6.49); low for a low-wind site | same |
| Rated power | P_rated = 3 822.19 W, held constant above rated | the reference blade's aerodynamic power at 11 m/s and λ = 6.5; a nameplate does not change with the blade | `operating.rated_power_w` and its comment |
| Rotor-speed ceiling | Ω_max = 300 rpm (62.8 m/s tip speed; V_c = 9.67 m/s) | residential-scale noise-conscious tip-speed limit (~60–65 m/s); nearest commercial analogues (Skystream 3.7 ≈ 330 rpm / 64 m/s; direct-drive 3 kW PMGs at 250–300 rpm) | `operating.max_rotor_speed_rpm` and its comment |
| Operating law | λ(V) = min(6.5, Ω_max R / V); P = min(P_aero, P_rated) | variable-speed fixed-TSR tracking with a rotor-speed ceiling and an ideal power hold above rated | `src/objective/power.py::tsr_schedule` |
| Airfoil | SG6043 throughout the span | low-Reynolds performance at this scale; S809 is used only for solver cross-validation | `data/airfoils/README.md`, `data/polars/sg6043/README.md` |

Each machine value is an assumption of the model with a stated basis. A
specified generator (rated rpm, nameplate) replaces the ceiling and rating
and every optimisation artefact is re-run, in the order given in
`verification/README.md`. The above-rated region is modelled as an ideal
power hold; loads are therefore constrained at and below rated only (§5).

## 4. Model

- **Solver.** Blade-element momentum in the single-residual Ning form, one
  state φ per station and operating point; Prandtl tip and hub loss on the
  momentum side; Glauert/Buhl high-induction relation with a C¹ blend at
  a = 0.4; bracketing root-find with guaranteed convergence; 25 strips
  (`src/bem/`; `docs/adjoint_derivation.md` §1–§5; `verification/phase_vi/`).
- **Polars.** XFOIL SG6043 cache, Re 40 k – 1 M, N_crit = 9, Viterna
  post-stall extension, C¹ (α, Re) interpolant with analytic partials; XFOIL
  never runs inside the loop (`data/polars/sg6043/`, `src/polars/`,
  `verification/polar_interpolant/`, `results/sg6043_uiuc_validation/`,
  `results/ncrit_sensitivity/`).
- **Energy functional.** AEP over 17 wind-speed bins of 1 m/s, Weibull bin
  masses from the CDF, the operating law above, capped bins as constants;
  J is C¹ with one C² defect where stations cross the Buhl threshold
  (`src/objective/`; `verification/smoothness_gate/README.md`).
- **Loads.** Per-blade normal load ½ρw²cC_n; root flapwise moment; spanwise
  moment; Euler–Bernoulli unit-load tip deflection with section stiffness
  ∝ c³; KS aggregate (ρ_KS = 100) over the nine operating points at or below
  rated (`src/objective/loads.py`; `docs/adjoint_derivation.md` §10–§11).
- **Material proxies.** Shell k_P ∫c dr (k_P = 2.048505, the SG6043
  perimeter per unit chord) is the objective; solid k_A ∫c² dr
  (k_A = 0.068502, the section area per unit chord²) is reported
  (`src/objective/mass.py`).
- **Laminate.** One stated construction: a 2 mm skin of the E-LT-5500/EP-3
  unidirectional E-glass/epoxy laminate (ρ = 1920 kg/m³, E_L = 41.8 GPa,
  UCS_L = 702 MPa; Griffith & Ashwill 2011, SAND2011-3779, Table 19) at
  constant thickness round the whole section, no spar cap, web or root
  insert; design allowable 702 / 3.5725 = 196.5 MPa (GL combined factor for
  a wet hand-laid, non-post-cured laminate). Thin-shell section
  coefficients k_I = 0.003293 (I = k_I c³t) and k_Z = 0.051643 (Z = k_Z c²t)
  exact from the SG6043 contour. Representative values for that
  construction, not measurements of a built blade
  (`docs/MATERIALS-STRUCTURAL-INPUTS.md`, `config/rotor_design.yaml::structure`).
- **Parameterisation.** Clamped cubic B-splines, 5 chord + 5 twist control
  points, d ∈ ℝ¹⁰, scaled to u ∈ [0, 1]¹⁰ by the box (chord 0.045–0.30 m,
  twist −2° to 35°); control-point count selected in
  `verification/representation_study/` (chord RMS 0.515 mm, twist RMS 0.106°
  against the analytic Schmitz blade at 5 + 5).
- **Reference blade x₀.** The Schmitz distribution with wake rotation at the
  SG6043 maximum-L/D point (Re 200 k, α = 5.36°, C_l = 1.282, L/D = 98.5),
  λ = 6.5, projected onto the spline (fit error 1.52 mm max / 0.52 mm RMS
  chord, 0.265° / 0.106° twist; −0.0154 % of AEP) —
  `verification/baseline/x0.json`, `verification/spline_fit_error/`.
- **Gradients.** Discrete adjoint on the diagonal ∂R/∂φ system for three
  right-hand sides (energy, root moment, tip deflection) sharing one forward
  state; the material gradient is analytic; verification in four tiers
  (`src/adjoint/`, `src/gradients/problem.py`; §6.4).

## 5. The design problem

    minimise    m(u) = k_P ∫ c dr / (k_P ∫ c dr)(x0)          shell material, geometric, exact gradient
    subject to  AEP(u) / AEP(x0) − (1 − δ)   ≥ 0             energy floor (δ = 0 is the headline), energy adjoint
                KS0 − KS(u)                  ≥ 0             root-moment cap, moment adjoint
                KS0 / c00² − KS(u) / c0(u)²  ≥ 0             root-stress proxy (Z ∝ c² t), moment adjoint
                D0 − D(u)                    ≥ 0             tip-deflection proxy (I ∝ c³ t), deflection adjoint
                c_i − c_{i+1} ≥ 0 (4), θ_i − θ_{i+1} ≥ 0 (4) monotone control points, constant Jacobian
                c_i − 0.060 ≥ 0 (5)                          buildable-tip floor, constant Jacobian
                polar envelope (50 rows), local solidity ≤ 0.5 (25 rows), u ∈ [0, 1]¹⁰

(`verification/mass_optimisation/README.md` "The problem";
`src/gradients/problem.py::constraints_for_mass_problem`; SLSQP,
ftol 1e-10, maxiter 400.)

Assumptions the rows carry:

- **The three structural rows are relative to x₀** and share one
  thin-shell assumption (mass ∝ c, section modulus ∝ c²t, stiffness ∝ c³t).
  They state that the design is no worse than the Schmitz blade on each
  measure, and no material input enters them. The laminate of §4 turns the
  same quantities into kilograms, megapascals and millimetres
  (`verification/absolute_material/`) and adds an absolute stress row
  (σ ≤ 196.5 MPa) that is verified and slack at every committed blade: at
  the operating loads of the load set it never comes within a factor of
  five of binding, so put in place of the relative row it would remove the
  stress constraint in effect (the optimum returns to the no-stress-row
  ablation, −4.53 %, a root at 1.70× the reference's stress per unit load).
  The sizing loads for a small blade's root — parked, gust, fatigue — are
  outside the model, so the relative rows are the design rows and the
  absolute figures are checks. The absolute deflection row needs a tip
  clearance, which is machine geometry and unresolved.
- **The buildable-tip floor of 60 mm** is a manufacturing judgement with a
  stated basis: a 60 mm SG6043 section is about 6 mm thick, which
  accommodates two 2 mm skins, a bond line and a web; the 45 mm box bound is
  the thickness of the whole section and is not a buildable limit. The floor
  lies below the reference tip chord (67 mm), so x₀ and every quantity
  normalised by it are unchanged. It is a row, not the box, because the box
  defines the scaling of u in which every committed gradient artefact was
  measured (`config/rotor_design.yaml` `manufacturing:` comment).
- **Loads are constrained at and below rated** (nine operating points; the
  rated point at 11 m/s, 300 rpm, λ = 5.71 carries KS weight 0.957). The
  cut-out condition is reported as a post-check only, because the model
  holds rated power above 11 m/s without a limiting mechanism
  (`verification/load_constraint/README.md`).
- **Monotone control points** imply a monotone spline (variation
  diminishing); the rows are on the control points, sufficient but not
  necessary.

## 6. Settled results

### 6.1 The reference blade x₀

AEP 10.2477 MWh/yr (band 8–12); rated point 11 m/s, 300 rpm, λ = 5.712,
C_P = 0.4516, C_T = 0.7029; root flapwise moment 177.38 N m per blade;
rotor thrust 517.5 N; chord control points 275.8 / 181.6 / 97.1 / 78.6 /
67.0 mm; twist 24.53 / 9.73 / 4.04 / 1.00 / 0.59°; shell area 0.43705 m² per
blade — `verification/baseline/baseline_reference.json`, `x0.json`;
`verification/mass_optimisation/reference_blades.json`.

### 6.2 Energy is nearly insensitive to blade shape (the design-basis fact)

| blade | AEP vs x₀ | shell vs x₀ | solid vs x₀ | artefact |
|---|---|---|---|---|
| x\* — unconstrained energy optimum | +0.147 % | | | `verification/adjoint_optimisation/result.json` (adjoint and FD drivers coincide; eight random starts converge to the same point, spread 0.0166 in u — `fd_optimisation_multistart/`) |
| x_c — energy optimum under the root-moment cap | +0.1465 % | +6.171 % | +14.196 % | `verification/load_constraint/result_eps0.json`; `mass_optimisation/reference_blades.json` |

The Schmitz reference is the analytic maximiser of C_P at λ_d for this polar,
and under the operating law the blade sees a single λ in the energy-carrying
bins; the verified optimum is therefore a tenth of a per cent above it while
carrying six per cent more shell material. Reducing the cap on the KS root
moment by 2 / 5 / 10 % costs 0.0206 / 0.1282 / 0.5785 % of AEP, with shadow
prices 0.0017–0.129 % AEP per 1 % KS (`load_constraint/pareto.json`,
`result_eps*.json`).

### 6.3 The minimum-material blade x_m

| | δ = 0 (at the reference energy) | δ = 0.5 % | artefact |
|---|---|---|---|
| shell material vs x₀ | **−3.379 %** | **−5.271 %** | `result_delta0.json`, `pareto.json` (`steps[2]`) |
| solid proxy vs x₀ | −7.278 % | −9.666 % | same |
| AEP vs x₀ | −1.4 × 10⁻¹⁰ % (floor active) | −0.500 % | same |
| KS / KS₀ (moment cap) | 0.9685, slack | | same |
| stress and deflection proxy ratios | 1.000, 1.000 (both active) | | same |
| shell mass per blade, recorded laminate | **1.622 kg** (x₀ 1.678 kg; −0.057 kg, −0.17 kg per rotor) | | `verification/absolute_material/result.json` |
| thin-shell root stress; tip deflection | 22.6 MPa (11.5 % of 196.5 MPa); 111.8 mm | | same |
| rated root moment; thrust | 171.74 N m (−3.18 %); 505.5 N (−2.3 %) | | `reference_blades.json` |
| chord control points [mm] | 271.4 / 149.7 / 135.2 / 60.0 / 60.0 | | `result_delta0.json` → `x_m` |
| twist control points [°] | 23.65 / 8.09 / 6.24 / −0.03 / −0.12 | | same |
| active set | AEP floor, stress, deflection, min-chord 3 and 4, monotone chord 3–4 | | `active_set`, `kkt` (NNLS) |
| exchange rate at the margin | 11.02 % material per 1 % energy | 2.31 | `kkt` multiplier of the AEP row |
| optimiser | SLSQP exit 0, 18 iterations, 31 objective solves, 19 energy adjoints, 22 load solves; no failed trial | | `result_delta0.json` |
| multi-start | ten starts (x₀, x_c, eight random) agree to 1.1 × 10⁻⁶ in u | | `multistart_delta0.json` |

Energy–material front (`pareto.json`): δ = 0 / 0.25 / 0.5 / 1 / 2 % → shell
saved 3.379 / 4.605 / 5.271 / 6.187 / 7.530 %, solid 7.278 / 8.557 / 9.666 /
11.490 / 14.294 %; KKT rates 11.02 / 3.11 / 2.31 / 1.56 / 1.22; every secant
lies between the KKT rates at its ends. Cold and warm starts agree to
≤ 5.7 × 10⁻⁶ in u.

Ablation at δ = 0 (`ablation.json`): full set 3.379 %; without deflection
4.303 % (deflection ratio 1.078); without stress 4.527 % (stress ratio
1.702); without both 6.532 % (1.981, 1.131); without the manufacturing block
4.238 % (tip on the 45 mm box bound); without the moment cap 3.379 %
(unchanged — the cap is slack).

### 6.4 Gradient verification

- Tier 1 (analytic partials vs complex step): ≤ 1 × 10⁻¹³ · max(1, |∂|) on
  the energy chain; ≤ 1.1 × 10⁻¹⁵ on the moment and deflection chains at x₀
  and x_m (`tests/test_adjoint_partials.py`, `mass_optimisation/checks.json`).
- Tier 2 (transpose identity): ≤ 1 × 10⁻¹⁴; 3.9 × 10⁻¹⁶ / 3.4 × 10⁻¹⁶ for
  the load systems.
- Tier 3 (assembled gradient vs central FD at h\* = 3.16 × 10⁻⁶): all
  components within 3 ε_j at the mid-run iterate and x\*; at x₀ one
  component (chord_4) exceeds the acceptance scale at 16.30. Tier 4 attributes
  its residual (8.8 × 10⁻¹⁰) to the round-off floor of J (4.0 × 10⁻¹⁰
  measured; 2.2 floors), not to the adjoint; the estimator ε_j is the fragile
  quantity. Recorded as found; the corresponding test is left failing so that
  the state is visible (`verification/gradient_verification/README.md`,
  `tier3.json`, `tier4.json`). Stress and deflection rows pass at x₀ and x_m
  (worst 1.014 ε_j).
- Cost: one adjoint gradient costs 1.09–1.12 objective evaluations
  independent of n (p = −0.003 over n = 10…160); central FD costs 2n;
  ratio 18 at n = 10, 297 at n = 160 (`verification/cost_scaling/README.md`).

### 6.5 Robustness of x_m (evaluation only, no re-optimisation)

`mass_optimisation/cross_evaluation.json`: over four resource corners
(k ∈ {1.5, 2.0}, c ∈ {7.049, 7.633} m/s) x_m stays within −0.006 … +0.011 %
of x₀'s energy. Under lower rotor-speed ceilings it loses energy relative to
x₀: −0.100 % at 60 m/s, −0.444 % at 55 m/s, −1.840 % at 50 m/s (+0.084 %
with no ceiling). The headline is therefore stated for the configured
300 rpm law. The resource band used here is an assumption of the
cross-evaluation.

### 6.6 Solver verification

Residual convergence to |R|/R₀ ≈ 10⁻¹⁵ at 247 stations across the NREL
Phase VI envelope (`verification/phase_vi/`); mean |C_P| deviation 0.51 %
against CCBlade and 1.83 % against pyBEMT on the Phase VI geometry
(`docs/validation/bem-cross-validation.md`). Measured Phase VI performance
data are not part of the verification set; this is cross-validation between
codes.

## 7. Scope and limitations (stated once)

- The material result is a shell-area saving in per cent; the kilogram
  figures are that saving under one stated laminate and a 2 mm skin, not a
  weighed blade.
- The structural design rows are relative proxies on a thin-shell
  assumption; no spar, root insert, buckling or fatigue model. The absolute
  stress and deflection are thin-shell checks at the operating loads, not
  a certification against IEC 61400-2's load cases; no tip clearance is
  specified.
- Loads are aerodynamic, steady, at and below rated; no gust, parked, yawed
  or fault case.
- One resource and one operating law; sensitivity is by evaluation (§6.5).
- Polars are XFOIL at N_crit = 9 validated against UIUC wind-tunnel data at
  Re 100–500 k; no CFD.
- The solver is cross-validated against other BEM codes, not against
  measured rotor performance.
- Rotor radius fixed; root cylinder drawn in renders only; no manufacture or
  field test.
- The parameterisation (5 + 5) was selected on the energy functional.

## 8. Artefact map

| quantity | file | key / section |
|---|---|---|
| site, atmosphere, resource | `config/site.yaml` | `location`, `atmosphere`, `wind_resource` |
| resource extrapolation | `verification/wind_resource/wind_resource_20m.json` | `extrapolation`, `cross_check_log_law` |
| rotor, operating law, bounds, rows | `config/rotor_design.yaml` | `geometry`, `operating`, `parameterisation`, `manufacturing`, `structure` |
| reference blade | `verification/baseline/baseline_reference.json`, `x0.json` | `design_point`, `aep_mwh_per_year`, `root_bending_moment_nm`, `peak_thrust_n` |
| spline fit error | `verification/spline_fit_error/fit_error.json`, README | table |
| control-point count | `verification/representation_study/README.md` | table |
| smoothness | `verification/smoothness_gate/README.md` | "Verdict" |
| FD step size | `verification/fd_step_size/sweep.json`, README | `h*_j`, `eps_j` |
| Tier 3/4 | `verification/gradient_verification/tier3.json`, `tier4.json`, README | "One failure, and what it is", "Tier 4" |
| cost scaling | `verification/cost_scaling/scaling.json`, README | table |
| energy optimum x\* | `verification/adjoint_optimisation/result.json`; `fd_optimisation_multistart/results.json` | `x_star`, `u_star`; spread |
| load-constrained optimum x_c, Pareto | `verification/load_constraint/result_eps0.json`, `pareto.json`, `checks.json` | `x_c`, `kkt` |
| x_m, front, ablation, checks, cross-evaluation | `verification/mass_optimisation/result_delta0.json`, `pareto.json`, `ablation.json`, `checks.json`, `cross_evaluation.json`, `reference_blades.json` | as named |
| kilograms, MPa, mm; the absolute rows | `verification/absolute_material/result.json`, README; `docs/MATERIALS-STRUCTURAL-INPUTS.md` | `blades[*].absolute`, `verification`, `runs`, `comparison` |
| adjoint derivation | `docs/adjoint_derivation.md` | §1–§8 energy, §10 moment, §11 deflection |
| cross-validation | `docs/validation/bem-cross-validation.md` | "Results" |
