# Discrete adjoint of the Ning-form BEM objective — derivation, implementation and verification

**Status: implemented and verified, B0–B6 complete (2026-09-13).** Written as a
plan at the end of the Stage A (finite-difference) run, per `the implementation plan`
§6, and completed in the Stage B run. Every partial below is derived by hand from
the code as it stands (`src/bem/station.py`, `src/bem/corrections.py`,
`src/bem/rotor.py`, `src/polars/polar.py`, `src/objective/*`), transcribed into
`src/adjoint/kernels.py`, and **verified** against a complex step of the code's own
residual and objective to `1e-13 · max(1, |partial|)` (Tier 1, §9.2). The four tiers'
measured numbers are in §9; the code is `src/adjoint/`, the tests
`tests/test_adjoint_*.py`, the artefacts `verification/gradient_verification/` and
`verification/adjoint_optimisation/`.

The bounds used in the measured runs of §9 are the ones recorded in each
artefact's README; the configured box is in `config/rotor_design.yaml`. Nothing
here reads a bound from `config/`.

**Revision 2026-09-19 — fixed generator rating.** The objective's cap above rated
was `P_aero(V_rated; d)`, floating with the design; `docs/AEP_GAIN_AUDIT.md` §3.2
found that indefensible (a nameplate does not grow because the blade improved) and
it is now the configured constant `operating.rated_power_w` (the reference
blade's own `P_aero(11 m/s; x0)` at λ = 6.5; a nameplate replaces it). For the adjoint this
removes the rated solve as an operating point — the system is **17 × 25 = 425**
states, not 450 — and turns the capped bins into a constant term in `J` with no
weight on any state (§7). §1 and §7 below are restated in the current form; the
measured numbers in §9 are the 2026-09-13 record of the 450-state system and are
left as written, with the re-run Tier 3/4 numbers in
`verification/gradient_verification/README.md`. Tier 1 and Tier 2 pass unchanged
on the 425-state system (`pytest tests/test_adjoint_*.py`).

**Revision 2026-09-19 (Phase 4).** §10 appended: the second right-hand side —
the per-blade flapwise root-moment integrand `m` and its partials in §6's
notation, the nine-point KS aggregate `KS_rho(M/M_ref)`, the same diagonal
adjoint on the 9 × 25 system, and the four tiers with the measured numbers from
`tests/test_loads.py` and `verification/load_constraint/`. The functional is
`src/objective/loads.py` (forward path) and `src/adjoint/loads.py`
(`RootMomentSystem`); the ScaledProblem constraint is
`src/gradients/problem.py::moment_constraint`.

---

## 1. State choice: `φ` only

**State:** `x = φ_{b,i}`, one inflow angle per (operating point `b`, station `i`).
17 operating points (the bin midpoints `V_b = 3.5 … 19.5 m/s`) × 25 stations =
**425 scalars**. (Until 2026-09-19 the rated solve at `V = 11 m/s` was an
eighteenth point, 450 scalars; see the revision note above.)

**Design:** `d ∈ ℝ¹⁰ = [c₀…c₄ (m), θ₀…θ₄ (rad)]`, the chord and twist control
points. Per-station chord and twist are the *constant* linear maps
`c = N_c d`, `θ = N_θ d` (`BladeParameterisation.dchord_dd()`, `dtwist_dd()`,
each (25, 10)). No function in the adjoint takes per-station values as a design
variable; `N_c`, `N_θ` are the only route from stations to `d`.

**Why `a`, `a'` are not states.** `station.residual(phi, station)` evaluates, at a
trial `φ`, the closed-form chain `φ → (Cl, Cd) → (Cn, Ct) → Y → a(Y, F)`. `a` is an
explicit function of `(φ, c, θ)`, not the result of an inner iteration; the same is
true of `a'` (which the residual does not even use — see the `residual` docstring on
the κ' substitution). Choosing `x = φ` therefore makes each `R_{b,i}` a scalar
equation in a scalar unknown with every other quantity explicit, so

    ∂R/∂x  is diagonal (425 × 425), and the adjoint "solve" is 425 scalar divisions.

Had `(a, a')` been taken as states alongside `φ`, the system would be 3 × 3 per
station with two trivially-satisfied rows, and the derivation would have to track
which of the two algebraically equal Buhl expressions the code selected. With `φ`
alone the Buhl branch is handled by implicit differentiation of the defining
equation (§4.4), which is branch-free. Tier 1 confirms the diagonal structure
directly: a complex step on a single `φ_{b,i}` leaves every other residual entry
exactly zero (`test_dR_dx_is_diagonal`).

**Operating points.** `Ω_b = λ V_b / R` with `λ = 6.5`, so `λ_{r,i} = λ r_i / R` is
bin-independent. Bins differ only through the Reynolds number
`Re_{b,i} = W_{b,i} c_i / ν`, `W_{b,i} = hypot(V_b, Ω_b r_i) = V_b √(1 + λ_{r,i}²)`.
`Re` is the **zero-induction estimate** in `rotor.solve_rotor`; it depends on `c_i`
only, not on `φ`. This is what makes the `Re(c)` path a design partial and not a
state partial.

---

## 2. The residual as coded

For one `(b, i)`; subscripts dropped. `s = sin φ`, `k = cos φ`.

    σ   = B c / (2π r)                              StationParams.solidity
    Re  = W c / ν                                   rotor.solve_rotor (constant W)
    α   = φ − θ
    Cl  = Cl(α, Re),  Cd = Cd(α, Re)                CachedPolar.cl/.cd (α in rad)
    Cn  = Cl k + Cd s
    Ct  = Cl s − Cd k
    F   = F_tip(φ) · F_hub(φ)                       corrections.combined_loss_factor
    Y   = σ Cn / s²
    a   = a(Y, F)                                   corrections.corrected_axial_induction
    T   = σ Ct / (4 F s)                            (the κ'-substituted tangential term)
    R   = s / (1 − a)  −  (k − T) / λ_r             station.residual

Power integrand at the solved `φ` (`rotor.solve_rotor` lines 326, 356–363;
`objective.power.aerodynamic_power`):

    w   = V (1 − a) / s
    q   = ½ ρ w² B c Ct r
    P_b = Ω_b · Σ_i t_i q_i                          trapezoid weights t_i (constant)
          t_0 = (r_1 − r_0)/2,  t_i = (r_{i+1} − r_{i−1})/2,  t_{n−1} = (r_{n−1} − r_{n−2})/2

(`solve_rotor` returns `Cp = P / (½ρV³πR²)` and `aerodynamic_power` multiplies it
back; the round trip is exact in the derivative.)

Objective (`objective.annual_energy_mwh`, `objective.objective`), with bin masses
`m_b = WeibullResource.probability_between(edges[:-1], edges[1:])`, `T_h = 8766`,
`limited_b = (P_b > P_rated)` exactly as `power_per_bin` flags it, and `P_rated`
the configured generator rating (`operating.rated_power_w`), a constant:

    P̃_b = P_rated  if limited_b  else  P_b
    J    = −(T_h / 10⁶) Σ_b m_b P̃_b                  [MWh/yr]

At `x0`: bins 11.5–19.5 m/s are limited (nine of them); `J = −10.270157`.

**As implemented.** `adjoint.kernels.station_partials` computes `R` and `q` with the
*same calls in the same order* — `corrected_axial_induction(Y, F)`,
`combined_loss_factor(...)`, the residual's own `axial − tangential` expression, the
integrand's `0.5·ρ·w²·B·c·Ct·r` — so the value of `R` matches
`bem.station.residual` to the bit (measured: max difference 0.0 over 450 stations at
both test designs) and `J(solve(d), d)` matches `annual_energy_mwh(d)` to 1.7e-16 and
3.5e-16 relative (assertion 1e-12). `BEMSystem.solve` is one call to
`aerodynamic_power` per operating point, `φ` read from the station records; no
second solve, no re-implemented root-find.

---

## 3. The polar partials the chain needs

`CachedPolar` exposes exactly the four first derivatives the chain needs, already
in the residual's units:

    Cl_α = dcl_dalpha(α)   per radian        Cl_Re = dcl_dre(α)   per unit Re
    Cd_α = dcd_dalpha(α)   per radian        Cd_Re = dcd_dre(α)   per unit Re

(`polars.interpolant` is verified to ~1e-14; it is not re-derived here. It is
complex-safe in both `α` and `Re`, which B0 relies on.)

    ∂Cl/∂φ = Cl_α        ∂Cl/∂θ = −Cl_α        ∂Cl/∂c = Cl_Re · W/ν
    ∂Cd/∂φ = Cd_α        ∂Cd/∂θ = −Cd_α        ∂Cd/∂c = Cd_Re · W/ν

The `∂/∂c` column is the **`Re(c)` path** (§4.3 of the brief) — the one most easily
forgotten because `c` also enters through `σ` and the explicit factor in `q`. In the
kernel `W/ν` is taken as `reynolds / chord`, which is exact because both are linear
in `c`. Tier 1's complex step on a chord control point sends a complex `c` into a
complex `Re` (B0's third edit makes `CachedPolar` keep it), so this column is on the
complex side of the comparison; `test_complex_chord_and_reynolds_propagate` checks
the polar-only path is non-zero on its own, so a dropped `Re(c)` term could not hide
behind the solidity term.

Force coefficients, using `dCn/dφ` picking up both the polar slope and the rotation
of the (Cl, Cd) frame:

    ∂Cn/∂φ = Cl_α k + Cd_α s − Ct              ∂Ct/∂φ = Cl_α s − Cd_α k + Cn
    ∂Cn/∂θ = −(Cl_α k + Cd_α s)                ∂Ct/∂θ = −(Cl_α s − Cd_α k)
    ∂Cn/∂c = (Cl_Re k + Cd_Re s) W/ν           ∂Ct/∂c = (Cl_Re s − Cd_Re k) W/ν

---

## 4. The loss factor and the induction

### 4.1 `dF/dφ`

`F_tip = (2/π) acos(e^{−f})`, `f = (B/2)(R − r)/(r s)`. On the momentum region
`φ ∈ (0, π/2]`, `s > 0`, so the `abs(sin φ)` in the code is inert, as are the
`min(1, e^{−f})` clamp (`f > 0` for `r < R`) and the `s < 1e-8` guard.

    df/dφ    = −f k / s
    dF_tip/df = (2/π) e^{−f} / √(1 − e^{−2f})
    dF_tip/dφ = −(2/π) f e^{−f} k / (s √(1 − e^{−2f}))

Hub identically with `f_h = (B/2)(r − r_hub)/(r_hub s)`. Then

    F = F_t F_h,   F_φ = F_t' F_h + F_t F_h'.

`F` has no `c` or `θ` dependence. Near the tip station (`r/R = 0.983`, `f ≈ 0.026/s`)
the factor `1/√(1 − e^{−2f})` is large but finite; Tier 1 shows the mixed tolerance
is met there (the worst `∂R/∂φ` error over all 450 stations is 4.9e-15). The
formula is also checked in isolation against a complex step of the code's own
`tip_loss_factor` / `hub_loss_factor` to 1e-13 relative
(`test_loss_factor_complex_step_matches_the_analytic_derivative`). The kernel takes
the *value* `F` from `combined_loss_factor` itself and mirrors its guards in the
derivative (derivative 0 wherever the code returns a constant).

### 4.2 Loading

    Y = σ Cn / s²
    ∂Y/∂φ = σ (∂Cn/∂φ) / s² − 2 σ Cn k / s³
    ∂Y/∂c = σ_c Cn / s² + σ (∂Cn/∂c) / s²,   σ_c = B / (2π r)
    ∂Y/∂θ = σ (∂Cn/∂θ) / s²

### 4.3 Momentum branch (`a_m = Y/(4F + Y) ≤ 0.4`)

    a = Y / (4F + Y)
    ∂a/∂Y = 4F / (4F + Y)²         ∂a/∂F = −4Y / (4F + Y)²

### 4.4 Buhl branch (`a_m > 0.4`), by implicit differentiation

The code computes `a` from Ning's γ-form and switches between two algebraically
equal expressions on conditioning (`corrected_axial_induction`). Both are roots
of the single defining equation

    G(a; Y, F) = Y (1 − a)² − Ct_B(a; F) = 0,
    Ct_B = 8/9 + (4F − 40/9) a + (50/9 − 4F) a²        (high_thrust_correction)

so differentiate `G` implicitly — independent of which expression the code chose:

    ∂G/∂a = −2 Y (1 − a) − (4F − 40/9) − 2 (50/9 − 4F) a
    ∂G/∂Y = (1 − a)²
    ∂G/∂F = −4 a (1 − a)                                (∂Ct_B/∂F = 4a − 4a²)

    ∂a/∂Y = −(∂G/∂Y) / (∂G/∂a)        ∂a/∂F = −(∂G/∂F) / (∂G/∂a)

At the blend `a = 0.4` the momentum and Buhl values of `∂a/∂Y`, `∂a/∂F` coincide
(`Ct_B` matches `4aF(1−a)` in value and slope there, asserted in
`tests/test_corrections.py`), so `a(Y, F)` is C¹ and the switch does not put a
kink in `∂R/∂φ`. It does put a C² break there, which is FD's problem, not the
adjoint's (Tier 4, §9.5).

The `γ₂ < 0` branch (`gamma1/gamma3`) is unreachable for positive loading and
never at a solved station; the adjoint does not model it, and
`induction_and_partials` raises if it is reached. Tier 1 asserts that no solved
station is on it. The kernel takes the *value* `a` from `corrected_axial_induction`
itself and decides the branch as the code does, on `a_m.real ≤ BUHL_AC`.

### 4.5 Chain into `a`

    a_φ = a_Y Y_φ + a_F F_φ         a_c = a_Y Y_c         a_θ = a_Y Y_θ

**69 of the 450 solved stations are above `a = 0.4` at `x0`** (4 of 25 at the rated
point), 47 at the perturbed Tier 1 design; the tests assert the count is `> 0`, so
the implicit partials are exercised, not just present.

---

## 5. Residual partials

With `T = σ Ct / (4 F s)`:

    T_φ = (σ/4) [ Ct_φ /(F s) − Ct F_φ /(F² s) − Ct k /(F s²) ]
    T_c = (σ_c Ct + σ Ct_c) / (4 F s)
    T_θ = σ Ct_θ / (4 F s)

    ∂R/∂φ = k/(1 − a) + s a_φ/(1 − a)²  +  (s + T_φ)/λ_r
    ∂R/∂c =              s a_c/(1 − a)²  +  T_c/λ_r
    ∂R/∂θ =              s a_θ/(1 − a)²  +  T_θ/λ_r

`λ_r`, `r`, `R`, `B`, `r_hub` are constants. `∂R/∂φ` is the diagonal of the state
Jacobian; it must be non-zero at every solved station (the bracket guarantees
`R` changes sign across the root). Measured `min |∂R/∂φ|` over all 450 stations:
0.752 at `x0`, 0.937 at the perturbed design (asserted `> 1e-2`).

---

## 6. Power-integrand partials

    w   = V (1 − a) / s
    w_φ = −V [ a_φ / s + (1 − a) k / s² ]
    w_c = −V a_c / s
    w_θ = −V a_θ / s

    q   = ½ ρ B r · c w² Ct
    q_φ = ½ ρ B r · c (2 w w_φ Ct + w² Ct_φ)
    q_c = ½ ρ B r · [ w² Ct + c (2 w w_c Ct + w² Ct_c) ]     ← explicit c, σ→a, Re→Cl,Cd
    q_θ = ½ ρ B r · c (2 w w_θ Ct + w² Ct_θ)

`P_b = Ω_b Σ_i t_i q_{b,i}`, so `∂P_b/∂φ_{b,i} = Ω_b t_i q_{φ,b,i}` and likewise for
`c_i`, `θ_i`. `Ω_b` is a constant (fixed `λ`, `V_b`).

---

## 7. Objective partials, rated-limit weighting, scaling

Let `L = {b : limited_b}` from the forward solve (`power_per_bin`; at an exact tie
`unlimited > rated` is False → not limited) and `M_L = Σ_{b∈L} m_b`. Define the
per-operating-point weight

    ω_b = −(T_h/10⁶) m_b     for unlimited b
    ω_b = 0                  for limited b        (its own solve does not enter J)

and the constant the capped bins contribute,

    J_L = −(T_h/10⁶) P_rated M_L                  (`BEMSystem.J_capped`; no state in it)

so `J = Σ_b ω_b P_b + J_L`. (Until 2026-09-19, with `P_rated = P_aero(V_rated; d)`,
the capped mass was instead a weight `ω_rated = −(T_h/10⁶) M_L` on an eighteenth
operating point, the rated solve. That weight, and that point, are gone.) Then,
for every operating point `b` and station `i`:

    ∂J/∂φ_{b,i} = ω_b Ω_b t_i q_{φ,b,i}
    ∂J/∂c_i     = Σ_b ω_b Ω_b t_i q_{c,b,i}          ∂J/∂θ_i = Σ_b ω_b Ω_b t_i q_{θ,b,i}
    ∂J/∂d       = N_cᵀ (∂J/∂c) + N_θᵀ (∂J/∂θ)                                (10,)

and for the residuals

    ∂R_{b,i}/∂d_j = (∂R/∂c)_{b,i} N_c[i, j] + (∂R/∂θ)_{b,i} N_θ[i, j]         (17, 25, 10)

Scaling (`DesignBounds`, `ScaledProblem`): `d = lo + u ⊙ span`, `fun = J/|J₀|`, so

    dfun/du = (dJ/dd ⊙ span) / |J₀|.

`J₀` is the value recorded by `ScaledProblem` at `u₀`; `ScaledProblem.jac_adjoint(u)`
is exactly `gradient(d) ⊙ span / |J₀|`, and every reported gradient in MWh/yr per
unit `u` is that with the `|J₀|` undone.

`limited` is a *fixed mask* taken from the forward solve; `J` is C¹ but not C²
where a below-rated bin crosses the limit as `d` moves (`objective.power` docstring),
and the adjoint at such a point is the one-sided derivative of whichever side the
forward solve landed on. This is a property of the objective, not of the gradient.
Tier 1 checks the weighting directly: `∂J/∂φ` is exactly zero on the nine limited
bins' rows and non-zero on every unlimited row, on the complex-step side as well
as the derived side (`test_dJ_dphi_matches_complex_step_at_every_station`).

---

## 8. Adjoint solve and assembly

`∂R/∂x` diagonal ⇒ the adjoint equation `(∂R/∂x)ᵀ ψ = −(∂J/∂x)ᵀ` is

    ψ_{b,i} = −(∂J/∂φ_{b,i}) / (∂R_{b,i}/∂φ_{b,i})                     (425 divisions)

    dJ/dd   = ∂J/∂d + Σ_{b,i} ψ_{b,i} ∂R_{b,i}/∂d                         (10,)

`BEMSystem.gradient(d)` does exactly this: one forward solve, the partials at
every station, `psi = -dJ_dx / dR_dx`, `dJ_dd = dJ_dd_explicit + apply_dR_dd_T(psi)`.
The forward-mode twin, `BEMSystem.tangent(v)` (`δφ = −(∂R/∂d v)/(∂R/∂x)`,
`δJ = ∂J/∂x·δφ + ∂J/∂d·v`), is a separate code path used only to test it (Tier 2).

**Honest cost note, with the measured numbers.** With a diagonal state Jacobian
and `n = 10` design variables the adjoint's advantage over central FD is
**accuracy, not solve count**, though the wall time is real:

| | per evaluation at `x0` | per SLSQP run (34 iterations) |
|---|---|---|
| `J` | 0.213 s | — |
| central FD gradient (`2n = 20` evaluations) | 4.27 s | 1,419 objective evaluations, 361 s (A4) |
| adjoint gradient (one solve + partials + assembly) | 0.244 s = 1.15 × `J` | 40 objective + 69 adjoint evaluations, 24.7 s (B5) |

A 14.6 × faster run for the same iterations; the ratio would grow linearly with
`n`. But at `n = 10` the run was already tolerable. What the adjoint actually
bought is a gradient free of the FD noise floor and of the O(h)·Δf″ contamination
FD picks up at the polar interpolant's C² breaks, and — the finding of Tier 3 and
B5 — an *independent* gradient that confirms the FD one was already good enough:
the two optima coincide to 8e-9 in `u`.

---

## 9. Implementation as executed (B0–B6), with the measured numbers

Cost basis from Stage A: `J` = 0.207 s per evaluation; central FD gradient at `x0`
= 20 evaluations ≈ 4.2 s; the 300-evaluation step-size sweep = 62 s; the FD-driven
SLSQP run (A4) = 34 iterations, 1,419 objective evaluations, 361 s, to the
floating-rating optimum of that system (the committed fixed-rating result is
+0.1467 %, `verification/adjoint_optimisation/`, `verification/load_constraint/`). Stage B, this run: `pytest -q` went from 391 to 450 passed (5 xfailed
throughout), 12 s to 29 s; the new tests, scripts and this document are ≈ 2,900
lines. Nothing needed a cache extension (§4.5 of the brief): every point visited
stayed inside Re 62.6 k … 496 k and α −0.44° … 6.13°.

### 9.1 B0 — complex-step safety, real path bit-identical (`46559f9`)

Three edits, the only ones under `src/` outside `adjoint/` and the one method added
to `gradients/problem.py` in B3:

1. `src/bem/station.py::_blade_element`: `sin`/`cos` via `_sin`, `_cos` helpers
   modelled on `corrections._sqrt` (`cmath` when `phi` is complex-typed).
2. `src/bem/corrections.py::tip_loss_factor`, `hub_loss_factor`: a complex branch
   (`_prandtl_factor_complex`) with every decision on `.real`; the real-path lines
   literally unchanged inside the `else`.
3. `src/polars/polar.py::CachedPolar.__init__`: range-check on `.real`, the
   Reynolds number stored as given (`float()` had dropped the imaginary part, and
   with it the `Re(c)` path); `__repr__` guarded.

`tests/test_complex_safety.py` (16 tests): the residual returns complex for a
complex `φ` at every rated-point station of `x0`; its complex-step `∂R/∂φ` agrees
with central FD at 1e-6 to 2.6e-10; complex twist and complex chord/Reynolds
propagate; `CachedPolar` accepts and keeps a complex Re and still range-checks; the
loss-factor complex branch matches the real one to 1e-14 and its complex-step
derivative matches §4.1 to 1e-13. The real path is bit-identical: the full suite
including `test_golden_regression`, `test_determinism` and `test_cost` passed
unchanged (407 passed after B0).

### 9.2 B1 — kernels, system, Tier 1 (`917658a`)

`src/adjoint/kernels.py::station_partials(phi, chord, twist, reynolds, r, lam_r, R,
B, r_hub, polar, v_inf, air_density, derivatives=True)` → `StationPartials` with
`residual, dR_dphi, dR_dc, dR_dtheta, q, dq_dphi, dq_dc, dq_dtheta` plus
`a, da_dphi, da_dc, da_dtheta, F, dF_dphi, alpha, cl, cd, cn, ct, buhl`. Two
departures from the plan's signature, both forced: `v_inf` and `air_density` are
needed for `q` (the plan's signature had no way to form it), and `derivatives=False`
is a values-only mode (the four polar-slope lookups skipped, every partial `None`)
so a complex step of `J` costs a value rather than a gradient — it is what made the
450-step per-station test of `∂J/∂φ` affordable (6 s per design).

`src/adjoint/system.py::BEMSystem(parameterisation, bounds, resource, polar_cache)`
as planned: `solve(d)` → `ForwardState(phi (18,25), power_w (18,), limited (17,), a,
reynolds)`; `residual(phi, d)` through `bem.station.residual`; `partials`, `dR_dx`,
`dR_dd` (via `N_c`, `N_θ`); `J`, `dJ_dx`, `dJ_dd` with §7's weights;
`apply_dR_dd` / `apply_dR_dd_T`; `tangent`; `gradient` → `GradientResult`;
`state_sensitivity` (for Tier 4). Rows are built as lists and `np.array` infers the
dtype: no real-dtype array on any complex path. `"adjoint"` was added to
`SOLVE_PATH_PACKAGES` in `tests/test_invariants.py` (the no-`xfoil` AST check).

**Tier 1, `tests/test_adjoint_partials.py` (26 tests)**, at `x0` and at
(chord × 1.2, twist + 3°), complex step `h = 1e-30` on the code's own residual and
on `J`, assertion `|err| ≤ 1e-13 · max(1, |partial|)` (`1e-12` for `J`). Measured
worst mixed errors:

| partial | how compared | `x0` | perturbed |
|---|---|---|---|
| `∂R/∂φ` (450 entries) | complex `φ` → `station.residual` | 4.9e-15 | 3.8e-15 |
| `∂R/∂d` (450 × 10) | complex `d` → complex `c, θ, Re` → `station.residual` | 6.5e-15 | 4.9e-15 |
| `∂J/∂φ` (450 entries, one step each) | complex `φ_{b,i}` → `J` | 3.4e-15 | 3.4e-15 |
| `∂J/∂d` (10) | complex `d` → `J` | 1.5e-15 | 4.0e-16 |
| `R` value vs `station.residual` | same calls | 0.0 | 0.0 |
| `J(solve(d), d)` vs `annual_energy_mwh` | relative | 1.7e-16 | 3.5e-16 |
| Buhl-branch stations in the set | count | 69 | 47 |
| `γ₂ < 0` stations | count | 0 | 0 |
| `min |∂R/∂φ|` | | 0.752 | 0.937 |

Two decades inside the tolerance everywhere, no plateau near 1e-6: no real dtype
on the path, and no partial wrong.

### 9.3 B2 — Tier 2 transpose identity (`a964d83`)

`tests/test_adjoint_transpose.py` (12 tests), eight random `(u, v)` draws at each
of the two designs, `⟨v, A u⟩ = ⟨Aᵀ v, u⟩` to `1e-14` for `A = ∂R/∂x` (diagonal),
`A = ∂R/∂d` (`apply_dR_dd` vs `apply_dR_dd_T`, separate code paths) and the full
operator `d ↦ dJ/dd` (`tangent` vs `gradient`). One definition was needed that the
plan left open: "relative" is relative to the Cauchy–Schwarz bound `‖v‖ ‖A u‖`,
the scale the operator's round-off lives on — an inner product of two random
vectors cancels, and relative-to-the-value would have tested the luck of the draw
(the first draft did, and failed at 1.06e-14 on a full-chain value of 0.074 whose
terms were O(10), a difference of 8e-16). Measured worst, on that scale:

| operator | `x0` | perturbed |
|---|---|---|
| `∂R/∂x` | 2.1e-17 | 1.8e-17 |
| `∂R/∂d` | 2.8e-17 | 1.8e-17 |
| full chain `d ↦ dJ/dd` | 7.5e-16 | 3.4e-16 |

Also: the matrix-free apply against the explicit (18, 25, 10) tensor with a
cancellation-aware scale; the support of a chord control point's column; and
`dJ/dd == explicit + ∂R/∂dᵀ ψ` with `ψ` solving the diagonal system exactly.

### 9.4 B3 — assembly and Tier 3 (`09f5ae3`)

`ScaledProblem.jac_adjoint(u) = gradient(d) ⊙ span / |J₀|` over a lazily built
`BEMSystem` (`adjoint_system()`); `PolarDomainError` logged and re-raised as in
`J`. The only edit to `src/gradients/`.

`verification/gradient_verification/run_tier3.py` → `tier3.json`, `adjoint_x0.json`,
`v_curve_vs_adjoint.png`, `tier3_agreement.png`, `README.md`. At `x0`, the FD-SLSQP
mid-run iterate `k = 17` and the FD optimum `u*`, all ten variables, adjoint vs
central FD at `h*_j` (A3's per-variable steps), acceptance `|adj − FD| ≤ 3 ε_j`
(`ε_j` A3's at `x0`, re-measured by the three-step local sweep at the other two
points) and `rel ≤ 1e-3`:

| point | worst variable | worst `|diff|/ε_j` | worst rel | Taylor ratios (min of 3 draws × 2 decades) |
|---|---|---|---|---|
| `x0` | `chord_4` (`h*` 3e-7) | 2.78 | 1.5e-7 | 98.4 |
| `k = 17` | `chord_3` (`h*` 1e-6) | 1.22 | 1.3e-5 | 98.6 |
| `u*` | `chord_4` (`h*` 3e-7) | 1.14 | 7.6e-5 | 98.6 |

All 30 pairs pass; 27 are within `1 ε_j`. The Taylor remainder falls 98.4–100.8 ×
per decade at every point and draw (assertion ≥ 30; a wrong gradient gives ≈ 10).
Wall time at `x0`: `J` 0.213 s, FD gradient 4.27 s, adjoint 0.244 s. A3's V-curves
re-plotted against the adjoint (its `--reference` machinery, output kept in
`gradient_verification/` so A3's figure is untouched) become the textbook V:
`h²` truncation down to a minimum of 1e-11 … 1e-10 fun units per `u` at
`h = 1e-5 … 3e-6`, `1/h` round-off below, the adjoint at the bottom for every
variable. `tests/test_adjoint_gradient.py` (4 tests) pins `jac_adjoint` to the
committed A3 reference within `3 ε_j` and asserts the Taylor ratios, so the suite
catches a drift on either side.

### 9.5 B4 — Tier 4, degradation at the polar interpolation (`7cb3ef6`)

`verification/gradient_verification/run_tier4.py` → `tier4.json`,
`tier4_attribution.png`, README section. From the adjoint system's own
linearisation (`dφ/du_j = −(∂R/∂d_j span_j)/(∂R/∂φ)`), for every `(b, i)` and
every variable at each of the three points: whether `α` crosses a 0.5° interpolant
node, `Re` a cached row, or `a` the Buhl blend inside `[u − h*_j, u + h*_j]`; the
round-off floor of `J` measured directly (`J` along a line of 1e-12 steps,
linear-fit residual std `δJ`, floor `δJ/h`); and for the worst Tier 3 variable
the whole 15-step V with counts, the smooth `C h²` law fitted on crossing-free
steps, and the excess over it.

**At `h*_j` no stencil crosses anything, at any point, for any variable** —
the nearest break is 2.4 stencil half-widths away in `α`, 53 in `Re`, 171 in `a`
— and the three headline disagreements are 0.92 ×, 1.37 × and 0.61 × the measured
round-off floor (`δJ` = 1.3e-15, 2.0e-15, 2.5e-15 MWh/yr). The attribution
sentence the brief asked for: **degradation at the polar interpolation at `h*_j`
is zero, attributable to zero knot crossings; Buhl crossings contribute zero
(first crossed at `h ≥ 3e-4` at `u*`, `h ≥ 3e-3` elsewhere, never inside an `h*_j`
stencil); elsewhere agreement is the round-off floor.** Across the whole grid the
only steps departing from `C h²` are steps whose stencil crosses a counted knot,
row or blend, all at `h ≥ 1e-4` — three to four decades above every `h*_j` — and
the truncation-side slopes are 1.69–1.96 (2 = smooth). At `u*` one Buhl crossing
plus four knots at `h = 1e-3` pull the FD error *below* the `h²` line (excess
−1.5e-4 MWh/yr per `u`) and the 3e-4 step above it (+4.7e-6): the non-monotone dip
the brief predicted for FD at a C² break, which the adjoint does not have. The
"shrink `h`" step of the plan was not needed; instead the largest crossing-free
step was found (1e-4, 3e-5, 1e-4) and the FD there shown to sit on the smooth law.
Nothing is left unexplained.

### 9.6 B5 — adjoint-driven SLSQP (`e6ac536`)

`verification/adjoint_optimisation/run_adjoint_slsqp.py` → `result.json`,
`iterates.json`, `optimised_blade.png`, `README.md`: A4's script with
`jac = problem.jac_adjoint`, everything else identical.

| | adjoint (B5) | FD (A4) |
|---|---|---|
| `nit` / `nfev` / `njev` | 34 / 37 / 34 | 34 / 37 / 34 |
| objective evaluations | 40 | 1,419 |
| wall time | 24.7 s | 360.9 s |
| AEP* | 10.292482 MWh/yr (floating-rating system) | 10.292482 MWh/yr (floating-rating system) |
| `‖u*_adj − u*_fd‖∞` | 8.2e-9 | tolerance 8.6e-3 (10 × the multi-start spread) |
| `AEP*_adj − AEP*_fd` | +1.4e-13 MWh/yr | tolerance ≈ 1e-6 (from `ftol`) |
| trajectory: `max_k ‖u_k − u_k^fd‖∞` | 6.8e-7 | — |

The optima coincide, and so do the 34 iterates on the way: the Tier 3 agreement
carried through every quasi-Newton update without the paths separating. No active
bound, no active envelope row, no `PolarDomainError`, α and Re inside the cache.
The optimum is the same as A4's, inside the defect gates and confirmed global by
the multi-start study; the committed fixed-rating figure is +0.1467 %.

### 9.7 B6 — this document (`Phase 3: derivation`)

Completed above: the state choice (§1), the residual as coded (§2), every partial
with its chain including the `Re(c)` path and both `a` branches (§3–§6), the
rated-limit weighting and the scaling chain (§7), the solve and assembly with the
honest cost note (§8), the four tiers with the measured numbers (§9.2–§9.5), the
optimisation comparison (§9.6).

### Failure modes the plan guarded against, and what happened

- **Silent wrong gradient** (R2): Tier 1/2 at machine precision as hard asserts
  (passed two decades inside tolerance); Tier 3 against `ε_j`, never looser (all
  30 pairs pass, worst 2.78 `ε`); Taylor remainder (98–101 × per decade); FD and
  adjoint optima coincide (8e-9 in `u`). No tolerance was adjusted to pass. The
  one scale definition made (Tier 2's Cauchy–Schwarz bound) is recorded in §9.3
  with the reason.
- **Forgetting `Re(c)`**: §3's `∂/∂c` column and `q_c`'s third term are in the
  kernel; Tier 1's complex `d → c → Re` step and B0's polar-only propagation test
  would have caught their absence.
- **Treating `a` as a state**: §1; `test_dR_dx_is_diagonal`.
- **Real-dtype arrays on a complex path**: lists and inferred dtypes throughout;
  no 1e-6 plateau anywhere in Tier 1.
- **Per-station DOF**: `N_c`, `N_θ` are the only route from stations to `d`;
  `dR_dd`, `dJ_dd`, `apply_dR_dd(_T)` all go through them.
- **Bounds**: read from each artefact's README; nothing in `config/` was changed
  by this work; no cache extension was needed.
- **Loosening a tier**: none loosened. Tier 3's residual disagreement was
  unexplained until Tier 4 explained it, and Tier 4 explained it as round-off with
  the polar interpolation contributing zero at `h*_j`.

---

## 10. A second right-hand side: the root moment and its KS aggregate

**Added 2026-09-19 (Phase 4, Steps 2a–2d).** The load constraint is a
*relative* flapwise root-moment cap on the BEM spanwise loading. It reuses the
whole Phase 3 machinery — the same residual `R`, the same diagonal `∂R/∂φ`,
the same `dR/dd` through `N_c`, `N_θ` — with a different scalar functional and
therefore a different right-hand side for the one adjoint solve. Nothing in
§1–§8 changes; this section gives the new functional and its tiers.

### 10.1 The integrand `m` and its partials

In §6's notation, at station `i` of operating point `b`, with
`arm = r_i − r_hub` a station constant (`r_hub = root_fraction · R`):

    m = 1/2 rho · arm · c · w^2 · Cn                              [N·m per station]

with `w = V (1−a)/sin φ` and `Cn = Cl cos φ + Cd sin φ` — the same `w`, `Cn`
the power integrand `q` uses. Two deliberate asymmetries against `q`: **no
blade count `B`** (this is one blade's root load; the rotor-summed figure would
overstate it by `B`) and **no rotor speed `Ω_b`** (the power assembly is
`P_b = Ω_b Σ_i t_i q`; the moment is `Σ_i t_i m` with no such factor).

The partials are the same product rule as `q_φ, q_c, q_θ`, with `c` entering
explicitly, through `σ → a → w`, and through `Re → Cl, Cd` in the third term of
`m_c`:

    m_φ     = 1/2 rho arm · c ( 2 w w_φ Cn + w^2 Cn_φ )
    m_c     = 1/2 rho arm · [ w^2 Cn + c ( 2 w w_c Cn + w^2 Cn_c ) ]
    m_θ     = 1/2 rho arm · c ( 2 w w_θ Cn + w^2 Cn_θ )

`Cn_φ, Cn_c, Cn_θ` already exist in the kernel (they feed `Y`), and `m_c`'s
third term is the `Re(c)` path `dre_dc = Re/c`. `r_hub is None` or `≤ 0` is
treated as `0` exactly as `loss_factor_and_derivative` does. The value code is
shared with `R` and `q`: with `derivatives=False` the kernel returns `m` and
`None` partials, so a complex step of the value path is available.

### 10.2 The load set `L`, the aggregate, the constraint

Per blade, `M_b = Σ_i t_i m_{b,i}` (`t_i` the trapezoid weights that make the
sum equal `bem.rotor._trapz`), and no `Ω_b`, no `B`.

The operating set is fixed at construction from `x0` and is **not** recomputed
per design (a set that moved with `d` would be a non-smooth constraint):

    L = { (V_b, λ_b) : bin midpoint b with P_b(x0) ≤ P_rated }
        ∪ { (11.0, λ(11.0) = 5.711986642890533) }
      = 3.5 … 9.5 m/s at λ = 6.5, 10.5 at λ = 5.984, 11.0 at 5.712   →  nine points

225 states, one scalar residual each. Normalising by
`M_ref = M(x0) = 177.3755406092969 N·m` (11 m/s on the 300 rpm ceiling) makes
`ρ` dimensionless and the constraint O(1):

    Mhat_b         = M_b / M_ref
    KS_ρ(Mhat)     = max_b Mhat_b + ln( Σ_b exp(ρ (Mhat_b − max)) ) / ρ
    ∂KS/∂Mhat_b    = exp(ρ (Mhat_b − max)) / Σ_b exp(...)          (softmax, sums to 1)

so, by the chain rule through the trapezoid and the normalisation,

    ∂KS/∂φ_{b,i}   = softmax_b · t_i · m_φ_{b,i} / M_ref

The adjoint is the same equation as §8 with this right-hand side:

    ψ_{b,i} = −(∂KS/∂φ_{b,i}) / (∂R_{b,i}/∂φ_{b,i})                (225 divisions)
    dKS/dd  = ∂KS/∂d|_explicit + Σ_{b,i} ψ_{b,i} ∂R_{b,i}/∂d

with the explicit part `Σ_b softmax_b t_i (m_c N_cᵀ + m_θ N_θᵀ)/M_ref`. The
constraint SLSQP sees, for reduction fraction `ε`, is

    g_ε(u) = (1 − ε) KS_ρ(x0) − KS_ρ(u)  ≥ 0,   KS_ρ(x0) = 1.0004357284444418

and the scaled Jacobian is `dg/du = −(dKS/dd) ⊙ span`. The limit is **`KS_ρ(x0)`
and not `Mhat_max(x0) = 1`**: the aggregate exceeds its own max by the
conservatism `KS − max` (0.00044 at `ρ = 100`), so a limit at 1 would make the
baseline infeasible at `ε = 0`; at `KS_ρ(x0)` the starting point has exactly
zero slack.

### 10.3 What the load model is, and is not

Out-of-plane moment, not rotated into the section flap axis; no centrifugal
relief, no gravity, no dynamic amplification, no in-plane component. The model
is the same integrand the baseline quotes, so `x0` vs optimised comparisons are
consistent; the absolute number is not a structural load. Because the cap is
relative no material allowable enters (stress `= M/Z` is a constant factor).
The design condition is the rated wind speed at the rotor-speed ceiling — the
highest B3-independent point — and the cut-out case
(20 m/s, λ = 3.14159, 300 rpm, α up to ≈ 30° on Viterna) is **reported as a
labelled post-check, never constrained on**.

### 10.4 Implementation

`src/objective/loads.py` is the forward path (`root_moment`, `root_moments`,
`load_operating_points`, `ks`, `ks_weights`), independent of the adjoint;
`design.baseline.root_bending_moment` is re-exported from it. `src/adjoint/loads.py`
is `RootMomentSystem(BEMSystem)`, which inherits `solve`, `residual`, `partials`,
`dR_dx`, `dR_dd`, the matrix-free operators and `state_sensitivity` unchanged
and adds `moments_from_m`, `moments`, `KS`, `dKS_dx`, `dKS_dd`, `tangent`,
`gradient` (`MomentGradientResult`). The parent's AEP methods (`J`, `weights`,
`J_capped`, `dJ_dx`, `dJ_dd`, `powers_from_q`) are not called: they carry `Ω_b`,
the bin masses and the rated mask. `ScaledProblem.moment_constraint` wraps it
for SLSQP and caches the last `gradient` keyed on the exact `u` bytes, so the
constraint's separate `fun`/`jac` at one `u` cost one 9-point solve
(`n_moment_solves` counts them).

### 10.5 The four tiers, measured

Tier 1 is a complex step (`h = 1e-30`) through the value code, reported as the
worst mixed error `|estimate − partial| / max(1, |partial|)`; Tier 2 is the
forward-mode tangent against the adjoint direction and the assembly identity;
Tier 3 is the constraint Jacobian against central FD at the global
`h* = 3.162277660168379e-06` and `h*/√10`, `h*√10`, with `ε_j` the three-step
local jitter and the round-off floor `δg/h*` measured from nine samples of `g`
along `u + t e_j`, `t = −4e-12 … 4e-12`. All numbers are at `src` commit
`a220e1b`, in `verification/load_constraint/checks.json`; Tier 1/2 also live in
`tests/test_loads.py`.

| point | Tier 1 worst mixed `m` | Tier 1 worst mixed KS | Tier 2 tangent | Tier 2 assembly | Tier 3 worst `|diff|/ε_j` | Taylor min ratio |
|---|---|---|---|---|---|---|
| `x0` | 3.19e-14 (`dm_dφ`) | 8.33e-16 (`dKS_dx`) | 8.9e-16 | 0.0 | 0.396 (`twist_0`) | 99.0 |
| `u*` (unconstrained optimum) | 3.62e-14 (`dm_dφ`) | 1.17e-15 (`dKS_dx`) | 8.8e-16 | 0.0 | 22.64 (`chord_1`), floor-limited | 94.8 |

At `x0` all ten variables are within `3 ε_j` (worst 0.396), and every
disagreement is 0.0–4.8 × the measured round-off floor. At `u*` the single
failure is `chord_1`: `ε_j = 3.42e-13`, two decades below every other variable,
while `|diff| = 7.74e-12` is **0.18 ×** the measured floor `4.41e-11` — the same
fragility §9.5 records for the objective at `x0`, where the flatness estimate
lands on a degenerate neighbour pair. The adjoint is not the suspect; the
acceptance scale is. No tolerance and no estimator was changed; `checks.json`
records `passes: false`, `floor_limited: true` and lists the point and variable
under `tier3_failures` (a failed check stays failed in the JSON; the reading of
it is MJ's, on the BLOCKING list). The production runs start at `x0`, where the
Jacobian passes with worst ratio 0.396.

The aggregate at `x0` is dominated by the rated point: softmax weight **0.957**
at 11 m/s (plus 0.043 at 10.5 m/s; the other seven sum below 1e-37), so at
`ρ = 100` the KS *is* the rated-point moment to within the conservatism. The
conservatism `KS − max` is 0.0115 / 0.00044 / 2.94e-7 at
`ρ = 30 / 100 / 300`.

### 10.6 What Step 2 measured

`verification/load_constraint/` (`result_eps{E}.json`, `pareto.json`). SLSQP
from `x0` at `ε = 0, 0.02, 0.05, 0.10`, exit 0 at every `ε`, constraint active
and `chord_0` on its 0.30 m upper bound at every `ε` (`chord_4` reaches its
0.045 m lower bound from `ε = 0.02`). At `ε = 0` — the Phase 5 production
optimum for now — `AEP = 10.262719 MWh/yr`, `+0.1465 %` over `x0` at a cost of
only **0.00017 %** against the unconstrained optimum: the unconstrained `u*`'s
extra +0.285 % KS bought essentially no energy. That is not an assertion but the
KKT multiplier: the moment row's least-squares multiplier (SLSQP exposes none;
`result_eps{E}.json::kkt`) is `1.67e-3` at `ε = 0`, i.e. a shadow price of
0.017 MWh/yr per unit KS — **0.0017 % AEP per 1 % KS** at the margin — rising
to `1.99e-2` at `ε = 0.02` and `0.129` at `ε = 0.10` (0.129 % per 1 % KS). The
multipliers are positive on every active row at every `ε`, so the points are
KKT points, and integrating the multiplier from 0 to 0.02 predicts a 0.021 %
cost against the measured 0.0206 %. At `ε = 0.10` the cost is **0.5785 %** for
a 10.000 % KS reduction, inside the 0.3–3 % band. The Pareto sweep is cold- and
warm-started at every `ε` and the two agree to the multi-start spread 0.0166;
cost is monotone in `ε`. The moment adjoint costs 0.127 s against the
objective's 0.205 s (9 points against 17).

Two things the `ε = 0` result says that a reader should not have to infer.
First, the cap is on `KS_ρ`, not on the rated-point moment: `KS_ρ(x0)` carries
0.00044 of conservatism from the 10.5 m/s point, and the optimiser spends part
of that budget — its rated moment is **+0.020 %** above `x0`'s while `KS`
equals `KS0` to `5e-8`. Second, for the same reason the Pareto steps are
`ε` in KS; the rated-point moment reductions they buy are 1.997 / 5.037 /
10.089 %.

## 11. A third right-hand side: static tip deflection and its KS aggregate

**Added 2026-09-20 (the mass problem, Step 2;
`docs/PLAN-mass-objective-2026-09-20.md`).** Under a material objective the
optimiser thins the blade; a root-moment cap does not stop that (a thinner
blade carries *less* moment) and a root-stress proxy holds only the root
section. The constraint that fights a thin outer blade is stiffness: for a
thin shell of constant laminate thickness the second moment falls as `c³`
while the load falls roughly as `c`. This section adds the static flapwise
tip deflection as a third functional over the Phase 4 load set `L`, on the
same residual, the same diagonal `∂R/∂φ`, the same `dR/dd`, and the same
station partials — nothing in §1–§10 changes, and the kernel is untouched.

### 11.1 From the station load to the tip deflection

At station `i` of operating point `b`, §10.1's integrand `m` carries the hub
arm `r_i − r_hub`. Dividing by that station constant gives the arm-free
normal load per blade,

    q_{b,i} = m_{b,i} / arm_i,      arm_i = r_i − r_hub  ≥ 0.034 m > 0

and, because `m` and each of `m_φ, m_c, m_θ` carry the same factor, its
partials are `m`'s divided by `arm_i` — exact, no new kernel code. The
flapwise moment at station `k` is the trapezoid rule on `[r_k, R]`:

    M_{b,k} = Σ_i W_{k,i} q_{b,i},     W_{k,i} = t_i (r_i − r_k) for i > k, 0 otherwise

`W` is a constant `25 × 25` matrix: the interior trapezoid weights of the
sub-grid from `r_k` outward are the full-grid `t_i` for `i > k`, and the
`i = k` term has zero arm. With `r_k → r_hub` the same rule is exactly
§10.2's root moment (`tests/test_deflection.py::test_w_with_the_hub_arm_is_the_root_moment`).

The unit-load (Euler–Bernoulli) tip deflection of a cantilever with
`EI(r) = E k_I t_shell c(r)³` — a thin shell of constant laminate thickness
`t_shell`, `k_I` the section's perimeter second moment per unit `c³ t` — is

    δ_b = ∫ M(r) (R − r) / (E I(r)) dr  →  δ_b = Σ_k G_k M_{b,k} / c_k³,   G_k = t_k (R − r_k)

per unit `E k_I t_shell`. The common factor is never applied: the constraint
is relative (`D(u) ≤ D(x0)`), so `E`, `k_I`, `t_shell` cancel exactly as the
material allowable does in §10, and `δ_b` is not a metre. The segment
`[r_hub, r_0]` (34 mm) is outside the station grid, as it is for the root
moment. Transposed onto the loads,

    δ_b = Σ_i A_i(c) q_{b,i},     A_i(c) = Σ_k G_k W_{k,i} / c_k³

`A_i` is a flexibility-weighted arm: the deflection the unit load at station
`i` produces through every inboard section's stiffness. The aggregate and
its normalisation are §10.2's, with the rated-point value at `x0` as the
reference:

    δ_ref = δ_rated(x0) = 30776.559398053847  (per unit E k_I t_shell)
    D     = KS_ρ(δ / δ_ref),   ρ = 100,   D(x0) = 1.0004420795243203

At `x0` the rated point carries softmax weight **0.957** (0.043 at 10.5 m/s),
the same split as the moment, and `δ_b/δ_ref` rises monotonically over `L`
from 0.107 at 3.5 m/s to 1 at 11 m/s.

### 11.2 The partials and the adjoint

With `w_b` the softmax weights,

    ∂D/∂φ_{b,i}      = (w_b / δ_ref) · A_i · m_φ_{b,i} / arm_i
    ∂D/∂c_j|explicit = Σ_b (w_b / δ_ref) [ A_j · m_c_{b,j} / arm_j  −  3 G_j M_{b,j} / c_j⁴ ]
    ∂D/∂θ_j|explicit = Σ_b (w_b / δ_ref) · A_j · m_θ_{b,j} / arm_j

The second term of `∂D/∂c_j` is the stiffness path — the chord at station
`j` enters `1/c_j³` in front of that station's own moment `M_{b,j}`; the
sum over `i` of `q_{b,i} ∂A_i/∂c_j` collapses to it because `Σ_i W_{j,i} q_{b,i}`
is `M_{b,j}`. Everything else is the load path §10 already has. The adjoint
is §8's equation with this right-hand side:

    ψ_{b,i} = −(∂D/∂φ_{b,i}) / (∂R_{b,i}/∂φ_{b,i})                (225 divisions)
    dD/dd   = N_cᵀ ∂D/∂c + N_θᵀ ∂D/∂θ + Σ_{b,i} ψ_{b,i} ∂R_{b,i}/∂d

and the constraint SLSQP sees is `g(u) = D(x0) − D(u) ≥ 0` with Jacobian
`−(dD/dd) ⊙ span`, zero-slack at `x0` exactly as the moment row is. At `x0`
the chord entries of `dD/dd` are negative except at the tip control point
(more chord inboard stiffens the blade more than it loads it; more chord at
the tip loads the whole cantilever) and the twist entries are negative
(more twist, less lift, less deflection) — the signs the physics requires.

### 11.3 Implementation and tiers

`src/adjoint/deflection.py` is `DeflectionSystem(RootMomentSystem)`: it adds
`W`, `G`, `arm`, `loads_from_m`, `spanwise_moments_from_q`,
`flexibility_arms`, `deflections`, `D`, `dD_dx`, `dD_dd`, and overrides
`tangent` and `gradient` (`DeflectionGradientResult`); the parent's moment
functional stays callable on the same instance and the same state, and
`gradient(d, state=…)` accepts the moment row's forward state so the two
rows share one 9-point solve. `src/objective/loads.py` is the forward twin
(`normal_load`, `spanwise_moments`, `tip_deflection`, `tip_deflection_at`)
from the solver's own station records through a plain sub-grid loop.

Measured at `src` commit of Step 2 (`tests/test_deflection.py`, both at
`x0` and at the ±2 % / ±0.5° perturbation of §10.5):

| check | at `x0` |
|---|---|
| forward twin vs system, every point of `L` | 2.2e-16 relative |
| `W` with the hub arm vs `RootMomentSystem.moments_from_m` | 1e-14 relative |
| Tier 1 worst mixed `∂D/∂φ` (225 complex steps of `D`) | 1.22e-15 |
| Tier 1 worst mixed `∂D/∂d` (10 complex steps of `D`) | 8.45e-16 |
| Tier 2 tangent − adjoint, unit `v` | 0.0 |

Tier 3 (the wrapped constraint's Jacobian against central FD at `h*`) is
Step 3's, with the stress row's, in `tests/test_mass_problem.py` and the
production artefact's `checks.json`.

### 11.4 What the deflection model is, and is not

One Euler–Bernoulli cantilever, one section, one shell of constant thickness,
static, out-of-plane loads at the B3-independent points, no centrifugal
stiffening, no gravity, no twist–bend coupling, no shear deformation, no
tower-clearance allowable (the absolute form needs `E`, `k_I`, `t_shell` and
a clearance, none of which the relative row requires). It is a *proxy* whose purpose is to stop a material
minimiser trading stiffness it does not see; the report states it as such
beside the root-stress proxy `KS_ρ(M̂)/c_0²`. Buckling and fatigue are not
modelled, for the reasons `docs/PROPOSAL-mass-objective-2026-09-19.md` §4.6
gives.
