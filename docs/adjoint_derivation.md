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

All bounds referred to are provisional: `chord_max_m = 0.45 m` is a placeholder
pending the hub-radius / root-attachment decision. Nothing here reads a bound from
`config/`.

**Revision 2026-09-19 — fixed generator rating.** The objective's cap above rated
was `P_aero(V_rated; d)`, floating with the design; `docs/AEP_GAIN_AUDIT.md` §3.2
found that indefensible (a nameplate does not grow because the blade improved) and
it is now the configured constant `operating.rated_power_w` (provisionally the
baseline's own `P_aero(11 m/s; x0)`, pending input B2). For the adjoint this
removes the rated solve as an operating point — the system is **17 × 25 = 425**
states, not 450 — and turns the capped bins into a constant term in `J` with no
weight on any state (§7). §1 and §7 below are restated in the current form; the
measured numbers in §9 are the 2026-09-13 record of the 450-state system and are
left as written, with the re-run Tier 3/4 numbers in
`verification/gradient_verification/README.md`. Tier 1 and Tier 2 pass unchanged
on the 425-state system (`pytest tests/test_adjoint_*.py`).

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
SLSQP run (A4) = 34 iterations, 1,419 objective evaluations, 361 s, for a +0.217 %
AEP gain. Stage B, this run: `pytest -q` went from 391 to 450 passed (5 xfailed
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
| AEP* | 10.292482 MWh/yr (+0.217 %) | 10.292482 MWh/yr (+0.217 %) |
| `‖u*_adj − u*_fd‖∞` | 8.2e-9 | tolerance 8.6e-3 (10 × the multi-start spread) |
| `AEP*_adj − AEP*_fd` | +1.4e-13 MWh/yr | tolerance ≈ 1e-6 (from `ftol`) |
| trajectory: `max_k ‖u_k − u_k^fd‖∞` | 6.8e-7 | — |

The optima coincide, and so do the 34 iterates on the way: the Tier 3 agreement
carried through every quasi-Newton update without the paths separating. No active
bound, no active envelope row, no `PolarDomainError`, α and Re inside the cache.
The gain is 0.22 % as in A4 (below the brief's 2–6 % expectation, above neither
defect gate; flagged for MJ there and confirmed global by the multi-start study).

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
- **Placeholders**: bounds provisional in every artefact; nothing in `config/`;
  no cache extension was needed.
- **Loosening a tier**: none loosened. Tier 3's residual disagreement was
  unexplained until Tier 4 explained it, and Tier 4 explained it as round-off with
  the polar interpolation contributing zero at `h*_j`.
