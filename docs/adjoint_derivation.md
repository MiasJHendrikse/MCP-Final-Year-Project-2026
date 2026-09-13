# Discrete adjoint of the Ning-form BEM objective — derivation and implementation plan

**Status: plan only.** Written 2026-09-13 at the end of the Stage A (finite-difference)
run, per `the implementation plan` §6. No code under `src/adjoint/` exists; none of
the B0–B6 checkpoints has been executed. Every partial below is derived by hand from
the code as it stands (`src/bem/station.py`, `src/bem/corrections.py`,
`src/bem/rotor.py`, `src/polars/polar.py`, `src/objective/*`) and is **unverified**
until Tier 1 (§8, B1) says otherwise. A derivation with an unverified partial is
still useful; it is not a result.

All bounds referred to are provisional: `chord_max_m = 0.45 m` is a placeholder
pending the hub-radius / root-attachment decision.

---

## 1. State choice: `φ` only

**State:** `x = φ_{b,i}`, one inflow angle per (operating point `b`, station `i`).
18 operating points (17 bin midpoints `V_b = 3.5 … 19.5 m/s` and the rated solve
`V = 11 m/s`) × 25 stations = **450 scalars**.

**Design:** `d ∈ ℝ¹⁰ = [c₀…c₄ (m), θ₀…θ₄ (rad)]`, the chord and twist control
points. Per-station chord and twist are the *constant* linear maps
`c = N_c d`, `θ = N_θ d` (`BladeParameterisation.dchord_dd()`, `dtwist_dd()`,
each (25, 10)). No function in the adjoint ever takes per-station values as a
design variable.

**Why `a`, `a'` are not states.** `station.residual(phi, station)` evaluates, at a
trial `φ`, the closed-form chain `φ → (Cl, Cd) → (Cn, Ct) → Y → a(Y, F)`. `a` is an
explicit function of `(φ, c, θ)`, not the result of an inner iteration; the same is
true of `a'` (which the residual does not even use — see the `residual` docstring on
the κ' substitution). Choosing `x = φ` therefore makes each `R_{b,i}` a scalar
equation in a scalar unknown with every other quantity explicit, so

    ∂R/∂x  is diagonal (450 × 450), and the adjoint "solve" is 450 scalar divisions.

Had `(a, a')` been taken as states alongside `φ`, the system would be 3 × 3 per
station with two trivially-satisfied rows, and the derivation would have to track
which of the two algebraically equal Buhl expressions the code selected. With `φ`
alone the Buhl branch is handled by implicit differentiation of the defining
equation (§4.4), which is branch-free.

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
`limited_b = (P_b > P_rated)` exactly as `power_per_bin` flags it:

    P̃_b = P_rated  if limited_b  else  P_b
    J    = −(T_h / 10⁶) Σ_b m_b P̃_b                  [MWh/yr]

At `x0`: bins 11.5–19.5 m/s are limited (nine of them); `J = −10.270157`.

---

## 3. The polar partials the chain needs

`CachedPolar` exposes exactly the four first derivatives the chain needs, already
in the residual's units:

    Cl_α = dcl_dalpha(α)   per radian        Cl_Re = dcl_dre(α)   per unit Re
    Cd_α = dcd_dalpha(α)   per radian        Cd_Re = dcd_dre(α)   per unit Re

(`polars.interpolant` is verified to ~1e-14; it is not re-derived here.)

    ∂Cl/∂φ = Cl_α        ∂Cl/∂θ = −Cl_α        ∂Cl/∂c = Cl_Re · W/ν
    ∂Cd/∂φ = Cd_α        ∂Cd/∂θ = −Cd_α        ∂Cd/∂c = Cd_Re · W/ν

The `∂/∂c` column is the **`Re(c)` path** (§4.3 of the brief) — the one most easily
forgotten because `c` also enters through `σ` and the explicit factor in `q`.

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
the factor `1/√(1 − e^{−2f})` is large but finite; complex step handles it, Tier 1
will show whether the mixed tolerance is met there.

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
adjoint's (Tier 4, §8).

The `γ₂ < 0` branch (`gamma1/gamma3`) is unreachable for positive loading and
never at a solved station; the adjoint does not model it. Tier 1 asserts that no
solved station is on it.

### 4.5 Chain into `a`

    a_φ = a_Y Y_φ + a_F F_φ         a_c = a_Y Y_c         a_θ = a_Y Y_θ

**4 of 25 stations are above `a = 0.4` at the rated point at `x0`** (brief §4.3);
the Tier 1 test set must contain them and assert the count `> 0`.

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
`R` changes sign across the root; Tier 1 also reports `min |∂R/∂φ|`).

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
    ω_rated = −(T_h/10⁶) M_L                      (the rated solve enters once per limited bin)

Then, for every operating point `b` (including `rated`) and station `i`:

    ∂J/∂φ_{b,i} = ω_b Ω_b t_i q_{φ,b,i}
    ∂J/∂c_i     = Σ_b ω_b Ω_b t_i q_{c,b,i}          ∂J/∂θ_i = Σ_b ω_b Ω_b t_i q_{θ,b,i}
    ∂J/∂d       = N_cᵀ (∂J/∂c) + N_θᵀ (∂J/∂θ)                                (10,)

and for the residuals

    ∂R_{b,i}/∂d_j = (∂R/∂c)_{b,i} N_c[i, j] + (∂R/∂θ)_{b,i} N_θ[i, j]         (18, 25, 10)

Scaling (`DesignBounds`, `ScaledProblem`): `d = lo + u ⊙ span`, `fun = J/|J₀|`, so

    dfun/du = (dJ/dd ⊙ span) / |J₀|.

`J₀` is the value recorded by `ScaledProblem` at `u₀`; the adjoint Jacobian is
reported in MWh/yr per unit `u` by undoing the `|J₀|`.

`limited` is a *fixed mask* taken from the forward solve; `J` is C¹ but not C²
where a below-rated bin crosses the limit as `d` moves (`objective.power` docstring),
and the adjoint at such a point is the one-sided derivative of whichever side the
forward solve landed on. This is a property of the objective, not of the gradient.

---

## 8. Adjoint solve and assembly

`∂R/∂x` diagonal ⇒ the adjoint equation `(∂R/∂x)ᵀ ψ = −(∂J/∂x)ᵀ` is

    ψ_{b,i} = −(∂J/∂φ_{b,i}) / (∂R_{b,i}/∂φ_{b,i})                     (450 divisions)

    dJ/dd   = ∂J/∂d + Σ_{b,i} ψ_{b,i} ∂R_{b,i}/∂d                         (10,)

Honest cost note: with a diagonal state Jacobian and `n = 10` design variables the
adjoint's advantage over central FD is **accuracy, not solve count**. FD costs
20 rotor-set solves ≈ 4.2 s at `x0` (A4 spent 1,419 objective evaluations, 361 s,
on 34 iterations); the adjoint costs one forward solve plus the partials
(≲ 2 × J ≈ 0.4 s, so the same run would be ≈ 35 × 0.6 s ≈ 20 s) but its real value is a gradient free of the
O(h)·Δf″ contamination FD picks up at the polar C² breaks, and a reference against
which the FD noise floor is measured, not merely a faster one.

---

## 9. Implementation plan (B0–B6), files, tests, artefacts, cost

Cost basis from Stage A (this run, 2026-09-13): `J` = 0.207 s per evaluation;
central FD gradient at `x0` = 20 evaluations ≈ 4.2 s; the 300-evaluation step-size
sweep = 62 s; the FD-driven SLSQP run (A4) = 34 iterations, 1,419 objective
evaluations, 361 s, for a +0.217 % AEP gain; `pytest -q` = 12–18 s (396 tests). Stage A (A1–A4: two source modules,
16 tests, two verification scripts with artefacts and READMEs, this document) cost
roughly the reading of ~2,300 lines plus ~1,300 lines written. The dollar figures
below are rough shares of a similar-sized budget, weighted by how much of each step
is derivation-and-debug (expensive) versus scripting (cheap).

### B0 — complex-step safety (real path bit-identical). ~10 %

Edits, three, in existing files (the *only* edits under `src/` outside `adjoint/`):

1. `src/bem/station.py::_blade_element`: `math.sin/cos(phi)` → complex-safe
   helpers `_sin`, `_cos` modelled on `corrections._sqrt`
   (`cmath` when `isinstance(phi, complex)`).
2. `src/bem/corrections.py::tip_loss_factor`, `hub_loss_factor`: complex-safe
   `sin`, `exp`, `acos`; every branch decision (`abs`, `min(1.0, …)`, the `< 1e-8`
   guard, `r >= R`) taken on `.real`; the real-path lines kept literally intact
   inside the `else`.
3. `src/polars/polar.py::CachedPolar.__init__`: range-check `reynolds.real`, store
   `reynolds` as given (drop the `float()`), guard `__repr__`.

Test: `tests/test_complex_safety.py` — `residual(complex(phi, 1e-30), station)`
returns complex with a finite imaginary part; `CachedPolar(interp, complex(Re, 1e-30))`
constructs; complex `twist` propagates. `pytest` unchanged including
`test_golden_regression.py`, `test_determinism.py`, `test_cost.py` (bit-identical
real path is the exit criterion, and the golden files are the proof).
⏹ `Phase 3: complex-step safety for the residual (real path unchanged)`

### B1 — kernels, system, Tier 1. ~35 % (the derivation lands here)

`src/adjoint/kernels.py::station_partials(phi, chord, twist, reynolds, r, lam_r, R,
B, r_hub, polar)` → `R, R_φ, R_c, R_θ, q, q_φ, q_c, q_θ` (scalars, complex-safe; §3–§6
transcribed; the same `corrections`/`CachedPolar` calls as `station.residual` so the
*value* of `R` matches to 1e-15 — tested). Uses `CachedPolar.dcl_dalpha/dcd_dalpha/
dcl_dre/dcd_dre` for the polar partials; the momentum/Buhl choice is made on
`a_m.real` exactly as the code does, and the Buhl partials are the implicit form
of §4.4.

`src/adjoint/system.py::BEMSystem(parameterisation, bounds, resource, polar_cache)`:
`solve(d)` — one `aerodynamic_power(geometry, V, λ, ρ, ν)` per operating point
(18), reading `phi` from `result["stations"][i]["phi"]`; no second solve, no
re-implemented root-find; `limited = P_b > P_rated`; returns `phi (18,25)`, `limited`,
`P_b`, `P_rated`. `residual(phi, d)`, `dR_dx(phi, d)` (18,25), `dR_dd(phi, d)`
(18,25,10) via `N_c`, `N_θ`; `J(phi, d)`, `dJ_dx(phi, d)` (18,25), `dJ_dd(phi, d)` (10,).
`J(solve(d), d) == annual_energy_mwh(d)` to 1e-12 relative (tested). No NumPy array
pre-allocated with a real dtype on any path a complex step passes through
(`np.result_type` or lists).

Tier 1, `tests/test_adjoint_partials.py`: at `x0` **and** at (chord × 1.2,
twist + 3°), every `(b, i)`: complex step `h = 1e-30` on `residual` w.r.t. `φ` vs
`dR_dx`; w.r.t. each `d_j` (complex `d` → complex `c, θ, Re`) vs `dR_dd`; the same
for `J` vs `dJ_dx`, `dJ_dd`. **Assert `|err| ≤ 1e-13 · max(1, |partial|)`** (1e-12
for `J`). Assert the set contains stations with `a > BUHL_AC` (count > 0) and none
on the `γ₂ < 0` branch. A plateau near 1e-6 means a real dtype, not mathematics.
Failure here is a bug: fix, never loosen.
⏹ `Phase 3: residual system and Tier 1 partials`

### B2 — Tier 2 transpose identity. ~10 %

`tests/test_adjoint_transpose.py`: random `u, v`, `⟨v, A u⟩ = ⟨Aᵀ v, u⟩` to 1e-14
relative for `A = ∂R/∂x` (diagonal), `A = ∂R/∂d` (450 × 10), and the assembled
operator `d ↦ dJ/dd` against its matrix-free adjoint application. "Apply `A`" and
"apply `Aᵀ`" are separate code paths (`system.apply_dR_dd(v)`,
`system.apply_dR_dd_T(w)`), not one matrix and its `.T`. Indexing errors across
(bin, station, design) show here.
⏹ `Phase 3: Tier 2 transpose identity`

### B3 — assembly and Tier 3. ~20 %

`BEMSystem.gradient(d)` per §8; `ScaledProblem.jac_adjoint(u) = gradient(d) ⊙ span / |J₀|`
(a new method on the existing class; no other edit to `src/gradients/`).

`verification/gradient_verification/run_tier3.py` → `tier3.json`, `README.md`: at
`x0`, at the mid-run iterate `k = 17` (`nit/2`, `nit = 34`) from
`verification/fd_optimisation/iterates.json`, and at the FD optimum `u*`; for all
10 variables: adjoint, FD at `h*_j` (`verification/fd_step_size/sweep.json`,
`h_star_per_variable`), abs error, rel error, `h*_j`, `ε_j`. Units MWh/yr per unit
`u`. **Acceptance, derived not invented:** `|adj_j − FD_j| ≤ 3 ε_j` with `ε_j` from
A3, and the gross-error tripwire `rel ≤ 1e-3`. Note that at `x0` `ε_j` is
1e-10–6e-9 MWh/yr per `u` (A3, relative 2e-9–1.6e-6), so the acceptance is tight;
at the two other points `ε_j` must be *re-measured* by a local three-step sweep
(`h*_j/√10, h*_j, h*_j√10`) since the noise floor is a property of the point.
Taylor-remainder test: random unit `v`, `|J(u + εv) − J(u) − ε gᵀv|` for
`ε = 1e-2, 1e-3, 1e-4`, ratio per decade ≥ 30 (asserted). Re-plot A3's V-curves with
`--reference adjoint.json` (`gradient_mwh_per_u`). Log wall time of `J`, FD gradient,
adjoint gradient at `x0` (expect ≈ 0.21 s / ≈ 4.2 s / ≲ 0.4 s).
⏹ `Phase 3: adjoint assembly and Tier 3`

### B4 — Tier 4: degradation at the polar interpolation. ~10 %

`verification/gradient_verification/run_tier4.py` + README section. For the
variable with the worst Tier 3 disagreement: (i) count `(b, i)` whose `α` crosses
a 0.5° knot inside `[u − h*, u + h*]` using `|dα/du_j| h*` (`dα/du_j = (∂φ/∂u_j −
N_θ[i,:] span)`, with `∂φ/∂u_j = −(∂R/∂d · span)_j / (∂R/∂φ)` from the adjoint
system itself) against the distance to the nearest knot; (ii) shrink `h` until no
station crosses an α-knot or a Re row, recompute FD, show the disagreement drops to
round-off; (iii) separately count stations crossing Buhl `a = 0.4` inside the
stencil and report their contribution on its own line. Nothing left unexplained.
A3's finding that the V-curve at `x0` is clean (no plateau) predicts small counts
at `x0`; the mid-run iterate and `u*` are where crossings are more likely.
⏹ `Phase 3: Tier 4 degradation at the polar interpolation`

### B5 — adjoint-driven SLSQP. ~10 %

`verification/adjoint_optimisation/run_adjoint_slsqp.py`: A4's script with
`jac=problem.jac_adjoint`, same options, same artefacts (`result.json`,
`iterates.json`, `optimised_blade.png`, `README.md`). Compare with A4
(`u*_fd`, `ΔAEP = +0.217 %`, `nit = 34`, `nfev = 37`, 361 s):
`‖u*_adj − u*_fd‖∞`, `ΔAEP`, `nit`, `nfev`, wall time. They must agree to SLSQP's
tolerance; a different optimum is a suspected silent gradient error to be
investigated before reporting.
⏹ `Phase 3: adjoint-driven SLSQP`

### B6 — this document, completed. ~5 %

Fill in the four tiers' measured numbers, the Tier 4 attribution sentence, the
cost table (J / FD / adjoint wall time; nfev of A4 vs B5), and the honest cost
note with numbers.
⏹ `Phase 3: derivation`

### Failure modes the plan guards against

- **Silent wrong gradient** (R2): Tier 1/2 at machine precision as hard asserts;
  Tier 3 against `ε_j`, never looser; Taylor remainder; FD and adjoint optima
  coincide. No tolerance is adjusted to pass; a failure is reported.
- **Forgetting `Re(c)`**: §3's `∂/∂c` column and `q_c`'s third term. Tier 1's
  complex `d → c → Re` step catches it.
- **Treating `a` as a state**: §1.
- **Real-dtype arrays on a complex path**: `np.result_type` / lists; the 1e-6
  plateau signature in Tier 1.
- **Per-station DOF**: `N_c`, `N_θ` are the only route from stations to `d`.
- **Placeholders**: bounds provisional in every sentence; nothing in `config/`.
