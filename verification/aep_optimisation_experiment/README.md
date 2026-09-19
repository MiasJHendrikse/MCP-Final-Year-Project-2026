# EXPERIMENT — direct AEP optimisation under a rotor-speed ceiling, with the discrete adjoint

> **Port note (2026-09-19).** This directory is a **frozen record, not a live
> deliverable.** It was written and run in a separate testing copy of the
> repository against commit `4c6feea` — i.e. *before* the configured bounds
> (`chord_max_m = 0.30 m`), *before* the 300 rpm operating law entered
> `config/rotor_design.yaml`, and *before* `PROVISIONAL_BOUNDS` was retired on
> 2026-09-19. It imports that retired `PROVISIONAL_BOUNDS` and treats the
> pre-ceiling `x*` as "the existing `x*`". **It is not expected to run in this
> repository as it stands, and it must not be re-run to refresh anything.** It
> is kept because it is the evidence behind decision **B1** (the rotor-speed
> ceiling is a machine fact and belongs in the config) and **O4** (the
> configured chord bound), both of which were adopted into
> `config/rotor_design.yaml` and then re-measured by the artefacts in
> `verification/fd_optimisation/`, `verification/adjoint_optimisation/` and
> `verification/fd_optimisation_multistart/`. Read §B/§E for the argument and
> `results/summary.json` for the numbers as they stood. Everywhere below, "the
> project" means the repository at `4c6feea`.

**Date:** 2026-09-19. **Status: experiment, not a project result.** Nothing in
`src/`, `config/`, any existing `verification/` artefact, the title, the
objective, the baseline, `x0`, `x*` or any documentation was changed. This
directory is self-contained and runs independently of the project workflow.

**Question.** If the same ten design variables are optimised directly for
Khomas Hochland AEP under a physically motivated variable-speed law with a
maximum rotor speed, instead of for the project's current objective (which
reduces to Cp at λ = 6.5 — `docs/AEP_GAIN_AUDIT.md` §1.1), does the optimiser
find a meaningfully better blade, and can the existing discrete adjoint
supply the gradient?

**Answer in four lines.**
1. **Technically: yes.** The adjoint extends to a per-bin λ with a 60-line
   subclass and no change to the derivation; it matches central FD to
   1e-8–1e-10 relative at `x0` and to the FD noise floor (|Δ| ≈ 1e-10 in
   scaled units) at every optimum. SLSQP converged (exit 0) in all four
   cases, in 27–49 iterations and 18–41 s each — 8–10× faster than the
   audit's FD-driven runs, to the same optima (‖Δu‖∞ ≤ 1.8e-6).
2. **Control:** with no ceiling the experiment reproduces the project's `x*`
   to 2e-15 in `u` — the objective *is* the project's, and its gain is
   +0.121 %.
3. **With a ceiling the gain over Schmitz is +0.34 % (60 m/s), +2.33 %
   (55 m/s), +7.84 % (50 m/s)** — the audit's FD numbers, now reached by the
   adjoint. The blades are wider (root +32/+60/+36 %, mid-span
   +66/+169/+107 %), more twisted (+1–7°), have a flatter Cp–λ curve, and
   cost 0.7–4.3 % more thrust at 8.5 m/s.
4. **The gain is real within the model but its size is set by the ceiling,
   which is an unsupplied machine fact (audit §5, B1), and at 50 m/s a
   large part of it is Schmitz being scored in post-stall (Viterna) territory.**
   The experiment produces evidence; it does not by itself justify changing
   the project's direction. See §E.

---

## A. Configuration

Everything below was read from the project's implementation at run time
(`config/site.yaml`, `config/rotor_design.yaml`, `src/objective`,
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS`), not assumed. The
brief's "approximately" figures were all confirmed.

| item | value | source |
|---|---|---|
| Wind resource | Weibull `k = 1.709226`, `c = 7.273759 m/s` at 20 m hub height (`V̄ = 6.49 m/s`) | `site.yaml`, `WeibullResource.from_config()` |
| Bins | 17 bins, 1 m/s, edges 3…20 m/s, midpoint rule, masses from the CDF (`Σ m_b = 0.79889 = F(20) − F(3)`) | `objective.power.wind_speed_bins`, `WeibullResource.probability_between` |
| Air | `ρ = 0.96839 kg/m³`, `ν = 1.8690e-5 m²/s` | `site.yaml` |
| Rotor | `R = 2.0 m`, `B = 3`, root cut-out 0.15 R, SG6043 polar cache, 25 strips | `rotor_design.yaml` |
| Annual hours | `T = 8766 h` | `objective.objective.HOURS_PER_YEAR` |
| Rating | `P_rated = 3822.189755 W`, fixed (a nameplate, not a function of the blade), `P = min(P_aero, P_rated)` | `rotor_design.yaml operating.rated_power_w` |
| BEM | `bem.rotor.solve_rotor`, unchanged; same convergence settings | `objective.power.aerodynamic_power` |
| Operating law | `λ(V) = min(6.5, V_tip,max / V)` — below `V_c = V_tip,max / 6.5` the project's λ = 6.5; above it `Ω = Ω_max` and λ falls as 1/V | `ceiling_objective.tsr_schedule` |
| Ceiling cases | `none` (control = the project's objective), **60 m/s (286 rpm, V_c = 9.2 m/s)**, **55 m/s (263 rpm, V_c = 8.5 m/s)**, **50 m/s (239 rpm, V_c = 7.7 m/s)** — the three cases the audit already tested (§3.1); none chosen here | `docs/AEP_GAIN_AUDIT.md` §3.1 |
| Design variables | the project's 10: 5 chord + 5 twist B-spline control points; no additions | `design.BladeParameterisation` |
| Bounds | provisional, as every project optimisation: chord 0.045–0.45 m, twist −2°…35°; scaled to `u ∈ [0,1]¹⁰` | `PROVISIONAL_BOUNDS` |
| Constraint | the polar-cache Reynolds envelope (linear, 50 rows, margin 5 %). Reused as-is: it assumes λ = 6.5 at every speed, which over-estimates `W` under a ceiling — conservative, the same choice `reoptimise.py` made. Never active (§B) | `ScaledProblem.envelope_constraint` |
| Optimiser | SLSQP, `ftol = 1e-8`, `maxiter = 200`, objective `fun(u) = J(u)/|J(u0)|`, `J = −AEP`, start `x0` (fitted Schmitz) | B5's set-up, `run_adjoint_slsqp.py` |
| Gradient | **discrete adjoint** — `CeilingBEMSystem` (below). FD is used only as a check | `ceiling_objective.CeilingProblem.jac_adjoint` |

### How the adjoint was extended (`ceiling_objective.py`)

`adjoint.system.BEMSystem` carries one `Ω_b = λ V_b / R` per operating point
in `self.omega`; the local speed ratio, the Reynolds estimate
`W = hypot(V_b, Ω_b r_i)` and the power assembly `P_b = Ω_b Σ_i t_i q_{b,i}`
all read that list. `CeilingBEMSystem` sets `Ω_b = λ_b V_b / R` with the
per-bin λ and re-states the forward solve (the parent passes the scalar
`self.tsr`). Nothing else changes: `∂R/∂x` stays diagonal, `∂R/∂d` keeps
its form, and the parent's `gradient()` already assembles

    dJ/dd = Σ_b ω_b dP_b/dd,   ω_b = −(T/1e6) m_b  (uncapped),  0  (capped)

which is the weighted sum of per-bin adjoint sensitivities the brief asks
for. **No limitation was found**; the audit's prediction ("a small change
to `power_per_bin` and `BEMSystem`, not a re-derivation") held.
`vtip_max = None` reproduces the parent bit-for-bit (checked: `J` and the
gradient at `x0` differ by exactly 0).

---

## B. Optimisation convergence

| case | AEP(x0) | AEP(x_AEP*) | nit | nfev | njev | status | ‖∇J‖ at optimum [MWh/yr per u] | envelope violation | active bounds | wall |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|---:|
| none (control) | 10.270157 | 10.282564 | 27 | 29 | 27 | 0, success | 3.6e-4 (free vars) | 0 | `twist_4` lower (−2°), as `x*` | 18 s |
| 60 m/s | 10.210635 | 10.245066 | 44 | 45 | 44 | 0, success | 1.3e-4 | 0 | none | 32 s |
| 55 m/s | 9.976818 | 10.208978 | 49 | 52 | 49 | 0, success | 3.6e-4 | 0 | none (root chord `u = 0.979`, 441 mm of 450) | 41 s |
| 50 m/s | 9.429005 | 10.168631 | 46 | 48 | 46 | 0, success | 1.7e-4 | 0 | none | 40 s |

`nfev`/`njev` are SciPy's counts (the iterate recorder adds one objective
and one adjoint per accepted step on top, as in B5). No `PolarDomainError`
in any case; the envelope's tightest row had 15–33 mm of slack. The
`u`-gradient at every optimum is ~1e-4 MWh/yr per unit of `u`, i.e. flat,
as the project's own optimum is. For comparison the audit's FD-driven runs
of the same problems took 570–1033 objective evaluations and 179–389 s.

---

## C. Main results

AEP in MWh/yr, all under the fixed 3822.19 W rating. "Cp @ λ = 6.5" is at
8.5 m/s (the audit's spanwise reference speed; the Cp–λ curves in
`results/cp_lambda.png` are at the same speed). `Ē[Cp]` is the mean Cp on
the case's schedule, weighted by each bin's unlimited energy `m_b V_b³` —
the Cp the site actually sees.

### No ceiling (control — the project's own objective)

| design | objective | AEP | gain vs Schmitz | gain vs x* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x0 | — | 10.270157 | — | | 0.4678 | 0.4707 | 0.8087 |
| Existing x* | Cp @ 6.5 (= project AEP) | 10.282564 | +0.121 % | — | 0.4687 | 0.4714 | 0.7938 |
| AEP optimum | AEP, no ceiling | 10.282564 | +0.121 % | +0.000 % | 0.4687 | 0.4714 | 0.7938 |

The experiment's optimum **is** `x*`: ‖u_AEP − u*‖∞ = 2e-15. The project's
objective and "AEP without a ceiling" are the same problem (audit §1.1),
and the pipeline reproduces the committed result exactly.

### V_tip,max = 60 m/s (286 rpm)

| design | objective | AEP | gain vs Schmitz | gain vs x* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x0 | — | 10.210635 | — | | 0.4678 | 0.3718 | 0.8087 |
| Existing x* | Cp @ 6.5 | 10.213119 | +0.024 % | — | 0.4687 | 0.3690 | 0.7938 |
| **AEP optimum** | AEP @ 60 m/s | **10.245066** | **+0.337 %** | +0.313 % | 0.4677 | 0.3961 | 0.8141 |

### V_tip,max = 55 m/s (263 rpm)

| design | objective | AEP | gain vs Schmitz | gain vs x* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x0 | — | 9.976818 | — | | 0.4678 | 0.3402 | 0.8087 |
| Existing x* | Cp @ 6.5 | 9.948790 | −0.281 % | — | 0.4687 | 0.3368 | 0.7938 |
| **AEP optimum** | AEP @ 55 m/s | **10.208978** | **+2.327 %** | +2.615 % | 0.4651 | 0.4040 | 0.8285 |

### V_tip,max = 50 m/s (239 rpm)

| design | objective | AEP | gain vs Schmitz | gain vs x* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x0 | — | 9.429005 | — | | 0.4678 | 0.3035 | 0.8087 |
| Existing x* | Cp @ 6.5 | 9.306357 | −1.301 % | — | 0.4687 | 0.2995 | 0.7938 |
| **AEP optimum** | AEP @ 50 m/s | **10.168631** | **+7.844 %** | +9.265 % | 0.4604 | 0.4013 | 0.8433 |

### Where the gain comes from (bin by bin, `gain_attribution_by_bin` in each case file)

The AEP optimum gives up 0.1–3.4 % of Cp in the bins that still run at
λ = 6.5 (3.5–7.5 m/s; together −23 / −13 / −6 % of the net gain at
60 / 55 / 50 m/s) and wins it back many times over in the bins between
`V_c` and rated where λ has dropped: at 55 m/s the 9.5, 10.5 and 11.5 m/s
bins contribute 113 % of the net gain (Cp 0.455→0.469, 0.419→0.462,
0.379→0.446). It is a multi-point design trading the design point for the
λ band; that is the physics the audit's §3.1 described.

### Robustness — every optimum under every schedule (`results/cross_evaluation.json`, % vs x0 under the same schedule)

| schedule ↓ / blade → | x0 | x* | x_AEP(60) | x_AEP(55) | x_AEP(50) |
|---|---:|---:|---:|---:|---:|
| none | 0 | +0.12 % | −0.07 % | −0.43 % | −0.98 % |
| 60 m/s | 0 | +0.02 % | **+0.34 %** | +0.20 % | −0.17 % |
| 55 m/s | 0 | −0.28 % | +1.88 % | **+2.33 %** | +2.15 % |
| 50 m/s | 0 | −1.30 % | +4.73 % | +7.57 % | **+7.84 %** |

A ceiling-designed blade loses little when the ceiling is wrong (the 55 m/s
blade is −0.43 % with no ceiling and +2.15 % at 50 m/s); the Schmitz and
`x*` blades lose steadily as the ceiling drops. `x*` is worse than Schmitz
under any ceiling ≤ 55 m/s: its Cp-at-one-point optimisation sharpened the
peak.

---

## D. Geometry comparison

Control points: chord in mm, twist in degrees (Δ chord in %, Δ twist in
degrees — twist percentages near 0° are meaningless). `results/geometry.png`
has the spanwise distributions.

| variable | x0 | x* | x_AEP(60) | Δ vs x0 | x_AEP(55) | Δ vs x0 | x_AEP(50) | Δ vs x0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| chord_0 (root) | 275.8 | 285.3 | 364.0 | +32.0 % | 441.4 | +60.0 % | 375.8 | +36.3 % |
| chord_1 | 181.6 | 173.4 | 179.6 | −1.1 % | 185.5 | +2.1 % | 293.5 | +61.6 % |
| chord_2 (mid) | 97.1 | 107.3 | 160.8 | +65.5 % | 261.3 | +169.0 % | 201.0 | +107.0 % |
| chord_3 | 78.6 | 77.9 | 73.7 | −6.3 % | 72.4 | −7.9 % | 127.9 | +62.6 % |
| chord_4 (tip) | 67.0 | 45.8 | 54.7 | −18.4 % | 66.2 | −1.3 % | 98.7 | +47.2 % |
| twist_0 (root) | 24.53 | 23.50 | 27.02 | +2.49° | 28.82 | +4.29° | 25.62 | +1.09° |
| twist_1 | 9.73 | 10.36 | 10.56 | +0.84° | 11.63 | +1.90° | 16.31 | +6.58° |
| twist_2 (mid) | 4.04 | 3.52 | 7.41 | +3.37° | 10.99 | +6.95° | 7.39 | +3.35° |
| twist_3 | 1.00 | 2.01 | 2.41 | +1.41° | 4.44 | +3.44° | 7.18 | +6.18° |
| twist_4 (tip) | 0.59 | −2.00 | −1.71 | −2.30° | −1.65 | −2.23° | 2.16 | +1.58° |
| ‖Δu‖∞ vs x0 | | 0.070 | 0.218 | | 0.409 | | 0.276 | |
| ‖Δu‖∞ vs x* | | — | 0.194 | | 0.385 | | 0.296 | |

For scale: `x*` moved 0.070 in `u` from Schmitz; the ceiling optima move
3–6× further, in a different direction (more solidity and more twist
inboard and mid-span, not the tip unloading `x*` did). The 55 m/s blade's
mid-span control point at 261 mm produces a visible bulge at r ≈ 1.0–1.3 m
(`geometry.png`); with five control points that is how the spline buys
solidity where the λ-band bins want it. Whether such a planform is
manufacturable is exactly the `chord_max_m` / solidity question
(OUTSTANDING-INPUTS §2) — the 55 m/s root is 9 mm inside the provisional
450 mm bound.

---

## E. Verification of the result

Every item below is in `results/case_<label>.json` under `checks`.

| check | none | 60 | 55 | 50 |
|---|---|---|---|---|
| 1. AEP at x_AEP* re-evaluated with a fresh parameterisation + resource | Δ = 0 | Δ = 0 | Δ = 0 | Δ = 0 |
| 1. AEP re-assembled from the station integrand `q` by the adjoint system (second code path) | Δ = −1.8e-15 | −3.6e-15 | 0 | −1.8e-15 |
| 1. AEP through `src/objective.annual_energy_mwh` (control only) | Δ = 0 | n/a | n/a | n/a |
| 2. adjoint vs central FD (h = 3e-6) at **x0**: max |Δ| / ‖∇‖ | 6.9e-10 | 5.5e-10 | 1.7e-10 | 2.3e-8 |
| 2. adjoint vs central FD (h = 3e-6) at **x_AEP\***: max |Δ| (scaled) / max |Δ| / ‖∇‖ | 1.2e-10 / 1.6e-7 | 7.9e-11 / 6.2e-6 | 6.6e-11 / 1.9e-6 | 4.1e-11 / 2.3e-6 |
| 2. the same at h = 1e-5 | 1.0e-9 | 7.2e-10 | 2.4e-10 | 2.1e-10 |
| 3. SLSQP status | 0 | 0 | 0 | 0 |
| 4. all `u ∈ [0,1]` | yes | yes | yes | yes |
| 5. BEM converged in all 17 bins at x_AEP* | yes | yes | yes | yes |
| 5. α on **uncapped** bins inside the XFOIL band (−8…18°) at x_AEP* | 0.9…6.5° | −0.8…6.7° | −2.0…6.5° | −0.9…6.2° |
| 6. bin masses = CDF(20) − CDF(3); adjoint weights = −(T/1e6)·m_b; Σ bin energies = total | yes / yes / Δ = 0 | same | same | same |
| 7. vs the audit's **FD-driven** optimum for the same ceiling: ‖Δu‖∞ / ΔAEP | 6e-10 / 2e-13 | 1.3e-8 / 5e-13 | 1.8e-6 / 3e-11 | 1.8e-8 / 5e-13 |

On item 2: at the optima the gradient is ~1e-4 in scaled units, so a
relative error of 1e-6 is an absolute error of 1e-10 — the central-FD
noise floor for this objective (it grows 10× from h = 3e-6 to 1e-5, the
signature of FD truncation error, not of an adjoint defect). At `x0`,
where the gradient is O(1), agreement is 1e-8 to 1e-10.

On item 7: two independent gradients (the audit's central FD, this
adjoint) driving the same optimiser on the same objective land on the same
point to within SLSQP's tolerance in all four cases. Together with item 2
that rules out a silent gradient error, and together with the control
reproducing `x*` it rules out a weighting or units error in the experiment's
objective.

**Nothing numerical is responsible for the improvement.** What *is*
responsible, and has to be said plainly:

- **The capped bins carry α up to 28–37° on every blade** (λ ≈ 2.5–3 at
  16–20 m/s under a ceiling): deep stall on the Viterna extrapolation. On
  the *optimised* blades those bins are all capped, so they add a constant
  to J, zero to the gradient, and cannot steer the optimiser. That is not
  true of the baselines:
  - at **50 m/s, Schmitz x0 does not reach rated power until 16.5 m/s**
    and is scored on 54 station-bins with α > 18° (up to 29.7°) in the
    bins 11.5–15.5 m/s; those bins carry **55 % of the +7.84 %**. `x*`
    never reaches rated at all under that ceiling (143 station-bins in
    post-stall, up to 36.7°), which is why "vs x*" reads +9.3 %.
  - at **55 m/s**, x0 is uncapped through 11.5 m/s with α ≤ 14.1° (inside
    the band); the 11.5 m/s bin carries 43 % of the +2.33 %.
  - at **60 m/s** every uncapped bin of every blade is inside the band.
  The 60 and 55 m/s gains are therefore within the polar cache's validated
  range; the 50 m/s gain is partly a measure of how badly a Schmitz blade
  stalls when forced to λ = 4 at 12.5 m/s, computed with post-stall polars
  BEM is least reliable for (audit §3.3).
- The envelope constraint is conservative under a ceiling (it assumes
  λ = 6.5); it was never active, so it did not shape any optimum.
- Bounds are provisional; the 55 m/s root chord is 9 mm from the box.

---

## F. Interpretation — the brief's eight questions

1. **Does direct AEP optimisation work with the current discrete-adjoint BEM framework?**
   Yes. A per-bin λ is a subclass of `BEMSystem` that sets `Ω_b` and
   re-states the forward solve; the derivation, kernels, transpose
   operators and assembly are untouched. FD-verified at every point
   checked. No limitation was encountered.

2. **Does it converge reliably?** Yes: exit 0 in all four cases, 27–49
   iterations, to the same optima the FD-driven runs found, at 1/8–1/10 of
   the wall time. The Cp/AEP objective is flat near the optimum in `u`
   (gradient ~1e-4), as the project's own is.

3. **How much AEP improvement over Schmitz?** +0.121 % (no ceiling — the
   existing result), **+0.34 % at 60 m/s, +2.33 % at 55 m/s, +7.84 % at
   50 m/s**, each over Schmitz evaluated under the *same* ceiling. Under
   no ceiling the improvement is not "meaningfully larger than 0.12 %" —
   it *is* 0.12 %.

4. **Over the existing x\*?** +0.31 % / +2.62 % / +9.27 % at 60/55/50 m/s.
   `x*` is *worse than Schmitz* under a 55 or 50 m/s ceiling (−0.28 %,
   −1.30 %): it sharpened the Cp peak at λ = 6.5 and pays for it off-peak.

5. **Is the blade meaningfully different from Schmitz?** Under a ceiling,
   yes: root chord +32–60 %, mid-span +66–169 %, +1–7° of twist, ‖Δu‖∞
   3–6× the move `x*` made, a visibly different planform and a flatter
   Cp–λ curve (Cp at λ = 4.5 and 8.5 m/s: 0.350 → 0.380 / 0.428 / 0.447).
   Thrust at 8.5 m/s rises 0.7 / 2.4 / 4.3 % where `x*` had cut it 1.8 %.
   Without a ceiling, no: it is `x*`.

6. **Does changing the objective from Cp @ 6.5 to AEP materially change the optimum?**
   Only if the operating law gives AEP a λ dimension. "AEP at fixed
   λ = 6.5" and "Cp @ 6.5" have the same optimum (the control proves it to
   machine precision). With a ceiling the optimum moves materially, and
   the further the ceiling bites below rated, the further it moves.

7. **Does the maximum-RPM constraint materially affect the optimum?** It is
   the *only* thing that does. 60 m/s (bites at 9.2 m/s): sub-1 %. 55 m/s
   (8.5 m/s): 2.3 %. 50 m/s (7.7 m/s): 7.8 %, of which roughly half is the
   baseline stalling in the Viterna region. The gain is a monotone
   function of a parameter this project has not been given.

8. **Is this strong enough to justify changing the project's direction?**
   The experiment establishes two things firmly: the framework can
   optimise site AEP under a variable-λ law with the adjoint, cheaply and
   verifiably; and a rotor-speed ceiling at or below ~55 m/s tip speed
   turns the +0.12 % problem into a multi-% one with a genuinely different
   blade. What it cannot establish is *which ceiling the machine has* —
   and the audit's §3.5 warning stands: picking 55 or 50 m/s because the
   number is bigger is the gimmick. The honest position is: **the
   direction is justified if and only if B1 (maximum rotor speed) is
   answered with a real ceiling** (generator rated rpm, a noise limit with
   its source), and `chord_max_m` / a solidity or load bound comes with it,
   because the ceiling-designed blades are wide and the 50 m/s case shows
   that without a fidelity guard the optimiser will happily be credited
   for the baseline's stall. This README makes no recommendation to change
   the project; it supplies the evidence for that decision.

---

## Files

| file | what |
|---|---|
| `ceiling_objective.py` | `tsr_schedule`, `CeilingBEMSystem` (the adjoint with per-bin λ), `evaluate_aep` (the forward AEP with per-bin record), `CeilingProblem` (the scaled SLSQP problem) |
| `run_experiment.py` | the four cases end to end: optimise, verify, cross-evaluate, plot |
| `results/case_<none,60,55,50>.json` | per case: configuration, convergence, AEP table, every check in §E, gradient checks at `x0` and the optimum, Cp metrics and Cp–λ sweeps, geometry deltas, per-bin rows, bin-by-bin gain attribution, every iterate |
| `results/summary.json` | the headline numbers of the four cases |
| `results/cross_evaluation.json` | every blade under every schedule |
| `results/geometry.png` | chord and twist: x0, x*, x_AEP* per case |
| `results/cp_lambda.png` | Cp–λ at 8.5 m/s of the three blades per case; shaded: the λ band the schedule sweeps below rated |

## Reproduce

    python verification/aep_optimisation_experiment/run_experiment.py            # all four cases, ~2.5 min
    python verification/aep_optimisation_experiment/run_experiment.py --cases 55  # one case

Deterministic: two runs gave identical JSON. Requires nothing beyond the
project's `requirements.txt`; reads `verification/baseline/x0.json`,
`verification/adjoint_optimisation/result.json` (for `x*`) and, for check 7,
`verification/aep_gain_audit/opt_<v>_fixed.json`. Writes only into
`results/`.
