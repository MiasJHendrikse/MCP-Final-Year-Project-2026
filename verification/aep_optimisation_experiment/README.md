# Experiment: maximising AEP under a rotor-speed ceiling

**This is a record of an earlier experiment, not a current result.** It was
run on an earlier version of the code, before the 300 rpm ceiling and the
0.30 m chord bound were added to `config/rotor_design.yaml`, and it imports a
set of design bounds that has since been removed. It won't run on the
repository as it is now, and it shouldn't be re-run to refresh anything. I've
kept it because it is the evidence behind two design choices: treating the
rotor-speed ceiling as a property of the machine (set in the config), and the
0.30 m maximum chord. Both were then re-measured by the current results in
`verification/fd_optimisation/`, `verification/adjoint_optimisation/` and
`verification/fd_optimisation_multistart/`. In this README, "the project"
means the code as it was then.

One later correction applies. The large gains below compare a
ceiling-optimised blade with a Schmitz blade designed for λ = 6.5 without a
ceiling. When the Schmitz blade is re-tuned for the same ceiling, the
optimiser is at most about 0.3 % ahead at every ceiling. That check was
exploratory and isn't committed here, but it is the reason the project's
final objective is material rather than energy (`docs/DESIGN-BASIS.md` §2).

**The question.** If the ten design variables are optimised directly for AEP
at the site under a variable-speed law with a maximum rotor speed, instead of
the original objective (which reduces to C_P at λ = 6.5), does the optimiser
find a much better blade? And can the existing adjoint provide the gradient?

**The answer.**

1. **Yes, technically.** The adjoint extends to a per-bin λ with a 60-line
   subclass and no change to the derivation. It matches central finite
   differences to 1e-8–1e-10 relative at x₀, and to the finite-difference
   noise floor (about 1e-10 in scaled units) at every optimum. SLSQP
   converged in all four cases, in 27–49 iterations and 18–41 s each. That is
   8–10 times faster than the earlier finite-difference runs, and it reaches
   the same optima (‖Δu‖∞ ≤ 1.8e-6).
2. **The control case checks out.** With no ceiling, the experiment
   reproduces the project's x\* to 2e-15 in u, with the same +0.121 % gain.
3. **With a ceiling, the gain over Schmitz is +0.34 % (60 m/s tip speed),
   +2.33 % (55 m/s) and +7.84 % (50 m/s).** The optimised blades are wider
   (root +32/+60/+36 %, mid-span +66/+169/+107 %) and more twisted (1–7°
   more). They have a flatter Cp–λ curve (Cp at λ = 4.5 and 8.5 m/s rises from
   0.350 to 0.380 / 0.428 / 0.447) and 0.7–4.3 % more thrust at 8.5 m/s,
   where x\* had cut thrust by 1.8 %.
4. **The gain is real within the model, but its size is set by the ceiling,
   and the ceiling is a machine property I didn't have.** At 50 m/s, much of
   the gain comes from the Schmitz blade being scored deep in stall (on the
   Viterna extrapolation). So the experiment provides evidence, but on its own
   it doesn't justify changing the project's direction.

## Configuration

Everything was read from the project's implementation at run time:

| item | value | source |
|---|---|---|
| Wind resource | Weibull `k = 1.709226`, `c = 7.273759 m/s` at 20 m (`V̄ = 6.49 m/s`) | `site.yaml`, `WeibullResource.from_config()` |
| Bins | 17 bins of 1 m/s from 3 to 20 m/s, midpoint rule, weights from the CDF (`Σ m_b = 0.79889 = F(20) − F(3)`) | `objective.power.wind_speed_bins`, `WeibullResource.probability_between` |
| Air | `ρ = 0.96839 kg/m³`, `ν = 1.8690e-5 m²/s` | `site.yaml` |
| Rotor | `R = 2.0 m`, `B = 3`, root cut-out 0.15 R, SG6043 polar table, 25 strips | `rotor_design.yaml` |
| Hours per year | `T = 8766 h` | `objective.objective.HOURS_PER_YEAR` |
| Rating | `P_rated = 3822.189755 W`, fixed, with `P = min(P_aero, P_rated)` | `rotor_design.yaml`, `operating.rated_power_w` |
| BEM | `bem.rotor.solve_rotor`, unchanged, same convergence settings | `objective.power.aerodynamic_power` |
| Operating law | `λ(V) = min(6.5, V_tip,max / V)`: λ = 6.5 below `V_c = V_tip,max / 6.5`, and constant rotor speed above it | `ceiling_objective.tsr_schedule` |
| Ceilings | none (the control), **60 m/s (286 rpm, V_c = 9.2 m/s)**, **55 m/s (263 rpm, V_c = 8.5 m/s)**, **50 m/s (239 rpm, V_c = 7.7 m/s)** | |
| Design variables | the project's 10: 5 chord and 5 twist control points | `design.BladeParameterisation` |
| Bounds | the bounds in use at the time: chord 0.045–0.45 m, twist −2° to 35°, scaled to `u ∈ [0,1]¹⁰` | the since-removed `PROVISIONAL_BOUNDS` |
| Constraint | the polar-table Reynolds envelope (linear, 50 rows, 5 % margin). It assumes λ = 6.5 at every speed, which overestimates `W` under a ceiling, so it's conservative. It was never active | `ScaledProblem.envelope_constraint` |
| Optimiser | SLSQP, `ftol = 1e-8`, `maxiter = 200`, objective `fun(u) = J(u)/|J(u0)|` with `J = −AEP`, starting from x₀ | as in `run_adjoint_slsqp.py` |
| Gradient | the discrete adjoint, `CeilingBEMSystem` (below); finite differences only as a check | `ceiling_objective.CeilingProblem.jac_adjoint` |

### How the adjoint was extended

`adjoint.system.BEMSystem` stores one rotor speed `Ω_b = λ V_b / R` per
operating point in `self.omega`. The local speed ratio, the Reynolds estimate
`W = hypot(V_b, Ω_b r_i)` and the power `P_b = Ω_b Σ_i t_i q_{b,i}` all read
from that list. `CeilingBEMSystem` sets `Ω_b = λ_b V_b / R` with the per-bin
λ and restates the forward solve (the parent class passes a single
`self.tsr`). Nothing else changes: `∂R/∂x` stays diagonal, `∂R/∂d` keeps its
form, and the parent's `gradient()` already assembles

    dJ/dd = Σ_b ω_b dP_b/dd,   ω_b = −(T/1e6) m_b  (below rated),  0  (above rated)

I found no limitation. Setting `vtip_max = None` reproduces the parent class
exactly (`J` and the gradient at x₀ differ by exactly 0).

## Convergence

| case | AEP(x0) | AEP(optimum) | nit | nfev | njev | status | ‖∇J‖ at optimum [MWh/yr per u] | envelope violation | active bounds | time |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|---:|
| none (control) | 10.270157 | 10.282564 | 27 | 29 | 27 | 0, success | 3.6e-4 (free variables) | 0 | `twist_4` lower (−2°), as for x\* | 18 s |
| 60 m/s | 10.210635 | 10.245066 | 44 | 45 | 44 | 0, success | 1.3e-4 | 0 | none | 32 s |
| 55 m/s | 9.976818 | 10.208978 | 49 | 52 | 49 | 0, success | 3.6e-4 | 0 | none (root chord `u = 0.979`, 441 of 450 mm) | 41 s |
| 50 m/s | 9.429005 | 10.168631 | 46 | 48 | 46 | 0, success | 1.7e-4 | 0 | none | 40 s |

`nfev` and `njev` are SciPy's counts; the iteration recorder adds one
objective and one adjoint per accepted step. No polar-range errors occurred,
and the tightest envelope row had 15–33 mm of slack. The gradient at every
optimum is about 1e-4 MWh/yr per unit of `u`, so the objective is flat there,
as it is at the project's own optimum. For comparison, the earlier
finite-difference runs of the same problems took 570–1033 objective
evaluations and 179–389 s.

## Results

AEP is in MWh/yr, with the rating fixed at 3822.19 W. "Cp @ λ = 6.5" is
measured at 8.5 m/s (the Cp–λ curves in `results/cp_lambda.png` are at the same
speed). `Ē[Cp]` is the mean Cp on each case's schedule, weighted by each bin's
unlimited energy `m_b V_b³`, so it's the Cp the site actually sees.

### No ceiling (the control)

| design | objective | AEP | gain vs Schmitz | gain vs x\* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x₀ | – | 10.270157 | – | | 0.4678 | 0.4707 | 0.8087 |
| existing x\* | Cp @ 6.5 (= the project's AEP) | 10.282564 | +0.121 % | – | 0.4687 | 0.4714 | 0.7938 |
| AEP optimum | AEP, no ceiling | 10.282564 | +0.121 % | +0.000 % | 0.4687 | 0.4714 | 0.7938 |

The experiment's optimum **is** x\* (‖u_AEP − u\*‖∞ = 2e-15). The original
objective and "AEP without a ceiling" are the same problem, and the pipeline
reproduces the committed result exactly.

### 60 m/s tip speed (286 rpm)

| design | objective | AEP | gain vs Schmitz | gain vs x\* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x₀ | – | 10.210635 | – | | 0.4678 | 0.3718 | 0.8087 |
| existing x\* | Cp @ 6.5 | 10.213119 | +0.024 % | – | 0.4687 | 0.3690 | 0.7938 |
| **AEP optimum** | AEP @ 60 m/s | **10.245066** | **+0.337 %** | +0.313 % | 0.4677 | 0.3961 | 0.8141 |

### 55 m/s tip speed (263 rpm)

| design | objective | AEP | gain vs Schmitz | gain vs x\* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x₀ | – | 9.976818 | – | | 0.4678 | 0.3402 | 0.8087 |
| existing x\* | Cp @ 6.5 | 9.948790 | −0.281 % | – | 0.4687 | 0.3368 | 0.7938 |
| **AEP optimum** | AEP @ 55 m/s | **10.208978** | **+2.327 %** | +2.615 % | 0.4651 | 0.4040 | 0.8285 |

### 50 m/s tip speed (239 rpm)

| design | objective | AEP | gain vs Schmitz | gain vs x\* | Cp @ 6.5 | Ē[Cp] | Ct @ 6.5 |
|---|---|---:|---:|---:|---:|---:|---:|
| Schmitz x₀ | – | 9.429005 | – | | 0.4678 | 0.3035 | 0.8087 |
| existing x\* | Cp @ 6.5 | 9.306357 | −1.301 % | – | 0.4687 | 0.2995 | 0.7938 |
| **AEP optimum** | AEP @ 50 m/s | **10.168631** | **+7.844 %** | +9.265 % | 0.4604 | 0.4013 | 0.8433 |

### Where the gain comes from

From `gain_attribution_by_bin` in each case file: the AEP optimum gives up
0.1–3.4 % of Cp in the bins that still run at λ = 6.5 (3.5–7.5 m/s; together
−23 / −13 / −6 % of the net gain at 60 / 55 / 50 m/s). It wins that back many
times over in the bins between `V_c` and rated, where λ has dropped. At
55 m/s, the 9.5, 10.5 and 11.5 m/s bins contribute 113 % of the net gain
(Cp 0.455 → 0.469, 0.419 → 0.462 and 0.379 → 0.446). In effect it is a
multi-point design, trading performance at the design point for the band of
λ below it.

### Robustness: every optimum under every schedule

`results/cross_evaluation.json`, as % change vs x₀ under the same schedule:

| schedule ↓ / blade → | x₀ | x\* | optimum (60) | optimum (55) | optimum (50) |
|---|---:|---:|---:|---:|---:|
| none | 0 | +0.12 % | −0.07 % | −0.43 % | −0.98 % |
| 60 m/s | 0 | +0.02 % | **+0.34 %** | +0.20 % | −0.17 % |
| 55 m/s | 0 | −0.28 % | +1.88 % | **+2.33 %** | +2.15 % |
| 50 m/s | 0 | −1.30 % | +4.73 % | +7.57 % | **+7.84 %** |

A blade designed for a ceiling loses little if the ceiling turns out
different (the 55 m/s blade loses 0.43 % with no ceiling and gains 2.15 % at
50 m/s). The Schmitz blade and x\* lose more and more as the ceiling drops.
x\* is worse than Schmitz under any ceiling of 55 m/s or below, because
optimising Cp at one point sharpened its peak.

## Geometry

Control points, with chord in mm and twist in degrees. Changes are in % for
chord and in degrees for twist, since percentages near 0° are meaningless.
`results/geometry.png` shows the distributions along the span.

| variable | x₀ | x\* | optimum (60) | Δ vs x₀ | optimum (55) | Δ vs x₀ | optimum (50) | Δ vs x₀ |
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
| ‖Δu‖∞ vs x₀ | | 0.070 | 0.218 | | 0.409 | | 0.276 | |
| ‖Δu‖∞ vs x\* | | – | 0.194 | | 0.385 | | 0.296 | |

For scale, x\* moved 0.070 in u away from Schmitz. The ceiling optima move
3–6 times further, and in a different direction: more solidity and twist
inboard and at mid-span, rather than the tip unloading that x\* did. The
55 m/s blade's mid-span control point at 261 mm makes a visible bulge at
r ≈ 1.0–1.3 m. With five control points, that is how the spline adds solidity
where the lower-λ bins want it. Whether such a planform could be built is the
question the maximum-chord and solidity limits answer
(`docs/OUTSTANDING-INPUTS.md` §2). The 55 m/s root chord was 9 mm inside the
450 mm bound in use at the time.

## Checks

Every item is recorded under `checks` in `results/case_<label>.json`.

| check | none | 60 | 55 | 50 |
|---|---|---|---|---|
| 1. AEP at the optimum re-evaluated with a fresh parameterisation and resource | Δ = 0 | Δ = 0 | Δ = 0 | Δ = 0 |
| 1. AEP re-assembled from the station integrand `q` by the adjoint system (a second code path) | Δ = −1.8e-15 | −3.6e-15 | 0 | −1.8e-15 |
| 1. AEP through `src/objective.annual_energy_mwh` (control only) | Δ = 0 | n/a | n/a | n/a |
| 2. adjoint vs central FD (h = 3e-6) at **x₀**: max |Δ| / ‖∇‖ | 6.9e-10 | 5.5e-10 | 1.7e-10 | 2.3e-8 |
| 2. adjoint vs central FD (h = 3e-6) at **the optimum**: max |Δ| (scaled) / max |Δ| / ‖∇‖ | 1.2e-10 / 1.6e-7 | 7.9e-11 / 6.2e-6 | 6.6e-11 / 1.9e-6 | 4.1e-11 / 2.3e-6 |
| 2. the same at h = 1e-5 | 1.0e-9 | 7.2e-10 | 2.4e-10 | 2.1e-10 |
| 3. SLSQP status | 0 | 0 | 0 | 0 |
| 4. all `u ∈ [0,1]` | yes | yes | yes | yes |
| 5. BEM converged in all 17 bins at the optimum | yes | yes | yes | yes |
| 5. α in bins below rated within the XFOIL range (−8° to 18°) at the optimum | 0.9° to 6.5° | −0.8° to 6.7° | −2.0° to 6.5° | −0.9° to 6.2° |
| 6. bin weights = CDF(20) − CDF(3); adjoint weights = −(T/1e6)·m_b; bin energies sum to the total | yes / yes / Δ = 0 | same | same | same |
| 7. vs the finite-difference-driven optimum for the same ceiling: ‖Δu‖∞ / ΔAEP | 6e-10 / 2e-13 | 1.3e-8 / 5e-13 | 1.8e-6 / 3e-11 | 1.8e-8 / 5e-13 |

On check 2: at the optima the gradient is about 1e-4 in scaled units, so a
relative error of 1e-6 is an absolute error of 1e-10, which is the
finite-difference noise floor for this objective. It grows ten times from
h = 3e-6 to 1e-5, which is the signature of truncation error in the finite
difference, not of an adjoint error. At x₀, where the gradient is of order
one, agreement is 1e-8 to 1e-10.

On check 7: two independent gradients driving the same optimiser on the same
objective land on the same point, within SLSQP's tolerance, in all four
cases. Together with check 2 that rules out a hidden gradient error, and
together with the control reproducing x\* it rules out a weighting or units
error.

**So nothing numerical is behind the improvement.** What is behind it:

- **In the bins above rated, the angle of attack reaches 28–37° on every
  blade** (λ ≈ 2.5–3 at 16–20 m/s under a ceiling), which is deep stall on the
  Viterna extrapolation. On the optimised blades those bins are all at rated
  power, so they add a constant to J and nothing to the gradient, and can't
  steer the optimiser. That isn't true of the reference blades:
  - At **50 m/s**, the Schmitz blade doesn't reach rated power until
    16.5 m/s. It is scored on 54 station-bins with α > 18° (up to 29.7°) in
    the 11.5–15.5 m/s bins, and those bins carry **55 % of the +7.84 %**.
    x\* never reaches rated under that ceiling (143 station-bins in post-stall,
    up to 36.7°), which is why its figure reads +9.3 %.
  - At **55 m/s**, the Schmitz blade stays below rated up to 11.5 m/s with
    α ≤ 14.1°, inside the valid range. The 11.5 m/s bin carries 43 % of the
    +2.33 %.
  - At **60 m/s**, every below-rated bin of every blade is within range.

  So the 60 and 55 m/s gains are within the polar table's validated range.
  The 50 m/s gain is partly a measure of how badly a Schmitz blade stalls
  when forced to λ = 4 at 12.5 m/s, computed with post-stall polars, which
  is where BEM is least reliable.
- The envelope constraint is conservative under a ceiling (it assumes
  λ = 6.5), but it was never active, so it didn't shape any optimum.
- The bounds were the ones in use at the time, and the 55 m/s root chord was
  9 mm from the upper bound.

## Conclusion

The experiment shows two things clearly. First, the framework can optimise
site AEP under a variable-λ law with the adjoint, cheaply and verifiably.
Second, a rotor-speed ceiling at or below about 55 m/s tip speed turns a
0.12 % problem into a multi-per-cent one, with a genuinely different blade.

What it can't show is which ceiling the real machine has, and picking 55 or
50 m/s because it gives a bigger number wouldn't be honest. The gain is a
steadily increasing function of a machine parameter I didn't have. Without a
fidelity guard, the 50 m/s case also shows that the optimiser will be
credited for the baseline blade stalling. The project went on to set a
300 rpm (62.8 m/s) ceiling on the basis in `docs/OUTSTANDING-INPUTS.md` §9,
where the energy gain is small, and to spend the design freedom on material
instead.

## Files

| file | contents |
|---|---|
| `ceiling_objective.py` | `tsr_schedule`, `CeilingBEMSystem` (the adjoint with per-bin λ), `evaluate_aep` (forward AEP with a per-bin record), `CeilingProblem` (the scaled SLSQP problem) |
| `run_experiment.py` | the four cases end to end: optimise, verify, cross-evaluate, plot |
| `results/case_<none,60,55,50>.json` | per case: configuration, convergence, the AEP table, every check above, gradient checks at x₀ and the optimum, Cp metrics and Cp–λ sweeps, geometry changes, per-bin rows, the gain by bin, and every iterate |
| `results/summary.json` | the headline numbers of the four cases |
| `results/cross_evaluation.json` | every blade under every schedule |
| `results/geometry.png` | chord and twist of x₀, x\* and each case's optimum |
| `results/cp_lambda.png` | Cp–λ at 8.5 m/s for the three blades in each case; the shaded region is the λ band the schedule covers below rated |

## Reproducing it (on the version of the code it was written for)

    python verification/aep_optimisation_experiment/run_experiment.py            # all four cases, ~2.5 min
    python verification/aep_optimisation_experiment/run_experiment.py --cases 55  # one case

The runs are deterministic (two runs gave identical JSON). They read
`verification/baseline/x0.json`, `verification/adjoint_optimisation/result.json`
(for x\*) and, for check 7, `verification/aep_gain_audit/opt_<v>_fixed.json`,
and write only into `results/`. The figures predate the shared figure style.
