# Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade

Final-year Mechanical Engineering project (DSP810S, NUST). An adjoint-BEM
aerodynamic optimisation framework for small wind turbine blades targeting
Namibian low-wind conditions (Khomas Hochland site, **V̄ = 6.49 m/s** at the
20 m hub height — measured, see `verification/wind_resource/`).

**Title:** *Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine
Blade for Low-Wind-Speed Conditions Using an Adjoint-BEM Framework*

**Scope note.** CFD validation is **out of scope** — see `PROJECT_PLAN.md`,
"Priority hierarchy": with ~11 weeks to 6 November and the WCE placement load,
CFD was cut to protect the core adjoint work. **Prof. van der Walt has signed
off on that scope boundary** (confirmed by MJ, 2026-09-13), so the title above
is shortened accordingly — it previously had to be quoted in its longer
registered form, *"…Adjoint-BEM Framework and CFD Validation"*, because
asserting an approval that had not happened would have been worse than a
self-contradictory title.

> If the title as *formally registered* with the department still carries "and
> CFD Validation", amending that registration is a separate administrative step
> from the scope sign-off. Worth checking which was agreed.

- **Author:** MJ Hendrikse
- **Supervisor:** Prof. Hannes van der Walt
- **Timeline:** 13 June – 6 November 2026

Progress, phase-by-phase, is tracked in
[`docs/journal/PROJECT_PLAN.md`](docs/journal/PROJECT_PLAN.md). The working
record of every session is in `docs/journal/Session Notes/`.

---

## Pipeline

XFOIL polars → BEM solver → objective function → finite-difference gradients →
discrete adjoint → structural constraints (three adjoint rows) → SLSQP
minimum-material optimisation at fixed energy.

**Re-pitched 2026-09-20.** Under the settled operating law the verified
energy optimum beats the polar-consistent Schmitz blade by **+0.147 %** and
carries +6.2 % material for it; annual energy is a plateau in the blade. The
production problem is therefore *minimise blade material subject to AEP ≥
Schmitz's, the root-moment cap, a root-stress proxy, a tip-deflection proxy
and monotone chord/twist* — the verified objective adjoint as the energy-floor
Jacobian, the moment adjoint as the cap and stress rows, and one new sibling
adjoint (deflection). Result: **−3.4 % shell material at exactly Schmitz's
energy, −5.3 % within 0.5 % of it**, with a 60 mm buildable-tip floor
(`verification/mass_optimisation/`); the plan and its cut list are
`docs/PLAN-mass-objective-2026-09-20.md`. Title unchanged; research question
rewritten in `docs/journal/PROJECT_PLAN.md`, Framing.

Phase numbering follows the plan rewrite of 2026-08-23:

| Phase | | Status |
|---|---|---|
| **0** | Tooling, solver, cross-validation | complete |
| **1** | Objective function (AEP, parameterisation, baseline, smoothness gate) | complete (2026-09-13); bounds decided 2026-09-19 (`chord_max_m = 0.30 m`, in `config/`) |
| **2** | Finite-difference gradient path | complete (2026-09-13) — `verification/fd_*` |
| **3** | Discrete adjoint | complete (2026-09-13) — Tiers 1–4 verified, `verification/gradient_verification/`, `docs/adjoint_derivation.md` |
| **4** | Structural constraint and cost scaling | complete (2026-09-19) — relative root-moment KS constraint, `verification/load_constraint/`; scaling law, `verification/cost_scaling/` |
| **5** | Production runs and results — **re-pitched to minimum-material design 2026-09-20** | **run (2026-09-20)** — `verification/mass_optimisation/`: `x_m` at −3.38 % shell / −7.28 % solid material at AEP(x₀) (−5.27 % / −9.67 % within 0.5 %), 60 mm tip floor, ten starts agree, δ ∈ {0 … 2 %} front with KKT exchange rates, ablation, Tiers 1–3 at `x₀` and `x_m`, rendered blades, cross-evaluation. The energy optimum `x_c` (`verification/load_constraint/result_eps0.json`, +0.147 %) is now the comparison blade |
| **6** | Report | next — the results are all committed; Phase 6 is the only remaining work |
| **7** | Optional (STEP export, DE, UI) | **cut 2026-09-20** |

The forward BEM solver is Phase **0**, not Phase 1, and the adjoint is Phase
**3**. Dates and exit criteria for each are in
[`docs/journal/PROJECT_PLAN.md`](docs/journal/PROJECT_PLAN.md).

---

## Repository layout

```
src/
  config/       The only reader of config/. Loader + frozen dataclasses.
    schema.py        frozen dataclasses, one per config file
    loader.py        YAML read, TODO handling, derived-value checks
    unresolved.py    the TODO sentinel that raises instead of defaulting
  bem/          Solver core. No I/O beyond the polar cache, no plotting.
    airfoil.py       LinearPolar, the analytic fixture (real polars: polars/)
    corrections.py   Prandtl tip/hub loss, Glauert/Buhl high-thrust
    station.py       single-station Ning (2014) residual + root-find
    rotor.py         spanwise loop, geometry, one operating point
    powercurve.py    wind-speed / TSR sweeps in dimensional units
  polars/       The cache-consuming layer. No XFOIL dependency, by design.
    viterna.py       post-stall extrapolation to +/-180 deg
    cache_format.py  the cache CSV schema and its provenance column
    envelope.py      the design rotor's Reynolds envelope, which sizes its cache
    cache.py         CSV load + ragged->rectangular reindexing into a PolarGrid
    interpolant.py   C1 (alpha, Re) surface with analytic partials
    polar.py         CachedPolar: one station's view of it, no clamping
  design/       The design vector and what it turns into. Phase 1.6/1.7.
    parameterisation.py  B-spline chord/twist; d -> RotorGeometry, linear in d
    bounds.py            box bounds, scaling, feasibility (bounds are TODO)
    schmitz.py           Schmitz optimum-rotor closed form, max-L/D point
    baseline.py          the Schmitz baseline blade and its x0
  objective/    What the optimiser will minimise. Phase 1.5.
    power.py         wind-speed bins, aerodynamic power, rated limiting
    weibull.py       WeibullResource; bin masses from the CDF
    height_extrapolation.py  Weibull k, c between heights (Justus & Mikhail)
    objective.py     AEP, the unit-weighted surrogate, J, the sanity band
    loads.py         the load set L, root moment, spanwise moments, tip
                     deflection (forward twins of the adjoint functionals), KS
    mass.py          the material proxies: shell k_P int c dr (the Phase 5
                     objective, exact gradient) and solid k_A int c^2 dr (reported)
  xfoil/        Polar generation and lookup.
    xfoil_runner.py       XFOIL subprocess wrapper
    build_polar_cache.py  sweep Re, write data/polars/<airfoil>/
    close_polar_gaps.py   repair a committed cache's holes, extend to +/-180
    polar_lookup.py       bilinear (alpha, Re) interpolation -- retired from
                          the solve path, kept as the Task 3 "before" baseline
  validation/   Cross-checks against external tools and measured data.
                Solver checks live in tests/ -- see below.
    validate_polars.py       5-check audit of a polar cache
    check_stitch_continuity.py  value/slope jumps where XFOIL meets Viterna
    ncrit_sensitivity.py     the {5,7,9,11,13} n_crit sweep behind the SG6043 cache
    uiuc_sg6043.py           digitised UIUC SG6043 wind-tunnel data
    compare_sg6043_uiuc.py   XFOIL vs that measured data
    compare_pybemt.py        cross-check vs pyBEMT
    compare_ccblade.py       cross-check vs CCBlade
    compare_qblade.py        cross-check vs QBlade CE
    export_qblade.py         write QBlade .bld / .plr files
    plot_bem_comparison.py   combine cross-check results into a plot
  demos/        Illustrative plotting scripts, not part of the solver.

config/         Versioned input configuration. Read only through src/config.
  site.yaml            Khomas Hochland site: location, atmosphere, wind resource
  rotor_design.yaml    the SG6043 design rotor being optimised
  rotor_phase_vi.yaml  the NREL Phase VI validation rotor's operating condition
  polars_s809.yaml     the committed S809 polar cache, as built
  polars_sg6043.yaml   the design rotor's cache, as built
data/           Inputs only.
  airfoils/       Selig-format coordinates (S809 validation, SG6043 design)
  polars/         the XFOIL polar cache, per airfoil, per Reynolds
                  s809/ validation, sg6043/ design, naca4412/ retained reference
  naca0012_validated/  Abbott & von Doenhoff + Ladson experimental data
  qblade/         QBlade .bld/.plr exports and their reference results
docs/
  journal/        Obsidian vault: project plan, overview, session notes
  validation/     cross-validation writeups, plots, reference case files
  OUTSTANDING-INPUTS.md  external data this project is waiting on — read this
                  before wondering why something raises
tests/            pytest suite: the machine-checkable Phase 1 exit criteria.
  golden/         Task 0 regression snapshot: Cp(lambda), spanwise a/a'/phi
  golden_reference.py, generate_golden.py   the snapshot's loader and writer
  test_invariants.py  AST-checked import rules (e.g. bem/ must not import xfoil)
  test_loads.py       Phase 4: the root-moment integrand, its adjoint, the KS
                      constraint on ScaledProblem (Tiers 1-3 at x0)
  test_mass.py, test_deflection.py, test_mass_problem.py   Phase 5: the
                      material proxies, the deflection adjoint (Tiers 1-3), and
                      the mass problem's rows, assembly, shared solve and guard
  test_operating_law_control.py  the pre-law numbers pinned bit-for-bit
verification/     Versioned report figures — committed evidence, not scratch.
                  verification/README.md is the index and the re-run order.
  polar_interpolant/     C1 interpolant vs the bilinear staircase
  phase_vi/              solver residual histories across the envelope
  wind_resource/         the 20 m Weibull fit behind AEP
  representation_study/  control-point count, justified against Schmitz
  baseline/              the Schmitz baseline blade and x0.json, the design
                         vector every later phase starts from
  smoothness_gate/       plan step 1.8: is J smooth enough to differentiate
  spline_fit_error/      the projection error of the baseline's spline
  fd_step_size/          the per-variable FD step h*_j and scale eps_j
  fd_optimisation/, fd_optimisation_multistart/   A4: FD-driven SLSQP, and its starts
  gradient_verification/ Tiers 3-4: adjoint vs FD, and what is left
  adjoint_optimisation/  B5: adjoint-driven SLSQP, agreement with A4
  load_constraint/       Phase 4: the root-moment KS constraint, Pareto 0-10 %
  cost_scaling/          Phase 4: gradient wall time vs n = 10..160
  mass_optimisation/     Phase 5: the minimum-material blade x_m, multi-start,
                         the energy-floor front, ablation, checks, cross-evaluation
  aep_gain_audit/, aep_optimisation_experiment/   frozen records behind
                         docs/AEP_GAIN_AUDIT.md and the B1/O4 decisions
results/          Generated plots and polars. results/_archive/ is scratch
                  (gitignored) — the XFOIL scripts write raw output there.
misc/             Gitignored. Scratch for things written for MJ rather than for
                  the record — session briefings and the like.
```

Each `verification/` subdirectory has a README stating what was run, what came
out, and how to reproduce it. `verification/smoothness_gate/README.md` carries
the Phase 1 exit verdict: **`J` is C¹ but not C²**, which is what the adjoint
needs, with the one non-smoothness traced to stations crossing the Buhl
`a = 0.4` threshold.

**Waiting on external input.** Several exit criteria are blocked on data that
cannot be inferred or defaulted: the manufacturability bounds, and three papers
needed for provenance. The wind resource **landed on 2026-09-13** — as a GASP
point extraction at 50 m, extrapolated to the 20 m hub height and recorded in
[`verification/wind_resource/`](verification/wind_resource/) — which in turn
exposed a defective exit criterion: plan §1.4's 4–6 MWh/yr AEP sanity band
compared an electrical capacity-factor estimate against a model that computes
aerodynamic shaft energy, and was never consistent with the plan's own resource
prior in the first place. Revised to 8–12 MWh/yr on 2026-09-13, with the
derivation recorded beside the numbers and a test asserting the band is still
narrow enough to catch a bug.
[`docs/OUTSTANDING-INPUTS.md`](docs/OUTSTANDING-INPUTS.md) is the single list,
with what each blocks and what happens in the meantime. Code that needs an
unresolved value raises `UnresolvedConfigError` naming the field, rather than
substituting anything.

**Rule of thumb:** `data/` and `config/` are inputs, `results/` is generated
output, `src/` is the only place Python lives.

**Two airfoils, two jobs.** **SG6043** is the **design** airfoil — the section
on the blade being optimised, chosen for its low-Reynolds performance at this
rotor's scale. **S809** is the **validation** airfoil: it is the NREL Phase VI
rotor's section, and it exists in this repository so the solver can be checked
against a rotor other people have also computed. Neither is a default for the
other, and there is no project-wide "primary" airfoil. `RotorGeometry` takes a
required `polar_cache` field with **no default anywhere**, so a rotor cannot be
built without saying which cache it reads: `phase_vi_geometry()` pins `"s809"`,
the design rotor pins `"sg6043"`. NACA 4412 was the original target before
either and is retained only as a secondary reference dataset.

**No site value is hard-coded.** Air density and kinematic viscosity are
**required keyword arguments** on `solve_rotor` and every `powercurve` entry
point, supplied from `config/` — omitting one is a `TypeError`, not a quietly
wrong number. The two rotors have two different atmospheres (Phase VI at sea
level because that is where the experiment ran; the design rotor at 1800 m),
so neither is a default for the other and there is no global fallback. See
[`config/README.md`](config/README.md).

---

## Setup

```bash
python -m pip install -r requirements.txt
```

Python 3.12+. XFOIL, pyBEMT, CCBlade and QBlade are **not** pip dependencies
and are not needed to run the solver — the polar cache is committed. See the
comments in `requirements.txt` for what each external tool is for.

---

## Running

Two working directories, and it matters which:

- **`pytest` and the `verification/` scripts run from the repo root.** No
  install step and no `cd`; `pyproject.toml` puts `src/` and `tests/` on the
  path.
- **The `python -m …` module commands run from `src/`**, since `bem`, `xfoil`,
  `validation`, `design`, `objective` and `demos` are sibling packages there.

### Validate the BEM solver

```bash
pytest          # from the repo root
```

**510 passed, 5 xfailed, 1 failed, ~85 s** (2026-09-19). The five `xfail`s
carry a documented physical reason (the `validate_polars` checks a cache is
known to fail — see `tests/test_polar_cache.py`, where each carries its
explanation). The one failure is deliberate and unchanged:
`tests/test_adjoint_gradient.py::test_adjoint_agrees_with_the_committed_fd_reference_at_x0`
on `chord_4`, ratio 16.302328995917193 — the adjoint is not the suspect, the
acceptance scale `eps_j` is; see `verification/gradient_verification/README.md`
("One failure, and what it is"). A different ratio, variable or second failure
is a regression.

**One test is designed to fail later, on purpose.** Where an external input
is still missing, the mechanism around it is built and tested, and a test
asserts that it *still raises*:

| test | file | fails when |
|---|---|---|
| `test_gwa_area_is_still_unresolved_and_behaves_like_it` | `test_config.py` | a Global Wind Atlas *area* extraction is done |

Two more were retired on 2026-09-19 when the bounds landed in `config/`
(`test_bounds_from_config_still_raise` is now
`test_bounds_from_config_carry_the_decided_values`;
`test_feasibility_reports_that_it_was_not_checked` is now
`test_x0_is_feasible_against_the_configured_bounds`), and three on 2026-09-13
when the wind resource landed.
**One of the 2026-09-13 three did not fire, and could not have** — it asserted against a
hard-coded string in `baseline.py` rather than against the thing that actually
changes state. A guard has to be attached to the mechanism, not to a
description of it; see the 2026-09-13 journal entry.

If one of these fails, nothing is broken — it means a value in
`docs/OUTSTANDING-INPUTS.md` arrived and the placeholder needs removing.
Documenting a known gap is easy; the point of these is guaranteeing somebody
notices when it closes.

The solver checks that used to be `validation/validate_stage1..4.py` and
`validate_powercurve.py` are now `tests/test_bem_stages.py` and
`tests/test_powercurve.py`. Same physics, same tolerances, same reference
values — what changed is that they are a suite rather than five scripts run
by hand, two of which used to exit non-zero in their known-good state.

### Audit the polar cache

```bash
python -m validation.validate_polars sg6043     # currently 3/5  (design)
python -m validation.validate_polars s809       # currently 4/5  (validation)
python -m validation.validate_polars naca4412   # currently 3/5  (reference)
```

These check the XFOIL-converged band only; the +/-180 deg Viterna extension the
cache files carry is checked separately, for continuity with the measured data
it is anchored to:

```bash
python -m validation.check_stitch_continuity s809
python -m validation.check_stitch_continuity sg6043
```

Every one of these audits has documented, deliberately-retained failures — see
the READMEs in `data/polars/*/` and the 2026-07-26 / 2026-07-28 journal
entries. **A non-zero exit here is the known state, not a new regression.**
`validate_polars` applies plausibility heuristics written for thick, mildly
cambered sections; SG6043 is a thin high-camber low-Reynolds section, and its
two failures (Cl(0) spread of 0.638, and Cd monotonicity) are the *physics* of
laminar-separation-bubble behaviour at 40k–100k rather than cache defects. That
argument is made in full, and checked against UIUC wind-tunnel measurements, in
[`data/polars/sg6043/README.md`](data/polars/sg6043/README.md).

### Use the solver

```python
from config import load_phase_vi_rotor
from bem.rotor import phase_vi_geometry, PHASE_VI_RATED_RPM
from bem.powercurve import power_curve, cp_lambda_curve, peak_cp

rotor = load_phase_vi_rotor()
geometry = phase_vi_geometry()          # pins polar_cache="s809"

# The atmosphere is required, not defaulted -- see "No site value is
# hard-coded" above. Phase VI ran at sea level, so it comes from that rotor's
# own config, never from the design rotor's.
air = dict(air_density=rotor.air_density,
           kinematic_viscosity=rotor.kinematic_viscosity)

# Fixed-speed machine (NREL Phase VI: 71.63 RPM synchronous)
curve = power_curve(geometry, [5, 7, 10, 13, 15, 20], rpm=PHASE_VI_RATED_RPM, **air)
for p in curve["points"]:
    print(f"V={p['v_inf']:5.1f} m/s  Cp={p['Cp']:.4f}  P={p['power']:9.1f} W")

# Cp-lambda curve
cl = cp_lambda_curve(geometry, [1.0 + 0.5 * i for i in range(14)], **air)
print(peak_cp(cl)["Cp"], peak_cp(cl)["tsr"])   # 0.4163 at TSR = 7.0
```

```
V=  5.0 m/s  Cp=0.4062  P=   2471.1 W
V=  7.0 m/s  Cp=0.3671  P=   6127.8 W
V= 10.0 m/s  Cp=0.2566  P=  12487.4 W
V= 13.0 m/s  Cp=0.1582  P=  16912.7 W
V= 15.0 m/s  Cp=0.1123  P=  18450.5 W
V= 20.0 m/s  Cp=0.0533  P=  20769.2 W
```

**Why the sweeps stop at 20 m/s and λ = 7.5.** Above those the blade's station
Reynolds numbers leave the committed S809 cache (which tops out at 1.3 M) and
`CachedPolar` raises `PolarDomainError` naming the station. That is deliberate:
the previous code silently clamped to the cache edge, and removing the clamp
moved Cp by up to −61 % at low λ. An out-of-range request is now an error
rather than a plausible wrong number. Sequence S also has a 25 m/s point; it is
recorded in `tests/test_powercurve.py` as a point this cache **cannot serve**,
not as a value to compare against.

### Regenerate the polar cache (needs XFOIL)

```bash
python -m xfoil.build_polar_cache s809
python -m xfoil.build_polar_cache sg6043
```

Set the XFOIL path via `_DEFAULT_XFOIL` in `src/xfoil/xfoil_runner.py`. This is
still a hard-coded absolute path into one machine's filesystem and should move
to an environment variable or `config/` (a known, unfixed chore — the repo-audit
work order raises it under plan item 1.2). It affects only cache regeneration,
never the solver, since the caches are committed.

**Neither cache needs regenerating to use this repository**, and the S809 cache
in particular is deliberately frozen: the golden regression and every
cross-validation figure are anchored to it exactly as committed.

### Reproduce the Phase 1 evidence

```bash
python verification/baseline/generate_baseline.py         # the Schmitz x0
python verification/representation_study/run_study.py     # control-point count
python verification/phase_vi/generate_residual_histories.py
python verification/smoothness_gate/run_gate.py           # ~11 min
```

---

## Validation status

The BEM solver is cross-validated against three independent BEM codes on the
real NREL Phase VI (Sequence S) rotor geometry:

| Reference | Mean \|Cp deviation\| |
|---|---|
| CCBlade (same Ning 2014 formulation) | **0.51 %** |
| pyBEMT (classic BEM, fewer corrections) | 1.83 % |
| QBlade CE (with matched Viterna extrapolation) | ~4–6 % in deep stall, ~1–2 % attached |

Full writeup, per-station diagnostics and exact reproduction steps:
[`docs/validation/bem-cross-validation.md`](docs/validation/bem-cross-validation.md).

**This is a cross-check between independent BEM codes, not a validation
against ground truth.** NREL Phase VI *experimental* performance data has not
been sourced (see the 2026-07-25 journal entry, Stage 5), so the Cp–λ curve
has never been checked against measurements. That gap is tracked as an open
item in `PROJECT_PLAN.md` and as item 6 of `docs/OUTSTANDING-INPUTS.md`.

---

## Phase 1 status

The objective function is built, and the question Phase 1 exists to answer has
been answered:

> **Is `J` smooth enough to differentiate?** Yes. `J` is **C¹ but not C²**,
> which is what the adjoint requires and what finite differences need to be
> interpretable.

3000 objective evaluations across all 10 design variables converged, with no
staircasing, no high-frequency noise and no discontinuous jumps. The single
non-smoothness is a curvature break where blade stations cross the Buhl
turbulent-wake threshold `a = 0.4` — inherent to a correction constructed to
match momentum theory in value and slope and nothing further, not an
implementation defect. Evidence and the full argument:
[`verification/smoothness_gate/README.md`](verification/smoothness_gate/README.md).

Status of the exit criteria:

| exit criterion | state |
|---|---|
| AEP of the baseline blade in the sanity band | **done** — 10.25 MWh/yr (10.2477 under the 300 rpm law, 2026-09-19), inside the 8–12 MWh/yr band. The band was revised from 4–6 on 2026-09-13; see `config/rotor_design.yaml` |
| Feasibility of `x0` against the bounds | **done** — checked against the configured bounds (`chord_max_m = 0.30 m`, `max_local_solidity 0.5`, decided 2026-09-19; the 0.45 m placeholder is retired), zero violations |
| Everything else in Phase 1 | done |

The mechanism around each hole is complete and tested; nothing anywhere
substitutes a default for a missing number, and where a computed result
disagrees with a stated expectation the disagreement is recorded rather than
tuned away. See
[`docs/OUTSTANDING-INPUTS.md`](docs/OUTSTANDING-INPUTS.md) for what is needed,
and the *"When the data lands"* section of
`docs/journal/Session Notes/2026-09-10.md` for what to do when it arrives.
