# Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade

Final-year Mechanical Engineering project (DSP810S, NUST). An adjoint-BEM
aerodynamic optimisation framework for small wind turbine blades targeting
Namibian low-wind conditions (~5–6 m/s mean, Windhoek reference site).

**Full title:** *Gradient-Based Aerodynamic Optimisation of a Small Wind
Turbine Blade for Low-Wind-Speed Conditions Using an Adjoint-BEM Framework
and CFD Validation*

- **Author:** MJ Hendrikse
- **Supervisor:** Prof. Hannes van der Walt
- **Timeline:** 13 June – 6 November 2026

Progress, phase-by-phase, is tracked in
[`docs/journal/PROJECT_PLAN.md`](docs/journal/PROJECT_PLAN.md). The working
record of every session is in `docs/journal/Session Notes/`.

---

## Pipeline

XFOIL polars → BEM solver → discrete adjoint → SLSQP optimisation →
Weibull/AEP → STEP export → (optional) Streamlit UI.

Phase 1 (the forward BEM solver) is built and cross-validated. Phase 2
(the discrete adjoint) is the next major piece of work.

---

## Repository layout

```
src/
  config/       The only reader of config/. Loader + frozen dataclasses.
    schema.py        frozen dataclasses, one per config file
    loader.py        YAML read, TODO handling, derived-value checks
    unresolved.py    the TODO sentinel that raises instead of defaulting
  bem/          Solver core. No I/O beyond the polar cache, no plotting.
    airfoil.py       polar adapters (LinearPolar, S809Polar)
    corrections.py   Prandtl tip/hub loss, Glauert/Buhl high-thrust
    station.py       single-station Ning (2014) residual + root-find
    rotor.py         spanwise loop, geometry, one operating point
    powercurve.py    wind-speed / TSR sweeps in dimensional units
  polars/       The cache-consuming layer. No XFOIL dependency, by design.
    viterna.py       post-stall extrapolation to +/-180 deg
    cache_format.py  the cache CSV schema and its provenance column
    envelope.py      the design rotor's Reynolds envelope, which sizes its cache
  xfoil/        Polar generation and lookup.
    xfoil_runner.py       XFOIL subprocess wrapper
    build_polar_cache.py  sweep Re, write data/polars/<airfoil>/
    close_polar_gaps.py   repair a committed cache's holes, extend to +/-180
    polar_lookup.py       bilinear (alpha, Re) interpolation
  validation/   Everything that checks or cross-checks the above.
    validate_stage1..4.py    BEM solver, stage by stage
    validate_powercurve.py   power curve / Cp-lambda sweep layer
    validate_polars.py       5-check audit of a polar cache
    check_stitch_continuity.py  value/slope jumps where XFOIL meets Viterna
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
  polars_sg6043.yaml   the design rotor's cache -- specified, not yet built
data/           Inputs only.
  airfoils/       Selig-format coordinates (S809)
  polars/         the XFOIL polar cache, per airfoil, per Reynolds
  naca0012_validated/  Abbott & von Doenhoff + Ladson experimental data
  qblade/         QBlade .bld/.plr exports and their reference results
docs/
  journal/        Obsidian vault: project plan, overview, session notes
  validation/     cross-validation writeups, plots, reference case files
results/          Generated plots and polars. results/_archive/ is scratch
                  (gitignored) — the XFOIL scripts write raw output there.
```

**Rule of thumb:** `data/` and `config/` are inputs, `results/` is generated
output, `src/` is the only place Python lives.

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

Every command below is run **from `src/`**, since `bem`, `xfoil`,
`validation` and `demos` are sibling packages there.

```bash
cd src
```

### Validate the BEM solver

Run in order; each is a regression check on everything before it.

```bash
python -m validation.validate_stage1      # single-station induction solve
python -m validation.validate_stage2      # + Prandtl tip/hub loss
python -m validation.validate_stage3      # + Glauert/Buhl high-thrust
python -m validation.validate_stage4      # + multi-station, real S809 polars
python -m validation.validate_powercurve  # + power curve / Cp-lambda sweep
```

Each ends with `ALL CHECKS PASSED`. Stage 3 also prints one intentionally
flagged, non-fatal note about a pre-existing residual pole — that is
expected, not a failure.

### Audit the polar cache

```bash
python -m validation.validate_polars s809       # currently 4/5
python -m validation.validate_polars naca4412   # currently 3/5
```

Both check the XFOIL-converged band only; the +/-180 deg Viterna extension the
cache files carry is checked separately, for continuity with the measured data
it is anchored to:

```bash
python -m validation.check_stitch_continuity s809
```

Both polar audits have one documented, deliberately-retained residual — see the READMEs
in `data/polars/*/` and the 2026-07-26 / 2026-07-28 journal entries. A
non-zero exit here is the known state, not a new regression.

### Use the solver

```python
from bem.rotor import phase_vi_geometry, PHASE_VI_RATED_RPM
from bem.powercurve import power_curve, cp_lambda_curve, peak_cp

geometry = phase_vi_geometry()

# Fixed-speed machine (NREL Phase VI: 71.63 RPM synchronous)
curve = power_curve(geometry, [5, 7, 10, 13, 15, 20, 25], rpm=PHASE_VI_RATED_RPM)
for p in curve["points"]:
    print(f"V={p['v_inf']:5.1f} m/s  Cp={p['Cp']:.4f}  P={p['power']:9.1f} W")

# Cp-lambda curve
cl = cp_lambda_curve(geometry, [1.0 + 0.5 * i for i in range(23)])
print(peak_cp(cl))   # peak Cp = 0.4159 at TSR = 7.0
```

### Regenerate the polar cache (needs XFOIL)

```bash
python -m xfoil.build_polar_cache s809
```

Set the XFOIL path via `_DEFAULT_XFOIL` in `src/xfoil/xfoil_runner.py`.

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
item in `PROJECT_PLAN.md`.
