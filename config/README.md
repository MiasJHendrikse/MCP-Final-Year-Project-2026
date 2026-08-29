# `config/` — versioned input configuration

Same rule as `data/`: **inputs only, committed, never generated.** `results/` is
where generated output goes.

Five files, one per thing that has its own basis:

| File | What it fixes |
|---|---|
| `site.yaml` | Khomas Hochland site: location, atmosphere (ρ, ν, μ, p, T), wind resource. Plan §1.1–1.3. |
| `rotor_design.yaml` | The SG6043 design rotor being optimised: R, B, λ, wind-speed envelope, parameterisation. Plan §1.4, §2.2, §4. |
| `rotor_phase_vi.yaml` | The NREL Phase VI validation rotor's operating condition, including its **sea-level** air. Plan §2.1. |
| `polars_s809.yaml` | The committed S809 polar cache, as built. Descriptive metadata, not a rebuild instruction. Plan §3.3. |
| `polars_sg6043.yaml` | The design rotor's cache — **not yet built**. A specification rather than a description: the Reynolds envelope it must cover, and the bounds set beyond it. Plan §3.3. |

The last two are the same shape for a reason. `load_polar_cache()` takes a
cache name and returns the same object either way, so code that reads a cache's
build settings does not need to know which of the two it has — and the
SG6043 file's unbuilt fields (`n_crit`) are `TODO` sentinels rather than
S809's numbers copied across.

## Read these only through `src/config`

```python
from config import load_site, load_design_rotor, load_phase_vi_rotor, load_polar_cache
```

`src/config` is the only module in the repo that opens a file in this
directory. That is deliberate and it is what makes plan item 1.1's *"nothing
downstream hard-codes a site value"* a checkable property rather than an
intention: there is exactly one place a site value can enter the code, and

```
grep -rn "1\.225\|1\.5e-5" --include="*.py" src/
```

returning nothing is the check. Before this directory existed, `air_density =
1.225` was a *default argument* on four solver signatures and a sea-level
kinematic viscosity was a module constant in two more — so an AEP call that
forgot the density argument came back 25 % high and still landed inside the
4–6 MWh/yr plausibility band. `solve_rotor` and every `powercurve` entry point
now take ρ and ν as **required keyword arguments**; omitting one is a
`TypeError`, not a wrong number.

## Two rotors, two atmospheres, no default

The Phase VI rotor runs at sea level because that is where the NASA Ames
experiment ran. The design rotor runs at ρ ≈ 0.968 kg/m³ because the site is at
1800 m. Neither is a default for the other, and there is no global fallback —
which is why the sea-level pair lives in `rotor_phase_vi.yaml` (a property of
*that rotor's* validation case) and the site pair lives in `site.yaml` (a
property of the *site*, which `rotor_design.yaml` deliberately does not
duplicate).

The same reasoning applies to the polar caches: S809 is built at `Ncrit = 5`
for a documented physical reason, SG6043 will be built at whatever its own
sensitivity study selects, so calibration is per-cache metadata and
`load_polar_cache()` takes the cache name rather than assuming one.

## `TODO:` fields are not placeholders

Any value written as a `TODO: …` string is **unresolved**. The loader turns it
into an `Unresolved` object that raises `UnresolvedConfigError` — naming the
field and quoting the note — on any attempt to use it as a number. It is never
substituted, defaulted, or filled with the "expected range" the plan records
for sanity-checking.

Outstanding at the time of writing:

- `site.yaml` → `location.gwa_area`, `latitude_deg`, `longitude_deg`, and the
  whole `wind_resource` block. Waiting on the Global Wind Atlas extraction. A
  fabricated wind resource propagates silently into every AEP figure and
  invalidates the results chapter, so these stay `TODO` (plan §1.3).
- `polars_sg6043.yaml` → `build.ncrit`. Plan §1.2's sensitivity study selects
  it. S809's `Ncrit = 5` is a result about a 21 % thick section at low Re and
  does not transfer; copying it across is exactly the per-cache-calibration
  mistake the two-file split exists to prevent.
- `rotor_design.yaml` → `parameterisation.bounds.*`. Plan §7.1 fixes that chord
  and twist bounds exist and are applied to the *scaled* variables; it commits
  to no numbers and no study has selected any.

Test with `config.is_resolved(value)` — truthiness deliberately raises, because
`if cfg.weibull_k:` is a silent fallback waiting to happen.

## Derived atmosphere values are checked, not trusted

`site.yaml`'s pressure, density and kinematic viscosity are *derived*
quantities:

```
p   = p0 · (1 − L·h/T0)^(g/(Rs·L))
ρ   = p / (Rs · T)
ν   = μ / ρ
```

The loader recomputes all three from `location.elevation_m`,
`atmosphere.temperature_k`, `atmosphere.dynamic_viscosity_pa_s` and
`atmosphere.standard_atmosphere`, and **rejects the file** if a recorded value
disagrees by more than 0.5 %. Editing the mean site temperature and forgetting
to update the density is otherwise a silent few-percent error in every AEP
figure downstream.

They are recorded as well as derived because they are reported numbers — the
plan's working conventions require a number that feeds a later step to be
written to the config file rather than left in prose.

> **Open discrepancy, deliberately left visible.** Plan §1.2 quotes ρ = 0.98
> kg/m³ and ν = 1.85 × 10⁻⁵ m²/s as its headline figures. Recomputing from the
> plan's own stated constants at T = 293.15 K gives **0.9684** and **1.869 ×
> 10⁻⁵** — 1.2 % and 1.1 % away. The recomputed values are what this directory
> records; the plan's pair appears to be rounded. Both sit inside plan item
> 1.1's 0.95–1.03 kg/m³ exit band. Flagged for the plan and report to
> reconcile, not silently averaged.

## Pointing somewhere else

`BLADE_CONFIG_DIR` overrides the directory. Plan §1.3 calls for the
optimisation to be re-run across a resource-uncertainty band; that is the
mechanism, so a sensitivity case is an alternative config directory rather than
an edit to a tracked file.
