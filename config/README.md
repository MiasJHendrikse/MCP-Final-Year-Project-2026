# Configuration

Like `data/`, this folder holds **inputs only**: committed by hand, never
generated. Generated output goes in `results/` and `verification/`.

There are five files, one for each thing that has its own basis:

| File | What it sets |
|---|---|
| `site.yaml` | the Khomas Hochland site: location, atmosphere (ρ, ν, μ, p, T) and wind resource |
| `rotor_design.yaml` | the SG6043 design rotor: R, B, λ, wind-speed range, operating law, parameterisation and bounds, the constraints of the design problem, and the laminate |
| `rotor_phase_vi.yaml` | the NREL Phase VI validation rotor's operating condition, including its **sea-level** air |
| `polars_s809.yaml` | the S809 polar table, as built |
| `polars_sg6043.yaml` | the SG6043 polar table, as built (`n_crit = 9`, selected by the sensitivity study in `results/ncrit_sensitivity/`) |

The two polar files have the same structure on purpose. `load_polar_cache()`
takes a table name and returns the same kind of object either way, so code
that reads build settings doesn't need to know which one it has.

## Read these only through `src/config`

```python
from config import load_site, load_design_rotor, load_phase_vi_rotor, load_polar_cache
```

`src/config` is the only module that opens files in this folder. That makes
"nothing downstream hard-codes a site value" something you can check rather
than just intend: there's exactly one place a site value can enter the code,
and

```
grep -rn "1\.225\|1\.5e-5" --include="*.py" src/
```

returning nothing is the check. Before this folder existed, `air_density =
1.225` was a *default argument* on four solver functions and a sea-level
viscosity was a module constant in two more. So an AEP call that forgot the
density argument came back 25 % high, and still landed inside the sanity band
in use at the time. The band has since been corrected to 8–12 MWh/yr and
checked against exactly that mistake: sea-level density now gives about 13.0
and **is** caught. A band centred on the right value is a stronger check than
a narrower one centred on the wrong value. `solve_rotor` and every
`powercurve` function now take ρ and ν as **required keyword arguments**, so
leaving one out is a `TypeError`, not a wrong number.

## Two rotors, two atmospheres, no default

The Phase VI rotor runs at sea level because that's where the NASA Ames
experiment was done. The design rotor runs at ρ ≈ 0.968 kg/m³ because the site
is at 1,800 m. Neither is a default for the other, and there's no global
fallback. That's why the sea-level values live in `rotor_phase_vi.yaml` (a
property of that rotor's validation case), and the site values live in
`site.yaml` (a property of the site, which `rotor_design.yaml` deliberately
doesn't duplicate).

The same idea applies to the polar tables. S809 is built at `Ncrit = 5` and
SG6043 at `Ncrit = 9`, each for its own documented reason, so calibration is
stored per table and `load_polar_cache()` always takes the table name.

## `TODO:` fields raise errors

A value written as a `TODO: …` string is **unresolved**. The loader turns it
into an `Unresolved` object that raises `UnresolvedConfigError`, naming the
field and quoting the note, as soon as anything tries to use it as a number.
It is never replaced by a default or a guess.

Two fields are currently unresolved, and no result depends on either:

- `site.yaml` → `location.gwa_area`. The wind data is a single-point
  extraction, so there's no area selection to record
  (`verification/wind_resource/README.md`).
- `rotor_design.yaml` → `structure.tip_clearance_m`. This depends on the tower
  and hub geometry, which aren't specified (`docs/OUTSTANDING-INPUTS.md` §11).

Test for them with `config.is_resolved(value)`. Truthiness deliberately raises,
because `if cfg.weibull_k:` is a silent fallback waiting to happen.

## Derived atmosphere values are checked

The pressure, density and kinematic viscosity in `site.yaml` are *derived*:

```
p   = p0 · (1 − L·h/T0)^(g/(Rs·L))
ρ   = p / (Rs · T)
ν   = μ / ρ
```

The loader recomputes all three from `location.elevation_m`,
`atmosphere.temperature_k`, `atmosphere.dynamic_viscosity_pa_s` and
`atmosphere.standard_atmosphere`, and **rejects the file** if a recorded value
is more than 0.5 % off. Otherwise, changing the site temperature and forgetting
to update the density would be a silent few-per-cent error in every AEP figure.

They're recorded as well as derived because they are reported numbers, and a
number that feeds later work should be written into the config rather than
left in prose.

One discrepancy is worth noting. My early project plan quoted ρ = 0.98 kg/m³
and ν = 1.85 × 10⁻⁵ m²/s. Recomputing from the stated constants at
T = 293.15 K gives **0.9684** and **1.869 × 10⁻⁵**, 1.2 % and 1.1 % different;
the plan's values appear to have been rounded. This folder records the
recomputed values. Both pairs are within the 0.95–1.03 kg/m³ range I set as a
sanity check.

## Using a different configuration

The `BLADE_CONFIG_DIR` environment variable points the loader at a different
folder. That's the intended way to run a sensitivity case, such as a different
wind resource, without editing the tracked files.
