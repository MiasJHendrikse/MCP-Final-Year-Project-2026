"""
The only reader of the repo-root `config/` directory.

Everything else in `src/` -- solver, validation scripts, demos -- goes through
the four `load_*` functions here and never opens a YAML file itself. That is
what makes plan item 1.1's "nothing downstream hard-codes a site value"
checkable rather than aspirational: there is exactly one place a site value can
enter the code, and `grep` over `src/` for the old literals is the test.

Two things this module refuses to do:

  * substitute a value for an unresolved `TODO` field (see `unresolved.py`);
  * accept an internally inconsistent atmosphere or wind resource. `rho`, `p`
    and `nu`, and likewise `mean_wind_speed_ms`, are derived quantities in
    site.yaml, so they are recomputed from the recorded inputs and the file is
    rejected if the recorded values disagree. Editing the mean site
    temperature and forgetting to update the density -- or editing Weibull `k`
    and forgetting the mean -- is otherwise a silent few-percent error in
    every AEP figure.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
import os
from functools import lru_cache

import yaml

from config.schema import (
    DesignRotorConfig,
    ParameterisationConfig,
    PolarCacheConfig,
    SiteConfig,
    StandardAtmosphere,
    ValidationRotorConfig,
)
from config.unresolved import Unresolved

_HERE = os.path.dirname(os.path.abspath(__file__))

#: Repo root, i.e. the parent of `src/`.
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))

#: Versioned input configuration. Overridable by `BLADE_CONFIG_DIR` so a
#: sensitivity study can point at an alternative set without editing code --
#: the mechanism plan section 1.3's "central case plus a sensitivity band"
#: will need.
CONFIG_DIR = os.environ.get("BLADE_CONFIG_DIR",
                            os.path.join(REPO_ROOT, "config"))

#: Recorded-vs-recomputed agreement required of site.yaml's derived
#: atmosphere. Tight enough to catch an edit that was not propagated, loose
#: enough to allow the recorded values to be rounded for reporting.
ATMOSPHERE_TOLERANCE = 5e-3

#: Recorded-vs-recomputed agreement required of `mean_wind_speed_ms` against
#: `c * Gamma(1 + 1/k)`. Tighter than the atmosphere's, because there is no
#: reason to record a rounded mean: the resource is written by
#: `verification/wind_resource/run_extrapolation.py` at full precision, and a
#: mean that disagrees at the fourth decimal means k or c was edited without
#: the mean being recomputed.
WIND_RESOURCE_TOLERANCE = 1e-4


class ConfigError(RuntimeError):
    """A config file is missing, malformed, or internally inconsistent."""


def _read_yaml(filename):
    path = os.path.join(CONFIG_DIR, filename)
    if not os.path.exists(path):
        raise ConfigError(
            f"config file not found: {path}. The repo-root `config/` "
            f"directory is versioned input data; it is not generated."
        )
    with open(path, encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise ConfigError(f"{path} did not parse to a mapping")
    return data


def _get(mapping, path, filename):
    """
    Fetch a dotted key, converting a `TODO:` string into `Unresolved`.

    A missing key is an error, not a `None`: the schema is fixed, and a config
    file that has quietly lost a field should say so here rather than three
    modules downstream.
    """

    node = mapping
    keys = path.split(".")
    for i, key in enumerate(keys):
        if not isinstance(node, dict) or key not in node:
            reached = ".".join(keys[:i + 1])
            raise ConfigError(
                f"{filename}: missing required field '{path}' "
                f"(stopped at '{reached}')"
            )
        node = node[key]

    if isinstance(node, str) and node.strip().startswith("TODO"):
        note = node.strip()[len("TODO"):].lstrip(": ").strip()
        return Unresolved(f"{filename}:{path}", note or "no note given")
    return node


def _check_atmosphere(site_yaml, filename):
    """
    Recompute pressure, density and kinematic viscosity from the recorded
    inputs and reject the file if the recorded values disagree.

        p   = p0 * (1 - L*h/T0) ** (g / (Rs * L))
        rho = p / (Rs * T)
        nu  = mu / rho
    """

    atm = site_yaml["atmosphere"]
    std = atm["standard_atmosphere"]
    h = site_yaml["location"]["elevation_m"]

    p0 = std["sea_level_pressure_pa"]
    t0 = std["sea_level_temperature_k"]
    lapse = std["lapse_rate_k_per_m"]
    r_s = std["specific_gas_constant"]
    g = std["gravity"]

    pressure = p0 * (1.0 - lapse * h / t0) ** (g / (r_s * lapse))
    density = pressure / (r_s * atm["temperature_k"])
    nu = atm["dynamic_viscosity_pa_s"] / density

    for label, recorded, computed in (
        ("pressure_pa", atm["pressure_pa"], pressure),
        ("air_density", atm["air_density"], density),
        ("kinematic_viscosity", atm["kinematic_viscosity"], nu),
    ):
        if abs(recorded - computed) > ATMOSPHERE_TOLERANCE * abs(computed):
            raise ConfigError(
                f"{filename}: atmosphere.{label} = {recorded!r} disagrees with "
                f"{computed!r}, recomputed from elevation_m, temperature_k, "
                f"dynamic_viscosity_pa_s and standard_atmosphere (tolerance "
                f"{ATMOSPHERE_TOLERANCE:.1%}). Recompute the derived "
                f"atmosphere rather than editing one field in isolation."
            )


def _check_wind_resource(site_yaml, filename):
    """
    Recompute the mean wind speed from `k` and `c` and reject a disagreement.

        V_bar = c * Gamma(1 + 1/k)

    The three wind-resource fields are not independent, and the mean is the
    derived one. It is recorded anyway, per the plan's working convention that
    a number feeding a later step is written down -- but recording a derived
    number creates the opportunity for it to drift from its definition, so the
    definition is enforced here.

    Skipped entirely while any of the three is still `TODO`: an unresolved
    resource is a legitimate state, and `WeibullResource.from_config()` is what
    reports it. Only a *partially* filled triple is an error, and that is
    caught below -- filling k and c but leaving the mean as TODO would
    otherwise sail through and hand a sentinel to whatever reads it.
    """

    wind = site_yaml["wind_resource"]
    fields = ("weibull_k", "weibull_c_ms", "mean_wind_speed_ms")
    resolved = {
        name: wind[name] for name in fields
        if not (isinstance(wind[name], str)
                and wind[name].strip().startswith("TODO"))
    }

    if not resolved:
        return
    if len(resolved) != len(fields):
        missing = [name for name in fields if name not in resolved]
        raise ConfigError(
            f"{filename}: wind_resource is partially resolved -- "
            f"{', '.join(sorted(resolved))} have values but "
            f"{', '.join(missing)} still TODO. The three are one quantity; "
            f"resolve them together or leave all three TODO."
        )

    k = float(resolved["weibull_k"])
    c = float(resolved["weibull_c_ms"])
    recorded = float(resolved["mean_wind_speed_ms"])
    computed = c * math.gamma(1.0 + 1.0 / k)

    if abs(recorded - computed) > WIND_RESOURCE_TOLERANCE * abs(computed):
        raise ConfigError(
            f"{filename}: wind_resource.mean_wind_speed_ms = {recorded!r} "
            f"disagrees with {computed!r}, recomputed as "
            f"c * Gamma(1 + 1/k) from weibull_c_ms = {c!r} and weibull_k = "
            f"{k!r} (tolerance {WIND_RESOURCE_TOLERANCE:.1e}). The mean is "
            f"derived, not independent -- recompute it rather than editing k "
            f"or c in isolation."
        )


@lru_cache(maxsize=None)
def load_site(filename="site.yaml"):
    """The site basis (plan sections 1.1-1.3). Cached; the file is read once."""

    data = _read_yaml(filename)
    _check_atmosphere(data, filename)
    _check_wind_resource(data, filename)

    def field(path):
        return _get(data, path, filename)

    std = data["atmosphere"]["standard_atmosphere"]
    return SiteConfig(
        name=field("name"),
        description=field("description"),
        gwa_area=field("location.gwa_area"),
        latitude_deg=field("location.latitude_deg"),
        longitude_deg=field("location.longitude_deg"),
        elevation_m=field("location.elevation_m"),
        hub_height_m=field("location.hub_height_m"),
        temperature_k=field("atmosphere.temperature_k"),
        pressure_pa=field("atmosphere.pressure_pa"),
        air_density=field("atmosphere.air_density"),
        dynamic_viscosity_pa_s=field("atmosphere.dynamic_viscosity_pa_s"),
        kinematic_viscosity=field("atmosphere.kinematic_viscosity"),
        standard_atmosphere=StandardAtmosphere(
            sea_level_pressure_pa=std["sea_level_pressure_pa"],
            sea_level_temperature_k=std["sea_level_temperature_k"],
            lapse_rate_k_per_m=std["lapse_rate_k_per_m"],
            specific_gas_constant=std["specific_gas_constant"],
            gravity=std["gravity"],
        ),
        weibull_k=field("wind_resource.weibull_k"),
        weibull_c_ms=field("wind_resource.weibull_c_ms"),
        mean_wind_speed_ms=field("wind_resource.mean_wind_speed_ms"),
    )


@lru_cache(maxsize=None)
def load_design_rotor(filename="rotor_design.yaml"):
    """The Khomas Hochland design rotor (plan sections 1.4, 2.2, 4)."""

    data = _read_yaml(filename)

    def field(path):
        return _get(data, path, filename)

    return DesignRotorConfig(
        name=field("name"),
        description=field("description"),
        radius_m=field("geometry.radius_m"),
        n_blades=field("geometry.n_blades"),
        root_fraction=field("geometry.root_fraction"),
        design_tsr=field("operating.design_tsr"),
        rated_wind_speed_ms=field("operating.rated_wind_speed_ms"),
        rated_power_w=field("operating.rated_power_w"),
        max_rotor_speed_rpm=field("operating.max_rotor_speed_rpm"),
        cut_in_wind_speed_ms=field("operating.cut_in_wind_speed_ms"),
        cut_out_wind_speed_ms=field("operating.cut_out_wind_speed_ms"),
        parameterisation=ParameterisationConfig(
            n_bem_strips=field("parameterisation.n_bem_strips"),
            n_control_points_chord=field("parameterisation.n_control_points_chord"),
            n_control_points_twist=field("parameterisation.n_control_points_twist"),
            chord_min_m=field("parameterisation.bounds.chord_min_m"),
            chord_max_m=field("parameterisation.bounds.chord_max_m"),
            twist_min_deg=field("parameterisation.bounds.twist_min_deg"),
            twist_max_deg=field("parameterisation.bounds.twist_max_deg"),
        ),
        max_local_solidity=field("constraints.max_local_solidity"),
        aep_mwh_per_year_min=field("sanity.aep_mwh_per_year_min"),
        aep_mwh_per_year_max=field("sanity.aep_mwh_per_year_max"),
    )


@lru_cache(maxsize=None)
def load_phase_vi_rotor(filename="rotor_phase_vi.yaml"):
    """The NREL Phase VI validation rotor's operating condition (plan 2.1)."""

    data = _read_yaml(filename)

    def field(path):
        return _get(data, path, filename)

    return ValidationRotorConfig(
        name=field("name"),
        description=field("description"),
        air_density=field("atmosphere.air_density"),
        kinematic_viscosity=field("atmosphere.kinematic_viscosity"),
        rated_rpm=field("operating.rated_rpm"),
        tip_pitch_deg=field("operating.tip_pitch_deg"),
        reference_wind_speed_ms=field("operating.reference_wind_speed_ms"),
    )


@lru_cache(maxsize=None)
def load_polar_cache(name="s809"):
    """
    One polar cache's as-built metadata (plan section 3.3).

    Per-cache by design: S809 is at Ncrit=5 for a documented physical reason
    and SG6043 will be at whatever its own sensitivity study selects, so there
    is no global default to fall back to and this function takes the cache
    name rather than assuming one.
    """

    filename = f"polars_{name}.yaml"
    data = _read_yaml(filename)

    def field(path):
        return _get(data, path, filename)

    return PolarCacheConfig(
        name=field("name"),
        airfoil=field("airfoil"),
        directory=field("directory"),
        description=field("description"),
        ncrit=field("build.ncrit"),
        n_panel=field("build.n_panel"),
        alpha_min_deg=field("build.alpha_min_deg"),
        alpha_max_deg=field("build.alpha_max_deg"),
        alpha_step_deg=field("build.alpha_step_deg"),
        reynolds=tuple(field("build.reynolds")),
        known_gaps=tuple(
            (gap["reynolds"], gap["alpha_deg"]) for gap in field("known_gaps")
        ),
        history=tuple(
            (entry["date"], entry["change"]) for entry in field("history")
        ),
        validation=dict(field("validation")),
    )
