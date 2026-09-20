"""
Frozen dataclasses for the four files under `config/`.

Frozen because a config object handed to a solver and then mutated by one
caller is indistinguishable from a correct run -- exactly the class of silent
wrongness this whole layer exists to remove. Field names match the YAML keys.

Fields that may still be `TODO` are annotated `float | Unresolved` (or
`str | Unresolved`); see `unresolved.py` for what happens if one is used.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
from dataclasses import dataclass

# Referenced by the string annotations below (never evaluated at runtime,
# but the name should resolve for anyone reading or type-checking them).
from config.unresolved import Unresolved  # noqa: F401


@dataclass(frozen=True)
class StandardAtmosphere:
    """ISA constants used by the barometric relation (plan section 1.2)."""

    sea_level_pressure_pa: float
    sea_level_temperature_k: float
    lapse_rate_k_per_m: float
    specific_gas_constant: float
    gravity: float


@dataclass(frozen=True)
class SiteConfig:
    """
    Everything downstream needs to know about the site.

    `air_density` and `kinematic_viscosity` are the two numbers that used to
    be hard-coded; they are the site's, not a default, and the design rotor
    takes both from here.
    """

    name: str
    description: str

    gwa_area: "str | Unresolved"
    latitude_deg: "float | Unresolved"
    longitude_deg: "float | Unresolved"
    elevation_m: float
    hub_height_m: float

    temperature_k: float
    pressure_pa: float
    air_density: float
    dynamic_viscosity_pa_s: float
    kinematic_viscosity: float
    standard_atmosphere: StandardAtmosphere

    weibull_k: "float | Unresolved"
    weibull_c_ms: "float | Unresolved"
    mean_wind_speed_ms: "float | Unresolved"


@dataclass(frozen=True)
class ParameterisationConfig:
    """
    Plan section 4. Strip count and design-variable count are independent
    quantities and are both carried here so they cannot be conflated.

    The bounds were `float | Unresolved` until 2026-09-19; they are resolved
    (config/rotor_design.yaml records the basis) and typed as such. The
    loader still converts a `TODO` string to `Unresolved`, so re-opening one
    is a matter of editing the YAML, and `DesignBounds.from_config()` would
    raise again.
    """

    n_bem_strips: int
    n_control_points_chord: int
    n_control_points_twist: int
    chord_min_m: float
    chord_max_m: float
    twist_min_deg: float
    twist_max_deg: float

    @property
    def n_design_variables(self):
        return self.n_control_points_chord + self.n_control_points_twist


@dataclass(frozen=True)
class DesignRotorConfig:
    """
    The optimisation subject (plan sections 1.4, 2.2). Deliberately carries no
    air properties -- those are the site's, from `SiteConfig`.

    `rated_power_w` is the generator rating the objective holds power at
    above rated. It is a property of the machine, not of the blade: the
    objective must never recompute it from the design being evaluated. Its
    current value is provisional (the baseline's own aerodynamic power at the
    rated wind speed, pending the nameplate -- outstanding input B2) and the
    YAML says so.

    `max_rotor_speed_rpm` is the rotor-speed ceiling (outstanding input B1,
    resolved provisionally 2026-09-19 at 300 rpm) that sets the per-bin
    tip-speed ratio `lambda(V) = min(design_tsr, Omega_max R / V)`; `None`
    is no ceiling, the pre-2026-09-19 objective. `max_local_solidity` is the
    station solidity cap (plan 7.4), resolved the same day.

    The mass problem (2026-09-20, `docs/PLAN-mass-objective-2026-09-20.md`):
    `mass_model` selects the material proxy (`"shell"`, the objective, or
    `"solid"`, reported only); `monotone_chord` / `monotone_twist` switch the
    manufacturability rows and `min_chord_m` is the buildable-tip floor on
    the chord control points (a row of the mass problem, not the box bound,
    so the Phase 1-4 scaling is untouched); `laminate_density_kg_m3` and `shell_thickness_m`
    are `Unresolved` until a laminate concept exists and are needed only to
    quote a mass in kg -- every constraint row is relative and cancels them.
    """

    name: str
    description: str
    radius_m: float
    n_blades: int
    root_fraction: float
    design_tsr: float
    rated_wind_speed_ms: float
    rated_power_w: float
    max_rotor_speed_rpm: "float | None"
    cut_in_wind_speed_ms: float
    cut_out_wind_speed_ms: float
    parameterisation: ParameterisationConfig
    max_local_solidity: float
    aep_mwh_per_year_min: float
    aep_mwh_per_year_max: float
    mass_model: str
    monotone_chord: bool
    monotone_twist: bool
    min_chord_m: float
    laminate_density_kg_m3: "float | Unresolved"
    shell_thickness_m: "float | Unresolved"

    @property
    def max_tip_speed_ms(self):
        """
        `Omega_max R`, m/s, or `None` with no ceiling. Derived, so the tip
        speed and the rpm cannot drift apart.
        """

        if self.max_rotor_speed_rpm is None:
            return None
        return float(self.max_rotor_speed_rpm) * 2.0 * math.pi / 60.0 * float(self.radius_m)


@dataclass(frozen=True)
class ValidationRotorConfig:
    """
    The NREL Phase VI solver-validation case (plan section 2.1), including its
    own sea-level air condition -- a choice for this rotor, not a default.
    """

    name: str
    description: str
    air_density: float
    kinematic_viscosity: float
    rated_rpm: float
    tip_pitch_deg: float
    reference_wind_speed_ms: float

    @property
    def dynamic_viscosity(self):
        """mu = rho * nu. Derived, so it cannot drift from the pair above."""

        return self.air_density * self.kinematic_viscosity


@dataclass(frozen=True)
class PolarCacheConfig:
    """
    One polar cache as built (plan section 3.3). Descriptive metadata, not a
    rebuild instruction, and per-cache by design: the two caches have
    different build settings and there is no global default.
    """

    name: str
    airfoil: str
    directory: str
    description: str
    ncrit: float
    n_panel: int
    alpha_min_deg: float
    alpha_max_deg: float
    alpha_step_deg: float
    reynolds: tuple
    known_gaps: tuple
    history: tuple
    validation: dict
