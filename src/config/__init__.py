"""
Versioned configuration: the single place site and rotor values enter the code.

Before this package existed, air density was a default argument in four solver
signatures and kinematic viscosity a module constant in two more, both at
sea-level values; the design site is at 1800 m. An AEP call that omitted the
density argument returned a figure 25 % high and still landed inside the
plausibility band, and every design-rotor Reynolds number would have come out
23 % high. A default there is worse than no default -- which is why
`solve_rotor` and the `powercurve` entry points now take rho and nu as required
arguments, threaded from here.

The corollary is that no air property appears as a literal anywhere under
`src/`. That is a checkable invariant, not a style preference: grep for one and
it should not be there.

Usage, from anywhere under `src/`:

    from config import load_phase_vi_rotor

    phase_vi = load_phase_vi_rotor()
    solve_rotor(geometry, tsr=6.0, v_inf=7.0,
                air_density=phase_vi.air_density,
                kinematic_viscosity=phase_vi.kinematic_viscosity)

Fields that are still `TODO` -- the Global Wind Atlas extraction, the
design-variable bounds -- come back as `Unresolved` objects that raise on any
use rather than as placeholder numbers. See `unresolved.py`.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from config.loader import (
    ATMOSPHERE_TOLERANCE,
    CONFIG_DIR,
    REPO_ROOT,
    ConfigError,
    load_design_rotor,
    load_phase_vi_rotor,
    load_polar_cache,
    load_site,
)
from config.schema import (
    DesignRotorConfig,
    ParameterisationConfig,
    PolarCacheConfig,
    SiteConfig,
    StandardAtmosphere,
    ValidationRotorConfig,
)
from config.unresolved import Unresolved, UnresolvedConfigError, is_resolved

__all__ = [
    "ATMOSPHERE_TOLERANCE",
    "CONFIG_DIR",
    "ConfigError",
    "DesignRotorConfig",
    "ParameterisationConfig",
    "PolarCacheConfig",
    "REPO_ROOT",
    "SiteConfig",
    "StandardAtmosphere",
    "Unresolved",
    "UnresolvedConfigError",
    "ValidationRotorConfig",
    "is_resolved",
    "load_design_rotor",
    "load_phase_vi_rotor",
    "load_polar_cache",
    "load_site",
]
