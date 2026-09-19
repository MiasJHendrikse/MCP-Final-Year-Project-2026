"""
The objective layer: from a design vector to AEP.

Plan step 1.5. The chain is
`design/` -> RotorGeometry -> `bem/` -> Cp -> P(V) -> weighted sum -> AEP,
and this package owns the last three links.

  `power`      wind-speed bins, the operating strategy, and P(V; d) with the
               rated-power limit applied.
  `weibull`    the wind-speed distribution, from `config/site.yaml`
               (`k = 1.709`, `c = 7.274 m/s` at 20 m since 2026-09-13).
  `objective`  AEP, J(d) = -AEP(d), and the unit-weighted surrogate the
               smoothness gate ran on while the resource was outstanding.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from objective.objective import (
    HOURS_PER_YEAR,
    annual_energy_mwh,
    bin_powers,
    energy_surrogate,
    objective,
    sanity_band,
)
from objective.power import BIN_WIDTH_MS, aerodynamic_power, power_per_bin, wind_speed_bins
from objective.weibull import WeibullResource

__all__ = [
    "BIN_WIDTH_MS",
    "HOURS_PER_YEAR",
    "WeibullResource",
    "aerodynamic_power",
    "annual_energy_mwh",
    "bin_powers",
    "energy_surrogate",
    "objective",
    "power_per_bin",
    "sanity_band",
    "wind_speed_bins",
]
