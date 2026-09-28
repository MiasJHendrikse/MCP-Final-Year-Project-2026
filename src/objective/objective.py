"""
The objective: AEP, and the surrogate the smoothness gate runs on.

Plan section 1.5 defines the optimisation objective as

    J(d) = -AEP(d) = -T * integral over V of P(V; d) * f_Weibull(V; k, c) dV

with T = 8766 hours (the mean Julian year, 365.25 days -- not 8760, which
drops a quarter-day per year and biases every AEP low by about 0.07 %).

Two entry points, and the difference between them matters
-----------------------------------------------------------
`annual_energy_mwh(d)` is the real thing. It needs the Weibull parameters,
which are `TODO`, so it raises today.

`energy_surrogate(d)` is the same sum with **unit weights instead of Weibull
weights**. It exists for one purpose: the smoothness gate, which
asks whether `J` is smooth enough to differentiate, must be able to run before
the wind resource arrives.

The substitution is defensible for that purpose and for no other. The Weibull
weights are a **fixed convex combination over the bins, independent of d** --
they do not depend on the blade at all. So the surrogate and the real objective
share the entire d-dependent chain:

    d -> control points -> chord/twist -> RotorGeometry -> BEM -> Cp -> P(V; d)

Any staircasing, kink, noise or discontinuous jump the gate is looking for
lives in that chain and appears identically in both. What differs is only how
the per-bin powers are weighted before being summed, and a positive fixed
weighting cannot create or remove a discontinuity.

What the surrogate is NOT is an estimate of AEP. It is not scaled to energy
units, it is not comparable to the 8-12 MWh/yr sanity band, and it must never be
quoted as a performance figure. `annual_energy_mwh` is the only thing that may
be, and it raises until the data lands.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import numpy as np

from config import load_design_rotor
from design.parameterisation import BladeParameterisation
from objective.power import power_per_bin, wind_speed_bins
from objective.weibull import WeibullResource

#: Hours in a mean Julian year (365.25 days). Not 8760: that drops a quarter
#: day per year and biases every AEP low by about 0.07 %, which is small but is
#: a bias rather than noise, and it is free to get right.
HOURS_PER_YEAR = 8766.0


def _geometry_for(design_vector, parameterisation, polar_cache):
    parameterisation = parameterisation or BladeParameterisation()
    return parameterisation.to_geometry(design_vector, polar_cache=polar_cache)


def bin_powers(design_vector, parameterisation=None, polar_cache="sg6043",
               tip_speed_ratio=None):
    """
    Per-bin limited power for a design vector -- the shared d-dependent part.

    Both `annual_energy_mwh` and `energy_surrogate` go through this, so they
    cannot drift apart on the physics.
    """

    geometry = _geometry_for(design_vector, parameterisation, polar_cache)
    return power_per_bin(geometry, tip_speed_ratio=tip_speed_ratio)


def annual_energy_mwh(design_vector, resource=None, parameterisation=None,
                      polar_cache="sg6043", tip_speed_ratio=None):
    """
    AEP in MWh per year.

    Parameters
    ----------
    resource : WeibullResource or None
        Defaults to `WeibullResource.from_config()`, which raises while the
        site parameters are `TODO`. Pass one explicitly for a study that states
        its own provisional values.

    Raises
    ------
    config.unresolved.UnresolvedConfigError
        If `resource` is None and the site wind resource is unresolved.
    """

    resource = resource or WeibullResource.from_config()
    edges, _midpoints, _width = wind_speed_bins()
    powers = bin_powers(design_vector, parameterisation, polar_cache,
                        tip_speed_ratio)

    # Exact bin masses from the CDF -- see WeibullResource.probability_between.
    mass = resource.probability_between(edges[:-1], edges[1:])
    energy_wh = float(np.sum(powers["power_w"] * mass) * HOURS_PER_YEAR)

    return energy_wh / 1e6


def energy_surrogate(design_vector, parameterisation=None,
                     polar_cache="sg6043", tip_speed_ratio=None):
    """
    Unit-weighted sum of per-bin power, in watts.

    The smoothness gate's stand-in for J while the wind resource is
    outstanding. See the module docstring for exactly why this is a valid
    substitute for *smoothness* testing and an invalid one for anything else.

    Returns
    -------
    float
        Sum over bins of the limited power. Not an energy, not an AEP, and not
        comparable to the sanity band.
    """

    powers = bin_powers(design_vector, parameterisation, polar_cache,
                        tip_speed_ratio)
    return float(np.sum(powers["power_w"]))


def objective(design_vector, resource=None, parameterisation=None,
              polar_cache="sg6043", tip_speed_ratio=None):
    """
    J(d) = -AEP(d), the quantity an optimiser minimises.

    Negated because the problem is stated as a minimisation and SciPy's
    optimisers minimise; the sign lives here, once, rather than at every call
    site.
    """

    return -annual_energy_mwh(design_vector, resource, parameterisation,
                              polar_cache, tip_speed_ratio)


def sanity_band():
    """The configured plausibility band for AEP, MWh/yr."""

    design = load_design_rotor()
    return design.aep_mwh_per_year_min, design.aep_mwh_per_year_max
