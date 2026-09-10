"""
Wind-speed bins and the operating strategy: P(V; d) (plan step 1.5).

Everything here is independent of the wind resource. The Weibull parameters
are still `TODO` (see `docs/OUTSTANDING-INPUTS.md`), but the bin scheme, the
operating strategy and the per-bin power are not blocked by that, and step
1.8's smoothness gate needs them.

Bin scheme, fixed and documented
---------------------------------
Cut-in 3 m/s to cut-out 20 m/s in **1 m/s bins** -- 17 bins, edges at the
integers, evaluated at the **midpoints** 3.5, 4.5 ... 19.5.

Midpoint rule rather than trapezoid, for a specific reason: the trapezoid rule
evaluates at the bin *edges*, so adjacent bins share an evaluation point and
the rated-power cut at exactly 11 m/s would land on a shared node. The midpoint
rule puts every evaluation strictly inside its own bin, so the cut falls
between nodes and no single evaluation has to be simultaneously below and above
rated. It is also second-order for a smooth integrand, which this is.

The bin count and the quadrature are properties of the *objective*, not tuning
knobs: changing either changes J, so they are constants here rather than
arguments, and any study that varies them has to say so.

Operating strategy (plan step 1.5)
-----------------------------------
Below rated: fixed tip-speed ratio at the design value, lambda = 6.5. Rotor
speed tracks the wind.

Above rated: **simple power limiting**. Power is held at the value the rotor
produces at the rated wind speed, i.e.

    P(V; d) = min( P_aero(V; d), P_aero(V_rated; d) )

Two things about this are worth stating because they matter to the smoothness
gate. It puts a kink in P against *V* -- deliberately, that is what limiting
is. It does **not** put a kink in P against *d* for the bins above rated,
because every one of them takes the value P_aero(V_rated; d), which is smooth
in d. A kink in d appears only if a *below*-rated bin crosses the limit as d
moves, and that is a real feature of the objective the gate should reveal
rather than something to be smoothed away in advance.

No generator, gearbox or electrical efficiency is applied. This is aerodynamic
rotor power, as `bem.powercurve` documents; a drivetrain model would multiply
through and belongs with the AEP write-up rather than here.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

import numpy as np

from bem.rotor import solve_rotor
from config import load_design_rotor, load_site

#: Bin width, m/s. See the module docstring on why this is a constant.
BIN_WIDTH_MS = 1.0


def wind_speed_bins():
    """
    (edges, midpoints, width) for the operating range.

    Returns
    -------
    (ndarray, ndarray, float)
        Edges from cut-in to cut-out inclusive, the midpoints power is
        evaluated at, and the bin width.
    """

    design = load_design_rotor()
    edges = np.arange(design.cut_in_wind_speed_ms,
                      design.cut_out_wind_speed_ms + 1e-9, BIN_WIDTH_MS)
    return edges, 0.5 * (edges[:-1] + edges[1:]), BIN_WIDTH_MS


def aerodynamic_power(geometry, v_inf, tip_speed_ratio, air_density,
                      kinematic_viscosity):
    """
    Unlimited aerodynamic rotor power at one wind speed, watts.

    Returns `(power, result)` so the caller keeps the solver status -- a
    non-converged operating point must stay visible rather than being reduced
    to a number.
    """

    result = solve_rotor(geometry, tsr=tip_speed_ratio, v_inf=v_inf,
                         air_density=air_density,
                         kinematic_viscosity=kinematic_viscosity)
    area = math.pi * geometry.R ** 2
    return result["Cp"] * 0.5 * air_density * v_inf ** 3 * area, result


def power_per_bin(geometry, tip_speed_ratio=None):
    """
    P(V) at every bin midpoint, with the rated-power limit applied.

    Returns
    -------
    dict
        `midpoints`, `power_w` (limited), `power_unlimited_w`, `rated_power_w`,
        `limited` (bool per bin), `converged` (bool per bin), and
        `all_converged`.

    Notes
    -----
    The rated power is computed by solving at the rated wind speed rather than
    by taking the largest below-rated bin: the limit is a property of the
    machine at its rated point, and reading it off the bin grid would make the
    objective depend on where the bins happen to fall.
    """

    design = load_design_rotor()
    site = load_site()
    tip_speed_ratio = design.design_tsr if tip_speed_ratio is None else tip_speed_ratio
    air = (site.air_density, site.kinematic_viscosity)

    rated_power, rated_result = aerodynamic_power(
        geometry, design.rated_wind_speed_ms, tip_speed_ratio, *air)

    _edges, midpoints, _width = wind_speed_bins()

    unlimited, converged = [], [rated_result["converged"]]
    for v_inf in midpoints:
        power, result = aerodynamic_power(geometry, float(v_inf),
                                          tip_speed_ratio, *air)
        unlimited.append(power)
        converged.append(result["converged"])

    unlimited = np.array(unlimited)
    limited = np.minimum(unlimited, rated_power)

    return {
        "midpoints": midpoints,
        "power_w": limited,
        "power_unlimited_w": unlimited,
        "rated_power_w": rated_power,
        "limited": unlimited > rated_power,
        "converged": converged,
        "all_converged": all(converged),
    }
