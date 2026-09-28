"""
The forward load path, independent of the adjoint: flapwise root bending
moment per blade, the load operating set `L`, and the KS aggregate over it.

The adjoint's own `RootMomentSystem` (`adjoint/loads.py`) is checked against
this module: `RootMomentSystem.moments(phi, d)` at the solved state must equal
`root_moments` here to round-off, because both read the same physics
(`station["w"]`, `station["Cn"]` in the solver are the kernel's `w`, `cn`).

The mass problem (2026-09-20) adds the deflection twin on the same station
records: `normal_load`, `spanwise_moments`, `tip_deflection`, checked to
round-off against `adjoint/deflection.py` in `tests/test_deflection.py`.

The load set `L` is the operating points at or below rated, fixed at
construction from the committed baseline `x0` and never recomputed per design
(a set that changed with `d` would be a non-smooth constraint):

    L = { (V_b, lambda_b) : bin midpoint b with P_b(x0) <= P_rated }
        union { (V_rated, lambda(V_rated)) }
      = 8 uncapped bins (3.5..10.5 m/s) + the rated point (11 m/s,
        lambda = 5.711986642890533 at 300 rpm)  ->  9 points.

Import direction: this module must not top-level import `objective.objective`
or `design.baseline` (that would close the cycle
`objective.objective -> design.parameterisation -> design.baseline ->
objective.loads -> objective.objective`). The baseline `x0` and its per-bin
powers are imported lazily inside `load_operating_points`.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import numpy as np

from bem.rotor import solve_rotor
from config import load_design_rotor
from objective.power import operating_points, tsr_schedule


def load_operating_points(design=None):
    """
    The load operating set `L`, in §3.1's order, rated point last.

    The uncapped bins are those whose *committed baseline* per-bin power is
    at or below the rating; the rated point is appended on the configured
    schedule. The set is a function of `x0` and the config, not of a design.
    """

    design = design or load_design_rotor()

    # Lazy: closing the import cycle with the design layer (see module docstring).
    from design.baseline import build_schmitz_baseline
    from objective.objective import bin_powers

    x0 = build_schmitz_baseline().design_vector
    powers = bin_powers(x0)
    limited = np.asarray(powers["limited"], dtype=bool)

    points = operating_points(design)
    selected = [points[b] for b in range(len(points)) if not limited[b]]

    rated_tsr = tsr_schedule(float(design.rated_wind_speed_ms),
                             design.design_tsr, design.max_tip_speed_ms)
    selected.append((float(design.rated_wind_speed_ms), rated_tsr))
    return selected


def root_bending_moment(stations, chords, air_density, n_blades, r_hub):
    """
    Flapwise root bending moment for one blade, N.m.

    M = integral of dT/dr * (r - r_hub) dr, with dT/dr the same normal-force
    integrand `solve_rotor` uses for thrust. Per *blade*, not per rotor: the
    root attachment carries one blade's load, and quoting a rotor-summed figure
    here would overstate it by a factor of B.
    """

    radii = [station["r"] for station in stations]
    load = [p * (r - r_hub) for p, r in zip(normal_load(stations, chords, air_density), radii)]
    return sum(0.5 * (load[i] + load[i + 1]) * (radii[i + 1] - radii[i])
               for i in range(len(radii) - 1))


def normal_load(stations, chords, air_density):
    """
    Flapwise (normal) aerodynamic load per unit span for one blade, N/m:
    `1/2 rho w^2 c Cn` at every station -- `root_bending_moment`'s integrand
    without its lever arm.
    """

    return [0.5 * air_density * station["w"] ** 2 * chord * station["Cn"]
            for station, chord in zip(stations, chords)]


def spanwise_moments(stations, chords, air_density):
    """
    Flapwise bending moment at every station for one blade, N.m:
    `M(r_k) = int_{r_k}^{R} q(r) (r - r_k) dr`, trapezoid on the stations
    from `k` outward (the `i = k` term has zero arm). The mass problem's
    deflection twin (2026-09-20); with the arm taken to `r_hub` instead the
    same rule is `root_bending_moment`.
    """

    radii = [station["r"] for station in stations]
    load = normal_load(stations, chords, air_density)
    moments = []
    for k in range(len(radii)):
        f = [q * (r - radii[k]) for q, r in zip(load, radii)]
        moments.append(sum(0.5 * (f[j] + f[j + 1]) * (radii[j + 1] - radii[j])
                           for j in range(k, len(radii) - 1)))
    return moments


def tip_deflection(stations, chords, air_density, radius_m):
    """
    Static flapwise tip deflection per unit `E k_I t_shell` for one blade:
    the unit-load (Euler-Bernoulli) integral `int M(r) (R - r) / c(r)^3 dr`
    over the stations, with `EI(r) = E k_I t_shell c(r)^3` for a thin shell
    of constant laminate thickness. Relative use only -- the common factor
    is never applied here; see `adjoint/deflection.py`.
    """

    radii = [station["r"] for station in stations]
    moments = spanwise_moments(stations, chords, air_density)
    f = [m * (radius_m - r) / chord ** 3 for m, r, chord in zip(moments, radii, chords)]
    return sum(0.5 * (f[i] + f[i + 1]) * (radii[i + 1] - radii[i])
               for i in range(len(radii) - 1))


def tip_deflection_at(geometry, v_inf, tsr, air_density, kinematic_viscosity):
    """
    `(delta, result)` at one solved operating point, as `root_moment` returns
    `(M, result)`.
    """

    result = solve_rotor(geometry, tsr=tsr, v_inf=v_inf,
                         air_density=air_density,
                         kinematic_viscosity=kinematic_viscosity)
    delta = tip_deflection(result["stations"], geometry.chord, air_density, geometry.R)
    return delta, result


def root_moment(geometry, v_inf, tsr, air_density, kinematic_viscosity):
    """
    Flapwise root bending moment per blade at one solved operating point.

    Returns `(M_nm, result)` where `result` is `solve_rotor`'s dict, so the
    caller keeps the solver status and every station field.
    """

    result = solve_rotor(geometry, tsr=tsr, v_inf=v_inf,
                         air_density=air_density,
                         kinematic_viscosity=kinematic_viscosity)
    moment = root_bending_moment(result["stations"], geometry.chord,
                                 air_density, geometry.n_blades, geometry.r_hub)
    return moment, result


def root_moments(geometry, points, air_density, kinematic_viscosity):
    """
    Per-point root bending moment over `points = [(V, lambda)]`.

    Returns `(np.ndarray of shape (n_points,), list[result])`.
    """

    moments, results = [], []
    for v, lam in points:
        moment, result = root_moment(geometry, v, lam, air_density,
                                     kinematic_viscosity)
        moments.append(moment)
        results.append(result)
    return np.array(moments), results


def ks(values, rho):
    """
    The KS aggregate `max + ln(sum exp(rho (v - max))) / rho` over `values`.

    Complex-safe: the max is taken on the real part, so a complex step of the
    underlying values propagates through the softmax (see `ks_weights`).
    """

    values = np.asarray(values)
    peak = np.max(values.real)
    shifted = values - peak
    return peak + np.log(np.sum(np.exp(rho * shifted))) / rho


def ks_weights(values, rho):
    """The analytic softmax weights `exp(rho (v - max)) / sum exp(...)`, sum to 1."""

    values = np.asarray(values)
    peak = np.max(values.real)
    shifted = values - peak
    exps = np.exp(rho * shifted)
    return exps / np.sum(exps)
