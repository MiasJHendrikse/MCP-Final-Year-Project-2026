"""
`BEMSystem`: the residual system `R(x; d) = 0` behind the objective, its
partials, and the discrete-adjoint gradient (Phase 3, B1-B3).

    state    x = phi_{b,i}     one inflow angle per (operating point, station)
                               17 x 25 = 425 scalars
    design   d in R^10         chord and twist control points
    R(x; d)  (17, 25)          bem.station.residual at every (b, i)
    J(x; d)                    -AEP [MWh/yr] assembled from the power
                               integrand q at the given state, plus the
                               constant the capped bins contribute

`dR/dx` is diagonal: each station's residual depends on its own `phi` only,
and `a`, `a'` are explicit functions of `(phi, c, theta)` rather than states
(`docs/adjoint_derivation.md` §1). The adjoint equation is therefore 425
scalar divisions,

    psi_{b,i} = -(dJ/dphi_{b,i}) / (dR_{b,i}/dphi_{b,i})
    dJ/dd     = dJ/dd|_explicit + sum_{b,i} psi_{b,i} dR_{b,i}/dd

Per-station chord and twist reach the design vector only through the
constant spline Jacobians `N_c`, `N_theta` (`BladeParameterisation`); nothing
here accepts per-station values as a design variable.

Operating points and the per-bin tip-speed ratio (2026-09-19)
--------------------------------------------------------------
The operating points are `objective.power.operating_points()`: the 17 bin
midpoints, each with its own `lambda_b = min(6.5, Omega_max R / V_b)` from
the configured rotor-speed ceiling. Everything lambda touches here reads
the per-point list `self.omega` (`Omega_b = lambda_b V_b / R`): the local
speed ratio `lam_r`, the zero-induction Reynolds estimate `W = hypot(V_b,
Omega_b r_i)`, and the power assembly `P_b = Omega_b sum_i t_i q_{b,i}`. The
schedule is a function of `V` alone, so `Omega_b` is a constant of the
station and nothing in the derivation changes: `dR/dx` stays diagonal,
`dR/dd` keeps its form, and `gradient()` assembles `dJ/dd = sum_b omega_b
dP_b/dd` as before. With no ceiling (`max_rotor_speed_rpm: null`) every
`lambda_b` is the design TSR and this is the pre-2026-09-19 system
bit-for-bit (`tests/test_operating_law_control.py`).

`points` may be passed explicitly as `[(V, lambda)]` for a system on other
operating points -- the load constraint's, for one -- with `mass=None`
meaning unit weights.

Forward solve
--------------
`solve(d)` calls `objective.power.aerodynamic_power` once per operating point
with that point's `lambda_b` and reads `phi` from the returned station
records. It does not re-implement the root-find and does not solve twice.
`limited = P_b > P_rated` exactly as `objective.power.power_per_bin` flags
it, with `P_rated` the configured generator rating
(`operating.rated_power_w`); that mask is a fixed input to every derivative
(§7 of the derivation).

The rating is fixed, so a capped bin contributes the constant `-(T/1e6) m_b
P_rated` to `J` and nothing to any derivative: its weight `omega_b` is zero
and its 25 states never enter the adjoint. Until 2026-09-19 the cap floated
with the design (`P_rated = P_aero(V_rated; d)`), the rated speed was an
eighteenth operating point and the capped mass was carried by its weight;
`docs/AEP_GAIN_AUDIT.md` §3.2 is why that changed.

Complex safety
---------------
`residual(phi, d)` and `J(phi, d)` accept complex `phi` and complex `d`, so
the Tier 1 tests can complex-step the code's own residual and the assembled
objective. No array is pre-allocated with a real dtype on those paths: rows
are built as lists and `np.array` infers the dtype.

`bounds` is taken as an argument so the scaling chain `span` is available
to `gradients.ScaledProblem.jac_adjoint`, and is read from nowhere else.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math
from dataclasses import dataclass

import numpy as np

from adjoint.kernels import station_partials
from bem.station import StationParams, residual as station_residual
from config import load_design_rotor, load_site
from objective.objective import HOURS_PER_YEAR
from objective.power import aerodynamic_power, operating_points, wind_speed_bins
from polars.polar import CachedPolar, interpolant_for


@dataclass(frozen=True)
class ForwardState:
    """
    The converged state of one forward solve at design `d`.

    `phi` is (n_points, n_stations), one row per bin midpoint; `power_w` is
    the unlimited power at every operating point in the same order;
    `limited` is the `power_per_bin` mask over the 17 bins.
    """

    d: np.ndarray
    phi: np.ndarray
    power_w: np.ndarray
    limited: np.ndarray
    a: np.ndarray
    reynolds: np.ndarray
    converged: bool


@dataclass(frozen=True)
class GradientResult:
    """`dJ/dd` and what it was assembled from."""

    dJ_dd: np.ndarray
    dJ_dd_explicit: np.ndarray
    psi: np.ndarray
    dR_dx: np.ndarray
    dJ_dx: np.ndarray
    state: ForwardState
    J: float


class BEMSystem:
    """
    The residual system and adjoint over the objective's 17 operating points.

    Parameters
    ----------
    parameterisation : design.BladeParameterisation
    bounds : design.DesignBounds
        Required; the scaling chain only (`DesignBounds.from_config()`).
    resource : objective.WeibullResource
    polar_cache : str
        The blade's polar cache; `"sg6043"`, the one default, matching
        `BladeParameterisation.to_geometry`.
    points : list of (float, float) or None
        Operating points `(V, lambda)`. `None` is the objective's,
        `objective.power.operating_points()`, with the bin masses from
        `resource`. Given explicitly, the masses are unit weights unless
        `mass` is passed, and `J` is then not the AEP -- subclasses on other
        operating points (the load constraint) define their own functional.
    mass : array-like or None
        Per-point weights to go with explicit `points`.
    """

    def __init__(self, parameterisation, bounds, resource, polar_cache="sg6043",
                 points=None, mass=None):
        if bounds.n_design_variables != parameterisation.n_design_variables:
            raise ValueError(
                f"bounds carry {bounds.n_design_variables} variables, the "
                f"parameterisation {parameterisation.n_design_variables}")

        self.parameterisation = parameterisation
        self.bounds = bounds
        self.resource = resource
        self.polar_cache = polar_cache

        self.design = load_design_rotor()
        self.site = load_site()
        self.interpolant = interpolant_for(polar_cache)

        self.tsr = float(self.design.design_tsr)
        self.air_density = float(self.site.air_density)
        self.kinematic_viscosity = float(self.site.kinematic_viscosity)
        self.R = float(parameterisation.radius_m)
        self.n_blades = int(self.design.n_blades)
        self.r_hub = float(parameterisation.root_fraction * parameterisation.radius_m)

        # Operating points: the 17 bin midpoints on the configured schedule,
        # each with its own tip-speed ratio. The rating is a number from
        # config, so nothing is solved at the rated wind speed.
        if points is None:
            edges, _midpoints, _width = wind_speed_bins()
            points = operating_points(self.design)
            mass = resource.probability_between(edges[:-1], edges[1:])
        elif mass is None:
            mass = np.ones(len(points))
        self.points = [(float(v), float(lam)) for v, lam in points]
        self.speeds = [v for v, _lam in self.points]
        self.tsr_per_bin = [lam for _v, lam in self.points]
        self.n_bins = len(self.points)
        self.n_points = len(self.speeds)
        self.p_rated = float(self.design.rated_power_w)
        self.mass = np.asarray(mass, dtype=float)
        if self.mass.shape != (self.n_points,):
            raise ValueError(f"mass has shape {self.mass.shape}, expected ({self.n_points},)")

        # Station constants: radii as `to_geometry` hands them to the solver
        # (float per station), the per-point rotor speed (the one list every
        # lambda-dependent quantity reads), the trapezoid weights.
        self.radii = [float(r) for r in parameterisation.radii]
        self.n_stations = len(self.radii)
        self.omega = [lam * v / self.R for v, lam in self.points]
        self.trapz_weights = _trapezoid_weights(self.radii)

        self.N_c = parameterisation.dchord_dd()
        self.N_theta = parameterisation.dtwist_dd()
        self.n_design = parameterisation.n_design_variables

    # -- per-station constants ----------------------------------------------

    def _w_approx(self, b, i):
        """`hypot(V, Omega r)`, in `solve_rotor`'s own operation order."""

        return math.hypot(self.speeds[b], self.omega[b] * self.radii[i])

    def reynolds(self, b, i, chord_i):
        """The zero-induction Reynolds estimate at `(b, i)` for chord `chord_i`."""

        return self._w_approx(b, i) * chord_i / self.kinematic_viscosity

    def _lam_r(self, b, i):
        return self.omega[b] * self.radii[i] / self.speeds[b]

    def _station(self, b, i, chord_i, twist_i):
        """`StationParams` for `(b, i)`, exactly as `solve_rotor` builds it."""

        polar = CachedPolar(self.interpolant, self.reynolds(b, i, chord_i))
        return StationParams(
            r=self.radii[i], chord=chord_i, twist=twist_i, airfoil=polar,
            tsr=self._lam_r(b, i), R=self.R, n_blades=self.n_blades, r_hub=self.r_hub)

    def chord_twist(self, d):
        """Per-station `(chord, twist)` from the design vector; dtype follows `d`."""

        d = np.asarray(d)
        chord = self.parameterisation.chord(d)
        twist = self.parameterisation.twist(d)
        return [x.item() for x in chord], [x.item() for x in twist]

    # -- forward solve --------------------------------------------------------

    def solve(self, d):
        """
        One forward solve: `phi` at every `(b, i)`, the unlimited powers, the
        `limited` mask. Raises if any station failed to converge -- the
        adjoint of an unconverged state means nothing.
        """

        d = np.asarray(d, dtype=float)
        geometry = self.parameterisation.to_geometry(d, polar_cache=self.polar_cache)

        phi, power, a, reynolds, converged = [], [], [], [], True
        for v, lam in self.points:
            p, result = aerodynamic_power(geometry, v, lam, self.air_density,
                                          self.kinematic_viscosity)
            stations = result["stations"]
            phi.append([s["phi"] for s in stations])
            a.append([s["a"] for s in stations])
            reynolds.append([s["reynolds"] for s in stations])
            power.append(float(p))
            converged = converged and bool(result["converged"])

        if not converged:
            raise RuntimeError(f"forward solve did not converge at d = {d.tolist()}")

        power = np.array(power)
        limited = power > self.p_rated
        return ForwardState(d=d, phi=np.array(phi), power_w=power, limited=limited,
                            a=np.array(a), reynolds=np.array(reynolds), converged=True)

    # -- the residual system ---------------------------------------------------

    def station_residual(self, b, i, phi_bi, d):
        """`bem.station.residual` at one `(b, i)`; complex-safe in `phi_bi`, `d`."""

        chord, twist = self.chord_twist(d)
        return station_residual(phi_bi, self._station(b, i, chord[i], twist[i]))

    def residual(self, phi, d):
        """
        `R(phi; d)`, shape (n_points, n_stations), through the code's own
        `bem.station.residual`. Complex-safe: the dtype follows the inputs.
        """

        phi = np.asarray(phi)
        chord, twist = self.chord_twist(d)
        rows = []
        for b in range(self.n_points):
            rows.append([station_residual(phi[b, i].item(),
                                          self._station(b, i, chord[i], twist[i]))
                         for i in range(self.n_stations)])
        return np.array(rows)

    def partials(self, phi, d, derivatives=True):
        """
        `adjoint.kernels.station_partials` at every `(b, i)`.

        Returns a dict of (n_points, n_stations) arrays keyed by the
        `StationPartials` field names, plus `buhl` as a boolean array. With
        `derivatives=False` only the value fields are present.
        """

        phi = np.asarray(phi)
        chord, twist = self.chord_twist(d)
        fields = ("residual", "q", "m", "a", "F", "dF_dphi", "alpha", "cl", "cd", "cn", "ct", "buhl")
        if derivatives:
            fields += ("dR_dphi", "dR_dc", "dR_dtheta", "dq_dphi", "dq_dc", "dq_dtheta",
                       "dm_dphi", "dm_dc", "dm_dtheta", "da_dphi", "da_dc", "da_dtheta")
        columns = {name: [] for name in fields}
        for b in range(self.n_points):
            row = {name: [] for name in fields}
            for i in range(self.n_stations):
                reynolds = self.reynolds(b, i, chord[i])
                polar = CachedPolar(self.interpolant, reynolds)
                out = station_partials(
                    phi[b, i].item(), chord[i], twist[i], reynolds, self.radii[i],
                    self._lam_r(b, i), self.R, self.n_blades, self.r_hub, polar,
                    self.speeds[b], self.air_density, derivatives=derivatives)
                for name in fields:
                    row[name].append(getattr(out, name))
            for name in fields:
                columns[name].append(row[name])
        return {name: np.array(values) for name, values in columns.items()}

    def dR_dx(self, phi, d):
        """The diagonal of `dR/dx`, shape (n_points, n_stations)."""

        return self.partials(phi, d)["dR_dphi"]

    def dR_dd(self, phi, d):
        """`dR/dd`, shape (n_points, n_stations, n_design), via `N_c`, `N_theta`."""

        parts = self.partials(phi, d)
        return self._assemble_dR_dd(parts)

    def _assemble_dR_dd(self, parts):
        return (parts["dR_dc"][:, :, None] * self.N_c[None, :, :]
                + parts["dR_dtheta"][:, :, None] * self.N_theta[None, :, :])

    # -- the objective from a state -------------------------------------------

    def powers_from_q(self, q):
        """`P_b = Omega_b sum_i t_i q_{b,i}` for every operating point; dtype follows `q`."""

        return np.array([self.omega[b] * sum(t * qi for t, qi in zip(self.trapz_weights, q[b]))
                         for b in range(self.n_points)])

    def limited_mask(self, power):
        """`P_b > P_rated` over the bins, decided on the real part."""

        return np.real(np.asarray(power)) > self.p_rated

    def weights(self, limited):
        """
        `omega_b`: the weight each operating point's power carries in `J`.

            omega_b = -(T/1e6) m_b        unlimited bin
            omega_b = 0                   limited bin (its own solve does not enter J)

        The limited bins' energy is the constant `J_capped(limited)`; it is
        not a weight on any state.
        """

        limited = np.asarray(limited, dtype=bool)
        factor = -HOURS_PER_YEAR / 1e6
        return np.where(limited, 0.0, factor * self.mass)

    def J_capped(self, limited):
        """`-(T/1e6) P_rated sum_{limited} m_b`: what the capped bins add to `J`."""

        limited = np.asarray(limited, dtype=bool)
        return -HOURS_PER_YEAR / 1e6 * self.p_rated * float(self.mass[limited].sum())

    def J(self, phi, d, limited=None):
        """
        `J = -AEP` [MWh/yr] assembled from `q` at the given state -- not
        re-solved. At the solved state this equals `objective.annual_energy_mwh`
        (tested to 1e-12 relative). Complex-safe.
        """

        parts = self.partials(phi, d, derivatives=False)
        power = self.powers_from_q(parts["q"])
        if limited is None:
            limited = self.limited_mask(power)
        omega = self.weights(limited)
        return (sum(omega[b] * power[b] for b in range(self.n_points))
                + self.J_capped(limited))

    def dJ_dx(self, phi, d, limited=None, parts=None):
        """`dJ/dphi_{b,i} = omega_b Omega_b t_i dq_{b,i}/dphi`, shape (n_points, n_stations)."""

        parts = self.partials(phi, d) if parts is None else parts
        omega = self.weights(self._mask(parts, limited))
        scale = omega[:, None] * np.asarray(self.omega)[:, None] * np.asarray(self.trapz_weights)[None, :]
        return scale * parts["dq_dphi"]

    def dJ_dd(self, phi, d, limited=None, parts=None):
        """The explicit `dJ/dd` (state held fixed), shape (n_design,)."""

        parts = self.partials(phi, d) if parts is None else parts
        omega = self.weights(self._mask(parts, limited))
        scale = omega[:, None] * np.asarray(self.omega)[:, None] * np.asarray(self.trapz_weights)[None, :]
        dJ_dc = (scale * parts["dq_dc"]).sum(axis=0)
        dJ_dtheta = (scale * parts["dq_dtheta"]).sum(axis=0)
        return self.N_c.T @ dJ_dc + self.N_theta.T @ dJ_dtheta

    def _mask(self, parts, limited):
        if limited is not None:
            return limited
        return self.limited_mask(self.powers_from_q(parts["q"]))

    def state_sensitivity(self, parts, j):
        """
        `dphi_{b,i}/dd_j` from the linearised system, shape (n_points, n_stations):
        `-(dR/dd_j) / (dR/dphi)`. The per-station alpha, Reynolds and
        induction sensitivities Tier 4 needs follow from it.
        """

        e = np.zeros(self.n_design)
        e[j] = 1.0
        return -self.apply_dR_dd(parts, e) / parts["dR_dphi"]

    # -- matrix-free operators (Tier 2) -----------------------------------------

    def apply_dR_dd(self, parts, v):
        """`(dR/dd) v` for `v` in R^n_design, shape (n_points, n_stations)."""

        v = np.asarray(v)
        dc = self.N_c @ v
        dtheta = self.N_theta @ v
        return parts["dR_dc"] * dc[None, :] + parts["dR_dtheta"] * dtheta[None, :]

    def apply_dR_dd_T(self, parts, w):
        """`(dR/dd)^T w` for `w` of shape (n_points, n_stations), shape (n_design,)."""

        w = np.asarray(w)
        by_c = (parts["dR_dc"] * w).sum(axis=0)
        by_theta = (parts["dR_dtheta"] * w).sum(axis=0)
        return self.N_c.T @ by_c + self.N_theta.T @ by_theta

    def tangent(self, phi, d, v, limited=None):
        """
        Forward-mode `dJ/dd . v`: perturb the design along `v`, propagate
        through the implicit state, sum into `J`. The independent code path
        the adjoint is tested against (Tier 2).
        """

        parts = self.partials(phi, d)
        dphi = -self.apply_dR_dd(parts, v) / parts["dR_dphi"]
        return (float(np.sum(self.dJ_dx(phi, d, limited, parts) * dphi))
                + float(self.dJ_dd(phi, d, limited, parts) @ np.asarray(v)))

    # -- the adjoint gradient (B3) -----------------------------------------------

    def gradient(self, d, state=None):
        """
        `dJ/dd` by the discrete adjoint at design `d`.

        One forward solve (or the one passed in as `state`), the partials at
        every station, 425 divisions for `psi`, one assembly.
        """

        state = self.solve(d) if state is None else state
        parts = self.partials(state.phi, state.d)
        limited = state.limited

        dR_dx = parts["dR_dphi"]
        dJ_dx = self.dJ_dx(state.phi, state.d, limited, parts)
        psi = -dJ_dx / dR_dx
        explicit = self.dJ_dd(state.phi, state.d, limited, parts)
        total = explicit + self.apply_dR_dd_T(parts, psi)

        power = self.powers_from_q(parts["q"])
        J = float(np.sum(self.weights(limited) * power)) + self.J_capped(limited)
        return GradientResult(dJ_dd=total, dJ_dd_explicit=explicit, psi=psi,
                              dR_dx=dR_dx, dJ_dx=dJ_dx, state=state, J=J)


def _trapezoid_weights(x):
    """
    `t_i` such that `sum_i t_i y_i == trapz(y, x)` for `bem.rotor._trapz`:
    `t_0 = (x_1 - x_0)/2`, `t_i = (x_{i+1} - x_{i-1})/2`, `t_{n-1} = (x_{n-1} - x_{n-2})/2`.
    """

    n = len(x)
    weights = [0.0] * n
    for i in range(n - 1):
        half = 0.5 * (x[i + 1] - x[i])
        weights[i] += half
        weights[i + 1] += half
    return weights
