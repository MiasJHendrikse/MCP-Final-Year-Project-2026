"""
`DeflectionSystem`: the third right-hand side -- the KS aggregate of the
static flapwise tip deflection over the load operating set `L`, and its
discrete adjoint (`docs/adjoint_derivation.md` §11).

It subclasses `RootMomentSystem` and inherits everything: the forward solve,
the residual, the station partials, the matrix-free operators, the load set
`L`, the KS machinery. The only new physics is a fixed linear map from the
station loads to the spanwise moment and a chord-dependent flexibility:

    q_{b,i}   = m_{b,i} / arm_i                    arm-free normal load per blade, N/m
                                                   (m is the kernel's root-moment
                                                   integrand; arm_i = r_i - r_hub > 0)
    M_{b,k}   = sum_i W_{k,i} q_{b,i}              spanwise flapwise moment at station k
                W_{k,i} = t_i (r_i - r_k), i > k;  0 otherwise      (constant)
    delta_b   = sum_k G_k M_{b,k} / c_k^3          unit-load (Euler-Bernoulli) tip
                G_k = t_k (R - r_k)                deflection with EI_k = c_k^3
    D         = KS_rho(delta / delta_ref)          the functional (the analogue of KS)

`EI_k = E k_I t_shell c_k^3` for a thin shell of constant laminate thickness;
`E`, `k_I`, `t_shell` are one common factor that cancels in the relative
constraint `D(u) <= D(x0)` and is set to one here -- `delta_b` is a
deflection per unit `E k_I t_shell`, never a metre. `delta_ref` is the rated
point's value at `x0`, passed in by the caller exactly as `m_ref_nm` is.

The transposed form `delta_b = sum_i A_i(c) q_{b,i}` with the
flexibility-weighted arm `A_i = sum_k G_k W_{k,i} / c_k^3` is what the state
derivative reads; the chord enters explicitly through `q` (as for `m`) and
through the stiffness, `d delta_b / d c_j |_stiffness = -3 G_j M_{b,j} / c_j^4`.
The adjoint is the parent's diagonal solve with this right-hand side:

    psi_{b,i} = -(dD/dphi_{b,i}) / (dR_{b,i}/dphi_{b,i})
    dD/dd     = dD/dd|_explicit + sum_{b,i} psi_{b,i} dR_{b,i}/dd

No kernel change: `q` and its partials are `m` and its partials divided by
the station constant `arm_i`, exact. The segment `[r_hub, r_0]` (34 mm) is
outside the station grid for the deflection exactly as it is for the root
moment. `gradient` and `tangent` here are the deflection's; the parent's
moment methods (`moments`, `KS`, `dKS_dx`, ...) remain callable on the same
instance and the same state.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

from dataclasses import dataclass

import numpy as np

from adjoint.loads import RootMomentSystem
from adjoint.system import ForwardState
from objective.loads import ks, ks_weights


@dataclass(frozen=True)
class DeflectionGradientResult:
    """`dD/dd` and what it was assembled from (the `MomentGradientResult` analogue)."""

    dD_dd: np.ndarray
    dD_dd_explicit: np.ndarray
    psi: np.ndarray
    dR_dx: np.ndarray
    dD_dx: np.ndarray
    state: ForwardState
    D: float
    deflections: np.ndarray
    weights: np.ndarray
    spanwise_moments_nm: np.ndarray


class DeflectionSystem(RootMomentSystem):
    """
    The tip-deflection KS functional and its adjoint over the load points.

    Parameters
    ----------
    parameterisation, bounds, resource, polar_cache, points, rho, m_ref_nm :
        as `RootMomentSystem` (`m_ref_nm` is still required: the parent's
        moment functional stays available on this instance).
    delta_ref : float
        The normalising deflection (per unit `E k_I t_shell`), the rated
        point's value at the committed `x0`; passed in by `ScaledProblem`
        so every artefact records one reference.
    """

    def __init__(self, parameterisation, bounds, resource, polar_cache="sg6043",
                 points=None, rho=100.0, m_ref_nm=None, delta_ref=None):
        super().__init__(parameterisation, bounds, resource, polar_cache=polar_cache,
                         points=points, rho=rho, m_ref_nm=m_ref_nm)
        if delta_ref is None:
            raise ValueError("delta_ref is required; pass the committed x0 rated-point deflection")
        self.delta_ref = float(delta_ref)

        radii = np.asarray(self.radii, dtype=float)
        weights = np.asarray(self.trapz_weights, dtype=float)
        self.arm = radii - self.r_hub
        if np.any(self.arm <= 0.0):
            raise ValueError("every station must lie outboard of the hub for q = m / arm")
        # W[k, i] = t_i (r_i - r_k) for i > k: the trapezoid rule on [r_k, R]
        # has the full-grid interior weights for i > k and a zero-arm term at
        # i = k, so the sub-grid and full-grid weights coincide.
        lever = radii[None, :] - radii[:, None]
        self.W = np.where(lever > 0.0, weights[None, :] * lever, 0.0)
        self.G = weights * (self.R - radii)

    # -- the functional --------------------------------------------------------

    def loads_from_m(self, m):
        """`q_{b,i} = m_{b,i} / arm_i`, shape (n_points, n_stations); dtype follows `m`."""

        return np.asarray(m) / self.arm[None, :]

    def spanwise_moments_from_q(self, q):
        """`M_{b,k} = sum_i W_{k,i} q_{b,i}`, shape (n_points, n_stations)."""

        return np.asarray(q) @ self.W.T

    def flexibility_arms(self, chord):
        """`A_i(c) = sum_k G_k W_{k,i} / c_k^3`, shape (n_stations,); dtype follows `chord`."""

        chord = np.asarray(chord)
        return (self.G / chord ** 3) @ self.W

    def deflections_from_m(self, m, d):
        """`delta_b = sum_k G_k M_{b,k} / c_k^3` for every point; dtype follows the inputs."""

        chord = self.parameterisation.chord(np.asarray(d))
        moments = self.spanwise_moments_from_q(self.loads_from_m(m))
        return moments @ (self.G / chord ** 3)

    def deflections(self, phi, d):
        """`delta_b` through the kernel's own `m` value code (not a re-solve)."""

        return self.deflections_from_m(self.partials(phi, d, derivatives=False)["m"], d)

    def D(self, phi, d):
        """`KS_rho(delta / delta_ref)` -- the functional."""

        return ks(self.deflections(phi, d) / self.delta_ref, self.rho)

    def _deflection_softmax(self, deflections):
        return ks_weights(deflections / self.delta_ref, self.rho)

    def _pieces(self, phi, d, parts):
        """Chord, `q`, `M`, `A`, `delta`, softmax weights at `(phi, d)` from `parts`."""

        chord = self.parameterisation.chord(np.asarray(d))
        q = self.loads_from_m(parts["m"])
        moments = self.spanwise_moments_from_q(q)
        arms = self.flexibility_arms(chord)
        deflections = moments @ (self.G / chord ** 3)
        weights = self._deflection_softmax(deflections)
        return chord, q, moments, arms, deflections, weights

    def dD_dx(self, phi, d, parts=None):
        """`dD/dphi_{b,i} = w_b A_i m_phi_{b,i} / (arm_i delta_ref)`, shape (n_points, n_stations)."""

        parts = self.partials(phi, d) if parts is None else parts
        _chord, _q, _moments, arms, _deflections, weights = self._pieces(phi, d, parts)
        scale = (weights[:, None] / self.delta_ref) * (arms / self.arm)[None, :]
        return scale * parts["dm_dphi"]

    def dD_dd(self, phi, d, parts=None):
        """
        The explicit `dD/dd` (state held fixed), shape (n_design,): the load's
        own chord and twist dependence through `A_i`, plus the stiffness term
        `-3 G_j M_{b,j} / c_j^4` on the chord block.
        """

        parts = self.partials(phi, d) if parts is None else parts
        chord, _q, moments, arms, _deflections, weights = self._pieces(phi, d, parts)
        scale = weights[:, None] / self.delta_ref
        through_load = scale * (arms / self.arm)[None, :]
        dD_dc = ((through_load * parts["dm_dc"]).sum(axis=0)
                 + (scale * moments).sum(axis=0) * (-3.0 * self.G / chord ** 4))
        dD_dtheta = (through_load * parts["dm_dtheta"]).sum(axis=0)
        return self.N_c.T @ dD_dc + self.N_theta.T @ dD_dtheta

    def tangent(self, phi, d, v):
        """Forward-mode `dD/dd . v` through the implicit state (the Tier 2 twin)."""

        parts = self.partials(phi, d)
        dphi = -self.apply_dR_dd(parts, v) / parts["dR_dphi"]
        return (float(np.sum(self.dD_dx(phi, d, parts) * dphi))
                + float(self.dD_dd(phi, d, parts) @ np.asarray(v)))

    def gradient(self, d, state=None):
        """
        `dD/dd` by the discrete adjoint at design `d`: one forward solve (or
        the one passed in as `state` -- the moment row's, for the same `d`),
        the partials, 225 divisions, one assembly.
        """

        state = self.solve(d) if state is None else state
        parts = self.partials(state.phi, state.d)

        dR_dx = parts["dR_dphi"]
        dD_dx = self.dD_dx(state.phi, state.d, parts)
        psi = -dD_dx / dR_dx
        explicit = self.dD_dd(state.phi, state.d, parts)
        total = explicit + self.apply_dR_dd_T(parts, psi)

        _chord, _q, moments, _arms, deflections, weights = self._pieces(state.phi, state.d, parts)
        D = ks(deflections / self.delta_ref, self.rho)
        return DeflectionGradientResult(dD_dd=total, dD_dd_explicit=explicit, psi=psi,
                                        dR_dx=dR_dx, dD_dx=dD_dx, state=state,
                                        D=float(D), deflections=deflections,
                                        weights=weights, spanwise_moments_nm=moments)
