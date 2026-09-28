"""
`RootMomentSystem`: the second right-hand side -- the KS aggregate
of the flapwise root-bending moment over the load operating set `L`, and its
discrete adjoint.

It subclasses `BEMSystem` and inherits the forward solve, the residual, the
partials, `dR/dx`, `dR/dd`, the matrix-free operators and `state_sensitivity`
unchanged. The only difference is the functional: instead of `J = -AEP` it
carries

    M_b   = sum_i t_i m_{b,i}                    (trapezoid; m is the kernel's
                                                  per-blade moment integrand)
    Mhat_b = M_b / M_ref                          (M_ref passed in explicitly)
    KS     = KS_rho(Mhat)                         (the smooth max)

The state `phi` is the same 9 x 25 system; the adjoint is the same diagonal
solve with a different right-hand side:

    psi_{b,i} = -(dKS/dphi_{b,i}) / (dR_{b,i}/dphi_{b,i})
    dKS/dd    = dKS/dd|_explicit + sum_{b,i} psi_{b,i} dR_{b,i}/dd

The parent's AEP methods (`J`, `weights`, `J_capped`, `dJ_dx`, `dJ_dd`,
`powers_from_q`) mean nothing on the load set and are not called here; `solve`
still returns a `limited` mask but every entry of `L` is below the rating at
`x0` (asserted in `tests/test_loads.py`).

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

from dataclasses import dataclass

import numpy as np

from adjoint.system import BEMSystem, ForwardState
from objective.loads import ks, ks_weights, load_operating_points


@dataclass(frozen=True)
class MomentGradientResult:
    """`dKS/dd` and what it was assembled from (the `GradientResult` analogue)."""

    dKS_dd: np.ndarray
    dKS_dd_explicit: np.ndarray
    psi: np.ndarray
    dR_dx: np.ndarray
    dKS_dx: np.ndarray
    state: ForwardState
    KS: float
    moments_nm: np.ndarray
    weights: np.ndarray


class RootMomentSystem(BEMSystem):
    """
    The root-moment KS functional and its adjoint over the 9 load points.

    Parameters
    ----------
    parameterisation, bounds, resource, polar_cache : as `BEMSystem`.
    points : list of (float, float) or None
        The load operating set `L`. `None` is
        `objective.loads.load_operating_points()`; given explicitly the
        masses are unit weights (the KS weights are the softmax, not a mass).
    rho : float
        The KS stiffness.
    m_ref_nm : float
        The normalising root moment [N m], one number in one place: the
        caller (`ScaledProblem`) reads it from the committed `x0` and passes
        it in, so every artefact records the same reference.
    """

    def __init__(self, parameterisation, bounds, resource, polar_cache="sg6043",
                 points=None, rho=100.0, m_ref_nm=None):
        if points is None:
            points = load_operating_points()
        super().__init__(parameterisation, bounds, resource, polar_cache=polar_cache,
                         points=points, mass=None)
        self.rho = float(rho)
        if m_ref_nm is None:
            raise ValueError("m_ref_nm is required; pass the committed x0 root moment")
        self.m_ref = float(m_ref_nm)

    # -- the functional --------------------------------------------------------

    def moments_from_m(self, m):
        """`M_b = sum_i t_i m_{b,i}`, shape (n_points,); dtype follows `m`."""

        return np.asarray(m) @ np.asarray(self.trapz_weights)

    def moments(self, phi, d):
        """`M_b` through the kernel's own `m` value code (not a re-solve)."""

        return self.moments_from_m(self.partials(phi, d, derivatives=False)["m"])

    def KS(self, phi, d):
        """`KS_rho(M/M_ref)` -- the functional (the analogue of `J`)."""

        return ks(self.moments(phi, d) / self.m_ref, self.rho)

    def _softmax(self, moments):
        """The KS softmax weights `w_b` over the normalised moments."""

        return ks_weights(moments / self.m_ref, self.rho)

    def dKS_dx(self, phi, d, parts=None):
        """`dKS/dphi_{b,i} = w_b t_i m_phi_{b,i} / M_ref`, shape (n_points, n_stations)."""

        parts = self.partials(phi, d) if parts is None else parts
        weights = self._softmax(self.moments_from_m(parts["m"]))
        scale = weights[:, None] * np.asarray(self.trapz_weights)[None, :] / self.m_ref
        return scale * parts["dm_dphi"]

    def dKS_dd(self, phi, d, parts=None):
        """The explicit `dKS/dd` (state held fixed), shape (n_design,)."""

        parts = self.partials(phi, d) if parts is None else parts
        weights = self._softmax(self.moments_from_m(parts["m"]))
        scale = weights[:, None] * np.asarray(self.trapz_weights)[None, :] / self.m_ref
        dKS_dc = (scale * parts["dm_dc"]).sum(axis=0)
        dKS_dtheta = (scale * parts["dm_dtheta"]).sum(axis=0)
        return self.N_c.T @ dKS_dc + self.N_theta.T @ dKS_dtheta

    def tangent(self, phi, d, v):
        """
        Forward-mode `dKS/dd . v`: perturb `d` along `v`, propagate through the
        implicit state, sum into `KS`. The Tier 2 twin of the adjoint.
        """

        parts = self.partials(phi, d)
        dphi = -self.apply_dR_dd(parts, v) / parts["dR_dphi"]
        return (float(np.sum(self.dKS_dx(phi, d, parts) * dphi))
                + float(self.dKS_dd(phi, d, parts) @ np.asarray(v)))

    def gradient(self, d, state=None):
        """
        `dKS/dd` by the discrete adjoint at design `d`.

        One forward solve (or the one passed in as `state`), the partials at
        every station, 225 divisions for `psi`, one assembly.
        """

        state = self.solve(d) if state is None else state
        parts = self.partials(state.phi, state.d)

        dR_dx = parts["dR_dphi"]
        dKS_dx = self.dKS_dx(state.phi, state.d, parts)
        psi = -dKS_dx / dR_dx
        explicit = self.dKS_dd(state.phi, state.d, parts)
        total = explicit + self.apply_dR_dd_T(parts, psi)

        moments = self.moments_from_m(parts["m"])
        weights = self._softmax(moments)
        KS = ks(moments / self.m_ref, self.rho)
        return MomentGradientResult(dKS_dd=total, dKS_dd_explicit=explicit,
                                    psi=psi, dR_dx=dR_dx, dKS_dx=dKS_dx,
                                    state=state, KS=float(KS), moments_nm=moments,
                                    weights=weights)
