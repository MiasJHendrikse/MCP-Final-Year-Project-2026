"""
The mass problem, Step 2 -- `DeflectionSystem`, the tip-deflection KS
functional and its adjoint over the load set `L`.

What is compared to what
-------------------------
  * the forward path `objective.loads.tip_deflection` -> `bem.rotor.solve_rotor`
    (the solver's own `w`, `Cn`; sub-grid trapezoid in a plain loop) against
    the system's `deflections` (the kernel's `m`, the constant map `W`);
  * a complex step (h = 1e-30) through the functional's *value* code
    (`D(phi, d)`) against `dD_dx`, `dD_dd` (Tier 1, at the functional level:
    TOL_J, as `dKS_dx` / `dKS_dd`);
  * the forward-mode tangent against the adjoint direction (Tier 2);
  * the constraint Jacobian against central FD at `u0` (Tier 3; the
    `ScaledProblem` wrapper arrives in Step 3 and the test is appended then).

Tolerances are the existing mixed ones; never loosen them.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os

import numpy as np
import pytest

from adjoint import BEMSystem, DeflectionSystem, RootMomentSystem
from config import load_site
from design import BladeParameterisation, DesignBounds
from objective import WeibullResource
from objective.loads import (
    load_operating_points,
    root_bending_moment,
    spanwise_moments,
    tip_deflection,
    tip_deflection_at,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(_HERE, ".."))
X0_PATH = os.path.join(ROOT, "verification", "baseline", "x0.json")

H = 1e-30
TOL_J = 1e-12
M_REFERENCE = 177.3755406092969
#: `tip_deflection` at `x0`, 11 m/s on the 300 rpm ceiling, per unit E k_I t_shell.
DELTA_REFERENCE = 30776.559398053847
RATED = (11.0, 5.711986642890533)


def _mixed(err, partial, tol):
    return np.abs(err) <= tol * np.maximum(1.0, np.abs(partial))


def _assert_mixed(estimate, partial, tol, label):
    err = np.asarray(estimate) - np.asarray(partial)
    ok = _mixed(err, partial, tol)
    if not np.all(ok):
        worst = np.unravel_index(np.argmax(np.abs(err) / np.maximum(1.0, np.abs(partial))),
                                 np.shape(err))
        pytest.fail(
            f"{label}: {int((~ok).sum())} of {ok.size} partials outside "
            f"{tol:g} max(1, |partial|); worst at {worst}: complex step "
            f"{np.asarray(estimate)[worst]!r}, derived {np.asarray(partial)[worst]!r}, "
            f"err {err[worst]:.3e}")


@pytest.fixture(scope="module")
def x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])


@pytest.fixture(scope="module")
def x_pert(x0):
    """The `test_loads.py` perturbation: +-2 % chord, +-0.5 deg twist, alternating."""

    d = x0.copy()
    for j in range(5):
        d[j] *= 1.02 if j % 2 == 0 else 0.98
    for j in range(5, len(d)):
        d[j] += math.radians(0.5 if j % 2 == 0 else -0.5)
    return d


@pytest.fixture(scope="module", params=["x0", "x_pert"])
def design(request, x0, x_pert):
    return x0 if request.param == "x0" else x_pert


@pytest.fixture(scope="module")
def system():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    return DeflectionSystem(parameterisation, bounds, WeibullResource.from_config(),
                            m_ref_nm=M_REFERENCE, delta_ref=DELTA_REFERENCE)


@pytest.fixture(scope="module")
def state(system, design):
    return system.solve(design)


@pytest.fixture(scope="module")
def parts(system, state):
    return system.partials(state.phi, state.d)


# ---------------------------------------------------------------------------
# the geometry of the functional
# ---------------------------------------------------------------------------

def test_delta_ref_is_required():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=5, n_twist=5)
    with pytest.raises(ValueError, match="delta_ref"):
        DeflectionSystem(parameterisation, bounds, WeibullResource.from_config(),
                         m_ref_nm=M_REFERENCE)


def test_w_is_strictly_outboard_and_constant(system):
    """`W[k, i] = t_i (r_i - r_k)` for `i > k`, zero on and below the diagonal."""

    W = system.W
    assert W.shape == (system.n_stations, system.n_stations)
    assert np.all(np.tril(W) == 0.0)
    radii = np.asarray(system.radii)
    t = np.asarray(system.trapz_weights)
    k, i = 3, 10
    assert W[k, i] == pytest.approx(t[i] * (radii[i] - radii[k]), rel=1e-15)
    assert np.all(system.arm > 0.0)


def test_w_with_the_hub_arm_is_the_root_moment(system, state, parts):
    """The same rule with `r_k -> r_hub` is `RootMomentSystem.moments_from_m`:
    `q` times the hub arm, trapezoid-summed, equals the parent's moments."""

    q = system.loads_from_m(parts["m"])
    via_arm = (q * system.arm[None, :] * np.asarray(system.trapz_weights)[None, :]).sum(axis=1)
    assert np.allclose(via_arm, system.moments_from_m(parts["m"]), rtol=1e-14)


def test_two_forms_of_the_deflection_agree(system, state, parts):
    """`sum_k G_k M_k / c_k^3` equals `sum_i A_i(c) q_i` (the transposed form)."""

    chord = system.parameterisation.chord(state.d)
    q = system.loads_from_m(parts["m"])
    transposed = q @ system.flexibility_arms(chord)
    assert np.allclose(transposed, system.deflections_from_m(parts["m"], state.d), rtol=1e-13)


# ---------------------------------------------------------------------------
# forward path
# ---------------------------------------------------------------------------

def test_spanwise_moment_at_the_first_station_is_below_the_root_moment(x0):
    """`M(r_0)` has a shorter arm than `M(r_hub)`; both from the same records."""

    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(x0)
    site = load_site()
    delta, result = tip_deflection_at(geometry, *RATED, site.air_density, site.kinematic_viscosity)
    moments = spanwise_moments(result["stations"], geometry.chord, site.air_density)
    root = root_bending_moment(result["stations"], geometry.chord, site.air_density,
                               geometry.n_blades, geometry.r_hub)
    assert root == pytest.approx(M_REFERENCE, rel=1e-12)
    assert 0.0 < moments[0] < root
    assert moments[-1] == 0.0
    assert all(moments[k] > moments[k + 1] for k in range(len(moments) - 1))
    assert delta == pytest.approx(DELTA_REFERENCE, rel=1e-12)


def test_system_deflections_equal_forward_path(system, state):
    """The kernel's `m` through `W` against the solver's records through the
    loop in `objective.loads` -- two code paths, one physics, at every point of `L`."""

    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(state.d)
    site = load_site()
    forward = np.array([tip_deflection_at(geometry, v, lam, site.air_density,
                                          site.kinematic_viscosity)[0]
                        for v, lam in load_operating_points()])
    assert np.allclose(system.deflections(state.phi, state.d), forward, rtol=1e-12)


def test_tip_deflection_is_the_unit_load_integral(x0):
    """`int M (R - r) / c^3 dr` on the trapezoid, spelled out."""

    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(x0)
    site = load_site()
    _delta, result = tip_deflection_at(geometry, *RATED, site.air_density, site.kinematic_viscosity)
    stations = result["stations"]
    radii = np.array([s["r"] for s in stations])
    moments = np.array(spanwise_moments(stations, geometry.chord, site.air_density))
    integrand = moments * (geometry.R - radii) / np.asarray(geometry.chord) ** 3
    expected = float(np.trapezoid(integrand, radii))
    assert tip_deflection(stations, geometry.chord, site.air_density, geometry.R) \
        == pytest.approx(expected, rel=1e-13)


def test_rated_point_dominates_the_aggregate_at_x0(system, x0):
    result = system.gradient(x0)
    assert result.deflections[-1] == pytest.approx(DELTA_REFERENCE, rel=1e-12)
    assert result.weights[-1] > 0.95
    assert result.D == pytest.approx(1.0004420795243203, rel=1e-10)
    assert np.all(np.diff(result.deflections) > 0.0)


# ---------------------------------------------------------------------------
# Tier 1: the functional's partials against a complex step of its value code
# ---------------------------------------------------------------------------

def test_dD_dx_matches_complex_step(system, state, parts):
    """All 225 entries: a complex step of `phi` through `D(phi, d)`."""

    dD_dx = system.dD_dx(state.phi, state.d, parts)
    estimate = np.empty_like(dD_dx)
    for b in range(system.n_points):
        for i in range(system.n_stations):
            phi = state.phi.astype(complex)
            phi[b, i] += 1j * H
            estimate[b, i] = system.D(phi, state.d).imag / H
    _assert_mixed(estimate, dD_dx, TOL_J, "dD/dphi")


def test_dD_dd_matches_complex_step(system, state, parts):
    """Every design variable: a complex `d` through `D(phi, d)` -- this is
    where the stiffness term `-3 G M / c^4` is checked."""

    dD_dd = system.dD_dd(state.phi, state.d, parts)
    for j in range(system.n_design):
        d = state.d.astype(complex)
        d[j] += 1j * H
        estimate = system.D(state.phi, d).imag / H
        _assert_mixed(estimate, dD_dd[j], TOL_J, f"dD/dd_{j}")


# ---------------------------------------------------------------------------
# Tier 2: the forward-mode tangent and the inherited transpose
# ---------------------------------------------------------------------------

def test_tangent_matches_adjoint_direction(system, state):
    rng = np.random.default_rng(42)
    gradient = system.gradient(state.d, state=state).dD_dd
    for _ in range(3):
        v = rng.normal(size=system.n_design)
        assert system.tangent(state.phi, state.d, v) == pytest.approx(float(gradient @ v), rel=1e-12)


def test_transpose_is_the_parent_operator():
    assert DeflectionSystem.apply_dR_dd_T is BEMSystem.apply_dR_dd_T
    assert DeflectionSystem.solve is BEMSystem.solve
    assert DeflectionSystem.partials is BEMSystem.partials


def test_gradient_reuses_a_passed_state(system, state):
    """`gradient(d, state=...)` does not re-solve: the result carries that state."""

    result = system.gradient(state.d, state=state)
    assert result.state is state


def test_moment_functional_still_available_on_the_same_instance(system, state):
    """The parent's KS is untouched by the subclass."""

    assert system.KS(state.phi, state.d) == pytest.approx(
        RootMomentSystem.KS(system, state.phi, state.d), rel=1e-15)
