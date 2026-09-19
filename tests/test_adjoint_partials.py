"""
Phase 3, B1 -- Tier 1: every hand-derived partial against a complex step of
the code's own residual and objective, at machine precision.

What is compared to what
-------------------------
`adjoint.kernels.station_partials` transcribes `docs/adjoint_derivation.md`
§2-§6. The reference on the other side is never the kernel: it is a complex
step (`h = 1e-30`, no subtraction, exact to round-off) through
`bem.station.residual` -- the function the forward solver actually
root-finds -- and through `BEMSystem.J`, which is pinned to
`objective.annual_energy_mwh` to 1e-12 relative below.

    dR/dphi     complex phi into station.residual, every (b, i)
    dR/dd_j     complex d -> complex chord, twist, Reynolds -> station.residual
    dJ/dphi     complex phi_{b,i} into J, every (b, i)   (425 evaluations)
    dJ/dd_j     complex d into J

at `x0` and at a perturbed design (chord x 1.2, twist + 3 deg), so the
partials are not verified at one lucky point.

The assertion is mixed, `|err| <= 1e-13 max(1, |partial|)` (`1e-12` for `J`,
a 425-term sum), so a partial that is legitimately ~0 is judged on absolute
error and one that is O(10) on relative. A plateau near 1e-6 would mean a
real dtype somewhere on the path, not mathematics. **Failure here is a bug:
fix it, never loosen the tolerance.**

The set must contain Buhl-branch stations (`a > BUHL_AC`) -- the implicit
differentiation of §4.4 is otherwise untested -- and none on the
`gamma2 < 0` branch, which the adjoint does not model.

Bounds enter only through `BEMSystem`'s constructor (the scaling chain),
from `DesignBounds.from_config()`.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os

import numpy as np
import pytest

from adjoint import BEMSystem
from bem.corrections import BUHL_AC, buhl_gammas
from design import BladeParameterisation, DesignBounds
from objective import WeibullResource, annual_energy_mwh

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.abspath(os.path.join(_HERE, "..", "verification", "baseline", "x0.json"))

H = 1e-30
TOL_R = 1e-13
TOL_J = 1e-12


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
def system():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist)
    return BEMSystem(parameterisation, bounds, WeibullResource.from_config())


@pytest.fixture(scope="module")
def x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])


@pytest.fixture(scope="module")
def perturbed(system, x0):
    """chord x 1.2, twist + 3 deg: a different loading, a different alpha band."""

    d = x0.copy()
    n_chord = system.parameterisation.n_chord
    d[:n_chord] *= 1.2
    d[n_chord:] += math.radians(3.0)
    return d


@pytest.fixture(scope="module", params=["x0", "perturbed"])
def design(request, x0, perturbed):
    return x0 if request.param == "x0" else perturbed


@pytest.fixture(scope="module")
def state(system, design):
    return system.solve(design)


@pytest.fixture(scope="module")
def parts(system, state):
    return system.partials(state.phi, state.d)


# ---------------------------------------------------------------------------
# the values: kernel R == station.residual, J == annual_energy_mwh
# ---------------------------------------------------------------------------

def test_kernel_residual_value_matches_station_residual(system, state, parts):
    """Same calls in the same order: the residual value is the code's own."""

    reference = system.residual(state.phi, state.d)
    assert reference.shape == (system.n_points, system.n_stations)
    assert np.max(np.abs(parts["residual"] - reference)) <= 1e-15


def test_state_is_a_root_of_the_residual(system, state):
    assert np.max(np.abs(system.residual(state.phi, state.d))) < 1e-12


def test_J_at_the_solved_state_matches_annual_energy_mwh(system, state):
    J = system.J(state.phi, state.d, state.limited)
    reference = -annual_energy_mwh(state.d, resource=system.resource,
                                   parameterisation=system.parameterisation)
    assert abs(J - reference) <= 1e-12 * abs(reference)


def test_the_limited_mask_is_the_forward_solves(system, state, parts):
    """The mask recomputed from `q` at the state equals `power_per_bin`'s."""

    power = system.powers_from_q(parts["q"])
    assert np.array_equal(system.limited_mask(power), state.limited)
    assert state.limited.any() and not state.limited.all()


# ---------------------------------------------------------------------------
# the set: Buhl-branch stations present, gamma2 < 0 absent, dR/dphi != 0
# ---------------------------------------------------------------------------

def test_set_contains_buhl_branch_stations_and_none_past_the_pole(system, state, parts):
    above = state.a > BUHL_AC
    assert above.sum() > 0, "no station on the Buhl branch: the implicit partials are untested"
    assert np.array_equal(parts["buhl"], above)

    chord, _twist = system.chord_twist(state.d)
    for b in range(system.n_points):
        for i in range(system.n_stations):
            if not parts["buhl"][b, i]:
                continue
            sigma = system.n_blades * chord[i] / (2.0 * math.pi * system.radii[i])
            Y = sigma * parts["cn"][b, i] / math.sin(state.phi[b, i]) ** 2
            F = parts["F"][b, i]
            _g1, gamma2, _g3 = buhl_gammas(Y / (4.0 * F), F)
            assert gamma2 > 0.0, (b, i, gamma2)


def test_dR_dphi_is_bounded_away_from_zero(parts):
    """The diagonal of dR/dx must be invertible; report its smallest entry."""

    smallest = np.min(np.abs(parts["dR_dphi"]))
    assert smallest > 1e-2, smallest


# ---------------------------------------------------------------------------
# dR/dphi and dR/dd against a complex step of station.residual
# ---------------------------------------------------------------------------

def test_dR_dphi_matches_complex_step(system, state, parts):
    """
    All 450 diagonal entries. Each station's residual sees only its own
    `phi`, so one complex step on every entry at once gives every diagonal
    entry; `test_dR_dx_is_diagonal` checks that structural claim separately.
    """

    estimate = system.residual(state.phi + 1j * H, state.d).imag / H
    _assert_mixed(estimate, parts["dR_dphi"], TOL_R, "dR/dphi")


def test_dR_dx_is_diagonal(system, state, parts):
    """Single-entry steps: the entry matches and every other entry is exactly 0."""

    rng = np.random.default_rng(3)
    for _ in range(12):
        b = int(rng.integers(system.n_points))
        i = int(rng.integers(system.n_stations))
        phi = state.phi.astype(complex)
        phi[b, i] += 1j * H
        estimate = system.residual(phi, state.d).imag / H
        assert _mixed(estimate[b, i] - parts["dR_dphi"][b, i], parts["dR_dphi"][b, i], TOL_R)
        estimate[b, i] = 0.0
        assert np.all(estimate == 0.0)


def test_dR_dd_matches_complex_step(system, state, parts):
    """
    Every design variable: a complex `d` becomes a complex chord, twist and
    Reynolds number, so the `Re(c)` path is on the complex side too.
    """

    dR_dd = system.dR_dd(state.phi, state.d)
    assert dR_dd.shape == (system.n_points, system.n_stations, system.n_design)
    for j in range(system.n_design):
        d = state.d.astype(complex)
        d[j] += 1j * H
        estimate = system.residual(state.phi, d).imag / H
        _assert_mixed(estimate, dR_dd[:, :, j], TOL_R, f"dR/dd_{j}")


def test_dR_dd_twist_columns_have_no_reynolds_path(system, state, parts):
    """Structural: twist control points do not move the Reynolds number."""

    n_chord = system.parameterisation.n_chord
    dR_dd = system.dR_dd(state.phi, state.d)
    expected = parts["dR_dtheta"][:, :, None] * system.N_theta[None, :, n_chord:]
    assert np.array_equal(dR_dd[:, :, n_chord:], expected)


# ---------------------------------------------------------------------------
# dJ/dphi and dJ/dd against a complex step of J
# ---------------------------------------------------------------------------

def test_dJ_dphi_matches_complex_step_at_every_station(system, state, parts):
    """425 complex steps of `J`, one per (b, i): the limited bins' entries
    must be exactly zero (the rating is a configured constant, so a capped
    bin's state does not enter J), every unlimited bin's non-zero."""

    dJ_dx = system.dJ_dx(state.phi, state.d, state.limited, parts)
    assert dJ_dx.shape == (system.n_points, system.n_stations)
    estimate = np.empty_like(dJ_dx)
    for b in range(system.n_points):
        for i in range(system.n_stations):
            phi = state.phi.astype(complex)
            phi[b, i] += 1j * H
            estimate[b, i] = system.J(phi, state.d, state.limited).imag / H
    _assert_mixed(estimate, dJ_dx, TOL_J, "dJ/dphi")

    limited_rows = np.flatnonzero(state.limited)
    assert np.all(dJ_dx[limited_rows] == 0.0)
    assert np.all(estimate[limited_rows] == 0.0)
    assert np.all(np.any(dJ_dx[~state.limited] != 0.0, axis=1))


def test_dJ_dd_matches_complex_step(system, state, parts):
    dJ_dd = system.dJ_dd(state.phi, state.d, state.limited, parts)
    assert dJ_dd.shape == (system.n_design,)
    for j in range(system.n_design):
        d = state.d.astype(complex)
        d[j] += 1j * H
        estimate = system.J(state.phi, d, state.limited).imag / H
        assert _mixed(estimate - dJ_dd[j], dJ_dd[j], TOL_J), (j, estimate, dJ_dd[j])


def test_no_plateau_at_the_real_dtype_signature(system, state, parts):
    """
    The 1e-6 plateau is the fingerprint of a real-dtype array on a complex
    path. The mixed tolerance above already excludes it; this names it.
    """

    estimate = system.residual(state.phi + 1j * H, state.d).imag / H
    worst = np.max(np.abs(estimate - parts["dR_dphi"]) / np.maximum(1.0, np.abs(parts["dR_dphi"])))
    assert worst < 1e-12, f"worst mixed error {worst:.3e}: a real dtype on the path?"
