"""
Phase 3, B2 -- Tier 2: the transpose identity `<v, A u> = <A^T v, u>`.

Tier 1 verified every *entry* of the partial arrays. What it cannot see is
an indexing error in how those entries are wired together: a bin index
swapped with a station index, a design column applied to the wrong
station, a transpose taken on the wrong axis. Those show up exactly here,
because "apply `A`" and "apply `A^T`" are written as separate code paths in
`adjoint.system` -- not one explicit matrix and its `.T`, which would test
NumPy -- and the identity holds only if both paths agree about which entry
multiplies which.

Three operators, at `x0` and at the perturbed design of Tier 1:

    A = dR/dx        (425 x 425, diagonal)          trivial; still done
    A = dR/dd        (425 x 10)                     apply_dR_dd / apply_dR_dd_T
    A = d -> dJ/dd   (1 x 10, the full chain)       tangent (forward mode)
                                                    vs the adjoint gradient

All to 1e-14 relative, where "relative" is to the Cauchy-Schwarz bound
`||v|| ||A u||` of the inner product -- the scale of the terms being summed.
An inner product of two random vectors cancels (its value can be any
fraction of that bound, including ~0), so relative-to-the-value would test
the luck of the draw rather than the operator; the bound is the scale on
which round-off of the operator application actually lives. Failure here is
a bug: fix, never loosen.

Bounds are provisional (`chord_max_m = 0.45 m` is a placeholder); they enter
only through `BEMSystem`'s constructor.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os

import numpy as np
import pytest

from adjoint import BEMSystem
from design import BladeParameterisation, DesignBounds
from objective import WeibullResource
from test_parameterisation import PROVISIONAL_BOUNDS

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.abspath(os.path.join(_HERE, "..", "verification", "baseline", "x0.json"))

RTOL = 1e-14
N_DRAWS = 8


@pytest.fixture(scope="module")
def system():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist, **PROVISIONAL_BOUNDS)
    return BEMSystem(parameterisation, bounds, WeibullResource.from_config())


@pytest.fixture(scope="module", params=["x0", "perturbed"])
def state(request, system):
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    d = np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])
    if request.param == "perturbed":
        n_chord = system.parameterisation.n_chord
        d[:n_chord] *= 1.2
        d[n_chord:] += math.radians(3.0)
    return system.solve(d)


@pytest.fixture(scope="module")
def parts(system, state):
    return system.partials(state.phi, state.d)


def _rng():
    return np.random.default_rng(2026)


def _assert_identity(left, right, scale, label):
    """`|left - right| <= RTOL * scale`, `scale` the Cauchy-Schwarz bound."""

    assert abs(left - right) <= RTOL * scale, (
        f"{label}: {left!r} vs {right!r}, diff {abs(left - right):.3e} on scale {scale:.3e}")


def _bound(v, Au):
    return float(np.linalg.norm(np.ravel(v)) * np.linalg.norm(np.ravel(Au)))


# ---------------------------------------------------------------------------
# A = dR/dx, diagonal
# ---------------------------------------------------------------------------

def _apply_dR_dx(parts, u):
    return parts["dR_dphi"] * u


def _apply_dR_dx_T(parts, v):
    # Written out rather than aliased to `_apply_dR_dx`: the point is that the
    # transpose of a diagonal operator is the same elementwise product, and
    # the identity below is what says so.
    return v * parts["dR_dphi"]


def test_dR_dx_transpose_identity(system, parts):
    rng = _rng()
    shape = (system.n_points, system.n_stations)
    for _ in range(N_DRAWS):
        u = rng.normal(size=shape)
        v = rng.normal(size=shape)
        Au = _apply_dR_dx(parts, u)
        _assert_identity(float(np.sum(v * Au)), float(np.sum(_apply_dR_dx_T(parts, v) * u)),
                         _bound(v, Au), "dR/dx")


# ---------------------------------------------------------------------------
# A = dR/dd, (18, 25) x 10
# ---------------------------------------------------------------------------

def test_dR_dd_transpose_identity(system, parts):
    rng = _rng()
    for _ in range(N_DRAWS):
        u = rng.normal(size=system.n_design)
        v = rng.normal(size=(system.n_points, system.n_stations))
        Au = system.apply_dR_dd(parts, u)
        _assert_identity(float(np.sum(v * Au)), float(system.apply_dR_dd_T(parts, v) @ u),
                         _bound(v, Au), "dR/dd")


def test_apply_dR_dd_agrees_with_the_explicit_tensor(system, state, parts):
    """The matrix-free apply against the (18, 25, 10) tensor `dR_dd` builds."""

    tensor = system._assemble_dR_dd(parts)
    rng = _rng()
    u = rng.normal(size=system.n_design)
    explicit = np.einsum("bij,j->bi", tensor, u)
    # Each entry is a sum of two terms (chord path, twist path) that can
    # cancel; the scale is the sum of their magnitudes, not the result.
    scale = (np.abs(parts["dR_dc"] * (system.N_c @ u)[None, :])
             + np.abs(parts["dR_dtheta"] * (system.N_theta @ u)[None, :]))
    assert np.all(np.abs(system.apply_dR_dd(parts, u) - explicit) <= RTOL * scale)
    v = rng.normal(size=(system.n_points, system.n_stations))
    explicit_T = np.einsum("bij,bi->j", tensor, v)
    scale_T = np.einsum("bij,bi->j", np.abs(tensor), np.abs(v))
    assert np.all(np.abs(system.apply_dR_dd_T(parts, v) - explicit_T) <= RTOL * scale_T)


def test_apply_dR_dd_T_is_not_a_reshape_of_apply_dR_dd(system, parts):
    """
    A wrong-axis sum would still pass a badly chosen identity if `u` and `v`
    were symmetric enough. Check the shapes are what the operator says.
    """

    u = np.zeros(system.n_design)
    u[0] = 1.0
    column = system.apply_dR_dd(parts, u)
    assert column.shape == (system.n_points, system.n_stations)
    # chord control point 0 reaches only the stations its basis function covers
    support = system.N_c[:, 0] != 0.0
    assert np.all(column[:, ~support] == 0.0)
    assert np.any(column[:, support] != 0.0)

    w = np.zeros((system.n_points, system.n_stations))
    w[3, 7] = 1.0
    row = system.apply_dR_dd_T(parts, w)
    assert row.shape == (system.n_design,)
    assert np.array_equal(row, system._assemble_dR_dd(parts)[3, 7])


# ---------------------------------------------------------------------------
# the full operator d -> dJ/dd: forward-mode tangent vs the adjoint
# ---------------------------------------------------------------------------

def test_full_chain_tangent_equals_adjoint_gradient(system, state):
    """
    `<w, A u> = <A^T w, u>` with `A: R^10 -> R` the whole implicit chain.
    Left side: `tangent(u)` (perturb `d`, solve the linearised state, sum
    into `J`). Right side: `gradient(d) . u` (adjoint). Different code paths
    through the same partials.
    """

    result = system.gradient(state.d, state=state)
    rng = _rng()
    for _ in range(N_DRAWS):
        u = rng.normal(size=system.n_design)
        w = float(rng.normal())
        Au = system.tangent(state.phi, state.d, u, state.limited)
        left = w * Au
        right = float((w * result.dJ_dd) @ u)
        _assert_identity(left, right, abs(w) * float(np.linalg.norm(result.dJ_dd)
                                                   * np.linalg.norm(u)), "full chain")


def test_adjoint_gradient_is_explicit_plus_implicit(system, state):
    """`dJ/dd = dJ/dd|explicit + dR/dd^T psi`, and psi solves the diagonal system."""

    result = system.gradient(state.d, state=state)
    parts = system.partials(state.phi, state.d)
    implicit = system.apply_dR_dd_T(parts, result.psi)
    assert np.array_equal(result.dJ_dd, result.dJ_dd_explicit + implicit)
    assert np.all(np.abs(result.dR_dx * result.psi + result.dJ_dx) <= RTOL * np.abs(result.dJ_dx))
    assert result.J == pytest.approx(float(np.sum(system.weights(state.limited) * state.power_w))
                                     + system.J_capped(state.limited), rel=1e-12)
