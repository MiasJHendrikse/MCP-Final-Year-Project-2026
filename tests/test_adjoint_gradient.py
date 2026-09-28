"""
The assembled adjoint gradient on the scaled problem.

Tier 1 and Tier 2 live in their own files. This one pins the two things
this step adds on top of them: `ScaledProblem.jac_adjoint` is the adjoint
gradient in `fun` units, and it agrees with the committed FD reference at
`x0` to within the step-size study's `eps_j` (a flatness estimate of the
FD's disagreement scale, not a measured noise floor) -- the Tier 3
acceptance, asserted here against `verification/fd_step_size/sweep.json` so
that a change to either side is caught by the suite, not only by re-running
the study. The
whole-chain Taylor remainder is asserted at `x0` as well.

Bounds are the configured ones (`DesignBounds.from_config()`).

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import os

import numpy as np
import pytest

from design import BladeParameterisation, DesignBounds
from gradients import ScaledProblem
from objective import WeibullResource

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.abspath(os.path.join(_HERE, "..", "verification", "baseline", "x0.json"))
SWEEP_PATH = os.path.abspath(os.path.join(_HERE, "..", "verification", "fd_step_size", "sweep.json"))


@pytest.fixture(scope="module")
def problem():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                          n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config())


@pytest.fixture(scope="module")
def u0(problem):
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    x0 = np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])
    u0 = problem.scaled(x0)
    problem.set_reference(u0)
    return u0


@pytest.fixture(scope="module")
def sweep():
    with open(SWEEP_PATH, encoding="utf-8") as handle:
        return json.load(handle)


@pytest.fixture(scope="module")
def adjoint_scaled(problem, u0):
    return problem.jac_adjoint(u0)


def test_jac_adjoint_requires_a_reference(problem):
    fresh = ScaledProblem(problem.parameterisation, problem.bounds, problem.resource)
    with pytest.raises(RuntimeError):
        fresh.jac_adjoint(np.full(fresh.n, 0.5))


def test_jac_adjoint_is_the_scaled_system_gradient(problem, u0, adjoint_scaled):
    assert adjoint_scaled.shape == (problem.n,)
    result = problem.adjoint_system().gradient(problem.physical(u0))
    expected = result.dJ_dd * problem.bounds.span() / abs(problem.J0)
    assert np.array_equal(adjoint_scaled, expected)
    assert problem.n_adjoint_evals >= 1


def test_adjoint_agrees_with_the_committed_fd_reference_at_x0(problem, sweep, adjoint_scaled):
    """
    Tier 3 at `x0`, against the committed step-size study: every variable
    within `3 eps_j` of the FD at `h*_j`, and within the 1e-3 tripwire.
    """

    assert sweep["J0_mwh_per_yr"] == pytest.approx(problem.J0, rel=1e-12)
    adjoint = problem.unscale(adjoint_scaled)
    fd = np.array(sweep["gradient_at_h_star_mwh_per_u"])
    eps = np.array(sweep["epsilon_mwh_per_u"])
    abs_err = np.abs(adjoint - fd)
    assert np.all(abs_err <= 3.0 * eps), (abs_err / eps).tolist()
    assert np.all(abs_err / np.abs(fd) <= 1e-3)


def test_taylor_remainder_falls_at_second_order(problem, u0, adjoint_scaled):
    """`|fun(u + eps v) - fun(u) - eps g.v|` drops ~100x per decade; assert >= 30."""

    rng = np.random.default_rng(11)
    v = rng.normal(size=problem.n)
    v /= np.linalg.norm(v)
    f0 = problem.fun(u0)
    slope = float(adjoint_scaled @ v)
    remainders = [abs(problem.fun(u0 + eps * v) - f0 - eps * slope) for eps in (1e-2, 1e-3, 1e-4)]
    ratios = [remainders[0] / remainders[1], remainders[1] / remainders[2]]
    assert min(ratios) >= 30.0, (remainders, ratios)
