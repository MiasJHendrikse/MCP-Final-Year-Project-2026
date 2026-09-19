"""
Phase 2, Stage A acceptance: the FD gradient and the scaled problem.

Everything the FD-driven SLSQP run rests on that can be checked without
running it: the central difference is exact where it should be exact, the
step-size sweep records an out-of-cache step rather than patching it, the
objective is normalised to exactly -1 at the start, and the polar-cache
envelope constraint says yes at `x0` and no where the box bounds would take
the optimiser outside the cache.

Bounds are the configured ones (`DesignBounds.from_config()`; `chord_max_m =
0.30 m` since 2026-09-19).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import os

import numpy as np
import pytest

from design import BladeParameterisation, DesignBounds
from gradients import ScaledProblem, central_difference, step_size_sweep
from gradients.finite_difference import OUT_OF_CACHE
from objective import WeibullResource
from polars.interpolant import PolarDomainError

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.abspath(os.path.join(_HERE, "..", "verification", "baseline", "x0.json"))


@pytest.fixture(scope="module")
def parameterisation():
    return BladeParameterisation()


@pytest.fixture(scope="module")
def bounds(parameterisation):
    return DesignBounds.from_config(n_chord=parameterisation.n_chord,
                        n_twist=parameterisation.n_twist)


@pytest.fixture(scope="module")
def x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"]
                    + artefact["twist_control_points_rad"])


@pytest.fixture(scope="module")
def problem(parameterisation, bounds):
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config())


@pytest.fixture(scope="module")
def u0(problem, x0):
    return problem.scaled(x0)


# ---------------------------------------------------------------------------
# central_difference: pure, exact on a quadratic
# ---------------------------------------------------------------------------

def _quadratic():
    rng = np.random.default_rng(7)
    matrix = rng.normal(size=(6, 6))
    matrix = matrix @ matrix.T
    linear = rng.normal(size=6)

    def fun(u):
        return 0.5 * u @ matrix @ u + linear @ u

    def grad(u):
        return matrix @ u + linear

    return fun, grad


def test_central_difference_is_exact_on_a_quadratic():
    """Central differences have no O(h^2) term on a quadratic: only round-off."""

    fun, grad = _quadratic()
    u = np.linspace(-1.0, 1.0, 6)
    estimate, n_evals = central_difference(fun, u, 1e-3)
    assert n_evals == 12
    assert np.max(np.abs(estimate - grad(u))) < 1e-12


def test_central_difference_accepts_a_per_variable_step():
    fun, grad = _quadratic()
    u = np.linspace(-1.0, 1.0, 6)
    steps = np.logspace(-2, -4, 6)
    estimate, _ = central_difference(fun, u, steps)
    assert np.max(np.abs(estimate - grad(u))) < 1e-10


def test_central_difference_rejects_a_non_positive_step():
    fun, _ = _quadratic()
    with pytest.raises(ValueError):
        central_difference(fun, np.zeros(6), 0.0)


def test_central_difference_is_central_not_forward():
    """A forward difference on a cubic has O(h) error; central has O(h^2)."""

    fun = lambda u: float(u[0] ** 3)  # noqa: E731
    u = np.array([1.0])
    h = 1e-3
    estimate, _ = central_difference(fun, u, h)
    # exact derivative 3; central error is exactly h^2 for u^3
    assert abs(estimate[0] - 3.0) == pytest.approx(h ** 2, rel=1e-6)


# ---------------------------------------------------------------------------
# step_size_sweep: every step kept, out-of-cache recorded not substituted
# ---------------------------------------------------------------------------

def test_step_size_sweep_keeps_every_step_and_records_out_of_cache():
    fun, grad = _quadratic()
    u = np.zeros(6)
    steps = np.logspace(-2, -4, 5)
    bad_h = float(steps[1])

    def fenced(v):
        # Only variable 2 stepping by exactly bad_h leaves the "cache".
        if abs(abs(v[2]) - bad_h) < 1e-15:
            raise PolarDomainError("simulated")
        return fun(v)

    results = step_size_sweep(fenced, u, steps)

    assert [h for h in results if h != OUT_OF_CACHE] == [float(h) for h in steps]
    assert results[OUT_OF_CACHE] == [(bad_h, 2)]
    assert np.isnan(results[bad_h][2])
    assert np.all(np.isfinite(np.delete(results[bad_h], 2)))
    for h in steps:
        h = float(h)
        finite = np.isfinite(results[h])
        assert np.max(np.abs(results[h][finite] - grad(u)[finite])) < 1e-10


# ---------------------------------------------------------------------------
# ScaledProblem
# ---------------------------------------------------------------------------

def test_fun_at_u0_is_exactly_minus_one(problem, u0):
    assert problem.fun(u0) == -1.0
    assert problem.J0 < 0.0
    assert problem.unscale(-1.0) == pytest.approx(problem.J0)


def test_envelope_constraint_is_satisfied_at_u0(problem, u0):
    constraint = problem.envelope_constraint()
    assert constraint["type"] == "ineq"
    g = constraint["fun"](u0)
    assert g.shape == (2 * problem.parameterisation.n_stations,)
    assert np.all(g > 0.0)


def test_envelope_constraint_is_violated_with_every_chord_at_its_bound(problem):
    """u = 1: all chord control points at chord_max_m -> tip Re > 1 M at 19.5 m/s."""

    g = problem.envelope_constraint()["fun"](np.ones(problem.n))
    ceiling_rows = g[problem.parameterisation.n_stations:]
    assert np.any(ceiling_rows < 0.0)


def test_envelope_jacobian_is_the_constant_derivative_of_its_fun(problem, u0):
    constraint = problem.envelope_constraint()
    jac = constraint["jac"](u0)
    assert jac.shape == (2 * problem.parameterisation.n_stations, problem.n)
    # twist columns are structurally zero: Reynolds depends on chord only
    assert np.all(jac[:, problem.parameterisation.n_chord:] == 0.0)
    fd, _ = central_difference(lambda u: float(constraint["fun"](u)[3]), u0, 1e-4)
    assert np.max(np.abs(fd - jac[3])) < 1e-10


def test_envelope_limits_match_the_documented_numbers(problem):
    """The chord ceiling/floor table in the implementation plan §4.4, to 3 s.f."""

    data = problem.envelope_data()
    r_over_R = data["radii"] / problem.parameterisation.radius_m
    tip = int(np.argmin(np.abs(r_over_R - 0.983)))
    root = int(np.argmin(np.abs(r_over_R - 0.167)))
    assert data["chord_max_m"][tip] == pytest.approx(0.148, abs=5e-4)
    assert data["chord_min_m"][root] == pytest.approx(0.145, abs=5e-4)
    assert data["reynolds_lo"] == 40_000.0
    assert data["reynolds_hi"] == 1_000_000.0


def test_solidity_constraint_requires_a_cap(problem):
    with pytest.raises(TypeError):
        problem.solidity_constraint()


def test_solidity_constraint_is_linear_when_a_cap_is_given(problem, u0):
    """Built, never run: only its shape and sign convention are checked."""

    constraint = problem.solidity_constraint(cap=1.0)
    g = constraint["fun"](u0)
    assert g.shape == (problem.parameterisation.n_stations,)
    assert np.all(g > 0.0)  # x0's solidity is far below 1
    assert constraint["jac"](u0).shape == (problem.parameterisation.n_stations, problem.n)


def test_jac_fd_returns_a_gradient_of_the_right_shape(problem, u0):
    grad = problem.jac_fd(u0, 1e-3)
    assert grad.shape == (problem.n,)
    assert np.all(np.isfinite(grad))


def test_alpha_check_at_x0_is_inside_the_cache(problem, u0):
    report = problem.alpha_check(u0)
    assert report["within"]
    assert report["all_converged"]
    assert report["alpha_limits_deg"] == [-8.0, 18.0]
    assert problem.reynolds_lo < report["reynolds_min"] < report["reynolds_max"] < problem.reynolds_hi


def test_fun_logs_and_reraises_a_polar_domain_error(parameterisation, bounds):
    fresh = ScaledProblem(parameterisation, bounds, WeibullResource.from_config())
    with pytest.raises(PolarDomainError):
        fresh.fun(np.ones(fresh.n))
    assert len(fresh.domain_errors) == 1
    assert fresh.domain_errors[0] == [1.0] * fresh.n
