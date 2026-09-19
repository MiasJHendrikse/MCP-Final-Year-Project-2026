"""
Phase 4, Step 2a -- the root-moment integrand `m` and the forward load path.

What is compared to what
-------------------------
`adjoint.kernels.station_partials` now also carries the flapwise root-moment
integrand `m = 1/2 rho (r - r_hub) c w^2 Cn` per blade and its partials
`dm_dphi, dm_dc, dm_dtheta` (docs/adjoint_derivation.md §10). The reference
on the other side is never the kernel:

  * the forward path `objective.loads.root_moment` -> `bem.rotor.solve_rotor`
    -> `root_bending_moment`, which reads the solver's own `w` and `Cn`;
  * a complex step (h = 1e-30, no subtraction) through the kernel's *value*
    code (`station_partials(..., derivatives=False)`), i.e. the same `m`
    expression the adjoint integrates.

Tolerances are the existing mixed ones (TOL_R = 1e-13, TOL_J = 1e-12);
do not invent looser ones. The per-station ``dm_dphi`` is checked at TOL_R.
The assembled ``dm_dd`` is checked at TOL_J because its twist columns pass
through ``polar.dcl_dalpha``, whose slope is only ~1e-13 relative, and the
``w^2 Cn_theta`` / ``2 w w_theta Cn`` cancellation in ``m_theta`` amplifies
that to ~6e-13 mixed -- the same w^2 structure that already places the
objective `q`/`J` partials at TOL_J. Failure here is a bug: fix it, never
loosen a tolerance.

The RootMomentSystem and the KS functional arrived in Step 2b. Step 2c adds
`ScaledProblem.moment_constraint` (the KS derivative with the adjoint
Jacobian) and the Tier 3 / sanity tests for it, at the bottom of this file.
The constraint is regrouped as `(KS0 - KS) - eps KS0`, algebraically the same
as `(1 - eps) KS0 - KS`, so the `eps = 0` slack at `u0` and the `eps`-linear
slack are exact to the bit.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os

import numpy as np
import pytest

from adjoint import BEMSystem, RootMomentSystem
from design import BladeParameterisation, DesignBounds
from gradients import ScaledProblem, central_difference
from objective import WeibullResource
from objective.loads import (
    ks,
    ks_weights,
    load_operating_points,
    root_bending_moment,
    root_moment,
    root_moments,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(_HERE, ".."))
X0_PATH = os.path.abspath(os.path.join(ROOT, "verification", "baseline", "x0.json"))
BASELINE_REFERENCE_PATH = os.path.abspath(
    os.path.join(ROOT, "verification", "baseline", "baseline_reference.json"))
ADJOINT_RESULT_PATH = os.path.abspath(
    os.path.join(ROOT, "verification", "adjoint_optimisation", "result.json"))

H = 1e-30
TOL_R = 1e-13
TOL_J = 1e-12
#: The committed global step of the A3 study (verification/fd_step_size).
H_STAR_GLOBAL = 3.162277660168379e-06

# The committed load operating set L (HANDOFF-2026-09-19 §3.1).
EXPECTED_LOAD_POINTS = [
    (3.5, 6.5),
    (4.5, 6.5),
    (5.5, 6.5),
    (6.5, 6.5),
    (7.5, 6.5),
    (8.5, 6.5),
    (9.5, 6.5),
    (10.5, 5.983986006837701),
    (11.0, 5.711986642890533),
]
RATED_TSR = 5.711986642890533
M_REFERENCE = 177.3755406092969


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
def x_pert(x0):
    """A fixed, recorded perturbation: +-2 % chord, +-0.5 deg twist on alternate
    control points, inside the box (HANDOFF-2026-09-19 §3.4)."""

    d = x0.copy()
    n_chord = 5
    for j in range(n_chord):
        d[j] *= 1.02 if j % 2 == 0 else 0.98
    for j in range(n_chord, len(d)):
        d[j] += math.radians(0.5 if j % 2 == 0 else -0.5)
    return d


@pytest.fixture(scope="module", params=["x0", "x_pert"])
def design(request, x0, x_pert):
    return x0 if request.param == "x0" else x_pert


@pytest.fixture(scope="module")
def state(system, design):
    return system.solve(design)


@pytest.fixture(scope="module")
def parts(system, state):
    return system.partials(state.phi, state.d)


@pytest.fixture(scope="module")
def problem():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    return ScaledProblem(parameterisation, bounds, WeibullResource.from_config())


@pytest.fixture(scope="module")
def u0(problem, x0):
    return problem.scaled(x0)


@pytest.fixture(scope="module")
def u_star():
    with open(ADJOINT_RESULT_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["u_star"], dtype=float)


@pytest.fixture(scope="module")
def load_system():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    return RootMomentSystem(parameterisation, bounds, WeibullResource.from_config(),
                            m_ref_nm=M_REFERENCE)


@pytest.fixture(scope="module")
def load_state(load_system, design):
    return load_system.solve(design)


@pytest.fixture(scope="module")
def load_parts(load_system, load_state):
    return load_system.partials(load_state.phi, load_state.d)


# ---------------------------------------------------------------------------
# forward path
# ---------------------------------------------------------------------------

def test_load_operating_points_is_the_nine_load_points():
    points = load_operating_points()
    assert points == EXPECTED_LOAD_POINTS


def test_root_moment_matches_baseline_reference(x0):
    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(x0)
    from config import load_design_rotor, load_site

    design = load_design_rotor()
    site = load_site()

    moment, result = root_moment(geometry, design.rated_wind_speed_ms, RATED_TSR,
                                 site.air_density, site.kinematic_viscosity)
    assert result["converged"]
    assert moment == pytest.approx(M_REFERENCE, abs=1e-9)


def test_root_moment_matches_adjoint_optimisation_reference(x0):
    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(x0)
    from config import load_design_rotor, load_site

    design = load_design_rotor()
    site = load_site()

    moment, _result = root_moment(geometry, design.rated_wind_speed_ms, RATED_TSR,
                                  site.air_density, site.kinematic_viscosity)

    with open(ADJOINT_RESULT_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    reference = artefact["loads_at_rated_ceiling"]["x0"]["root_bending_moment_nm"]
    assert moment == pytest.approx(reference, abs=1e-9)


def test_root_moments_over_the_load_set_run(x0):
    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(x0)
    from config import load_site

    site = load_site()
    moments, results = root_moments(geometry, EXPECTED_LOAD_POINTS,
                                    site.air_density, site.kinematic_viscosity)
    assert moments.shape == (len(EXPECTED_LOAD_POINTS),)
    assert all(result["converged"] for result in results)
    # The rated point (last) is the design condition and must reproduce M_ref.
    assert moments[-1] == pytest.approx(M_REFERENCE, abs=1e-9)


def test_kernel_r_hub_matches_geometry_r_hub(system, x0):
    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(x0)
    assert system.r_hub == pytest.approx(geometry.r_hub, abs=1e-12)


def test_kernel_m_integrates_to_the_forward_root_moment(x0):
    """The kernel's `m` value, trapezoid-integrated at the rated point, equals
    the forward path's `root_moment` -- a second code path for `w`, `Cn`."""

    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    from config import load_site

    site = load_site()
    rated_system = BEMSystem(parameterisation, bounds, WeibullResource.from_config(),
                             points=[(11.0, RATED_TSR)])
    state = rated_system.solve(x0)
    m = rated_system.partials(state.phi, state.d, derivatives=False)["m"][0]
    kernel_moment = float(sum(t * mi for t, mi in zip(rated_system.trapz_weights, m)))

    geometry = parameterisation.to_geometry(x0)
    forward_moment, _result = root_moment(geometry, 11.0, RATED_TSR,
                                          site.air_density, site.kinematic_viscosity)
    assert kernel_moment == pytest.approx(forward_moment, rel=1e-12)


def test_root_bending_moment_is_per_blade(x0):
    """The forward integrand has no blade count: it quotes one blade's load."""

    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(x0)
    from config import load_design_rotor, load_site

    design = load_design_rotor()
    site = load_site()
    moment, result = root_moment(geometry, design.rated_wind_speed_ms, RATED_TSR,
                                 site.air_density, site.kinematic_viscosity)
    stations = result["stations"]
    manual = root_bending_moment(stations, geometry.chord, site.air_density,
                                 geometry.n_blades, geometry.r_hub)
    assert moment == pytest.approx(manual, rel=1e-14)


def test_load_system_moments_equal_forward_path(load_system, load_state):
    """`RootMomentSystem.moments` (the kernel's `m`, trapezoid-integrated) is
    the forward path's `root_moments` (the solver's `w`, `Cn`) at every point
    of `L` -- a second, independent code path for the same physics."""

    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(load_state.d)
    from config import load_site

    site = load_site()
    forward_moments, results = root_moments(geometry, EXPECTED_LOAD_POINTS,
                                            site.air_density, site.kinematic_viscosity)
    assert all(result["converged"] for result in results)
    adjoint_moments = load_system.moments(load_state.phi, load_state.d)
    assert adjoint_moments.shape == forward_moments.shape
    assert np.allclose(adjoint_moments, forward_moments, rtol=1e-12)


def test_load_system_omega_matches_the_schedule(load_system):
    """The load system's operating points are `L`; `omega = lambda V / R`."""

    assert load_system.points == EXPECTED_LOAD_POINTS
    assert load_system.tsr_per_bin == [lam for _v, lam in EXPECTED_LOAD_POINTS]
    for b, (v, lam) in enumerate(EXPECTED_LOAD_POINTS):
        assert load_system.omega[b] == pytest.approx(lam * v / load_system.R, rel=1e-15)


def test_load_set_is_below_the_rating_at_x0(load_system, load_state):
    """The parent's `limited` mask is all-False on `L` at `x0` -- the KS
    functional carries every point, none of the AEP cap logic applies."""

    assert load_state.limited.shape == (len(EXPECTED_LOAD_POINTS),)
    assert not np.any(load_state.limited)


# ---------------------------------------------------------------------------
# Tier 1: m and its partials against complex steps of the value code
# ---------------------------------------------------------------------------

def test_dm_dphi_matches_complex_step(system, state, parts):
    """All 425 entries: a complex step of phi through the value path's `m`."""

    estimate = np.empty_like(parts["dm_dphi"])
    for b in range(system.n_points):
        for i in range(system.n_stations):
            phi = state.phi.astype(complex)
            phi[b, i] += 1j * H
            estimate[b, i] = system.partials(phi, state.d, derivatives=False)["m"].imag[b, i] / H
    _assert_mixed(estimate, parts["dm_dphi"], TOL_R, "dm/dphi")


def test_dm_dd_matches_complex_step(system, state, parts):
    """Every design variable: a complex `d` becomes complex chord, twist and
    Reynolds, so the Re(c) path is on the complex side too.  TOL_J (not
    TOL_R): the twist columns ride the polar's analytic slope and the
    w^2 Cn_theta cancellation, ~6e-13 mixed at x0 -- see the module docstring."""

    dm_dd = parts["dm_dc"][:, :, None] * system.N_c[None, :, :] \
        + parts["dm_dtheta"][:, :, None] * system.N_theta[None, :, :]
    assert dm_dd.shape == (system.n_points, system.n_stations, system.n_design)
    for j in range(system.n_design):
        d = state.d.astype(complex)
        d[j] += 1j * H
        estimate = system.partials(state.phi, d, derivatives=False)["m"].imag / H
        _assert_mixed(estimate, dm_dd[:, :, j], TOL_J, f"dm/dd_{j}")


def test_no_plateau_in_m_at_the_real_dtype_signature(system, state, parts):
    """The 1e-6 plateau is the fingerprint of a real-dtype array on a complex
    path; the mixed tolerance above already excludes it -- this names it."""

    estimate = np.empty_like(parts["dm_dphi"])
    for b in range(system.n_points):
        for i in range(system.n_stations):
            phi = state.phi.astype(complex)
            phi[b, i] += 1j * H
            estimate[b, i] = system.partials(phi, state.d, derivatives=False)["m"].imag[b, i] / H
    worst = np.max(np.abs(estimate - parts["dm_dphi"])
                   / np.maximum(1.0, np.abs(parts["dm_dphi"])))
    assert worst < 1e-12, f"worst mixed error {worst:.3e}: a real dtype on the path?"


# ---------------------------------------------------------------------------
# Tier 1: the KS functional's partials (RootMomentSystem) against complex
# steps of the value code, at x0 and x_pert
# ---------------------------------------------------------------------------

def test_dKS_dx_matches_complex_step(load_system, load_state, load_parts):
    """All 225 entries: a complex step of `phi` through `KS(phi, d)`."""

    dKS_dx = load_system.dKS_dx(load_state.phi, load_state.d, load_parts)
    estimate = np.empty_like(dKS_dx)
    for b in range(load_system.n_points):
        for i in range(load_system.n_stations):
            phi = load_state.phi.astype(complex)
            phi[b, i] += 1j * H
            estimate[b, i] = load_system.KS(phi, load_state.d).imag / H
    _assert_mixed(estimate, dKS_dx, TOL_J, "dKS/dphi")


def test_dKS_dd_matches_complex_step(load_system, load_state, load_parts):
    """Every design variable: a complex `d` through `KS(phi, d)` (TOL_J, the
    functional level, exactly as `dJ_dd`)."""

    dKS_dd = load_system.dKS_dd(load_state.phi, load_state.d, load_parts)
    for j in range(load_system.n_design):
        d = load_state.d.astype(complex)
        d[j] += 1j * H
        estimate = load_system.KS(load_state.phi, d).imag / H
        _assert_mixed(estimate, dKS_dd[j], TOL_J, f"dKS/dd_{j}")


# ---------------------------------------------------------------------------
# Tier 2: the forward-mode tangent and the inherited transpose
# ---------------------------------------------------------------------------

def test_tangent_matches_adjoint_direction(load_system, load_state):
    """`tangent(phi, d, v) == gradient(d).dKS_dd @ v` for three seeded `v`."""

    rng = np.random.default_rng(42)
    for _ in range(3):
        v = rng.normal(size=load_system.n_design)
        tangent = load_system.tangent(load_state.phi, load_state.d, v)
        adjoint = float(load_system.gradient(load_state.d).dKS_dd @ v)
        assert tangent == pytest.approx(adjoint, rel=1e-12)


def test_moment_transpose_is_the_parent_operator(load_system):
    """The matrix-free transpose is the parent's; its Tier 2 identity carries
    over unchanged."""

    assert RootMomentSystem.apply_dR_dd_T is BEMSystem.apply_dR_dd_T


# ---------------------------------------------------------------------------
# KS aggregate (the functional is Step 2b; the scalar is here because the
# forward path defines it)
# ---------------------------------------------------------------------------

def test_ks_bounds_and_weights():
    values = np.array([0.1, 0.3, 0.9, 0.5])
    for rho in (30.0, 100.0, 300.0):
        value = ks(values, rho)
        assert float(value) >= float(np.max(values))
        assert float(value) - float(np.max(values)) <= math.log(len(values)) / rho
        weights = ks_weights(values, rho)
        assert weights.shape == values.shape
        assert np.allclose(np.sum(weights), 1.0, rtol=1e-15)
        assert np.allclose(weights, np.exp(rho * (values - np.max(values)))
                           / np.sum(np.exp(rho * (values - np.max(values)))), rtol=1e-15)


def test_ks_is_complex_safe():
    """A complex step of `ks` is exactly the analytic softmax weights."""

    values = np.array([0.1, 0.3, 0.9, 0.5], dtype=float)
    rho = 100.0
    estimate = np.array([ks(values + 1j * H * (np.arange(len(values)) == j), rho).imag / H
                         for j in range(len(values))])
    assert np.allclose(estimate, ks_weights(values, rho), atol=1e-12)


def test_ks_decreases_toward_max_as_rho_increases():
    values = np.array([0.1, 0.3, 0.9, 0.5])
    at_30 = ks(values, 30.0)
    at_100 = ks(values, 100.0)
    at_300 = ks(values, 300.0)
    assert at_30 >= at_100 >= at_300
    assert at_300 >= float(np.max(values))


# ---------------------------------------------------------------------------
# Step 2c: the ScaledProblem moment constraint -- sanity
# ---------------------------------------------------------------------------

def test_moment_constraint_is_zero_slack_at_the_scaled_reference(problem, u0):
    """`eps = 0` caps the KS at the baseline's own value, so `g(u0) = 0`
    exactly -- the sanity gate the sign of the whole constraint rests on."""

    assert float(problem.moment_constraint(0.0)["fun"](u0)[0]) == 0.0


def test_moment_constraint_scales_exactly_with_eps_at_u0(problem, u0):
    """`g_eps(u0) = -eps KS0` to the bit, the regrouped form's whole point."""

    assert float(problem.moment_constraint(0.05)["fun"](u0)[0]) == -0.05 * problem.KS0


def test_moment_constraint_makes_the_unconstrained_optimum_infeasible(problem, u_star):
    """The unconstrained optimum carries +0.31 % moment, so the `eps = 0` cap
    is violated there, by exactly `KS(u*) - KS(u0)`."""

    g = float(problem.moment_constraint(0.0)["fun"](u_star)[0])
    assert g < 0.0
    assert abs(g) == problem.moment_ks(u_star) - problem.KS0


def test_constraints_with_moment_appends_the_moment_row(problem, u0):
    cons = problem.constraints_with_moment(0.0)
    assert len(cons) == 3
    assert all(c["type"] == "ineq" for c in cons)
    assert cons[2]["fun"](u0).shape == (1,)
    assert cons[2]["jac"](u0).shape == (1, problem.n)
    active = problem.active_set(u0, moment_eps=0.0)
    assert "moment" in active and active["moment"]["slack"] == 0.0


def test_moment_state_is_cached_between_fun_and_jac(problem, u0):
    """SLSQP calls `fun` and `jac` at the same `u`: one solve, not two."""

    constraint = problem.moment_constraint(0.0)
    constraint["fun"](u0)  # warm the reference into the cache
    before = problem.n_moment_solves
    constraint["fun"](u0)
    constraint["jac"](u0)
    assert problem.n_moment_solves == before

    rng = np.random.default_rng(3)
    u = u0 + 1e-6 * rng.normal(size=problem.n)
    fresh = problem.moment_constraint(0.0)
    before = problem.n_moment_solves
    fresh["fun"](u)
    fresh["jac"](u)
    assert problem.n_moment_solves == before + 1


def test_moment_report_quotes_the_rated_weight_and_conservatism(problem, u0):
    report = problem.moment_report(u0)
    assert report["softmax_weight_on_rated"] == pytest.approx(0.957, abs=5e-4)
    assert report["KS_by_rho"]["100"] == pytest.approx(problem.KS0, rel=1e-12)
    assert (report["KS_by_rho"]["30"] >= report["KS_by_rho"]["100"]
            >= report["KS_by_rho"]["300"])
    for key in ("30", "100", "300"):
        assert report["conservatism_by_rho"][key] >= 0.0
    assert report["m_ref_nm"] == pytest.approx(M_REFERENCE, rel=1e-14)


def test_cut_out_post_check_reports_alpha_above_the_xfoil_band(problem, u0):
    """The cut-out point is reported, labelled B3-dependent, never constrained."""

    cut = problem.moment_report(u0)["cut_out"]
    assert cut["B3_dependent"] is True
    assert cut["v_ms"] == 20.0
    assert cut["rpm"] == pytest.approx(300.0, rel=1e-9)
    assert cut["alpha_max_deg"] > 18.0
    assert cut["converged"]


# ---------------------------------------------------------------------------
# Step 2c: Tier 3 -- the constraint Jacobian against the FD noise floor
# ---------------------------------------------------------------------------

def test_moment_constraint_jacobian_matches_fd_at_x0(problem, u0):
    """`moment_constraint(0)["jac"](u0)` vs central FD of its `fun` at the
    committed global `h*` and `h* / sqrt(10)`, `h* sqrt(10)`, with `eps_j`
    re-measured as run_tier3.py does. Acceptance `|adj - fd| <= 3 eps_j`; the
    measured round-off floor (the hand-off's `delta_g / h*` at `h*/1000`) is
    reported in any failure, for the §3.4 xfail route."""

    constraint = problem.moment_constraint(0.0)
    jac = constraint["jac"](u0)[0]
    g = lambda u: float(constraint["fun"](u)[0])  # noqa: E731

    def fd(h):
        return central_difference(g, u0, h)[0]

    g_mid = fd(H_STAR_GLOBAL)
    g_lo = fd(H_STAR_GLOBAL / math.sqrt(10.0))
    g_hi = fd(H_STAR_GLOBAL * math.sqrt(10.0))
    eps = np.maximum(np.abs(g_mid - g_lo), np.abs(g_mid - g_hi))
    abs_err = np.abs(jac - g_mid)

    small = H_STAR_GLOBAL / 1000.0
    e = np.eye(problem.n)
    floor = np.array([abs(g(u0 + small * e[j]) - g(u0 - small * e[j])) / H_STAR_GLOBAL
                      for j in range(problem.n)])

    if not np.all(abs_err <= 3.0 * eps):
        worst = int(np.argmax(abs_err / eps))
        pytest.fail(
            f"moment constraint Jacobian outside 3 eps at variable {worst}: "
            f"|adj - fd| = {abs_err[worst]:.3e}, eps = {eps[worst]:.3e}, "
            f"ratio = {abs_err[worst] / eps[worst]:.3f}, measured floor = "
            f"{floor[worst]:.3e}, eps_below_floor = {bool(eps[worst] < floor[worst])}")


def test_moment_constraint_taylor_remainder_is_second_order(problem, u0):
    """Random unit `v`: the constraint's remainder falls ~100x per decade over
    1e-2, 1e-3, 1e-4; assert >= 30, exactly as run_tier3.py does."""

    constraint = problem.moment_constraint(0.0)
    jac = constraint["jac"](u0)[0]
    g = lambda u: float(constraint["fun"](u)[0])  # noqa: E731
    steps = (1e-2, 1e-3, 1e-4)

    rng = np.random.default_rng(11)
    envelope = problem.envelope_constraint()["fun"]
    for _attempt in range(50):
        v = rng.normal(size=problem.n)
        v /= np.linalg.norm(v)
        if all(np.all(envelope(u0 + h * v) > 0.0) for h in steps):
            break
    else:
        pytest.fail("no envelope-feasible direction found for the moment Taylor test")

    g0 = g(u0)
    slope = float(jac @ v)
    remainders = [abs(g(u0 + h * v) - g0 - h * slope) for h in steps]
    ratios = [remainders[0] / remainders[1], remainders[1] / remainders[2]]
    assert min(ratios) >= 30.0, (remainders, ratios)
