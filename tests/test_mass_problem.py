"""
The mass problem, Step 3 -- the rows `ScaledProblem` hands SLSQP
(`docs/PLAN-mass-objective-2026-09-20.md`): the material objective, the AEP
floor, the root-stress proxy, the manufacturability rows, the assembly, the
shared load solve and the failed-trial guard. The deflection row's own
Tier 3 is in `tests/test_deflection.py`.

Every relative row is zero-slack at `u0` to the bit; the stress row's
Jacobian is checked against central FD at the committed `h*` on the Phase 4
pattern (floor reported only on failure); the guard returns `-FAILED_SLACK`
and logs, and never substitutes a Jacobian.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import dataclasses
import json
import math
import os

import numpy as np
import pytest

import config
from design import BladeParameterisation, DesignBounds
from gradients import ScaledProblem, central_difference
from gradients.problem import (ABSOLUTE_ROWS, ALL_MASS_PROBLEM_ROWS, FAILED_SLACK,
                               MASS_PROBLEM_ROWS)
from objective import WeibullResource
from objective.mass import section_coefficients
from polars.interpolant import PolarDomainError

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(_HERE, ".."))
XC_PATH = os.path.join(ROOT, "verification", "load_constraint", "result_eps0.json")
XM_PATH = os.path.join(ROOT, "verification", "mass_optimisation", "result_delta0.json")

H_STAR_GLOBAL = 3.162277660168379e-06


@pytest.fixture(scope="module")
def problem():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    problem = ScaledProblem(parameterisation, bounds, WeibullResource.from_config())
    problem.set_reference(problem.scaled(problem.reference_design()))
    return problem


@pytest.fixture(scope="module")
def u0(problem):
    return np.asarray(problem.u0, dtype=float)


@pytest.fixture(scope="module")
def u_c(problem):
    with open(XC_PATH, encoding="utf-8") as handle:
        return problem.scaled(np.array(json.load(handle)["x_c"], dtype=float))


@pytest.fixture(scope="module")
def u_m(problem):
    with open(XM_PATH, encoding="utf-8") as handle:
        return problem.scaled(np.array(json.load(handle)["x_m"], dtype=float))


def _tier3(problem, constraint, u, label):
    """The `test_moment_constraint_jacobian_matches_fd_at_x0` acceptance, verbatim."""

    jac = constraint["jac"](u)[0]
    g = lambda v: float(constraint["fun"](v)[0])  # noqa: E731

    def fd(h):
        return central_difference(g, u, h)[0]

    g_mid = fd(H_STAR_GLOBAL)
    g_lo = fd(H_STAR_GLOBAL / math.sqrt(10.0))
    g_hi = fd(H_STAR_GLOBAL * math.sqrt(10.0))
    eps = np.maximum(np.abs(g_mid - g_lo), np.abs(g_mid - g_hi))
    abs_err = np.abs(jac - g_mid)

    if not np.all(abs_err <= 3.0 * eps):
        worst = int(np.argmax(abs_err / eps))
        ts = np.linspace(-4e-12, 4e-12, 9)
        e = np.zeros(problem.n)
        e[worst] = 1.0
        values = np.array([g(u + t * e) for t in ts])
        residual = values - np.polyval(np.polyfit(ts, values, 1), ts)
        floor = float(np.sqrt(np.mean(residual ** 2))) / H_STAR_GLOBAL
        pytest.fail(
            f"{label} Jacobian outside 3 eps at variable {worst}: "
            f"|adj - fd| = {abs_err[worst]:.3e}, eps = {eps[worst]:.3e}, "
            f"ratio = {abs_err[worst] / eps[worst]:.3f}, measured round-off floor = "
            f"{floor:.3e}, eps_below_floor = {bool(eps[worst] < floor)}")


def _taylor(problem, constraint, u):
    jac = constraint["jac"](u)[0]
    g = lambda v: float(constraint["fun"](v)[0])  # noqa: E731
    steps = (1e-2, 1e-3, 1e-4)
    rng = np.random.default_rng(11)
    envelope = problem.envelope_constraint()["fun"]
    for _attempt in range(50):
        v = rng.normal(size=problem.n)
        v /= np.linalg.norm(v)
        if all(np.all(envelope(u + h * v) > 0.0) for h in steps):
            break
    else:
        pytest.fail("no envelope-feasible direction found")
    g0 = g(u)
    slope = float(jac @ v)
    remainders = [abs(g(u + h * v) - g0 - h * slope) for h in steps]
    ratios = [remainders[0] / remainders[1], remainders[1] / remainders[2]]
    assert min(ratios) >= 30.0, (remainders, ratios)


# ---------------------------------------------------------------------------
# the objective
# ---------------------------------------------------------------------------

def test_mass_is_exactly_one_at_the_reference(problem, u0):
    assert problem.mass(u0) == 1.0
    assert problem.material_ref == pytest.approx(0.43705266291483974, rel=1e-12)


def test_mass_jac_matches_central_differences(problem, u0, u_c):
    for u in (u0, u_c):
        fd, _n = central_difference(problem.mass, u, 1e-4)
        assert np.allclose(problem.mass_jac(u), fd, rtol=1e-9, atol=1e-12)
    assert np.array_equal(problem.mass_jac(u0), problem.mass_jac(u_c))  # linear objective


def test_material_report_is_relative_to_the_reference(problem, u0, u_c):
    assert problem.material_report(u0)["shell_over_reference"] == 1.0
    report = problem.material_report(u_c)
    assert report["shell_over_reference"] == pytest.approx(1.06171, abs=5e-5)
    assert report["solid_over_reference"] == pytest.approx(1.14196, abs=5e-5)


# ---------------------------------------------------------------------------
# the AEP floor
# ---------------------------------------------------------------------------

def test_aep_floor_is_zero_slack_at_the_reference(problem, u0):
    assert float(problem.aep_floor_constraint(0.0)["fun"](u0)[0]) == 0.0
    assert float(problem.aep_floor_constraint(0.005)["fun"](u0)[0]) == pytest.approx(0.005, abs=1e-15)


def test_aep_floor_jacobian_is_minus_the_objective_adjoint(problem, u0):
    jac = problem.aep_floor_constraint(0.0)["jac"](u0)
    assert jac.shape == (1, problem.n)
    assert np.array_equal(jac[0], -problem.jac_adjoint(u0))


def test_aep_floor_needs_the_reference_set():
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=5, n_twist=5)
    fresh = ScaledProblem(parameterisation, bounds, WeibullResource.from_config())
    with pytest.raises(RuntimeError, match="J0 not set"):
        fresh.aep_floor_constraint(0.0)
    with pytest.raises(ValueError, match="delta"):
        fresh.set_reference(fresh.scaled(fresh.reference_design()))
        fresh.aep_floor_constraint(1.0)


def test_the_energy_optimum_has_the_committed_slack(problem, u_c):
    """`x_c` sits +0.14649 % above `x0`: the floor's slack at `x_c` is that number."""

    assert float(problem.aep_floor_constraint(0.0)["fun"](u_c)[0]) == pytest.approx(0.0014649, abs=2e-7)


# ---------------------------------------------------------------------------
# the root-stress proxy
# ---------------------------------------------------------------------------

def test_stress_is_zero_slack_at_the_reference(problem, u0):
    assert float(problem.stress_constraint()["fun"](u0)[0]) == 0.0
    assert problem.stress_ratio(u0) == 1.0


def test_stress_ratio_falls_with_a_fatter_root(problem, u_c):
    """`x_c` holds `KS = KS0` on a 300 mm root against Schmitz's 276 mm: `(276/300)^2`."""

    ratio = problem.stress_ratio(u_c)
    c00 = float(problem.physical(problem.u0)[0])
    c0 = float(problem.physical(u_c)[0])
    ks = problem.moment_ks(u_c) / problem.KS0
    assert ratio == pytest.approx(ks * (c00 / c0) ** 2, rel=1e-12)
    assert 0.84 < ratio < 0.85


def test_stress_jacobian_matches_fd_at_x0(problem, u0):
    _tier3(problem, problem.stress_constraint(), u0, "stress constraint")


def test_stress_taylor_remainder_is_second_order(problem, u0):
    _taylor(problem, problem.stress_constraint(), u0)


# ---------------------------------------------------------------------------
# the absolute rows (2026-09-20, evening)
# ---------------------------------------------------------------------------

def test_absolute_stress_is_the_ks_moment_over_the_thin_shell_section_modulus(problem, u0, u_c):
    """
    `sigma = KS M_ref / (k_Z c0^2 t)` with the recorded 2 mm skin: 22.6 MPa
    at `x0`, against a 196.5 MPa design allowable -- slack by a factor of
    nine at the operating loads (problem.py, "which rows are the design
    rows"). `x_c`'s fatter root carries less.
    """

    design = config.load_design_rotor()
    k_z = section_coefficients().section_modulus
    for u in (u0, u_c):
        c0 = float(problem.physical(u)[0])
        expected = problem.moment_ks(u) * problem.m_ref / (k_z * c0 ** 2 * design.shell_thickness_m)
        assert problem.stress_absolute_pa(u) == pytest.approx(expected, rel=1e-12)
    assert problem.stress_absolute_pa(u0) / 1e6 == pytest.approx(22.588, abs=2e-3)
    assert problem.stress_absolute_pa(u_c) < problem.stress_absolute_pa(u0)
    assert problem.design_allowable_stress_pa() / 1e6 == pytest.approx(196.50, abs=1e-2)
    assert problem.stress_absolute_pa(u_c) / problem.stress_absolute_pa(u0) == pytest.approx(
        problem.stress_ratio(u_c), rel=1e-12)     # the same functional, constants carried


def test_absolute_stress_row_is_the_inequality_over_its_allowable(problem, u0):
    g = float(problem.stress_constraint_absolute()["fun"](u0)[0])
    sigma = problem.stress_absolute_pa(u0)
    allowable = problem.design_allowable_stress_pa()
    assert g == pytest.approx(1.0 - sigma / allowable, rel=1e-12)
    assert 0.88 < g < 0.89
    slacks = problem.mass_problem_slacks(u0, 0.0)["stress_absolute"]
    assert slacks["available"] is True and not slacks["active"]
    assert slacks["stress_mpa"] == pytest.approx(sigma / 1e6, rel=1e-12)
    assert slacks["slack_mpa"] == pytest.approx((allowable - sigma) / 1e6, rel=1e-9)


def test_absolute_stress_jacobian_matches_fd_at_x0_and_x_m(problem, u0, u_m):
    _tier3(problem, problem.stress_constraint_absolute(), u0, "absolute stress at x0")
    _tier3(problem, problem.stress_constraint_absolute(), u_m, "absolute stress at x_m")


def test_absolute_stress_taylor_remainder_is_second_order(problem, u0):
    _taylor(problem, problem.stress_constraint_absolute(), u0)


def test_absolute_stress_jacobian_is_the_relative_one_scaled(problem, u_c):
    """One chain rule, two constants: `M_ref / (k_Z t sigma_design)`."""

    design = config.load_design_rotor()
    scale = problem.m_ref / (section_coefficients().section_modulus * design.shell_thickness_m
                             * problem.design_allowable_stress_pa())
    relative = problem.stress_constraint()["jac"](u_c)
    absolute = problem.stress_constraint_absolute()["jac"](u_c)
    assert np.allclose(absolute, scale * relative, rtol=1e-12, atol=0.0)


def test_absolute_deflection_is_the_unit_deflection_over_e_k_i_t(problem, u0, u_c):
    """
    `delta_tip = delta_ref D / (E k_I t)`: 111.8 mm at `x0`'s rated point
    with the 2 mm E-glass skin (the KS aggregate; the rated point alone is
    `delta_ref / (E k_I t)` = 111.79 mm), 85 mm for `x_c`.
    """

    design = config.load_design_rotor()
    factor = design.youngs_modulus_pa * section_coefficients().second_moment * design.shell_thickness_m
    assert problem.deflection_absolute_m(u0) == pytest.approx(
        problem.delta_ref * problem.D0 / factor, rel=1e-12)
    assert problem.deflection_absolute_m(u0) * 1e3 == pytest.approx(111.84, abs=1e-2)
    assert problem.delta_ref / factor * 1e3 == pytest.approx(111.79, abs=1e-2)
    assert problem.deflection_absolute_m(u_c) * 1e3 == pytest.approx(85.2, abs=0.1)
    slacks = problem.mass_problem_slacks(u0, 0.0)["deflection_absolute"]
    assert slacks["tip_deflection_mm"] == pytest.approx(111.84, abs=1e-2)


def test_absolute_deflection_row_refuses_to_assemble_while_the_clearance_is_todo(problem):
    available = problem.absolute_rows_available()
    assert available["stress_absolute"] is True
    assert available["deflection_absolute"].startswith("TODO: tip_clearance_m")
    with pytest.raises(config.UnresolvedConfigError, match="tip_clearance_m"):
        problem.deflection_constraint_absolute()
    with pytest.raises(config.UnresolvedConfigError, match="tip_clearance_m"):
        problem.mass_problem_rows(0.0, include=("deflection_absolute",))
    slacks = problem.mass_problem_slacks(problem.u0, 0.0)["deflection_absolute"]
    assert slacks["available"].startswith("TODO") and "slack" not in slacks


@pytest.fixture(scope="module")
def problem_with_clearance():
    """A problem whose design carries an INJECTED 150 mm clearance -- a test value, not a decision."""

    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    problem = ScaledProblem(parameterisation, bounds, WeibullResource.from_config())
    problem.design = dataclasses.replace(problem.design, tip_clearance_m=0.150)
    problem.set_reference(problem.scaled(problem.reference_design()))
    return problem


def test_absolute_deflection_row_with_an_injected_clearance(problem_with_clearance, u0, u_m):
    """
    With a clearance supplied the row assembles, its slack is
    `1 - delta_tip / clearance`, its Jacobian passes Tier 3 at `x0` and
    `x_m`, and `ALL_MASS_PROBLEM_ROWS` assembles in order.
    """

    p = problem_with_clearance
    assert p.absolute_rows_available()["deflection_absolute"] is True
    row = p.deflection_constraint_absolute()
    g = float(row["fun"](u0)[0])
    assert g == pytest.approx(1.0 - p.deflection_absolute_m(u0) / 0.150, rel=1e-12)
    assert 0.25 < g < 0.26
    _tier3(p, row, u0, "absolute deflection at x0")
    _tier3(p, row, u_m, "absolute deflection at x_m")
    _taylor(p, row, u0)
    scale = p.delta_ref / (p.material_model().stiffness_factor_pa_m() * 0.150)
    assert np.allclose(row["jac"](u_m), scale * p.deflection_constraint()["jac"](u_m),
                       rtol=1e-12, atol=0.0)

    rows = p.mass_problem_rows(0.0, include=ALL_MASS_PROBLEM_ROWS)
    assert [label for label, _c in rows] == list(ALL_MASS_PROBLEM_ROWS)
    slacks = p.mass_problem_slacks(u0, 0.0)["deflection_absolute"]
    assert slacks["clearance_mm"] == 150.0 and not slacks["active"]
    assert slacks["slack_mm"] == pytest.approx(150.0 - slacks["tip_deflection_mm"], rel=1e-9)


def test_absolute_rows_are_additional_and_slack_at_the_committed_optimum(problem, u_m):
    """
    `MASS_PROBLEM_ROWS` is the committed relative set, unchanged; the
    absolute stress row joins by name and is slack at `x_m` (the operating
    loads do not size this root), so the committed optimum is unchanged by
    it -- the reason the relative rows stay the design rows.
    """

    assert MASS_PROBLEM_ROWS == ("envelope", "solidity", "manufacturing", "aep_floor",
                                 "moment", "stress", "deflection")
    assert ABSOLUTE_ROWS == ("stress_absolute", "deflection_absolute")
    assert ALL_MASS_PROBLEM_ROWS == MASS_PROBLEM_ROWS + ABSOLUTE_ROWS
    rows = problem.mass_problem_rows(0.0, include=MASS_PROBLEM_ROWS + ("stress_absolute",))
    assert [label for label, _c in rows][-1] == "stress_absolute"
    assert len(problem.constraints_for_mass_problem(0.0)) == 7      # the default is the committed set
    g = float(dict(rows)["stress_absolute"]["fun"](u_m)[0])
    assert g > 0.88
    assert problem.mass_kg(u_m) == pytest.approx(1.6216, abs=5e-4)


# ---------------------------------------------------------------------------
# manufacturability
# ---------------------------------------------------------------------------

def test_manufacturing_rows_are_the_control_point_differences(problem, u0):
    matrix, labels = problem.manufacturing_rows()
    assert matrix.shape == (8, problem.n)
    assert labels[0] == "monotone chord_0-chord_1" and labels[4] == "monotone twist_0-twist_1"
    for k in range(4):
        assert matrix[k, k] == 1.0 and matrix[k, k + 1] == -1.0
        assert matrix[4 + k, 5 + k] == 1.0 and matrix[4 + k, 6 + k] == -1.0
    assert np.count_nonzero(matrix) == 16

    constraint = problem.manufacturing_constraints()
    d = problem.physical(u0)
    expected = np.concatenate([d[:5][:-1] - d[:5][1:], d[5:][:-1] - d[5:][1:]])
    values = constraint["fun"](u0)
    assert np.allclose(values[:8], expected, rtol=1e-14, atol=1e-16)
    assert np.array_equal(constraint["jac"](u0)[:8], matrix * problem.bounds.span()[None, :])


def test_min_chord_floor_is_a_row_on_every_chord_control_point(problem, u0):
    """`manufacturing.min_chord_m` (0.060, 2026-09-20) as rows `c_i - 0.060 >= 0`:
    Schmitz clears it by 7 mm at the tip; a 50 mm tip violates it; the box
    bound (0.045) is untouched so the Phase 1-4 scaling is the committed one."""

    matrix, rhs, labels = problem.min_chord_rows()
    assert matrix.shape == (5, problem.n) and np.array_equal(rhs, np.full(5, 0.060))
    assert labels == [f"min chord chord_{i}" for i in range(5)]
    assert problem.bounds.chord_min_m == pytest.approx(0.045)
    assert problem.manufacturing_row_labels()[8:] == labels

    d = problem.physical(u0)
    values = problem.manufacturing_constraints()["fun"](u0)
    assert values.shape == (13,)
    assert np.allclose(values[8:], d[:5] - 0.060, rtol=1e-14, atol=1e-16)
    assert values[12] == pytest.approx(0.067 - 0.060, abs=1e-3)      # Schmitz tip 67 mm

    thin = d.copy()
    thin[4] = 0.050
    assert problem.manufacturing_constraints()["fun"](problem.scaled(thin))[12] < 0.0
    slacks = problem.mass_problem_slacks(problem.scaled(thin), 0.0)
    assert "min chord chord_4" in slacks["manufacturing"]["rows"]


def test_schmitz_is_strictly_monotone_and_a_wavy_blade_is_not(problem, u0):
    constraint = problem.manufacturing_constraints()
    assert np.all(constraint["fun"](u0) > 1e-3)
    wavy = problem.physical(u0).copy()
    wavy[1], wavy[2] = wavy[2], wavy[1]
    assert np.any(constraint["fun"](problem.scaled(wavy)) < 0.0)


# ---------------------------------------------------------------------------
# the assembly, the shared solve, the guard
# ---------------------------------------------------------------------------

def test_mass_problem_rows_are_in_order_and_zero_slack_at_x0(problem, u0):
    rows = problem.mass_problem_rows(0.0)
    assert [label for label, _c in rows] == list(MASS_PROBLEM_ROWS)
    assert len(problem.constraints_for_mass_problem(0.0)) == 7
    for label, constraint in rows:
        g = constraint["fun"](u0)
        assert constraint["jac"](u0).shape == (len(g), problem.n)
        if label in ("aep_floor", "moment", "stress", "deflection"):
            assert float(g[0]) == 0.0, label
        else:
            assert np.all(g > 0.0), label


def test_ablation_drops_rows_without_reordering(problem):
    rows = problem.mass_problem_rows(0.0, include=("aep_floor", "envelope", "manufacturing"))
    assert [label for label, _c in rows] == ["envelope", "manufacturing", "aep_floor"]
    with pytest.raises(ValueError, match="unknown"):
        problem.mass_problem_rows(0.0, include=("buckling",))


def test_the_energy_optimum_is_feasible_for_the_mass_problem(problem, u_c):
    """`x_c` (the counter-example) satisfies every state row at `delta = 0`: the
    moment cap is active (Phase 4), the stress and deflection rows are slack,
    the monotone rows are slack. The one row it violates is the 60 mm
    min-chord floor of 2026-09-20 (its tip is 47.7 mm, from the 45 mm box it
    was optimised under): recorded, not patched -- it is the energy
    reference, evaluated as it is."""

    slacks = problem.mass_problem_slacks(u_c, 0.0)
    assert not slacks["aep_floor"]["active"] and slacks["aep_floor"]["slack"] > 0.0
    assert slacks["moment"]["active"]
    assert slacks["stress"]["slack"] > 0.0 and slacks["stress"]["stress_ratio"] < 0.9
    assert slacks["deflection"]["slack"] > 0.0 and 0.7 < slacks["deflection"]["deflection_ratio"] < 0.8
    assert slacks["manufacturing"]["rows"] == ["min chord chord_4"]
    assert slacks["manufacturing"]["slack_min"] == pytest.approx(0.0477 - 0.060, abs=1e-3)
    report = problem.active_set(u_c, mass_delta=0.0)
    assert report["mass_problem"]["moment"]["active"]


def test_one_load_solve_serves_the_moment_stress_and_deflection_rows(problem, u_c):
    rng = np.random.default_rng(5)
    u = u_c + 1e-6 * rng.normal(size=problem.n)
    rows = dict(problem.mass_problem_rows(0.0))
    before = (problem.n_moment_solves, problem.n_deflection_solves)
    rows["moment"]["fun"](u)
    rows["stress"]["fun"](u)
    rows["stress"]["jac"](u)
    rows["deflection"]["fun"](u)
    rows["deflection"]["jac"](u)
    assert (problem.n_moment_solves, problem.n_deflection_solves) == (before[0] + 1, before[1])

    # And the other way round: a deflection evaluation first is reused by the moment row.
    u = u_c + 2e-6 * rng.normal(size=problem.n)
    before = (problem.n_moment_solves, problem.n_deflection_solves)
    rows["deflection"]["fun"](u)
    rows["moment"]["fun"](u)
    assert (problem.n_moment_solves, problem.n_deflection_solves) == (before[0], before[1] + 1)


def test_guard_reports_a_violation_and_logs_but_never_a_stale_jacobian(problem):
    """`u = 1` is outside the polar cache (`test_fun_logs_and_reraises_a_polar_domain_error`)."""

    bad = np.ones(problem.n)
    rows = dict(problem.mass_problem_rows(0.0))
    n_before = len(problem.evaluation_failures)
    for label in ("aep_floor", "moment", "stress", "deflection"):
        g = rows[label]["fun"](bad)
        assert g.shape == (1,) and float(g[0]) == -FAILED_SLACK, label
        assert problem.evaluation_failures[-1]["row"] == label
        assert problem.evaluation_failures[-1]["u"] == [1.0] * problem.n
        with pytest.raises(PolarDomainError):
            rows[label]["jac"](bad)
    assert len(problem.evaluation_failures) == n_before + 4
    with pytest.raises(PolarDomainError):
        problem.mass_problem_slacks(bad, 0.0)  # unguarded, by design


def test_a_failed_load_solve_is_remembered_not_repeated_by_the_next_row(problem, u_c):
    """Review 2026-09-20: the three state rows share one forward solve; when
    it fails at a trial point, the second and third row re-raise the remembered
    failure instead of running the same failing solve again."""

    bad = np.ones(problem.n)
    problem._load_failure = None  # the guard test above already failed at this point
    rows = dict(problem.mass_problem_rows(0.0))
    before = (problem.n_moment_solves, problem.n_deflection_solves)
    for label in ("moment", "stress", "deflection"):
        g = rows[label]["fun"](bad)
        assert float(g[0]) == -FAILED_SLACK, label
    with pytest.raises(PolarDomainError):
        rows["deflection"]["jac"](bad)
    assert (problem.n_moment_solves, problem.n_deflection_solves) == (before[0] + 1, before[1])

    # A different point is not poisoned by the remembered failure.
    good = u_c + 1e-6 * np.random.default_rng(7).normal(size=problem.n)
    assert float(rows["moment"]["fun"](good)[0]) > -FAILED_SLACK
    assert problem.n_moment_solves == before[0] + 2


def test_run_mass_slsqp_rejects_an_unknown_row_name():
    """Review 2026-09-20: a mistyped `--rows` token used to be dropped silently."""

    import importlib
    import sys
    here = os.path.join(ROOT, "verification", "mass_optimisation")
    sys.path.insert(0, here)
    try:
        runner = importlib.import_module("run_mass_slsqp")
    finally:
        sys.path.remove(here)
    assert runner.parse_rows("envelope, solidity,aep_floor") == ("envelope", "solidity", "aep_floor")
    assert runner.parse_rows(",".join(MASS_PROBLEM_ROWS)) == tuple(MASS_PROBLEM_ROWS)
    with pytest.raises(ValueError, match="solidty"):
        runner.parse_rows("envelope,solidty,aep_floor")
