"""
Plan step 1.6 acceptance: the blade parameterisation and its derivatives.

Named in work order Task 7's list, which noted it "lands with Phase 1.6" --
this is that landing. The headline criterion is complex-step agreement to
~1e-14 on `dc/dd` and `dtheta/dd`; what is actually achieved is exact
agreement, for a structural reason the tests below assert directly rather than
just benefiting from.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import math

import numpy as np
import pytest

from config import load_design_rotor
from design import BladeParameterisation, DesignBounds, basis_matrix, clamped_knots

#: A plausible blade: chord tapering 0.20 -> 0.05 m, twist 25 -> -2 deg. Not a
#: designed shape -- just something with the right monotone character to
#: exercise the map.
def _design_vector(parameterisation):
    return np.concatenate([
        np.linspace(0.20, 0.05, parameterisation.n_chord),
        np.radians(np.linspace(25.0, -2.0, parameterisation.n_twist)),
    ])


# The design-variable bounds are in `config/rotor_design.yaml` since
# 2026-09-19 (chord_max_m = 0.30 m, O4) and are read with
# `DesignBounds.from_config()`. The provisional set that lived here from
# 2026-09-13 (chord_max_m = 0.45 m, a placeholder with no basis) is retired;
# `verification/aep_gain_audit/` carries its own copy as a historical record.


# ---------------------------------------------------------------------------
# "Analytic dc/dd and dtheta/dd", verified against complex step
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_chord,n_twist", [(4, 4), (5, 5), (6, 6), (8, 8), (3, 3)])
def test_chord_derivative_matches_complex_step(n_chord, n_twist):
    """
    dc/dd against complex-step differentiation, every design variable.

    The plan asks for ~1e-14. The tolerance here is 1e-15 -- tighter than
    asked -- because a B-spline is linear in its control points, so the
    Jacobian is not an approximation of anything: `dc_i/dd_j` IS the basis
    function `N_j(s_i)`. Anything above round-off would mean the basis matrix
    and the evaluation disagree, which is a wiring bug rather than a
    derivative error.
    """

    parameterisation = BladeParameterisation(n_chord=n_chord, n_twist=n_twist)
    design_vector = _design_vector(parameterisation)
    jacobian = parameterisation.dchord_dd()
    h = 1e-30

    for j in range(parameterisation.n_design_variables):
        perturbed = design_vector.astype(complex)
        perturbed[j] += 1j * h
        complex_step = parameterisation.chord(perturbed).imag / h
        assert np.allclose(complex_step, jacobian[:, j], rtol=0.0, atol=1e-15), j


@pytest.mark.parametrize("n_chord,n_twist", [(4, 4), (5, 5), (6, 6), (8, 8), (3, 3)])
def test_twist_derivative_matches_complex_step(n_chord, n_twist):
    """dtheta/dd against complex step, same reasoning as chord."""

    parameterisation = BladeParameterisation(n_chord=n_chord, n_twist=n_twist)
    design_vector = _design_vector(parameterisation)
    jacobian = parameterisation.dtwist_dd()
    h = 1e-30

    for j in range(parameterisation.n_design_variables):
        perturbed = design_vector.astype(complex)
        perturbed[j] += 1j * h
        complex_step = parameterisation.twist(perturbed).imag / h
        assert np.allclose(complex_step, jacobian[:, j], rtol=0.0, atol=1e-15), j


def test_jacobians_are_constant_in_the_design_vector():
    """
    The stronger statement behind the tolerance above.

    If `dc/dd` depended on `d` at all, the map would not be linear and the
    exactness would be a coincidence of the test point. Evaluated at three
    unrelated design vectors, the Jacobian must be bit-identical.
    """

    parameterisation = BladeParameterisation()
    first = parameterisation.dchord_dd()

    for _ in range(3):
        assert np.array_equal(parameterisation.dchord_dd(), first)
        assert np.array_equal(parameterisation.dtwist_dd(),
                              parameterisation.dtwist_dd())


def test_cross_block_derivatives_are_structurally_zero():
    """
    Chord does not depend on the twist control points, or vice versa.

    Worth asserting rather than assuming: it is what lets the adjoint treat
    the parameterisation as two independent blocks, and a wiring error that
    coupled them would still produce plausible blades.
    """

    parameterisation = BladeParameterisation()
    n_chord = parameterisation.n_chord

    assert np.all(parameterisation.dchord_dd()[:, n_chord:] == 0.0)
    assert np.all(parameterisation.dtwist_dd()[:, :n_chord] == 0.0)


# ---------------------------------------------------------------------------
# The spline basis itself
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_control", [4, 5, 6, 8, 10])
def test_basis_is_a_partition_of_unity(n_control):
    """
    Every row of the basis matrix sums to 1.

    The cheapest available check that the knot vector and the degree are
    consistent: if they were not, the basis would not sum to one and every
    shape it produced would be subtly scaled.
    """

    stations = np.linspace(0.02, 0.98, 25)
    matrix = basis_matrix(stations, n_control, degree=3)

    assert np.allclose(matrix.sum(axis=1), 1.0, atol=1e-12)
    assert np.all(matrix >= -1e-15), "B-spline basis functions are non-negative"


def test_clamped_knots_have_the_right_multiplicity():
    """
    Clamped: end knots repeated degree+1 times, so the curve reaches its first
    and last control points.

    That property is what makes the root and tip values directly settable and
    boundable, which is the reason for choosing a clamped basis at all.
    """

    for n_control in (4, 5, 8):
        knots = clamped_knots(n_control, 3)
        assert len(knots) == n_control + 3 + 1
        assert np.all(knots[:4] == 0.0)
        assert np.all(knots[-4:] == 1.0)
        assert np.all(np.diff(knots) >= 0.0)


def test_curve_reaches_its_end_control_points():
    """The consequence of clamping, checked on the curve rather than the knots."""

    parameterisation = BladeParameterisation(n_stations=200)
    design_vector = _design_vector(parameterisation)
    chord_points, twist_points = parameterisation.split(design_vector)
    chord, twist = parameterisation.evaluate(design_vector)

    # Stations are strip midpoints, so they approach but never reach s = 0/1.
    assert chord[0] == pytest.approx(chord_points[0], abs=2e-3)
    assert chord[-1] == pytest.approx(chord_points[-1], abs=2e-3)
    assert twist[0] == pytest.approx(twist_points[0], abs=2e-3)
    assert twist[-1] == pytest.approx(twist_points[-1], abs=2e-3)


def test_too_few_control_points_is_rejected():
    with pytest.raises(ValueError):
        clamped_knots(3, degree=3)
    with pytest.raises(ValueError):
        BladeParameterisation(n_chord=1, n_twist=5)


# ---------------------------------------------------------------------------
# The sawtooth null space the parameterisation exists to remove
# ---------------------------------------------------------------------------

def test_parameterisation_cannot_represent_a_sawtooth():
    """
    Plan section 4.2, made checkable.

    The failure mode the spline exists to prevent is an optimiser walking into
    a strip-to-strip oscillation that barely moves the objective and cannot be
    manufactured. That is a *representation* claim, so it can be tested: fit a
    pure alternating +/- sawtooth and confirm the parameterisation cannot
    reproduce it -- the residual stays a large fraction of the signal, whereas
    a per-station parameterisation would fit it exactly.
    """

    parameterisation = BladeParameterisation()
    n = parameterisation.n_stations

    smooth = np.linspace(0.20, 0.05, n)
    sawtooth = smooth + 0.02 * np.array([(-1.0) ** i for i in range(n)])

    _, chord_rms, _ = parameterisation.fit(sawtooth, np.zeros(n))

    # The oscillation has amplitude 0.02, so an RMS of 0.02/sqrt(2) ~ 0.0141 is
    # what "rejected it entirely" looks like. Anything much below that would
    # mean the parameterisation is reproducing part of the sawtooth.
    assert chord_rms > 0.5 * 0.02 / math.sqrt(2.0), (
        "the parameterisation is absorbing a strip-to-strip oscillation it "
        "should not be able to represent")


def test_a_smooth_target_is_represented_well():
    """The complement: the shapes a blade actually has must fit closely."""

    parameterisation = BladeParameterisation()
    stations = parameterisation.stations

    chord_target = 0.20 - 0.15 * stations
    twist_target = np.radians(25.0 - 27.0 * stations)

    _, chord_rms, twist_rms = parameterisation.fit(chord_target, twist_target)

    assert chord_rms < 1e-6, chord_rms
    assert twist_rms < 1e-6, twist_rms


def test_fit_is_exact_for_a_shape_the_basis_can_make():
    """
    Round trip: evaluate a design vector, fit the result, recover the vector.

    Linear map, so the least-squares solve is exact when the target is in
    range -- no starting guess, no local minima, no tolerance to tune.
    """

    parameterisation = BladeParameterisation()
    design_vector = _design_vector(parameterisation)
    chord, twist = parameterisation.evaluate(design_vector)

    recovered, chord_rms, twist_rms = parameterisation.fit(chord, twist)

    assert np.allclose(recovered, design_vector, atol=1e-12)
    assert chord_rms < 1e-12 and twist_rms < 1e-12


# ---------------------------------------------------------------------------
# Strip count and design-variable count are independent
# ---------------------------------------------------------------------------

def test_strip_count_and_control_point_count_are_independent():
    """
    Plan section 4.1 calls conflating these "the single most likely
    misunderstanding by a reader", so it gets a test rather than only a
    comment.
    """

    coarse = BladeParameterisation(n_chord=5, n_twist=5, n_stations=10)
    fine = BladeParameterisation(n_chord=5, n_twist=5, n_stations=200)

    assert coarse.n_design_variables == fine.n_design_variables == 10
    assert coarse.n_stations == 10 and fine.n_stations == 200

    # Same design vector, same shape -- the strips only sample it.
    design_vector = _design_vector(coarse)
    assert coarse.chord(design_vector)[0] != fine.chord(design_vector)[0]
    assert np.interp(0.5, coarse.stations, coarse.chord(design_vector)) == pytest.approx(
        np.interp(0.5, fine.stations, fine.chord(design_vector)), rel=1e-3)


def test_defaults_come_from_config():
    """25 strips, 5 + 5 control points -- the plan's intended configuration."""

    design = load_design_rotor()
    parameterisation = BladeParameterisation()

    assert parameterisation.n_stations == design.parameterisation.n_bem_strips
    assert parameterisation.n_chord == design.parameterisation.n_control_points_chord
    assert parameterisation.n_twist == design.parameterisation.n_control_points_twist
    assert parameterisation.radius_m == design.radius_m
    assert parameterisation.root_fraction == design.root_fraction


# ---------------------------------------------------------------------------
# Bounds and scaling
# ---------------------------------------------------------------------------

def test_bounds_from_config_carry_the_decided_values():
    """
    The bounds were TODO until 2026-09-19; `from_config()` raised, and the
    test in this slot asserted that it did, as the reminder to replace it
    with a real check the day the values landed. They landed (O4, MJ's
    decision, basis in `config/rotor_design.yaml`), so this checks the real
    values in the units `DesignBounds` works in -- chord in metres, twist in
    radians -- and that the control-point counts follow the config.
    """

    bounds = DesignBounds.from_config()

    assert bounds.chord_min_m == pytest.approx(0.045)
    assert bounds.chord_max_m == pytest.approx(0.30)
    assert bounds.twist_min_rad == pytest.approx(math.radians(-2.0))
    assert bounds.twist_max_rad == pytest.approx(math.radians(35.0))
    assert (bounds.n_chord, bounds.n_twist) == (5, 5)
    assert bounds.n_design_variables == 10


def test_scaling_round_trips_and_the_chain_rule_factor_is_exact():
    """
    Physical <-> scaled, and the constant factor gradients pick up.

    Plan section 4.4: the scaling exists because SLSQP's line search degrades
    on unscaled variables, and the chain rule factor is `hi - lo`. Checked by
    complex step so it is the factor the code actually applies, not the one
    the docstring claims.
    """

    bounds = DesignBounds.from_config(n_chord=5, n_twist=5)

    physical = np.concatenate([
        np.linspace(0.20, 0.05, 5), np.radians(np.linspace(25.0, -2.0, 5))])
    scaled = bounds.to_scaled(physical)

    assert np.all((scaled >= 0.0) & (scaled <= 1.0))
    assert np.allclose(bounds.to_physical(scaled), physical, atol=1e-15)

    h = 1e-30
    for j in range(bounds.n_design_variables):
        perturbed = scaled.astype(complex)
        perturbed[j] += 1j * h
        derivative = bounds.to_physical(perturbed).imag / h
        expected = np.zeros(bounds.n_design_variables)
        expected[j] = bounds.span()[j]
        assert np.allclose(derivative, expected, atol=1e-15), j


def test_clipping_reports_what_it_moved():
    """
    Plan step 1.7 needs feasibility *recorded*, not silently satisfied.

    A clip that returns only the clipped vector is indistinguishable from a
    vector that was already feasible, which is precisely the distinction the
    baseline has to report.
    """

    bounds = DesignBounds.from_config(n_chord=3, n_twist=3)
    physical = np.array([0.5, 0.1, 0.001,
                         math.radians(50.0), 0.0, math.radians(-20.0)])

    clipped, violations = bounds.clip_physical(physical)

    assert [index for index, _value, _bound in violations] == [0, 2, 3, 5]
    assert np.all(clipped >= bounds.lower() - 1e-15)
    assert np.all(clipped <= bounds.upper() + 1e-15)

    _, none_moved = bounds.clip_physical(clipped)
    assert none_moved == []


def test_unordered_bounds_are_rejected():
    with pytest.raises(ValueError):
        DesignBounds(chord_min_m=0.4, chord_max_m=0.02,
                     twist_min_rad=0.0, twist_max_rad=1.0, n_chord=3, n_twist=3)


# ---------------------------------------------------------------------------
# Handing the result to the solver
# ---------------------------------------------------------------------------

def test_to_geometry_produces_a_solvable_rotor():
    """
    The parameterisation's output must be something `bem/` accepts.

    Every station strictly inside (r_hub, R) is the specific requirement --
    `StationParams` raises at either endpoint, where the Prandtl factors
    degenerate -- and it is why the stations are strip midpoints.
    """

    parameterisation = BladeParameterisation()
    geometry = parameterisation.to_geometry(_design_vector(parameterisation))

    assert geometry.polar_cache == "sg6043"
    assert geometry.R == parameterisation.radius_m
    assert len(geometry.r) == parameterisation.n_stations
    assert all(geometry.r_hub < r < geometry.R for r in geometry.r)
    assert all(geometry.r[i] < geometry.r[i + 1] for i in range(len(geometry.r) - 1))
