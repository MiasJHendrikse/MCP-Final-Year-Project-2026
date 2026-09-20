"""
The material proxy (`objective/mass.py`): the section coefficients are the
airfoil's (including the thin-shell `k_I`, `k_Z`), the span rule is the
adjoint's, the gradient is exact, the mass in kilograms is `rho t A_shell`
with the recorded laminate and nothing else, and it goes back to `None` the
moment a structural input is re-opened.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import dataclasses
import json
import os

import numpy as np
import pytest

import config
from adjoint.system import _trapezoid_weights
from design import BladeParameterisation
from gradients import central_difference
from objective.mass import MODEL_POWER, MaterialModel, planform_weights, section_coefficients

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.join(_HERE, "..", "verification", "baseline", "x0.json")
XC_PATH = os.path.join(_HERE, "..", "verification", "load_constraint", "result_eps0.json")


@pytest.fixture(scope="module")
def parameterisation():
    return BladeParameterisation()


@pytest.fixture(scope="module")
def model(parameterisation):
    return MaterialModel(parameterisation)


@pytest.fixture(scope="module")
def x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])


@pytest.fixture(scope="module")
def x_c():
    with open(XC_PATH, encoding="utf-8") as handle:
        return np.array(json.load(handle)["x_c"])


# -- the section ---------------------------------------------------------------

def test_sg6043_section_coefficients_are_read_from_the_coordinates():
    """`k_A`, `k_P` from `data/airfoils/sg6043.dat` (81 points), pinned."""

    section = section_coefficients("sg6043")
    assert section.n_points == 81
    assert section.area == pytest.approx(0.06850241963050001, rel=1e-12)
    assert section.perimeter == pytest.approx(2.048504996756389, rel=1e-12)
    # A 10-12 % thick section: area well inside (0.05, 0.09) c^2, perimeter
    # just over 2 c.
    assert 0.05 < section.area < 0.09
    assert 2.0 < section.perimeter < 2.1


def test_sg6043_thin_shell_coefficients_are_pinned_and_consistent():
    """
    `k_I`, `k_Z` of the constant-thickness skin about the chord-parallel axis
    through the contour centroid (mass.py, "the thin-shell section
    coefficients"): pinned, `k_Z = k_I / extreme fibre`, the extreme fibre is
    the upper surface, and the principal-axis rotation neglected is < 0.1 deg.
    """

    section = section_coefficients("sg6043")
    assert section.second_moment == pytest.approx(0.003293118852977795, rel=1e-12)
    assert section.section_modulus == pytest.approx(0.051642701138200756, rel=1e-12)
    assert section.section_modulus == pytest.approx(section.second_moment / section.extreme_fibre,
                                                    rel=1e-14)
    assert section.shell_centroid_y == pytest.approx(0.038922637325357994, rel=1e-12)
    assert section.extreme_fibre == pytest.approx(0.10269 - section.shell_centroid_y, rel=1e-12)
    assert abs(section.principal_rotation_deg) < 0.1
    # Order of magnitude: a 10 % thick contour has y ~ +-0.05 c over 2.05 c of
    # arc, so k_I ~ 2.05 * 0.05^2 / 2 ~ 0.0025 per unit thickness.
    assert 0.002 < section.second_moment < 0.005
    assert 0.03 < section.section_modulus < 0.08


def test_thin_shell_second_moment_matches_a_direct_arc_length_quadrature():
    """The per-segment formula against a brute-force subdivision of the contour."""

    section = section_coefficients("sg6043")
    points = np.loadtxt(section.path, skiprows=1)
    x, y = points[:, 0], points[:, 1]
    x_next, y_next = np.roll(x, -1), np.roll(y, -1)
    n_sub = 400
    s = (np.arange(n_sub) + 0.5) / n_sub
    ys = (y[:, None] + (y_next - y)[:, None] * s[None, :]).ravel()
    ds = np.repeat(np.hypot(x_next - x, y_next - y) / n_sub, n_sub)
    y_bar = float(np.sum(ys * ds) / np.sum(ds))
    i_xx = float(np.sum(ds * (ys - y_bar) ** 2))
    assert y_bar == pytest.approx(section.shell_centroid_y, rel=1e-12)
    assert i_xx == pytest.approx(section.second_moment, rel=1e-6)


def test_the_configured_model_is_the_shell(model):
    assert model.model == "shell"
    assert model.power == 1
    assert MODEL_POWER == {"shell": 1, "solid": 2}


# -- the span rule ---------------------------------------------------------------

def test_planform_weights_are_the_adjoints_trapezoid_rule(parameterisation):
    """One integration rule for material, energy and load."""

    radii = parameterisation.radii
    assert np.array_equal(planform_weights(radii), np.asarray(_trapezoid_weights(list(radii))))


def test_planform_area_matches_numpy_trapezoid_at_x0(model, parameterisation, x0):
    chord = parameterisation.chord(x0)
    expected = float(np.trapezoid(chord, parameterisation.radii))
    assert model.planform_integral(x0, 1) == pytest.approx(expected, rel=1e-14)


# -- the values at the committed blades ---------------------------------------------

def test_reference_material_at_x0_is_pinned(model, x0):
    """The normalisers every relative figure divides by."""

    assert model.shell_area(x0) == pytest.approx(0.43705266291483974, rel=1e-12)
    assert model.solid_volume(x0) == pytest.approx(0.0022341394994333107, rel=1e-12)
    assert model.value(x0) == model.shell_area(x0)


def test_the_energy_optimum_carries_more_material_than_schmitz(model, x0, x_c):
    """The proposal's +14 % (solid proxy) and the shell figure beside it."""

    report = model.report(x_c, reference=x0)
    assert report["solid_over_reference"] == pytest.approx(1.14196, abs=5e-5)
    assert report["shell_over_reference"] == pytest.approx(1.06171, abs=5e-5)
    assert report["mass_kg"] == pytest.approx(1.7819, abs=5e-4)


# -- the gradient -----------------------------------------------------------------

@pytest.mark.parametrize("power", [1, 2])
def test_planform_gradient_matches_central_differences_to_round_off(model, x0, power):
    """Polynomial in `d`: central FD is exact up to round-off, and so is the gradient."""

    f = lambda d: float(model.planform_integral(d, power))  # noqa: E731
    fd, _n = central_difference(f, x0, 1e-4)
    analytic = model.planform_integral_gradient(x0, power)
    assert np.allclose(analytic, fd, rtol=1e-9, atol=1e-12)
    assert np.all(analytic[model.n_chord:] == 0.0)


def test_shell_gradient_is_constant_and_positive(model, x0, x_c):
    """Linear objective: the same gradient at every design, every chord entry positive."""

    g0, gc = model.gradient(x0), model.gradient(x_c)
    assert np.array_equal(g0, gc)
    assert np.all(g0[:model.n_chord] > 0.0)
    assert np.all(g0[model.n_chord:] == 0.0)


def test_solid_gradient_matches_fd(parameterisation, x0):
    solid = MaterialModel(parameterisation, model="solid")
    fd, _n = central_difference(lambda d: float(solid.value(d)), x0, 1e-4)
    assert np.allclose(solid.gradient(x0), fd, rtol=1e-9, atol=1e-12)


def test_planform_integral_is_complex_safe(model, x0):
    d = x0.astype(complex)
    d[0] += 1e-20j
    value = model.planform_integral(d, 1)
    assert value.imag / 1e-20 == pytest.approx(model.planform_integral_gradient(x0, 1)[0], rel=1e-12)


# -- kilograms: the recorded laminate, and nothing without it -----------------------------

def test_mass_in_kg_is_rho_t_times_the_shell_area_with_the_recorded_laminate(model, x0, x_c):
    """
    `1920 kg/m^3 * 2 mm * k_P int c dr` per blade (docs/MATERIALS-STRUCTURAL-INPUTS.md):
    the reference blade's skin is 1.678 kg, the energy optimum's 1.782 kg;
    the solid figure beside it is the laminate's density on the section
    volume. No other factor -- no web, no root, no gelcoat.
    """

    design = config.load_design_rotor()
    rho, t = design.laminate_density_kg_m3, design.shell_thickness_m
    assert model.mass_kg(x0) == pytest.approx(rho * t * model.shell_area(x0), rel=1e-14)
    assert model.mass_kg(x0) == pytest.approx(1.6782822255929846, rel=1e-12)
    assert model.mass_kg(x_c) == pytest.approx(1.7819, abs=5e-4)
    assert model.solid_mass_kg(x0) == pytest.approx(rho * model.solid_volume(x0), rel=1e-14)
    assert model.solid_mass_kg(x0) == pytest.approx(4.2895, abs=5e-4)
    report = model.report(x0)
    assert report["mass_kg"] == report["shell_mass_kg"] == model.mass_kg(x0)
    assert report["solid_mass_kg"] == model.solid_mass_kg(x0)
    assert report["section_coefficients"]["k_I"] == section_coefficients().second_moment
    assert report["section_coefficients"]["k_Z"] == section_coefficients().section_modulus


def test_structural_factors_carry_the_recorded_inputs(model):
    """`E k_I t` and `k_Z t`: the constants the absolute rows carry instead of cancelling."""

    design = config.load_design_rotor()
    section = section_coefficients()
    assert model.stiffness_factor_pa_m() == pytest.approx(
        design.youngs_modulus_pa * section.second_moment * design.shell_thickness_m, rel=1e-14)
    assert model.section_modulus_factor_m() == pytest.approx(
        section.section_modulus * design.shell_thickness_m, rel=1e-14)
    inputs = model.structural_inputs()
    assert inputs["laminate_density_kg_m3"] == 1920.0
    assert inputs["tip_clearance_m"].startswith("TODO:")


def test_mass_in_kg_is_none_again_if_a_structural_input_is_reopened(parameterisation, x0):
    """Re-opening a field to TODO takes the number away; nothing is substituted."""

    todo = config.Unresolved("rotor_design.yaml:structure.laminate_density_kg_m3", "re-opened")
    design = dataclasses.replace(config.load_design_rotor(), laminate_density_kg_m3=todo)
    shell = MaterialModel(parameterisation, design=design)
    assert shell.mass_kg(x0) is None
    assert shell.shell_mass_kg(x0) is None
    assert shell.solid_mass_kg(x0) is None
    assert shell.report(x0)["mass_kg"] is None
    assert shell.stiffness_factor_pa_m() is not None      # needs E and t, not rho

    todo_t = config.Unresolved("rotor_design.yaml:structure.shell_thickness_m", "re-opened")
    design = dataclasses.replace(config.load_design_rotor(), shell_thickness_m=todo_t)
    shell = MaterialModel(parameterisation, design=design)
    assert shell.mass_kg(x0) is None
    assert shell.solid_mass_kg(x0) is not None            # the solid needs only rho
    assert shell.stiffness_factor_pa_m() is None
    assert shell.section_modulus_factor_m() is None
    assert shell.structural_inputs()["shell_thickness_m"] == "TODO: re-opened"


def test_mass_in_kg_when_a_laminate_is_supplied(parameterisation, x0):
    """`rho t A_shell` for the shell, `rho V_solid` for the solid -- no other factor."""

    design = dataclasses.replace(config.load_design_rotor(),
                                 laminate_density_kg_m3=1800.0, shell_thickness_m=0.003)
    shell = MaterialModel(parameterisation, design=design)
    solid = MaterialModel(parameterisation, design=design, model="solid")
    assert shell.mass_kg(x0) == pytest.approx(1800.0 * 0.003 * shell.shell_area(x0), rel=1e-14)
    assert solid.mass_kg(x0) == pytest.approx(1800.0 * solid.solid_volume(x0), rel=1e-14)


def test_unknown_model_is_rejected(parameterisation):
    with pytest.raises(ValueError, match="unknown mass model"):
        MaterialModel(parameterisation, model="hollow")
