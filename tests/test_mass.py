"""
The material proxy (`objective/mass.py`): the section coefficients are the
airfoil's, the span rule is the adjoint's, the gradient is exact, and no mass
in kilograms exists until the structural inputs do.

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
    assert report["mass_kg"] is None


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


# -- kilograms only when the inputs exist -----------------------------------------------

def test_mass_in_kg_is_none_while_the_structure_is_todo(model, x0):
    design = config.load_design_rotor()
    assert not config.is_resolved(design.laminate_density_kg_m3)
    assert not config.is_resolved(design.shell_thickness_m)
    assert model.mass_kg(x0) is None


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
