"""
The control case: with `max_rotor_speed_rpm: null` the operating law is the
pre-2026-09-19 objective, bit for bit.

On 2026-09-19 the per-bin tip-speed ratio `lambda_b = min(6.5, Omega_max R /
V_b)` went into `objective.power` and `adjoint.system`. The change is meant
to be exactly nothing when there is no ceiling: same schedule, same solves
in the same order, same assembly. This module pins that against numbers
captured from the repository at `4c6feea` (the last commit before the law
landed) -- the AEP of the Schmitz baseline `x0`, the AEP of the optimum `x*`
of that objective, and the discrete-adjoint gradient at `x0` in physical
units -- with `==` and `np.array_equal`, not tolerances. If any of these
moves, the control case is no longer the old objective and the reason has
to be found, not the tolerance widened.

The numbers are literals here rather than read from an artefact because the
artefacts have since been regenerated under the ceiling.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import os

import numpy as np
import pytest

from adjoint.system import BEMSystem
from design import BladeParameterisation, DesignBounds
from objective import WeibullResource
from objective.objective import annual_energy_mwh
from objective.power import operating_points, power_per_bin
from test_config import _load_raw, _temporary_config

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.abspath(os.path.join(_HERE, "..", "verification", "baseline", "x0.json"))

#: Captured at `4c6feea` (2026-09-19, before the operating law), from
#: `objective.annual_energy_mwh` and `BEMSystem.gradient` on the committed
#: `x0` and the committed `x*` of `verification/adjoint_optimisation/result.json`
#: at that commit.
REFERENCE_COMMIT = "4c6feea"
AEP_X0_MWH = 10.2701574120755
AEP_X_STAR_MWH = 10.28256415375132
X_STAR_PRE_CEILING = [
    0.285331831248933, 0.17343287278160208, 0.10734485804097259,
    0.07790724468032187, 0.04582107770409465,
    0.41018455558122086, 0.18086208658505795, 0.06139525981804147,
    0.03501802390844747, -0.03490658503988659,
]
#: `BEMSystem.gradient(x0).J` -- assembled from the station integrand `q`,
#: a second code path that agrees with `annual_energy_mwh` to 2e-15.
J_X0_ADJOINT = -10.270157412075502
DJ_DD_X0 = [
    -0.002675751593509079, -0.00032978586365173257, 0.1830579013228233,
    1.3775977735205505, 0.9542324246417548,
    0.006639062788918015, -0.022729859823558574, -0.08447789592385968,
    -0.38968449063159305, -0.14976245456271164,
]


@pytest.fixture(scope="module")
def x0():
    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    return np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])


@pytest.fixture
def no_ceiling():
    """The configured rotor with `max_rotor_speed_rpm: null`."""

    rotor_yaml = _load_raw("rotor_design.yaml")
    rotor_yaml["operating"]["max_rotor_speed_rpm"] = None
    with _temporary_config({"rotor_design.yaml": rotor_yaml}):
        yield


def test_the_schedule_is_flat_without_a_ceiling(no_ceiling):
    points = operating_points()
    assert len(points) == 17
    assert all(lam == 6.5 for _v, lam in points)


def test_aep_at_x0_is_the_pre_ceiling_number_exactly(no_ceiling, x0):
    assert annual_energy_mwh(x0) == AEP_X0_MWH


def test_aep_at_the_pre_ceiling_optimum_is_exact(no_ceiling):
    assert annual_energy_mwh(np.array(X_STAR_PRE_CEILING)) == AEP_X_STAR_MWH


def test_the_adjoint_gradient_at_x0_is_the_pre_ceiling_gradient_exactly(no_ceiling, x0):
    parameterisation = BladeParameterisation()
    bounds = DesignBounds.from_config(n_chord=parameterisation.n_chord,
                                      n_twist=parameterisation.n_twist)
    system = BEMSystem(parameterisation, bounds, WeibullResource.from_config())
    assert all(lam == 6.5 for lam in system.tsr_per_bin)
    result = system.gradient(x0)
    assert result.J == J_X0_ADJOINT
    assert np.array_equal(result.dJ_dd, np.array(DJ_DD_X0))


def test_a_fixed_tsr_override_equals_no_ceiling(no_ceiling, x0):
    """`power_per_bin(..., tip_speed_ratio=6.5)` is the same solve as the control."""

    geometry = BladeParameterisation().to_geometry(x0)
    scheduled = power_per_bin(geometry)
    fixed = power_per_bin(geometry, tip_speed_ratio=6.5)
    assert np.array_equal(scheduled["power_w"], fixed["power_w"])
    assert np.array_equal(scheduled["tsr"], fixed["tsr"])


def test_the_ceiling_moves_the_numbers(x0):
    """
    And with the configured 300 rpm ceiling they move: the sanity check that
    the control fixture is testing something. Only the 10.5 m/s bin is
    uncapped and below lambda = 6.5, so AEP(x0) drops by ~0.2 %.
    """

    points = operating_points()
    assert [lam for _v, lam in points][:7] == [6.5] * 7
    assert points[7][1] == pytest.approx(300.0 * 2 * np.pi / 60 * 2.0 / 10.5)
    aep = annual_energy_mwh(x0)
    assert aep != AEP_X0_MWH
    assert 0.001 < (AEP_X0_MWH - aep) / AEP_X0_MWH < 0.005
