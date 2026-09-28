"""
Complex-step safety of the residual, with the real path
bit-identical.

The discrete adjoint's Tier 1 check (`test_adjoint_partials.py`) verifies
every hand-derived partial against a complex-step derivative of the *code's
own* residual, `bem.station.residual`, with `h = 1e-30`. That only works if
a complex perturbation to `phi`, to the twist, or to the chord (which
arrives at the polar as a complex Reynolds number) passes through the whole
residual with its imaginary part intact. Three edits made that true:

  1. `station._blade_element`: `sin`/`cos` dispatch on type (`_sin`, `_cos`).
  2. `corrections.tip_loss_factor` / `hub_loss_factor`: a complex branch with
     every decision taken on `.real`; the real-path lines untouched.
  3. `polars.polar.CachedPolar.__init__`: range-check on `.real`, store the
     Reynolds number as given rather than through `float()`.

The exit criterion for "real path unchanged" is the rest of the suite --
`test_golden_regression.py`, `test_determinism.py`, `test_cost.py` -- which
this file does not duplicate. What it checks is the other half: that the
complex path exists, returns complex, and carries a derivative that agrees
with a central difference to the precision a central difference has.

Bounds do not enter here. `x0` is the committed baseline design.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import math
import os

import numpy as np
import pytest

from bem.corrections import BUHL_AC, combined_loss_factor, hub_loss_factor, tip_loss_factor
from bem.station import StationParams, residual, solve_station
from config import load_design_rotor, load_site
from design import BladeParameterisation
from polars.interpolant import PolarDomainError
from polars.polar import CachedPolar, interpolant_for

_HERE = os.path.dirname(os.path.abspath(__file__))
X0_PATH = os.path.abspath(os.path.join(_HERE, "..", "verification", "baseline", "x0.json"))

#: Complex-step size. Anything below ~1e-150 risks underflow in squares; at
#: 1e-30 the imaginary part is exact to round-off and no subtraction occurs.
H = 1e-30


@pytest.fixture(scope="module")
def stations():
    """Every station of `x0` at the rated point, solved, with its polar."""

    parameterisation = BladeParameterisation()
    design = load_design_rotor()
    site = load_site()
    interpolant = interpolant_for("sg6043")

    with open(X0_PATH, encoding="utf-8") as handle:
        artefact = json.load(handle)
    d = np.array(artefact["chord_control_points_m"] + artefact["twist_control_points_rad"])
    chord, twist = parameterisation.evaluate(d)

    v_inf = float(design.rated_wind_speed_ms)
    omega = design.design_tsr * v_inf / parameterisation.radius_m
    r_hub = parameterisation.root_fraction * parameterisation.radius_m

    out = []
    for r, c, theta in zip(parameterisation.radii, chord, twist):
        reynolds = math.hypot(v_inf, omega * r) * c / site.kinematic_viscosity
        station = StationParams(
            r=float(r), chord=float(c), twist=float(theta),
            airfoil=CachedPolar(interpolant, reynolds), tsr=omega * r / v_inf,
            R=parameterisation.radius_m, n_blades=design.n_blades, r_hub=r_hub)
        solved = solve_station(station)
        assert solved["converged"]
        out.append((station, solved, reynolds, interpolant))
    return out


# ---------------------------------------------------------------------------
# the residual accepts a complex phi and carries the derivative
# ---------------------------------------------------------------------------

def test_residual_returns_complex_for_complex_phi(stations):
    for station, solved, _re, _ip in stations:
        value = residual(complex(solved["phi"], H), station)
        assert isinstance(value, complex)
        assert math.isfinite(value.imag) and value.imag != 0.0


def test_complex_step_dR_dphi_agrees_with_a_central_difference(stations):
    """
    Complex step at 1e-30 against central FD at 1e-6: they should agree to
    the FD's own truncation, ~1e-9 relative, at every station. This is not
    the Tier 1 assertion (that is against the hand-derived partial, at
    1e-13); it establishes that the imaginary part is a derivative at all.
    """

    for station, solved, _re, _ip in stations:
        phi = solved["phi"]
        cs = residual(complex(phi, H), station).imag / H
        h = 1e-6
        fd = (residual(phi + h, station) - residual(phi - h, station)) / (2.0 * h)
        assert abs(cs - fd) <= 1e-7 * max(1.0, abs(fd)), (station.r, cs, fd)


def test_the_complex_path_real_part_matches_the_real_path(stations):
    """
    The complex path is a different set of `cmath` calls, so its real part
    may differ from the real path in the last bit -- but not more. (The
    real path itself is unchanged: that is what the golden tests assert.)
    """

    for station, solved, _re, _ip in stations:
        phi = solved["phi"]
        assert abs(residual(complex(phi, 0.0), station).real - residual(phi, station)) < 1e-13


def test_rated_point_includes_buhl_branch_stations(stations):
    """The test set must exercise the turbulent-wake branch, not just momentum."""

    assert sum(solved["a"] > BUHL_AC for _s, solved, _re, _ip in stations) > 0


# ---------------------------------------------------------------------------
# twist and chord (via the Reynolds number) perturbations propagate
# ---------------------------------------------------------------------------

def test_complex_twist_propagates(stations):
    for station, solved, reynolds, interpolant in stations:
        perturbed = StationParams(
            r=station.r, chord=station.chord, twist=complex(station.twist, H),
            airfoil=station.airfoil, tsr=station.tsr, R=station.R,
            n_blades=station.n_blades, r_hub=station.r_hub)
        value = residual(solved["phi"], perturbed)
        assert isinstance(value, complex) and value.imag != 0.0


def test_complex_chord_and_reynolds_propagate(stations):
    """
    A chord perturbation reaches the residual twice: through the solidity
    and through the Reynolds number the polar is read at. Both must carry.
    The polar-only path is checked on its own so a zero there (the `Re(c)`
    path silently dropped) cannot hide behind the solidity term.
    """

    for station, solved, reynolds, interpolant in stations:
        d_re = complex(reynolds, H * reynolds / station.chord)
        both = StationParams(
            r=station.r, chord=complex(station.chord, H), twist=station.twist,
            airfoil=CachedPolar(interpolant, d_re), tsr=station.tsr, R=station.R,
            n_blades=station.n_blades, r_hub=station.r_hub)
        polar_only = StationParams(
            r=station.r, chord=station.chord, twist=station.twist,
            airfoil=CachedPolar(interpolant, d_re), tsr=station.tsr, R=station.R,
            n_blades=station.n_blades, r_hub=station.r_hub)
        value = residual(solved["phi"], both)
        assert isinstance(value, complex) and value.imag != 0.0
        assert residual(solved["phi"], polar_only).imag != 0.0


# ---------------------------------------------------------------------------
# CachedPolar
# ---------------------------------------------------------------------------

def test_cached_polar_accepts_a_complex_reynolds_and_keeps_it():
    interpolant = interpolant_for("sg6043")
    reynolds = complex(2.0e5, H)
    polar = CachedPolar(interpolant, reynolds)
    assert polar.reynolds == reynolds
    assert isinstance(polar.cl(0.05), complex)
    assert polar.cl(0.05).imag / H == pytest.approx(polar.dcl_dre(0.05), rel=1e-12)
    assert "CachedPolar" in repr(polar)


def test_cached_polar_still_range_checks_on_the_real_part():
    interpolant = interpolant_for("sg6043")
    with pytest.raises(PolarDomainError):
        CachedPolar(interpolant, complex(2.0e6, H))
    with pytest.raises(PolarDomainError):
        CachedPolar(interpolant, 2.0e6)


def test_cached_polar_real_reynolds_is_stored_as_a_float():
    interpolant = interpolant_for("sg6043")
    polar = CachedPolar(interpolant, 200_000)
    assert isinstance(polar.reynolds, float)
    assert repr(polar) == "CachedPolar(reynolds=200,000)"


# ---------------------------------------------------------------------------
# loss factors: complex branch mirrors the real one, decisions on .real
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("phi", [0.02, 0.1, 0.4, 1.0, math.pi / 2])
def test_loss_factor_complex_branch_matches_real_branch(phi):
    real = combined_loss_factor(1.9, 2.0, 3, phi, 0.3)
    cplx = combined_loss_factor(1.9, 2.0, 3, complex(phi, 0.0), 0.3)
    assert isinstance(cplx, complex)
    assert cplx.real == pytest.approx(real, rel=1e-14, abs=1e-15)
    assert cplx.imag == 0.0


def test_loss_factor_complex_step_matches_the_analytic_derivative():
    """dF_tip/dphi = -(2/pi) f e^{-f} cos(phi) / (sin(phi) sqrt(1 - e^{-2f}))."""

    r, R, B, phi = 1.7, 2.0, 3, 0.35
    f = (B / 2.0) * (R - r) / (r * math.sin(phi))
    analytic = -(2.0 / math.pi) * f * math.exp(-f) * math.cos(phi) / (
        math.sin(phi) * math.sqrt(1.0 - math.exp(-2.0 * f)))
    cs = tip_loss_factor(r, R, B, complex(phi, H)).imag / H
    assert cs == pytest.approx(analytic, rel=1e-13)

    f_h = (B / 2.0) * (r - 0.3) / (0.3 * math.sin(phi))
    analytic_h = -(2.0 / math.pi) * f_h * math.exp(-f_h) * math.cos(phi) / (
        math.sin(phi) * math.sqrt(1.0 - math.exp(-2.0 * f_h)))
    cs_h = hub_loss_factor(r, 0.3, B, complex(phi, H)).imag / H
    assert cs_h == pytest.approx(analytic_h, rel=1e-13)


def test_loss_factor_guards_are_taken_on_the_real_part():
    assert tip_loss_factor(1.9, 2.0, 3, complex(1e-9, H)) == 1.0
    assert tip_loss_factor(2.0, 2.0, 3, complex(0.3, H)) == 0.0
    assert hub_loss_factor(0.3, 0.3, 3, complex(0.3, H)) == 0.0
    assert hub_loss_factor(1.0, None, 3, complex(0.3, H)) == 1.0
