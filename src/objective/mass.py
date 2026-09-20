"""
The material proxy: blade shell area and section volume from the planform,
and their exact gradients through the spline basis (the mass problem, Step 1;
`docs/PLAN-mass-objective-2026-09-20.md` section 1).

What is measured
-----------------
For a blade of one airfoil section scaled by the local chord `c(r)`, two
geometric quantities are exact functions of the planform alone:

    shell area     A_shell = k_P * int c(r)   dr     [m^2]   (skin of constant thickness)
    section volume V_solid = k_A * int c(r)^2 dr     [m^3]   (solid, or t/c-scaled laminate)

with `k_P` the section perimeter per unit chord and `k_A` the section area
per unit chord squared, both read once from the airfoil coordinates
(`data/airfoils/<cache>.dat`: SG6043 gives 2.0485 and 0.0685). They are
geometry, not assumptions. What turns either into a mass is a laminate
concept -- `rho_lam * t_shell` for the shell, `rho` for the solid -- and
those are `TODO` in `config/rotor_design.yaml::structure` until one exists.
`mass_kg` therefore returns `None` rather than a number until they are
resolved; the production headline is a percentage of the reference blade's
material and needs neither.

The objective is the shell (`objective.mass_model: shell`): a constant-
thickness hand-laid skin is how a small blade of this size is built, and the
same shell assumption is what the root-stress (`Z ~ c^2 t`) and deflection
(`I ~ c^3 t`) proxies rest on, so the three are one model. The solid proxy
is reported beside it and never optimised.

The integral and its gradient
------------------------------
The span integral uses the *same* trapezoid weights `t_i` over the 25
strip-midpoint stations as the adjoint's power and moment assembly
(`adjoint.system._trapezoid_weights`), so material, energy and load are
integrated on one rule. The half-strips beyond the first and last station
are outside the rule for all three; stated, not fixed.

`c = N_c d` is linear in the design vector, so

    d/dd int c^p dr = N_c^T (p * t (.) c^(p-1))

exactly, with `N_c = parameterisation.dchord_dd()` (constant) and zero on
the twist block. For `p = 1` the objective is linear in `d` and its gradient
is a constant vector; for `p = 2` it is a convex quadratic. No BEM state is
involved: `tests/test_mass.py` checks the gradient against central
differences to round-off.

Import direction: `adjoint` imports `objective`, so the trapezoid rule is
imported lazily inside `planform_weights` (the `objective.loads` precedent).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import os
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from config import REPO_ROOT, is_resolved, load_design_rotor

#: Where the airfoil coordinate files live; the file is `<polar_cache>.dat`.
AIRFOIL_DIR = os.path.join(REPO_ROOT, "data", "airfoils")

#: The exponent each model integrates the chord at.
MODEL_POWER = {"shell": 1, "solid": 2}


@dataclass(frozen=True)
class SectionCoefficients:
    """
    `k_A = A / c^2` and `k_P = P / c` of the unit-chord section, from the
    closed coordinate loop, plus where they came from.
    """

    area: float
    perimeter: float
    path: str
    n_points: int


@lru_cache(maxsize=None)
def section_coefficients(polar_cache="sg6043"):
    """
    Section area and perimeter per unit chord from the coordinate file.

    Shoelace over the closed loop for the area, polyline length for the
    perimeter (the loop closes on itself at the trailing edge; the closing
    segment is included in both). The file's first line is the name.
    """

    path = os.path.join(AIRFOIL_DIR, f"{polar_cache}.dat")
    points = np.loadtxt(path, skiprows=1)
    x, y = points[:, 0], points[:, 1]
    x_next, y_next = np.roll(x, -1), np.roll(y, -1)
    area = 0.5 * abs(float(np.sum(x * y_next - x_next * y)))
    perimeter = float(np.sum(np.hypot(x_next - x, y_next - y)))
    return SectionCoefficients(area=area, perimeter=perimeter, path=path,
                               n_points=int(len(x)))


def planform_weights(radii):
    """The adjoint's trapezoid weights `t_i` over the station radii."""

    # Lazy: `adjoint.system` imports `objective` at module level.
    from adjoint.system import _trapezoid_weights

    return np.asarray(_trapezoid_weights([float(r) for r in radii]), dtype=float)


class MaterialModel:
    """
    The material proxy and its gradient over one parameterisation.

    Parameters
    ----------
    parameterisation : design.BladeParameterisation
    polar_cache : str
        Which airfoil's section coefficients to read; `"sg6043"`, the one
        default, as everywhere else.
    model : str or None
        `"shell"` or `"solid"`; `None` reads `objective.mass_model` from the
        design config.
    design : config.DesignRotorConfig or None
        Where the model name and the structural inputs are read from;
        `None` loads the config.
    """

    def __init__(self, parameterisation, polar_cache="sg6043", model=None, design=None):
        self.parameterisation = parameterisation
        self.polar_cache = polar_cache
        self.design = load_design_rotor() if design is None else design
        self.model = str(self.design.mass_model if model is None else model)
        if self.model not in MODEL_POWER:
            raise ValueError(f"unknown mass model {self.model!r}; one of {tuple(MODEL_POWER)}")
        self.power = MODEL_POWER[self.model]
        self.coefficients = section_coefficients(polar_cache)
        self.radii = np.asarray(parameterisation.radii, dtype=float)
        self.weights = planform_weights(self.radii)
        self.N_c = parameterisation.dchord_dd()
        self.n_chord = parameterisation.n_chord

    # -- the integrals ---------------------------------------------------------

    def chord(self, d):
        """Station chords `N_c d`; dtype follows `d`."""

        return self.parameterisation.chord(np.asarray(d))

    def planform_integral(self, d, power):
        """`sum_i t_i c_i^p` -- `int c^p dr` on the trapezoid rule; complex-safe."""

        c = self.chord(d)
        return (self.weights * c ** power).sum()

    def planform_integral_gradient(self, d, power):
        """`N_c^T (p t c^(p-1))`, zero on the twist block; shape (n_design,)."""

        c = self.chord(d)
        return self.N_c.T @ (power * self.weights * c ** (power - 1))

    def shell_area(self, d):
        """`k_P int c dr` per blade, m^2."""

        return self.coefficients.perimeter * self.planform_integral(d, 1)

    def solid_volume(self, d):
        """`k_A int c^2 dr` per blade, m^3."""

        return self.coefficients.area * self.planform_integral(d, 2)

    # -- the objective -----------------------------------------------------------

    def value(self, d):
        """The configured proxy at `d`: shell area or solid volume."""

        return self.shell_area(d) if self.model == "shell" else self.solid_volume(d)

    def gradient(self, d):
        """`d value / dd`, exact."""

        k = self.coefficients.perimeter if self.model == "shell" else self.coefficients.area
        return k * self.planform_integral_gradient(d, self.power)

    # -- reporting ---------------------------------------------------------------

    def mass_kg(self, d):
        """
        The proxy as a mass, or `None` while the structural inputs are TODO.

        shell: `rho_lam * t_shell * A_shell`; solid: `rho_lam * V_solid`
        (a solid section of the laminate's density -- labelled as such by
        the caller). Never a substitute number.
        """

        rho = self.design.laminate_density_kg_m3
        if not is_resolved(rho):
            return None
        if self.model == "shell":
            t = self.design.shell_thickness_m
            if not is_resolved(t):
                return None
            return float(rho) * float(t) * float(self.shell_area(d))
        return float(rho) * float(self.solid_volume(d))

    def report(self, d, reference=None):
        """
        Both proxies at `d`, and as a fraction of `reference` if given: what
        every artefact records for a blade.
        """

        d = np.asarray(d, dtype=float)
        out = {
            "model": self.model,
            "shell_area_m2": float(self.shell_area(d)),
            "solid_volume_m3": float(self.solid_volume(d)),
            "planform_area_m2": float(self.planform_integral(d, 1)),
            "mass_kg": self.mass_kg(d),
            "section_coefficients": {"k_A": self.coefficients.area,
                                     "k_P": self.coefficients.perimeter},
        }
        if reference is not None:
            reference = np.asarray(reference, dtype=float)
            out["shell_over_reference"] = float(self.shell_area(d) / self.shell_area(reference))
            out["solid_over_reference"] = float(self.solid_volume(d) / self.solid_volume(reference))
        return out
