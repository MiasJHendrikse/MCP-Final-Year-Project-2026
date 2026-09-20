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
geometry, not assumptions. What turns either into a mass is a laminate --
`rho_lam * t_shell` for the shell, `rho` for the solid -- read from
`config/rotor_design.yaml::structure` (resolved 2026-09-20 evening:
1920 kg/m^3 and 2 mm, provenance in `docs/MATERIALS-STRUCTURAL-INPUTS.md`).
`mass_kg` returns a mass per blade in kg while they are resolved and `None`
if either is re-opened to `TODO` -- never a substitute number. The
production headline is still the percentage of the reference blade's
material, which needs neither.

The objective is the shell (`objective.mass_model: shell`): a constant-
thickness hand-laid skin is how a small blade of this size is built, and the
same shell assumption is what the root-stress (`Z ~ c^2 t`) and deflection
(`I ~ c^3 t`) proxies rest on, so the three are one model. The solid proxy
is reported beside it and never optimised.

The thin-shell section coefficients (2026-09-20, evening)
---------------------------------------------------------
The same idealisation -- one skin of constant thickness `t` following the
section contour, and nothing else -- has two more exact coefficients, so
that the relative proxies' `c^2 t` and `c^3 t` carry a constant and become
absolute:

    second moment    I = k_I c^3 t,   k_I = oint (y - y_bar)^2 ds   [per unit chord]
    section modulus  Z = k_Z c^2 t,   k_Z = k_I / max |y - y_bar|

with `s` the arc length round the unit-chord contour and `y_bar` the
centroid of the contour (perimeter-weighted, not the area centroid: a thin
shell's material sits on the line). The integral is exact per straight
segment (`L (y_m^2 + dy^2 / 12)`), the axis is chord-parallel through that
centroid -- flapwise bending -- and the extreme fibre is the upper surface
(0.0638 c above the centroid for the SG6043 against 0.0535 c below). The
principal-axis rotation of the cambered contour is 0.086 deg and is
neglected; the section's twist relative to the rotor plane is neglected
exactly as the relative proxies neglect it. SG6043: `k_I = 0.003293`,
`k_Z = 0.05164`. What these encode is a blade whose *skin* is the whole
structure: no spar cap, no web contribution, no root insert. The absolute
stress `KS M_ref / Z` and deflection `delta_unit / (E k_I t)` are therefore
the thin-shell blade's, and the module says so wherever it quotes them.

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
    closed coordinate loop; the thin-shell `k_I = I / (c^3 t)` and
    `k_Z = Z / (c^2 t)` about the chord-parallel axis through the contour's
    centroid (module docstring); and where they came from.
    """

    area: float
    perimeter: float
    second_moment: float
    section_modulus: float
    shell_centroid_y: float
    extreme_fibre: float
    principal_rotation_deg: float
    path: str
    n_points: int


@lru_cache(maxsize=None)
def section_coefficients(polar_cache="sg6043"):
    """
    Section area, perimeter, thin-shell second moment and section modulus
    per unit chord from the coordinate file.

    Shoelace over the closed loop for the area, polyline length for the
    perimeter (the loop closes on itself at the trailing edge; the closing
    segment is included in both). The shell moments treat each segment as a
    straight line of material: `I = sum L (y_m^2 + dy^2 / 12)` about the
    perimeter-weighted centroid, exact for the polyline. The file's first
    line is the name.
    """

    path = os.path.join(AIRFOIL_DIR, f"{polar_cache}.dat")
    points = np.loadtxt(path, skiprows=1)
    x, y = points[:, 0], points[:, 1]
    x_next, y_next = np.roll(x, -1), np.roll(y, -1)
    area = 0.5 * abs(float(np.sum(x * y_next - x_next * y)))
    lengths = np.hypot(x_next - x, y_next - y)
    perimeter = float(np.sum(lengths))

    # The thin shell: material on the contour, weighted by arc length.
    x_mid, y_mid = 0.5 * (x + x_next), 0.5 * (y + y_next)
    dx, dy = x_next - x, y_next - y
    x_bar = float(np.sum(lengths * x_mid) / perimeter)
    y_bar = float(np.sum(lengths * y_mid) / perimeter)
    i_xx = float(np.sum(lengths * ((y_mid - y_bar) ** 2 + dy ** 2 / 12.0)))
    i_yy = float(np.sum(lengths * ((x_mid - x_bar) ** 2 + dx ** 2 / 12.0)))
    i_xy = float(np.sum(lengths * ((x_mid - x_bar) * (y_mid - y_bar) + dx * dy / 12.0)))
    extreme = float(np.max(np.abs(y - y_bar)))
    rotation = 0.5 * np.degrees(np.arctan2(-2.0 * i_xy, i_yy - i_xx))
    return SectionCoefficients(area=area, perimeter=perimeter,
                               second_moment=i_xx, section_modulus=i_xx / extreme,
                               shell_centroid_y=y_bar, extreme_fibre=extreme,
                               principal_rotation_deg=float(rotation),
                               path=path, n_points=int(len(x)))


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

    def shell_mass_kg(self, d):
        """`rho_lam * t_shell * A_shell` per blade, kg, or `None` while either input is TODO."""

        rho, t = self.design.laminate_density_kg_m3, self.design.shell_thickness_m
        if not (is_resolved(rho) and is_resolved(t)):
            return None
        return float(rho) * float(t) * float(self.shell_area(d))

    def solid_mass_kg(self, d):
        """
        `rho_lam * V_solid` per blade, kg -- a solid section of the laminate's
        density, labelled as such by the caller -- or `None` while the
        density is TODO.
        """

        rho = self.design.laminate_density_kg_m3
        if not is_resolved(rho):
            return None
        return float(rho) * float(self.solid_volume(d))

    def mass_kg(self, d):
        """
        The configured proxy as a mass per blade, kg (`shell_mass_kg` or
        `solid_mass_kg`), or `None` while a structural input it needs is
        TODO. Never a substitute number.
        """

        return self.shell_mass_kg(d) if self.model == "shell" else self.solid_mass_kg(d)

    # -- the absolute structural factors -----------------------------------------

    def stiffness_factor_pa_m(self):
        """
        `E k_I t_shell` [Pa m]: what divides the unit-stiffness deflection
        `int M (R - r) / c^3 dr` to make metres. `None` while an input is TODO.
        """

        e, t = self.design.youngs_modulus_pa, self.design.shell_thickness_m
        if not (is_resolved(e) and is_resolved(t)):
            return None
        return float(e) * self.coefficients.second_moment * float(t)

    def section_modulus_factor_m(self):
        """`k_Z t_shell` [m]: `Z = (k_Z t) c^2`. `None` while the thickness is TODO."""

        t = self.design.shell_thickness_m
        if not is_resolved(t):
            return None
        return self.coefficients.section_modulus * float(t)

    def structural_inputs(self):
        """
        The structural inputs as recorded, for the artefacts: each field's
        value or the `TODO` note, plus the section coefficients. Reports, never
        substitutes.
        """

        out = {}
        for name in ("laminate_density_kg_m3", "shell_thickness_m", "youngs_modulus_pa",
                     "allowable_stress_pa", "safety_factor", "tip_clearance_m"):
            value = getattr(self.design, name)
            out[name] = float(value) if is_resolved(value) else f"TODO: {value.note}"
        out["section_coefficients"] = self.section_coefficient_record()
        return out

    def section_coefficient_record(self):
        c = self.coefficients
        return {"k_A": c.area, "k_P": c.perimeter, "k_I": c.second_moment,
                "k_Z": c.section_modulus, "shell_centroid_y_over_c": c.shell_centroid_y,
                "extreme_fibre_over_c": c.extreme_fibre,
                "principal_rotation_deg": c.principal_rotation_deg}

    def report(self, d, reference=None):
        """
        Both proxies at `d`, both masses in kg per blade (or `None`), and as a
        fraction of `reference` if given: what every artefact records for a
        blade.
        """

        d = np.asarray(d, dtype=float)
        out = {
            "model": self.model,
            "shell_area_m2": float(self.shell_area(d)),
            "solid_volume_m3": float(self.solid_volume(d)),
            "planform_area_m2": float(self.planform_integral(d, 1)),
            "mass_kg": self.mass_kg(d),
            "shell_mass_kg": self.shell_mass_kg(d),
            "solid_mass_kg": self.solid_mass_kg(d),
            "section_coefficients": self.section_coefficient_record(),
        }
        if reference is not None:
            reference = np.asarray(reference, dtype=float)
            out["shell_over_reference"] = float(self.shell_area(d) / self.shell_area(reference))
            out["solid_over_reference"] = float(self.solid_volume(d) / self.solid_volume(reference))
        return out
