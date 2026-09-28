"""
Design-variable bounds and the scaling to O(1).

The bounds live in `config/rotor_design.yaml` (`parameterisation.bounds`)
with their basis written beside them, and `DesignBounds.from_config()` reads
them. Resolved 2026-09-19: three of the four were grounded on 2026-09-13
(`chord_min_m` from the laminate minimum, the twist range from the Schmitz
baseline with margin) and `chord_max_m = 0.30 m` was decided on the local
solidity at the root cut-out (sigma = 0.48), the commercial c/R range and the
measured insensitivity of the result to it (docs/OUTSTANDING-INPUTS.md
section 2). Until then the file carried `TODO` and `from_config()` raised
`UnresolvedConfigError`; a provisional set (`chord_max_m = 0.45` as a
placeholder with no basis) lived in `tests/test_parameterisation.py` so
bounds-dependent work could proceed without the placeholder leaking into
config. That mechanism is kept -- a `TODO` in the YAML still raises here --
and the placeholder is retired.

Scaling, and why it is not cosmetic
------------------------------------
Chord is O(0.1) m and twist is O(0.1) rad, but their *sensitivities* differ by
far more, and SLSQP's line search degrades badly on
unscaled variables. Every design variable is therefore mapped to [0, 1] across
its own bound interval before it reaches an optimiser, and the chain rule is
applied to gradients accordingly:

    d_physical = lo + u * (hi - lo)          u in [0, 1]
    dJ/du      = dJ/d_physical * (hi - lo)

That factor is exact and constant, which matters for the adjoint: it composes into
the adjoint as a diagonal scaling rather than as anything needing its own
verification.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math
from dataclasses import dataclass

import numpy as np

from config import is_resolved, load_design_rotor
from config.unresolved import UnresolvedConfigError


@dataclass(frozen=True)
class DesignBounds:
    """
    Lower and upper bounds for the chord and twist control points.

    Frozen, and stored in the units the parameterisation works in: chord in
    metres, twist in radians. The config file records twist bounds in degrees
    because that is what a human writes; the conversion happens once, here, on
    the way in.

    Parameters
    ----------
    chord_min_m, chord_max_m : float
    twist_min_rad, twist_max_rad : float
    n_chord, n_twist : int
        Control-point counts, so `lower()` and `upper()` can build the full
        vectors without the caller restating them.
    """

    chord_min_m: float
    chord_max_m: float
    twist_min_rad: float
    twist_max_rad: float
    n_chord: int
    n_twist: int

    def __post_init__(self):
        if not self.chord_min_m < self.chord_max_m:
            raise ValueError(
                f"chord bounds are not ordered: {self.chord_min_m} !< {self.chord_max_m}")
        if not self.twist_min_rad < self.twist_max_rad:
            raise ValueError(
                f"twist bounds are not ordered: {self.twist_min_rad} !< {self.twist_max_rad}")

    @classmethod
    def from_config(cls, n_chord=None, n_twist=None):
        """
        Build from `config/rotor_design.yaml`.

        Raises
        ------
        config.unresolved.UnresolvedConfigError
            If any bound is `TODO` in the YAML (none is, since 2026-09-19).
            The error names the field, so a caller that reaches here gets
            told what is missing rather than a plausible number.
        """

        design = load_design_rotor()
        parameterisation = design.parameterisation

        fields = {
            "chord_min_m": parameterisation.chord_min_m,
            "chord_max_m": parameterisation.chord_max_m,
            "twist_min_deg": parameterisation.twist_min_deg,
            "twist_max_deg": parameterisation.twist_max_deg,
        }
        missing = [name for name, value in fields.items() if not is_resolved(value)]
        if missing:
            raise UnresolvedConfigError(
                "design-variable bounds are still TODO in "
                f"config/rotor_design.yaml: {', '.join(missing)}. Plan section "
                "7.1 requires them to come from a manufacturability study, not "
                "from a default. See docs/OUTSTANDING-INPUTS.md. Pass explicit "
                "values to DesignBounds(...) for a study that states its own "
                "provisional range."
            )

        return cls(
            chord_min_m=float(fields["chord_min_m"]),
            chord_max_m=float(fields["chord_max_m"]),
            twist_min_rad=math.radians(float(fields["twist_min_deg"])),
            twist_max_rad=math.radians(float(fields["twist_max_deg"])),
            n_chord=n_chord if n_chord is not None else parameterisation.n_control_points_chord,
            n_twist=n_twist if n_twist is not None else parameterisation.n_control_points_twist,
        )

    @property
    def n_design_variables(self):
        return self.n_chord + self.n_twist

    def lower(self):
        """Lower bound vector, chord block then twist block."""

        return np.concatenate([
            np.full(self.n_chord, self.chord_min_m),
            np.full(self.n_twist, self.twist_min_rad),
        ])

    def upper(self):
        """Upper bound vector, same ordering."""

        return np.concatenate([
            np.full(self.n_chord, self.chord_max_m),
            np.full(self.n_twist, self.twist_max_rad),
        ])

    def as_record(self):
        """The four bounds as a JSON-ready dict, twist in both rad and deg."""

        return {
            "chord_min_m": float(self.chord_min_m),
            "chord_max_m": float(self.chord_max_m),
            "twist_min_rad": float(self.twist_min_rad),
            "twist_max_rad": float(self.twist_max_rad),
            "twist_min_deg": math.degrees(self.twist_min_rad),
            "twist_max_deg": math.degrees(self.twist_max_rad),
        }

    def span(self):
        """`upper - lower`, the constant chain-rule factor for the scaling."""

        return self.upper() - self.lower()

    def to_scaled(self, physical):
        """Physical design vector -> scaled to [0, 1] across the bounds."""

        return (np.asarray(physical, dtype=float) - self.lower()) / self.span()

    def to_physical(self, scaled):
        """
        Scaled design vector -> physical units.

        Complex-safe: the arithmetic is a multiply and an add against real
        constants, so a complex-step perturbation to `scaled` propagates
        through unchanged (`np.asarray` without a dtype preserves complex).
        """

        return self.lower() + np.asarray(scaled) * self.span()

    def clip_physical(self, physical):
        """
        Clip a physical design vector into the bounds, and report what moved.

        Returns `(clipped, violations)` where `violations` is a list of
        `(index, value, bound)`. Plan step 1.7 requires the Schmitz baseline's
        feasibility to be *recorded* -- clipped and noted if violated -- rather
        than silently satisfied, so the report is part of the return value and
        not something the caller has to reconstruct.
        """

        physical = np.asarray(physical, dtype=float)
        lower, upper = self.lower(), self.upper()

        violations = []
        for index, (value, lo, hi) in enumerate(zip(physical, lower, upper)):
            if value < lo:
                violations.append((index, float(value), float(lo)))
            elif value > hi:
                violations.append((index, float(value), float(hi)))

        return np.clip(physical, lower, upper), violations
