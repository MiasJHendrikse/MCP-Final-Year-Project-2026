"""
Weibull parameter extrapolation between heights.

The site's wind resource was extracted at **50 m**. The rotor's hub height is
**20 m** (`config/site.yaml:hub_height_m`). The two Weibull parameters must
therefore be moved down 30 m, and the project log (2026-09-10) is explicit
about how:

    "not at GWA's default 10 m or 50 m, and if only those are available the
     height extrapolation must be done and recorded, not fudged."

This module is that extrapolation, as mechanism. The numbers it produced for
this site are recorded in `verification/wind_resource/` and written into
`config/site.yaml`; nothing downstream calls this at run time.

The method
----------
Justus & Mikhail (1976). Both Weibull parameters vary with height, and the
correlation gives each from the reference-level pair alone:

    A(z) = A_ref * (z / z_ref) ** n

    n = (0.37 - 0.0881 * ln A_ref) / (1 - 0.0881 * ln(z_ref / 10))

    k(z) = k_ref * (1 - 0.0881 * ln(z_ref / 10)) / (1 - 0.0881 * ln(z / 10))

Two properties are worth stating because they are the reason this method was
chosen over a bare power law:

1. **k is not held constant.** A power law applied to A alone implicitly
   asserts the distribution keeps its shape with height, which is false --
   shear narrows the distribution upward and widens it downward. Holding k
   fixed at a value read 30 m above hub height would bias the AEP integral in
   a direction that no later step could detect.

2. **It needs no surface roughness.** The exponent comes from A_ref, not from
   a z0 that would have to be assumed for this terrain. That matters here:
   no roughness survey exists for the site, so a log-law central case would
   have required inventing one, which I won't do. z0 appears in
   this project only in the *cross-check* band, where a range is honest.

The 10 m in both denominators is the correlation's own anchor height, not a
property of this site, and is not configurable.

PROVENANCE -- READ BEFORE CITING
--------------------------------
The formulae above are transcribed from general knowledge of the method, NOT
from a retrieved copy of Justus & Mikhail (1976). No primary copy was
obtainable; the paper is not in `docs/references/`. This is exactly the
position `src/bem/corrections.py` and `src/bem/station.py` are in for Buhl and
Ning respectively, and it is recorded the same way -- as an open provenance
item in `docs/OUTSTANDING-INPUTS.md`, not as a settled citation.

What *is* established meanwhile is internal consistency and physical
behaviour, as committed tests in `tests/test_height_extrapolation.py`: the
identity at z = z_ref, monotonicity in z, the sign of dk/dz, and -- the check
with the most power -- agreement with an independent log-law calculation over
the plausible roughness range for this terrain. A transcription error in
either constant would break that agreement. What those tests cannot catch is
a faithful transcription of the wrong correlation, which is what the
provenance item is for.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math
from dataclasses import dataclass

#: The correlation's anchor height, m. Part of the published method, not a
#: site property -- see the module docstring.
JUSTUS_ANCHOR_HEIGHT_M = 10.0

#: The two published constants of Justus & Mikhail (1976).
JUSTUS_C0 = 0.37
JUSTUS_C1 = 0.0881


def _justus_denominator(height_m):
    """`1 - 0.0881 * ln(z / 10)`, which appears in both relations."""

    return 1.0 - JUSTUS_C1 * math.log(height_m / JUSTUS_ANCHOR_HEIGHT_M)


def shear_exponent(scale_ms, reference_height_m):
    """
    The Justus & Mikhail power-law exponent `n` for `A(z)`.

    Depends on the reference *scale parameter* and the reference *height*, and
    on nothing else -- notably not on surface roughness.
    """

    if scale_ms <= 0.0:
        raise ValueError(f"Weibull scale must be positive, got {scale_ms}")
    if reference_height_m <= 0.0:
        raise ValueError(
            f"reference height must be positive, got {reference_height_m}")

    return ((JUSTUS_C0 - JUSTUS_C1 * math.log(scale_ms))
            / _justus_denominator(reference_height_m))


@dataclass(frozen=True)
class HeightExtrapolation:
    """
    A recorded 2-parameter Weibull extrapolation from one height to another.

    Carries the inputs alongside the outputs deliberately: this object is
    serialised into the verification artefact, and an extrapolated pair that
    does not say what it came from is not auditable.
    """

    reference_height_m: float
    reference_k: float
    reference_scale_ms: float
    target_height_m: float
    shear_exponent: float
    k: float
    scale_ms: float

    @property
    def mean_speed_ms(self):
        """V-bar = A * Gamma(1 + 1/k) at the target height, m/s."""

        return self.scale_ms * math.gamma(1.0 + 1.0 / self.k)

    @property
    def reference_mean_speed_ms(self):
        """V-bar at the reference height -- the cross-check on the input pair."""

        return self.reference_scale_ms * math.gamma(
            1.0 + 1.0 / self.reference_k)

    def as_dict(self):
        return {
            "method": "Justus & Mikhail (1976)",
            "reference_height_m": self.reference_height_m,
            "reference_k": self.reference_k,
            "reference_scale_ms": self.reference_scale_ms,
            "reference_mean_speed_ms": self.reference_mean_speed_ms,
            "target_height_m": self.target_height_m,
            "shear_exponent": self.shear_exponent,
            "k": self.k,
            "scale_ms": self.scale_ms,
            "mean_speed_ms": self.mean_speed_ms,
        }


def extrapolate_weibull(k, scale_ms, reference_height_m, target_height_m):
    """
    Move a Weibull pair from `reference_height_m` to `target_height_m`.

    Parameters
    ----------
    k : float
        Shape parameter at the reference height, dimensionless.
    scale_ms : float
        Scale parameter `A` (equivalently `c`) at the reference height, m/s.
    reference_height_m, target_height_m : float
        Heights above ground, m. Both must be positive. Extrapolation may run
        in either direction; for this project it runs *downward*, 50 -> 20 m.

    Returns
    -------
    HeightExtrapolation
    """

    if k <= 0.0:
        raise ValueError(f"Weibull k must be positive, got {k}")
    if target_height_m <= 0.0:
        raise ValueError(
            f"target height must be positive, got {target_height_m}")

    n = shear_exponent(scale_ms, reference_height_m)

    return HeightExtrapolation(
        reference_height_m=reference_height_m,
        reference_k=k,
        reference_scale_ms=scale_ms,
        target_height_m=target_height_m,
        shear_exponent=n,
        k=k * (_justus_denominator(reference_height_m)
               / _justus_denominator(target_height_m)),
        scale_ms=scale_ms * (target_height_m / reference_height_m) ** n,
    )


def log_law_scale_ms(scale_ms, reference_height_m, target_height_m,
                     roughness_length_m):
    """
    The same scale shift under the logarithmic wind profile, for cross-checking.

        A(z) / A(z_ref) = ln(z / z0) / ln(z_ref / z0)

    This is deliberately NOT the project's central case. It requires a surface
    roughness length, and no roughness survey exists for this site -- so it is
    evaluated across a *range* of plausible z0 to bound the Justus & Mikhail
    result, never to supply a single number of its own.
    """

    if roughness_length_m <= 0.0:
        raise ValueError(
            f"roughness length must be positive, got {roughness_length_m}")
    if min(reference_height_m, target_height_m) <= roughness_length_m:
        raise ValueError(
            "the log law is only defined well above the roughness length; got "
            f"z0 = {roughness_length_m} against heights "
            f"{reference_height_m} and {target_height_m}"
        )

    return scale_ms * (math.log(target_height_m / roughness_length_m)
                       / math.log(reference_height_m / roughness_length_m))
