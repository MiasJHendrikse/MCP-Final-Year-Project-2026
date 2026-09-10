"""
The design layer: the map from a design vector to a blade.

Plan section 4. Everything here is about *shape* -- how a handful of numbers an
optimiser can move become a chord and twist distribution, what their bounds
are, and how the derivatives of that map are obtained. No aerodynamics: this
package hands `bem/` a `RotorGeometry` and knows nothing about what happens to
it afterwards.

  `parameterisation`  B-spline chord and twist over the normalised span, with
                      the constant analytic Jacobians dc/dd and dtheta/dd.
  `bounds`            design-variable bounds and the scaling to O(1). The
                      bounds themselves are still TODO in config -- this is
                      the mechanism, not the numbers.

Still to come: `schmitz.py` and `baseline.py` (plan step 1.7), which construct
the x0 the whole results chapter is measured against.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from design.bounds import DesignBounds
from design.parameterisation import (
    DEFAULT_DEGREE,
    BladeParameterisation,
    basis_matrix,
    clamped_knots,
)

__all__ = [
    "DEFAULT_DEGREE",
    "BladeParameterisation",
    "DesignBounds",
    "basis_matrix",
    "clamped_knots",
]
