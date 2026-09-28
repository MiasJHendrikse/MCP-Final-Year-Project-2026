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

  `schmitz`           the analytic Schmitz optimum blade.
  `baseline`          x0: Schmitz fitted onto the parameterisation, its
                      feasibility status, and the reference numbers the whole
                      results chapter is measured against.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

from design.baseline import (
    BaselineBlade,
    build_schmitz_baseline,
    evaluate_baseline,
    root_bending_moment,
)
from design.bounds import DesignBounds
from design.parameterisation import (
    DEFAULT_DEGREE,
    BladeParameterisation,
    basis_matrix,
    clamped_knots,
)

__all__ = [
    "BaselineBlade",
    "DEFAULT_DEGREE",
    "BladeParameterisation",
    "DesignBounds",
    "basis_matrix",
    "build_schmitz_baseline",
    "clamped_knots",
    "evaluate_baseline",
    "root_bending_moment",
]
