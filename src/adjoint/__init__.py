"""
The discrete adjoint of the Ning-form BEM objective (Phase 3).

    state  phi_{b,i}  (18 operating points x 25 stations)
    design d in R^10  (chord and twist control points)

  `kernels`  one station: the residual `R`, the power integrand `q`, and
             every first partial in `(phi, c, theta)` -- the derivation of
             `docs/adjoint_derivation.md` §2-§6 transcribed, complex-safe.
  `system`   `BEMSystem`: the forward solve (through the objective's own
             solver, once), `R(phi; d)` through the code's own residual, the
             partial arrays, the matrix-free `dR/dd` and its transpose, the
             forward-mode tangent, and the adjoint gradient.

Verified in four tiers (`docs/adjoint_derivation.md` §9): Tier 1 every
partial against a complex step of the code's own residual and objective;
Tier 2 the transpose identity; Tier 3 the assembled gradient against the FD
noise floor; Tier 4 the residual disagreement attributed to the polar
interpolation.

Bounds are provisional (`chord_max_m = 0.45 m` is a placeholder); nothing
here reads a bound from `config/`.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from adjoint.deflection import DeflectionGradientResult, DeflectionSystem
from adjoint.kernels import StationPartials, station_partials
from adjoint.loads import MomentGradientResult, RootMomentSystem
from adjoint.system import BEMSystem, ForwardState, GradientResult

__all__ = ["BEMSystem", "DeflectionGradientResult", "DeflectionSystem", "ForwardState",
           "GradientResult", "MomentGradientResult", "RootMomentSystem", "StationPartials",
           "station_partials"]
