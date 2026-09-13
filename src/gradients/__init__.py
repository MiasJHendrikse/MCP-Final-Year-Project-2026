"""
The gradient layer: finite-difference gradients of the scaled objective, and
the scaled optimisation problem SLSQP is handed (Phase 2, Stage A).

Plan sections 9-10 (`docs/PROJECT_DIRECTION_v2.md`). The chain is

    u in [0,1]^10 -> d = lo + u * span -> chord/twist -> BEM -> J = -AEP

and this package owns the two things that sit *outside* that chain: how `J`
is differentiated numerically, and how the problem is presented to an
optimiser (scaling, the polar-cache envelope constraint, post-checks).

  `finite_difference`  central differences with a scalar or per-variable
                       step, and the step-size sweep the noise floor is read
                       off. Pure functions: no config, no rotor.
  `problem`            `ScaledProblem`: J(u)/|J(u0)|, its FD Jacobian, the
                       mandatory Reynolds-envelope linear constraint, and the
                       angle-of-attack post-check.

The discrete adjoint (Phase 3) is a separate package, `adjoint/`, not yet
created; its plan is `docs/adjoint_derivation.md`.

Bounds are provisional: `chord_max_m = 0.45 m` is a placeholder pending the
hub-radius / root-attachment decision, and every optimisation result produced
through this package is "under provisional bounds". Nothing here reads a bound
from `config/`; a `DesignBounds` instance is always passed in explicitly.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from gradients.finite_difference import central_difference, step_size_sweep
from gradients.problem import ScaledProblem

__all__ = ["ScaledProblem", "central_difference", "step_size_sweep"]
