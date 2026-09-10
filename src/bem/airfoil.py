"""
`LinearPolar`: a synthetic airfoil polar with exactly-known derivatives.

Real cached polars live in `polars.polar.CachedPolar` since work order Task 4.
This module used to hold `S809Polar` alongside, adapting the XFOIL cache
through `xfoil.polar_lookup` with alpha and Reynolds clamped to the table
bounds and a silent inward-stepping fallback over cache gaps. All three are
gone -- see `polars/polar.py`'s module docstring for what each one was hiding
and why a gap-free +/-180 deg cache (Task 2) plus a C1 interpolant (Task 3)
makes them unnecessary. Nothing under `bem/` imports `xfoil` any more, which
`tests/test_invariants.py` asserts.

What is left here is deliberately kept, not left behind. An analytic polar
whose Cl(alpha) and its derivative are known in closed form is the right
fixture for complex-step-verifying the BEM partials in Phase 3: it puts no
interpolation noise between the residual and the number being checked, so a
disagreement there is the solver's, not the polar's.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from dataclasses import dataclass


@dataclass
class LinearPolar:
    """
    Cl(alpha) = cl_slope * alpha, clamped to +/- cl_max beyond +/- stall_alpha.
    Cd(alpha) = cd0, constant (no stall drag rise modelled at this stage).

    Parameters
    ----------
    cl_slope : float
        Lift-curve slope, per radian (thin-airfoil default: 2*pi).
    stall_alpha : float
        Angle of attack, radians, at which Cl clamps to cl_max.
    cd0 : float
        Constant drag coefficient.
    """

    cl_slope: float = 6.2831853071795865  # 2*pi
    stall_alpha: float = 0.20943951023931953  # 12 deg, radians
    cd0: float = 0.01

    @property
    def cl_max(self):
        return self.cl_slope * self.stall_alpha

    def cl(self, alpha):
        """Lift coefficient at angle of attack alpha (radians)."""
        if alpha > self.stall_alpha:
            return self.cl_max
        if alpha < -self.stall_alpha:
            return -self.cl_max
        return self.cl_slope * alpha

    def cd(self, alpha):
        """Drag coefficient at angle of attack alpha (radians)."""
        return self.cd0

    def __call__(self, alpha):
        return self.cl(alpha), self.cd(alpha)
