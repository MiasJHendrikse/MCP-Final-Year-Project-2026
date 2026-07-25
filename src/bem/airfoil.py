"""
Synthetic airfoil polar for Stage 1 BEM development.

Linear Cl-alpha up to a stall angle (then clamped flat) and constant Cd.
This stands in for a real (e.g. S809) polar until the BEM solver is wired up
to polar_lookup.get_polar in a later session — Stage 1 only needs a smooth,
differentiable, closed-form Cl(alpha)/Cd(alpha) to validate the induction
solve itself.

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
