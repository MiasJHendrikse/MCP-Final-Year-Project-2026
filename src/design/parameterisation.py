"""
B-spline chord and twist over the normalised span, with analytic derivatives.

Why a reduced smooth parameterisation is a correctness requirement, not a
convenience
-------------------------------------------------------------------------------
Plan section 4.2: independent per-station chord and twist admit a **sawtooth
null space**. BEM strips are aerodynamically uncoupled and the objective is a
sum over strips, so adjacent stations can oscillate against one another -- one
strip up, the next down -- with almost no net change in the objective. An
optimiser has no reason to avoid that and will walk into it, producing serrated
geometry that is unmanufacturable and yet a legitimate local optimum of the
discrete problem. A spline removes the null space by construction, which is why
no separate smoothness or monotonicity constraint is needed (section 7.2).

The property that makes the derivatives exact
-----------------------------------------------
A B-spline curve is **linear in its control points**:

    c(s) = sum_j  N_j(s) * d_j

so

    dc(s)/dd_j = N_j(s)

exactly, and independently of `d`. The Jacobian is therefore a *constant
matrix* -- the basis functions evaluated once at the strip stations -- not
something to be differentiated per design vector. Plan section 4.5 asks for
`dc/dd` and `dtheta/dd` analytically; here they are not merely analytic but
constant, which is the strongest form that requirement can take. It also means
the complex-step verification in `tests/test_parameterisation.py` is checking
the basis evaluation rather than a derivative derivation, and that the adjoint's
chain rule through this layer is a single matrix multiply.

Chord and twist are kept as separate blocks of one design vector rather than
two vectors, because that is the shape an optimiser and an adjoint both want:
`d = [chord control points, twist control points]`, with the split recorded on
the object so nothing downstream has to guess it.

Knots and degree
-----------------
Cubic (degree 3), with a *clamped* uniform knot vector: the first and last
knots are repeated `degree + 1` times, so the curve interpolates its first and
last control points exactly. That matters for a blade -- the root and tip
values are the two a designer most wants to be able to set and bound directly,
and under an unclamped basis neither control point would be reached.

The minimum control-point count for degree 3 is 4. The representation study
(`verification/representation_study/`) sweeps 6/8/10/12/16 *total*
across chord and twist, so per-block counts of 3 are reachable; the degree is
lowered automatically in that case rather than failing, and `degree` records
what was actually used.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import numpy as np
from scipy.interpolate import BSpline

from config import load_design_rotor

#: Spline degree. Cubic gives C2 continuity in the shape, which is more than
#: the objective needs but costs nothing and keeps curvature continuous for
#: any structural post-processing.
DEFAULT_DEGREE = 3


def clamped_knots(n_control, degree):
    """
    Clamped uniform knot vector for `n_control` control points of `degree`.

    Length is `n_control + degree + 1`, with the end knots repeated
    `degree + 1` times so the curve starts at the first control point and ends
    at the last.
    """

    if n_control < degree + 1:
        raise ValueError(
            f"a degree-{degree} B-spline needs at least {degree + 1} control "
            f"points, got {n_control}")

    interior = n_control - degree - 1
    return np.concatenate([
        np.zeros(degree + 1),
        np.arange(1, interior + 1) / (interior + 1),
        np.ones(degree + 1),
    ])


def basis_matrix(stations, n_control, degree):
    """
    The (n_stations, n_control) matrix `N` with `N[i, j] = N_j(s_i)`.

    Built by evaluating each basis function once, via a unit control vector --
    scipy exposes the basis only through the curve, and a curve whose control
    points are the j-th unit vector *is* the j-th basis function. Done once at
    construction, never in a hot path.

    Every row sums to 1 (partition of unity); `test_parameterisation.py`
    asserts that, since it is the cheapest available check that the knot
    vector and degree are consistent.
    """

    stations = np.asarray(stations, dtype=float)
    knots = clamped_knots(n_control, degree)

    matrix = np.empty((len(stations), n_control))
    for j in range(n_control):
        coefficients = np.zeros(n_control)
        coefficients[j] = 1.0
        matrix[:, j] = BSpline(knots, coefficients, degree, extrapolate=False)(stations)

    return matrix


class BladeParameterisation:
    """
    Maps a design vector to chord and twist distributions over the span.

    Parameters
    ----------
    n_chord, n_twist : int
        Control points per block. Default from
        `config/rotor_design.yaml` (5 and 5, the 10-parameter
        configuration -- a starting point for the representation study, not
        its conclusion).
    n_stations : int
        BEM strips. Default from config (25). Deliberately independent of the
        control-point counts: conflating the two is the single most likely
        misunderstanding for a reader, so they are separate
        arguments read from separate config fields.
    degree : int
        Spline degree, lowered automatically if a block has too few control
        points for it.
    radius_m, root_fraction : float
        Span, from config.

    Attributes
    ----------
    stations : ndarray
        (n_stations,) normalised span coordinates in [0, 1].
    radii : ndarray
        (n_stations,) station radii in metres.
    chord_basis, twist_basis : ndarray
        The constant Jacobian blocks: `dchord/dd_chord` and `dtwist/dd_twist`.
    """

    def __init__(self, n_chord=None, n_twist=None, n_stations=None,
                 degree=DEFAULT_DEGREE, radius_m=None, root_fraction=None):
        design = load_design_rotor()
        parameterisation = design.parameterisation

        self.n_chord = int(n_chord if n_chord is not None
                           else parameterisation.n_control_points_chord)
        self.n_twist = int(n_twist if n_twist is not None
                           else parameterisation.n_control_points_twist)
        self.n_stations = int(n_stations if n_stations is not None
                              else parameterisation.n_bem_strips)
        self.radius_m = float(radius_m if radius_m is not None else design.radius_m)
        self.root_fraction = float(root_fraction if root_fraction is not None
                                   else design.root_fraction)

        # Degree is capped by the smaller block: a 3-control-point block cannot
        # carry a cubic. Lowering is preferable to raising an error here --
        # the representation study sweeps counts down to 3 per block, and it is
        # the study's job to show what the lower counts cost.
        self.degree = int(min(degree, self.n_chord - 1, self.n_twist - 1))
        if self.degree < 1:
            raise ValueError(
                f"need at least 2 control points per block, got "
                f"n_chord={self.n_chord}, n_twist={self.n_twist}")

        # Station midpoints, not edges: each BEM strip is represented by its
        # centre, and placing a station at s = 0 or s = 1 would put it exactly
        # at the hub or the tip, where the Prandtl factors degenerate and
        # StationParams raises (see bem/station.py's r -> R discussion).
        edges = np.linspace(0.0, 1.0, self.n_stations + 1)
        self.stations = 0.5 * (edges[:-1] + edges[1:])
        self.radii = self.radius_m * (
            self.root_fraction + (1.0 - self.root_fraction) * self.stations)

        self.chord_basis = basis_matrix(self.stations, self.n_chord, self.degree)
        self.twist_basis = basis_matrix(self.stations, self.n_twist, self.degree)

    # -- shape --------------------------------------------------------------

    @property
    def n_design_variables(self):
        return self.n_chord + self.n_twist

    def split(self, design_vector):
        """Design vector -> (chord control points, twist control points)."""

        design_vector = np.asarray(design_vector)
        if design_vector.shape[-1] != self.n_design_variables:
            raise ValueError(
                f"expected {self.n_design_variables} design variables "
                f"({self.n_chord} chord + {self.n_twist} twist), got "
                f"{design_vector.shape[-1]}")
        return design_vector[:self.n_chord], design_vector[self.n_chord:]

    # -- evaluation ---------------------------------------------------------

    def chord(self, design_vector):
        """Chord at every station, metres. Linear in the design vector."""

        chord_points, _ = self.split(design_vector)
        return self.chord_basis @ chord_points

    def twist(self, design_vector):
        """Twist at every station, radians. Linear in the design vector."""

        _, twist_points = self.split(design_vector)
        return self.twist_basis @ twist_points

    def evaluate(self, design_vector):
        """(chord, twist) at every station."""

        return self.chord(design_vector), self.twist(design_vector)

    # -- derivatives --------------------------------------------------------

    def dchord_dd(self):
        """
        `dchord/dd`, shape (n_stations, n_design_variables).

        Constant -- see the module docstring. The twist block is structurally
        zero: chord does not depend on the twist control points, and saying so
        explicitly is cheaper and clearer than letting a caller discover it.
        """

        jacobian = np.zeros((self.n_stations, self.n_design_variables))
        jacobian[:, :self.n_chord] = self.chord_basis
        return jacobian

    def dtwist_dd(self):
        """`dtwist/dd`, shape (n_stations, n_design_variables). Constant."""

        jacobian = np.zeros((self.n_stations, self.n_design_variables))
        jacobian[:, self.n_chord:] = self.twist_basis
        return jacobian

    # -- fitting ------------------------------------------------------------

    def fit(self, chord_target, twist_target):
        """
        Least-squares control points reproducing given station distributions.

        This is how the analytic Schmitz blade is projected onto the
        parameterisation the optimiser actually uses, and how the
        representation study measures what each control-point count can
        represent. Because the map is linear, the fit is a single linear
        least-squares solve with no starting guess and no local minima -- the
        answer is the answer.

        Returns
        -------
        (design_vector, rms_error_chord_m, rms_error_twist_rad)
        """

        chord_target = np.asarray(chord_target, dtype=float)
        twist_target = np.asarray(twist_target, dtype=float)

        chord_points, *_ = np.linalg.lstsq(self.chord_basis, chord_target, rcond=None)
        twist_points, *_ = np.linalg.lstsq(self.twist_basis, twist_target, rcond=None)

        design_vector = np.concatenate([chord_points, twist_points])
        chord_error = self.chord_basis @ chord_points - chord_target
        twist_error = self.twist_basis @ twist_points - twist_target

        return (design_vector,
                float(np.sqrt(np.mean(chord_error ** 2))),
                float(np.sqrt(np.mean(twist_error ** 2))))

    # -- to a solvable rotor ------------------------------------------------

    def to_geometry(self, design_vector, polar_cache="sg6043"):
        """
        A `bem.rotor.RotorGeometry` for this design vector.

        Imported locally so `design/` does not pull `bem/` in at module import
        time -- the dependency runs one way, and keeping it lazy means the
        parameterisation and its tests stay usable without the solver.
        """

        from bem.rotor import RotorGeometry

        design = load_design_rotor()
        chord, twist = self.evaluate(design_vector)

        return RotorGeometry(
            r=[float(value) for value in self.radii],
            chord=[float(value) for value in chord],
            twist=[float(value) for value in twist],
            R=self.radius_m,
            polar_cache=polar_cache,
            n_blades=design.n_blades,
            r_hub=self.root_fraction * self.radius_m,
        )
