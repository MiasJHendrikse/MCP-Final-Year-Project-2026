"""
The design rotor's per-station Reynolds envelope -- what the SG6043 cache has
to cover.

This is load-bearing: the polar layer has no Reynolds clamp, and an
out-of-range lookup raises; if the cache floor
sits above the real envelope that turns a silent wrong answer into a hard crash
across the low wind-speed bins, which at this site carry most of the AEP. So
the envelope is established first, the cache bounds are set beyond it with
margin, and only then does the clamp go.

Everything here is derived from `config/` plus one geometric model. Nothing is
read from a polar cache -- the point of the module is to size a cache that does
not exist yet.

The chord model
---------------
The design rotor's chord distribution is not fixed in advance: the Schmitz
design sets a starting point and the optimiser moves it. What is fixed is the family it
comes from, so the envelope is computed over the Schmitz optimum chord

    c(r) = (16 pi r) / (B C_L) * sin^2( (1/3) arctan( R / (lambda r) ) )

which is the same model used to choose R = 2.0 m, scaled by a
band of perturbation factors standing in for everything that will move it: the
design-lift lever (chord scales as 1/C_L, so C_L in [0.71, 1.43] about
C_L = 1.0 is the same thing as a 0.7-1.4 chord
factor), the optimiser's chord control points, and the manufacturability
bounds, which were set later. The band is an assumption, stated here and
recorded in `config/polars_sg6043.yaml`, not a config value -- and the envelope
is reported as a function of it so a later, narrower band can be read straight
off rather than re-derived.

Control strategy
----------------
Reynolds number at a station is `Re = W c / nu` with `W = hypot(V, Omega r)`,
so it depends on how rotor speed is scheduled against wind speed. Three cases
are evaluated, because the machine's above-rated behaviour is not yet decided
and they do not agree:

  below-rated     variable speed at the design tip-speed ratio, and at the
                  neighbouring ratios a Cp-lambda sweep or a retuned design
                  would visit.
  above-rated     rotor speed held at its rated value -- the conventional
                  small-machine assumption, and the one that keeps tip speed
                  physical at 20 m/s.
  unregulated     constant tip-speed ratio all the way to cut-out. Not a
                  sensible machine (a 2 m rotor at lambda = 6.5 in 20 m/s wind
                  turns at ~620 rpm), but it is exactly what
                  `bem.powercurve.power_curve(tsr=...)` computes today, so a
                  cache that does not cover it will raise rather than answer.

Run directly (from `src/`):

    python -m polars.envelope

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math

import numpy as np

from config import load_design_rotor, load_site

#: Radial band the stations are placed over, as a fraction of R -- now read
#: from `config/rotor_design.yaml` rather than declared here. It used to be a
#: local `0.15`; the parameterisation needs the same span, and two
#: copies of a number that must agree is how they stop agreeing. The value is
#: unchanged, so the envelope this module computes is unchanged.
_ROOT_FRACTION = load_design_rotor().root_fraction

#: Chord perturbation band (see the module docstring). Reported alongside the
#: envelope so the assumption travels with the number.
_CHORD_FACTORS = (0.7, 1.0, 1.4)

#: Tip-speed ratios evaluated below rated, about the design value.
_TSR_OFFSETS = (-1.5, 0.0, 1.5)

#: Design lift coefficient the envelope is drawn at -- the same basis used
#: for the R = 2.0 m choice. The chord factors above cover
#: moving it.
_DESIGN_CL = 1.0

#: Wind-speed bin width for the sweep, matching the AEP bins.
_BIN_WIDTH_MS = 1.0


def schmitz_chord(r, radius, n_blades, tsr, design_cl=_DESIGN_CL):
    """
    Schmitz optimum chord at radius `r`, metres.

    Parameters
    ----------
    r : float or array_like
        Station radius, m.
    radius : float
        Rotor tip radius, m.
    n_blades : int
    tsr : float
        Design tip-speed ratio the chord distribution is drawn for.
    design_cl : float, optional
        Section lift coefficient the blade is designed at. Chord scales as
        1/C_L; see the module docstring.
    """

    r = np.asarray(r, dtype=float)
    phi = np.arctan2(radius, tsr * r) / 3.0
    return (16.0 * math.pi * r) / (n_blades * design_cl) * np.sin(phi) ** 2


def _cases(design):
    """
    (label, wind speed, rotor speed) triples spanning the operating envelope.

    See the module docstring for why three control strategies rather than one.
    """

    v_in = design.cut_in_wind_speed_ms
    v_rated = design.rated_wind_speed_ms
    v_out = design.cut_out_wind_speed_ms
    lam = design.design_tsr
    radius = design.radius_m

    n_bins = int(round((v_out - v_in) / _BIN_WIDTH_MS)) + 1
    winds = v_in + _BIN_WIDTH_MS * np.arange(n_bins)
    omega_rated = lam * v_rated / radius

    for v in winds:
        if v <= v_rated:
            for offset in _TSR_OFFSETS:
                yield "below-rated", float(v), (lam + offset) * v / radius
        else:
            yield "above-rated", float(v), omega_rated
        yield "unregulated", float(v), lam * v / radius


def reynolds_envelope(n_stations=None, chord_factors=_CHORD_FACTORS):
    """
    Per-station Reynolds extremes over the whole design operating envelope.

    Returns
    -------
    dict
        Per-case and overall (min, max) Reynolds numbers, each with the
        (case, wind speed, r/R, chord factor) that produced it, plus the
        inputs the numbers depend on.
    """

    design = load_design_rotor()
    site = load_site()
    radius = design.radius_m
    nu = site.kinematic_viscosity

    if n_stations is None:
        n_stations = design.parameterisation.n_bem_strips
    # Station midpoints, plus the tip itself: the tip has the smallest chord
    # and the highest speed, so it sets both ends of the envelope and must not
    # be missed by a midpoint grid.
    edges = np.linspace(_ROOT_FRACTION * radius, radius, n_stations + 1)
    r = np.append(0.5 * (edges[:-1] + edges[1:]), radius)

    chord_base = schmitz_chord(r, radius, design.n_blades, design.design_tsr)

    extremes = {}
    for case, v, omega in _cases(design):
        w = np.hypot(v, omega * r)
        for factor in chord_factors:
            reynolds = w * chord_base * factor / nu
            for label, index, value in (
                ("min", int(np.argmin(reynolds)), float(reynolds.min())),
                ("max", int(np.argmax(reynolds)), float(reynolds.max())),
            ):
                point = {
                    "reynolds": value,
                    "case": case,
                    "v_inf_ms": v,
                    "r_over_R": float(r[index] / radius),
                    "chord_m": float(chord_base[index] * factor),
                    "chord_factor": factor,
                }
                for key in (f"overall_{label}", f"{case}_{label}"):
                    held = extremes.get(key)
                    if held is None or (value < held["reynolds"] if label == "min"
                                        else value > held["reynolds"]):
                        extremes[key] = point

    extremes["inputs"] = {
        "radius_m": radius,
        "n_blades": design.n_blades,
        "design_tsr": design.design_tsr,
        "tsr_offsets": list(_TSR_OFFSETS),
        "cut_in_wind_speed_ms": design.cut_in_wind_speed_ms,
        "rated_wind_speed_ms": design.rated_wind_speed_ms,
        "cut_out_wind_speed_ms": design.cut_out_wind_speed_ms,
        "kinematic_viscosity": nu,
        "design_cl": _DESIGN_CL,
        "chord_factors": list(chord_factors),
        "root_fraction": _ROOT_FRACTION,
        "n_stations": int(n_stations),
    }
    return extremes


def blade_aspect_ratio(n_stations=None):
    """
    Schmitz-baseline blade aspect ratio, R / mean chord.

    Wanted by `polars.viterna.cd_max_finite_blade`, which turns it into a
    derived post-stall CD_MAX for the design rotor's cache instead of the
    fitted 1.8 the S809 cache carries.
    """

    design = load_design_rotor()
    radius = design.radius_m
    if n_stations is None:
        n_stations = design.parameterisation.n_bem_strips
    r = np.linspace(_ROOT_FRACTION * radius, radius, n_stations)
    chord = schmitz_chord(r, radius, design.n_blades, design.design_tsr)
    return float(radius / chord.mean())


def _format(point):
    return (f"Re = {point['reynolds']:>9,.0f}  ({point['case']}, "
            f"V = {point['v_inf_ms']:>4.1f} m/s, r/R = {point['r_over_R']:.2f}, "
            f"c = {point['chord_m']:.3f} m at factor {point['chord_factor']:.1f})")


def main():
    envelope = reynolds_envelope()
    inputs = envelope["inputs"]

    print("=" * 78)
    print("SG6043 design-rotor Reynolds envelope")
    print("=" * 78)
    print(f"R = {inputs['radius_m']} m, B = {inputs['n_blades']}, "
          f"lambda_design = {inputs['design_tsr']} "
          f"(+/- {max(inputs['tsr_offsets'])} below rated), "
          f"V = {inputs['cut_in_wind_speed_ms']}..{inputs['cut_out_wind_speed_ms']} m/s "
          f"in {_BIN_WIDTH_MS:g} m/s bins, rated at "
          f"{inputs['rated_wind_speed_ms']} m/s")
    print(f"Schmitz chord at C_L = {inputs['design_cl']}, perturbed by "
          f"{inputs['chord_factors']}; {inputs['n_stations']} stations over "
          f"r/R = {inputs['root_fraction']}..1.0")
    print(f"nu = {inputs['kinematic_viscosity']:.4e} m^2/s (site config, 1800 m)")
    print(f"Schmitz-baseline blade aspect ratio = {blade_aspect_ratio():.1f}")
    print()

    for case in ("below-rated", "above-rated", "unregulated"):
        print(f"{case}:")
        print(f"  min  {_format(envelope[f'{case}_min'])}")
        print(f"  max  {_format(envelope[f'{case}_max'])}")
    print()
    print("overall:")
    print(f"  min  {_format(envelope['overall_min'])}")
    print(f"  max  {_format(envelope['overall_max'])}")
    return envelope


if __name__ == "__main__":
    main()
