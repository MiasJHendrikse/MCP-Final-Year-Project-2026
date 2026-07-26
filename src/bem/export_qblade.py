"""
Export the NREL Phase VI blade and its S809 polar to QBlade's file formats,
for manual cross-checking of this project's BEM solver against QBlade CE.

Writes two files:

  <name>.plr  -- QBlade "Multi RE Polar File": one Reynolds number's S809
                 polar over the full -180..180 deg range QBlade requires,
                 built from the XFOIL cache in data/polars/s809/ and
                 Viterna-extrapolated outside the converged range.
  <name>.bld  -- QBlade "Blade Definition File": the Phase VI Sequence S
                 station table (rotor.phase_vi_geometry), every station
                 pointing at the .plr above.

This regenerates, reproducibly, what was previously done by an uncommitted
ad hoc script -- which is why the original pair carried the corrupted
pre-2026-07-26 Re=100k polar with no record of how it had been produced.

Extrapolation
-------------
QBlade needs a polar defined over the whole -180..180 deg circle; XFOIL only
converges a narrow band. The scheme below reproduces exactly what the original
export used, so a regenerated file differs from its predecessor only where the
underlying XFOIL data differs:

  alpha_n <= a <= alpha_s   the cached XFOIL data verbatim (alpha_n, alpha_s
                            are the most negative / most positive converged
                            alphas).
  alpha_s < a <= 90         classical Viterna, with Cd_max = 1.8 and the
                            A1/A2/B1/B2 coefficients matched to the airfoil's
                            state at alpha_s.
  90 < a <= 180             flat-plate reversed-flow branch,
                            Cl = -0.63 sin(2(180-a)), Cd mirrored from the
                            Viterna branch and floored at 0.020.
  -alpha_s <= a < alpha_n   linear ramp in both Cl and Cd, from the XFOIL
                            endpoint at alpha_n onto the negative branch's
                            value at -alpha_s. (XFOIL rarely reaches negative
                            stall, so there is no measured anchor to Viterna
                            against on this side.)
  a < -alpha_s              negative branch: Cl = -0.7 * Cl_viterna(|a|),
                            Cd = Cd_viterna(|a|); beyond -90 deg it is the
                            negation of the positive branch.

  Cm is held constant outside the converged band at whichever endpoint value
  is nearer -- QBlade's BEM does not use Cm for power/thrust, so this is a
  placeholder rather than a model.

Known asymmetry with our own solver, deliberately preserved
-----------------------------------------------------------
bem.airfoil.S809Polar *clamps* alpha to the cached range (flat extrapolation)
rather than extrapolating, so past the converged band the two solvers are
fed genuinely different aerodynamics. At a normal operating point neither
solver converges to an alpha out there, so the comparison stays valid -- but
if a station does stall in QBlade, that is the first thing to suspect.

Run directly: `python -m bem.export_qblade` (from src/).

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import datetime
import math
import os

import numpy as np

from bem.rotor import (
    PHASE_VI_SEQUENCE_S_TIP_PITCH_DEG,
    PHASE_VI_N_BLADES,
    phase_vi_geometry,
)
from xfoil.polar_lookup import DATA_DIR

#: Post-stall drag maximum used by the Viterna fit. 1.8 is the conventional
#: value for a finite blade and is what the original export used.
CD_MAX = 1.8

#: Reversed-flow (a > 90 deg) lift amplitude, and the Cd floor applied there --
#: a reversed airfoil aligned with the flow still has some profile drag.
_REVERSED_CL_AMPLITUDE = 0.63
_REVERSED_CD_FLOOR = 0.020

#: Scaling applied to the mirrored Viterna curve on the negative-alpha side.
_NEGATIVE_CL_SCALE = 0.7

#: S809 maximum thickness, percent chord (matches data/airfoils/s809.dat,
#: measured 20.99% -- see the 2026-07-26 journal entry).
S809_THICKNESS_PCT = 21.0

_COL = 20  # QBlade's fixed column width in these files


def _viterna_coefficients(alpha_s_deg, cl_s, cd_s):
    """Classical Viterna A1/A2/B1/B2, matched to the airfoil state at alpha_s."""
    a = math.radians(alpha_s_deg)
    sin_a, cos_a = math.sin(a), math.cos(a)
    b1 = CD_MAX
    b2 = (cd_s - CD_MAX * sin_a ** 2) / cos_a
    a1 = b1 / 2.0
    a2 = (cl_s - CD_MAX * sin_a * cos_a) * sin_a / cos_a ** 2
    return a1, a2, b1, b2


def _viterna(alpha_deg, coeffs):
    """Viterna Cl, Cd at alpha (deg), valid for 0 < alpha <= 90."""
    a1, a2, b1, b2 = coeffs
    r = math.radians(alpha_deg)
    sin_r, cos_r = math.sin(r), math.cos(r)
    cl = a1 * math.sin(2 * r) + a2 * cos_r ** 2 / sin_r
    cd = b1 * sin_r ** 2 + b2 * cos_r
    return cl, cd


def build_full_range_polar(alpha_xf, cl_xf, cd_xf, cm_xf, step=0.5):
    """
    Extend a converged XFOIL polar to the full -180..180 deg range QBlade needs.

    Parameters
    ----------
    alpha_xf, cl_xf, cd_xf, cm_xf : array_like
        The cached XFOIL polar, ascending in alpha.
    step : float, optional
        Output alpha resolution in degrees (default 0.5).

    Returns
    -------
    numpy.ndarray
        Columns [alpha, cl, cd, cm] over -180..180 deg.
    """

    alpha_xf = np.asarray(alpha_xf, dtype=float)
    cl_xf = np.asarray(cl_xf, dtype=float)
    cd_xf = np.asarray(cd_xf, dtype=float)
    cm_xf = np.asarray(cm_xf, dtype=float)

    alpha_s, cl_s, cd_s = alpha_xf[-1], cl_xf[-1], cd_xf[-1]
    alpha_n, cl_n, cd_n = alpha_xf[0], cl_xf[0], cd_xf[0]
    cm_hi, cm_lo = cm_xf[-1], cm_xf[0]

    coeffs = _viterna_coefficients(alpha_s, cl_s, cd_s)

    # Anchor for the negative-side ramp: the negative branch evaluated at
    # -alpha_s, which is where that branch takes over.
    cl_anchor_pos, cd_anchor = _viterna(alpha_s, coeffs)
    cl_anchor = -_NEGATIVE_CL_SCALE * cl_anchor_pos

    n = int(round(360.0 / step)) + 1
    out = np.empty((n, 4))

    for i in range(n):
        a = -180.0 + i * step
        a = round(a, 6)

        if alpha_n - 1e-9 <= a <= alpha_s + 1e-9:
            # Inside the converged band: use the XFOIL data as-is. np.interp is
            # exact at the sample points and covers the odd non-converged alpha
            # the cache may be missing internally.
            cl = float(np.interp(a, alpha_xf, cl_xf))
            cd = float(np.interp(a, alpha_xf, cd_xf))
            cm = float(np.interp(a, alpha_xf, cm_xf))

        elif a > alpha_s:
            cm = cm_hi
            if a <= 90.0:
                cl, cd = _viterna(a, coeffs)
            else:
                mirror = 180.0 - a
                cl = -_REVERSED_CL_AMPLITUDE * math.sin(2 * math.radians(mirror))
                if mirror <= 1e-9:
                    cd = _REVERSED_CD_FLOOR
                else:
                    cd = max(_viterna(mirror, coeffs)[1], _REVERSED_CD_FLOOR)

        else:  # a < alpha_n
            cm = cm_lo
            if a >= -alpha_s:
                # Linear ramp from the XFOIL endpoint onto the negative branch.
                t = (a - alpha_n) / (-alpha_s - alpha_n)
                cl = cl_n + t * (cl_anchor - cl_n)
                cd = cd_n + t * (cd_anchor - cd_n)
            else:
                mag = -a
                if mag <= 90.0:
                    cl_p, cd = _viterna(mag, coeffs)
                    cl = -_NEGATIVE_CL_SCALE * cl_p
                else:
                    mirror = 180.0 - mag
                    cl = _REVERSED_CL_AMPLITUDE * math.sin(2 * math.radians(mirror))
                    if mirror <= 1e-9:
                        cd = _REVERSED_CD_FLOOR
                    else:
                        cd = max(_viterna(mirror, coeffs)[1], _REVERSED_CD_FLOOR)

        out[i] = (a, cl, cd, cm)

    return out


def _stamp():
    now = datetime.datetime.now()
    return now.strftime("%H:%M:%S"), now.strftime("%d.%m.%Y")


def _field(value, width=_COL):
    return f"{value:<{width}}"


def write_plr(path, polar, polar_name, foil_name, reynolds,
              thickness_pct=S809_THICKNESS_PCT):
    """Write a QBlade Multi RE Polar File."""
    time_s, date_s = _stamp()
    lines = [
        "----------------------------------------QBlade Multi RE Polar File"
        "--------------------------------------------------",
        "Generated with : bem/export_qblade.py",
        "Archive Format: 310001",
        f"Time : {time_s}",
        f"Date : {date_s}",
        "",
        "----------------------------------------Object Names"
        "----------------------------------------------------------------",
        f"{_field(polar_name, 42)}POLARNAME",
        f"{_field(foil_name, 42)}FOILNAME",
        "",
        "----------------------------------------Parameters"
        "------------------------------------------------------------------",
        f"{_field(f'{thickness_pct:.1f}', 42)}THICKNESS",
        f"{_field('0', 42)}ISDECOMPOSED",
        f"REYNOLDS            {reynolds:.4E}",
        "",
        "----------------------------------------Polar Data"
        "------------------------------------------------------------------",
        "".join(_field(h) for h in ("AOA", "CL", "CD", "CM")),
    ]
    for a, cl, cd, cm in polar:
        lines.append("".join(_field(f"{v:.6f}") for v in (a, cl, cd, cm)))

    with open(path, "w", newline="\r\n") as f:
        f.write("\n".join(lines))
    return path


def write_bld(path, geometry, object_name, polar_filename, taxis=0.25):
    """Write a QBlade Blade Definition File for the given rotor geometry."""
    time_s, date_s = _stamp()
    lines = [
        "----------------------------------------QBlade Blade Definition File"
        "------------------------------------------------------",
        "Generated with : bem/export_qblade.py",
        "Archive Format: 310002",
        f"Time : {time_s}",
        f"Date : {date_s}",
        "",
        "----------------------------------------Object Name"
        "-----------------------------------------------------------------",
        f"{_field(object_name, 42)}OBJECTNAME         - the name of the blade object",
        "",
        "----------------------------------------Parameters"
        "------------------------------------------------------------------",
        f"{_field('HAWT', 42)}ROTORTYPE          - the rotor type",
        f"{_field(str(PHASE_VI_N_BLADES), 42)}NUMBLADES          - number of blades",
        "",
        "----------------------------------------Blade Data"
        "-------------------------------------------------------------------",
        "".join(_field(h) for h in ("POS [m]", "CHORD [m]", "TWIST [deg]",
                                    "OFFSET_X [m]", "OFFSET_Y [m]", "TAXIS [-]"))
        + "POLAR_FILE",
    ]
    for r, chord, twist_rad in zip(geometry.r, geometry.chord, geometry.twist):
        row = [f"{r:.4f}", f"{chord:.4f}", f"{math.degrees(twist_rad):.4f}",
               "0.0000", "0.0000", f"{taxis:.4f}"]
        lines.append("".join(_field(v) for v in row) + polar_filename)

    with open(path, "w", newline="\r\n") as f:
        f.write("\n".join(lines) + "\n")
    return path


def export(out_dir, reynolds=100_000, tip_pitch_deg=PHASE_VI_SEQUENCE_S_TIP_PITCH_DEG,
           object_name="PhaseVI_SequenceS", stem="phase_vi_blade", plr_name="S809.plr"):
    """
    Write <stem>.bld and <plr_name> into out_dir from the current polar cache.

    Parameters
    ----------
    out_dir : str
        Destination directory.
    reynolds : int
        Which cached Reynolds number to embed in the polar file.
    tip_pitch_deg : float
        Collective pitch passed to phase_vi_geometry.
    object_name : str
        QBlade blade-object name written into the .bld.
    stem : str
        Basename (no extension) for the .bld.
    plr_name : str
        Filename for the .plr, and the value written into the .bld's
        POLAR_FILE column — the two must agree or QBlade cannot resolve it.

    Returns
    -------
    tuple of str
        (bld_path, plr_path)
    """

    csv_path = os.path.join(DATA_DIR, "polars", "s809", f"S809_Re{reynolds}.csv")
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    polar = build_full_range_polar(data[:, 0], data[:, 1], data[:, 2], data[:, 3])

    plr_path = write_plr(os.path.join(out_dir, plr_name), polar,
                         polar_name=f"S809_Re{reynolds}", foil_name="S809",
                         reynolds=float(reynolds))

    geometry = phase_vi_geometry(tip_pitch_deg=tip_pitch_deg)
    bld_path = write_bld(os.path.join(out_dir, f"{stem}.bld"), geometry,
                         object_name=object_name, polar_filename=plr_name)

    return bld_path, plr_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export Phase VI blade + S809 polar for QBlade.")
    parser.add_argument("--out-dir", default=os.path.join(os.path.expanduser("~"), "Downloads"),
                        help="Directory to write the .bld and .plr into (default: ~/Downloads).")
    parser.add_argument("--reynolds", type=int, default=100_000,
                        help="Which cached Reynolds number to export (default: 100000).")
    parser.add_argument("--stem", default="phase_vi_blade",
                        help="Basename for the .bld (default: phase_vi_blade).")
    parser.add_argument("--plr-name", default="S809.plr",
                        help="Filename for the .plr (default: S809.plr).")
    parser.add_argument("--object-name", default="PhaseVI_SequenceS",
                        help="QBlade blade-object name (default: PhaseVI_SequenceS).")
    args = parser.parse_args()

    bld, plr = export(args.out_dir, reynolds=args.reynolds, stem=args.stem,
                      plr_name=args.plr_name, object_name=args.object_name)
    print(f"Wrote {bld}")
    print(f"Wrote {plr}")
