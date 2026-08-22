"""
Run the project's own BEM solver (bem.rotor.solve_rotor) over the same NREL
Phase VI Sequence S blade that is exported to QBlade as phase_vi_blade*.bld /
S809*.plr (see bem/export_qblade.py), so the two solvers' results can be
compared directly.

To keep this an apples-to-apples comparison with the QBlade run, every
station is forced onto the exact same single-Reynolds S809 polar the matching
.plr contains (rotor.solve_rotor normally estimates a *per-station* Reynolds
number and interpolates the real XFOIL cache continuously in Re -- see
bem/rotor.py's solve_rotor docstring, "airfoils" parameter). Forcing a single
fixed Re here means both solvers are querying the same 2D polar curve; any
difference in Ct/Cp or per-station loads then comes from the solvers
themselves (BEM implementation, induction/tip-loss treatment, extrapolation
used past the XFOIL-converged band), not from using different input
aerodynamic data.

Two Reynolds numbers are provided as matched pairs with data/qblade/*.bld:

  100,000  (phase_vi_blade.bld / S809.plr) -- the original comparison point.
  500,000  (phase_vi_blade_Re500k.bld / S809_Re500k.plr) -- matches this
           rotor's actual per-station Reynolds range (577k-945k estimated at
           the Sequence S rated operating point), which all clamp to the
           cache's Re=500k ceiling -- so this pair reproduces solve_rotor's
           own default (continuous per-station Re) behaviour exactly, and is
           the more physically meaningful of the two for comparing against
           QBlade.

Edit V_INF_MS and RPM below to match whatever operating point you run in
QBlade -- tip-speed ratio is derived from them the same way QBlade derives
it internally.

Run from src/: `python -m validation.compare_qblade [--reynolds 100000|500000]`
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

from bem.airfoil import S809Polar  # noqa: E402
from bem.rotor import (  # noqa: E402
    PHASE_VI_RATED_RPM,
    PHASE_VI_SEQUENCE_S_TIP_PITCH_DEG,
    phase_vi_geometry,
    solve_rotor,
)

# Results are written next to the QBlade .bld/.plr exports they pair with, so
# a comparison run and the files it was run against stay together.
_HERE = os.path.dirname(os.path.abspath(__file__))
QBLADE_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "data", "qblade"))

# ---------------------------------------------------------------------------
# Operating point -- match these to your QBlade simulation settings.
# ---------------------------------------------------------------------------
V_INF_MS = 7.0                          # freestream wind speed, m/s
RPM = PHASE_VI_RATED_RPM                # 71.63 RPM, Sequence S rated speed
TIP_PITCH_DEG = PHASE_VI_SEQUENCE_S_TIP_PITCH_DEG  # 3.0 deg, matches the .bld export
AIR_DENSITY = 1.225                     # kg/m^3, QBlade's default too

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reynolds", type=int, default=100_000, choices=[100_000, 500_000],
                    help="Which fixed Reynolds number to pin every station to, matching "
                         "the .plr of the same value (default: 100000).")
args = parser.parse_args()
FIXED_REYNOLDS = args.reynolds

geometry = phase_vi_geometry(tip_pitch_deg=TIP_PITCH_DEG)
omega_rad_s = RPM * 2.0 * 3.141592653589793 / 60.0
tsr = omega_rad_s * geometry.R / V_INF_MS

# One S809Polar per station, all pinned to the same Re so both solvers query
# the identical XFOIL curve (see module docstring).
airfoils = [S809Polar(FIXED_REYNOLDS) for _ in geometry.r]

result = solve_rotor(geometry, tsr=tsr, v_inf=V_INF_MS, air_density=AIR_DENSITY, airfoils=airfoils)

print(f"Fixed Reynolds = {FIXED_REYNOLDS:,} (matches "
      f"{'S809.plr' if FIXED_REYNOLDS == 100_000 else 'S809_Re500k.plr'})")
print(f"V_inf = {V_INF_MS} m/s, RPM = {RPM}, TSR = {tsr:.3f}, tip pitch = {TIP_PITCH_DEG} deg")
print(f"Rotor Ct = {result['Ct']:.4f}")
print(f"Rotor Cp = {result['Cp']:.4f}")
print()
print(f"{'r [m]':>8} {'phi [deg]':>10} {'alpha [deg]':>12} {'a [-]':>8} {'a_prime [-]':>12} "
      f"{'Cl [-]':>8} {'Cd [-]':>8} {'Cn [-]':>8} {'F [-]':>8}")
for s in result["stations"]:
    print(f"{s['r']:8.4f} {s['phi']*57.29577951308232:10.3f} {s['alpha']*57.29577951308232:12.3f} "
          f"{s['a']:8.4f} {s['a_prime']:12.4f} {s['Cl']:8.4f} {s['Cd']:8.4f} {s['Cn']:8.4f} {s['F']:8.4f}")

# CSV alongside the .bld/.plr exports for easy overlay against QBlade's
# per-station BEM output (QBlade: Simulation -> BEM -> export blade results),
# one file per Reynolds so the 100k and 500k runs don't overwrite each other.
out_csv = os.path.join(QBLADE_DIR, f"bem_solver_phase_vi_result_Re{FIXED_REYNOLDS}.csv")
with open(out_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["r_m", "phi_deg", "alpha_deg", "a", "a_prime", "Cl", "Cd", "Cn", "Ct_tangential", "F", "reynolds"])
    for s in result["stations"]:
        writer.writerow([
            s["r"], s["phi"] * 57.29577951308232, s["alpha"] * 57.29577951308232,
            s["a"], s["a_prime"], s["Cl"], s["Cd"], s["Cn"], s["Ct"], s["F"], s["reynolds"],
        ])

print(f"\nPer-station results written to {out_csv}")
