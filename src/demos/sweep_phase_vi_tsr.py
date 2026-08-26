"""
TSR sweep of the project's own BEM solver (bem.rotor.solve_rotor) over the
NREL Phase VI Sequence S blade -- companion to validation/compare_qblade.py,
which runs a single operating point. This sweeps TSR = 1.0 .. 12.0 in 0.2
steps at fixed V_inf (i.e. varying rotor speed, same convention as a
QBlade Cp-lambda / Ct-lambda sweep) and plots:

  1. Cp vs TSR
  2. Ct vs TSR
  3. Angle of attack vs radius, one curve per TSR in TSR_FAMILY
  4. Axial induction factor a vs radius, one curve per TSR in TSR_FAMILY

("Axial induction" is the factor `a` in station.py's result dict, not the
Prandtl tip/hub-loss factor `F` -- station.py returns both, but only `a` is
what's conventionally called the axial induction factor.)

Plots 3-4 use radius as the x-axis; a single TSR sweep produces 56 spanwise
profiles (one per TSR step), too many to show at once, so only a handful
of representative TSRs (TSR_FAMILY below) are drawn as separate curves --
edit that list to change which operating points are shown.

As in validation/compare_qblade.py, every station is forced onto a single fixed
S809 polar (default Re=100,000, matching S809.plr; pass --reynolds 500000 to
match S809_Re500k.plr instead), so this stays comparable to a QBlade run
station-for-station.

Run from src/: `python -m demos.sweep_phase_vi_tsr [--reynolds 100000|500000]`
"""

import argparse
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

from bem.airfoil import S809Polar  # noqa: E402
from bem.rotor import phase_vi_geometry, solve_rotor  # noqa: E402
from config import load_phase_vi_rotor  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Operating point
# ---------------------------------------------------------------------------
# The Phase VI rotor's operating condition, including its sea-level air, comes
# from config/rotor_phase_vi.yaml -- not from constants declared here. Four
# scripts used to carry their own copy of the density and viscosity.
PHASE_VI = load_phase_vi_rotor()
V_INF_MS = PHASE_VI.reference_wind_speed_ms
TIP_PITCH_DEG = PHASE_VI.tip_pitch_deg

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reynolds", type=int, default=100_000, choices=[100_000, 500_000],
                    help="Fixed Reynolds number for every station (default: 100000).")
args = parser.parse_args()
FIXED_REYNOLDS = args.reynolds

TSR_MIN, TSR_MAX, TSR_STEP = 1.0, 12.0, 0.2
TSR_FAMILY = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]  # spanwise curves drawn for these

# Okabe-Ito colorblind-safe categorical palette, fixed order (never cycled/reassigned).
CATEGORICAL = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7"]

# Generated sweep data and plots go to results/, not next to the source.
OUT_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "results", "phase_vi_sweep"))
os.makedirs(OUT_DIR, exist_ok=True)

geometry = phase_vi_geometry(tip_pitch_deg=TIP_PITCH_DEG)
airfoils = [S809Polar(FIXED_REYNOLDS) for _ in geometry.r]

n_steps = round((TSR_MAX - TSR_MIN) / TSR_STEP) + 1
tsr_values = [round(TSR_MIN + i * TSR_STEP, 4) for i in range(n_steps)]

ct_rotor, cp_rotor = [], []
spanwise = {}  # tsr -> result dict, only kept for TSR_FAMILY

for tsr in tsr_values:
    result = solve_rotor(geometry, tsr=tsr, v_inf=V_INF_MS,
                         air_density=PHASE_VI.air_density,
                         kinematic_viscosity=PHASE_VI.kinematic_viscosity,
                         airfoils=airfoils)
    ct_rotor.append(result["Ct"])
    cp_rotor.append(result["Cp"])
    if any(abs(tsr - fam) < 1e-6 for fam in TSR_FAMILY):
        spanwise[tsr] = result

# ---------------------------------------------------------------------------
# Write the raw sweep + spanwise data alongside the plot, for anything not
# covered by the four panels below.
# ---------------------------------------------------------------------------
sweep_csv = os.path.join(OUT_DIR, "phase_vi_tsr_sweep_Re" + str(FIXED_REYNOLDS) + ".csv")
with open(sweep_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["tsr", "Ct", "Cp"])
    writer.writerows(zip(tsr_values, ct_rotor, cp_rotor))

spanwise_csv = os.path.join(OUT_DIR, "phase_vi_tsr_family_spanwise_Re" + str(FIXED_REYNOLDS) + ".csv")
with open(spanwise_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["tsr", "r_m", "alpha_deg", "a"])
    for tsr, result in spanwise.items():
        for s in result["stations"]:
            writer.writerow([tsr, s["r"], s["alpha"] * 180.0 / np.pi, s["a"]])

# ---------------------------------------------------------------------------
# Plot: 2x2 grid, Cp/Ct vs TSR on top, AoA/a vs radius (TSR family) on bottom.
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 10,
    "axes.edgecolor": "#888888",
    "axes.labelcolor": "#333333",
    "xtick.color": "#555555",
    "ytick.color": "#555555",
    "axes.grid": True,
    "grid.color": "#e0e0e0",
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
ax_cp, ax_ct, ax_aoa, ax_a = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

ax_cp.plot(tsr_values, cp_rotor, color=CATEGORICAL[0], linewidth=2)
ax_cp.set_xlabel("TSR [-]")
ax_cp.set_ylabel("Cp [-]")
ax_cp.set_title("Power coefficient vs tip-speed ratio")

ax_ct.plot(tsr_values, ct_rotor, color=CATEGORICAL[1], linewidth=2)
ax_ct.set_xlabel("TSR [-]")
ax_ct.set_ylabel("Ct [-]")
ax_ct.set_title("Thrust coefficient vs tip-speed ratio")

for i, tsr in enumerate(TSR_FAMILY):
    if tsr not in spanwise:
        continue
    result = spanwise[tsr]
    r = [s["r"] for s in result["stations"]]
    alpha_deg = [s["alpha"] * 180.0 / np.pi for s in result["stations"]]
    a = [s["a"] for s in result["stations"]]
    color = CATEGORICAL[i % len(CATEGORICAL)]
    ax_aoa.plot(r, alpha_deg, color=color, linewidth=2, marker="o", markersize=4, label=f"TSR={tsr:g}")
    ax_a.plot(r, a, color=color, linewidth=2, marker="o", markersize=4, label=f"TSR={tsr:g}")

ax_aoa.set_xlabel("r [m]")
ax_aoa.set_ylabel("Angle of attack [deg]")
ax_aoa.set_title("Spanwise AoA")
ax_aoa.legend(frameon=False, fontsize=8)

ax_a.set_xlabel("r [m]")
ax_a.set_ylabel("Axial induction factor a [-]")
ax_a.set_title("Spanwise axial induction")
ax_a.legend(frameon=False, fontsize=8)

fig.suptitle("NREL Phase VI (Sequence S) -- bem.rotor.solve_rotor sweep", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])

png_path = os.path.join(OUT_DIR, "phase_vi_tsr_sweep_Re" + str(FIXED_REYNOLDS) + ".png")
fig.savefig(png_path, dpi=150)
print(f"Sweep data: {sweep_csv}")
print(f"Spanwise (TSR family) data: {spanwise_csv}")
print(f"Plot: {png_path}")
