"""
Airfoil comparison at a fixed Reynolds number.

Runs XFOIL for several NACA airfoils at Re = 200,000 and plots a 2x2 grid:
  - CL vs alpha
  - CD vs alpha
  - CL vs CD (drag polar)
  - L/D vs alpha
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from xfoil.xfoil_runner import run_xfoil_polar, RESULTS_DIR


# ---------------------------------------------------------------------------
# Plotting helper
# ---------------------------------------------------------------------------
def plot_airfoil_comparison(polars_by_airfoil, reynolds, save_path=None):
    """
    Compare multiple airfoils at a fixed Reynolds number on a 2x2 grid.

    Parameters
    ----------
    polars_by_airfoil : dict
        Dictionary mapping airfoil name (str) to polar array.
    reynolds : float
        The Reynolds number shared by all polars (used in the title).
    save_path : str or None
        If provided, save the figure to this path.
    """

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Airfoil Comparison at Re = {reynolds:,}", fontsize=14, fontweight="bold")

    names = [n for n, p in polars_by_airfoil.items() if p is not None]
    colours = plt.cm.tab10(np.linspace(0, 1, len(names)))

    for name, colour in zip(names, colours):
        polar = polars_by_airfoil[name]
        alpha = polar[:, 0]
        cl = polar[:, 1]
        cd = polar[:, 2]
        ld = np.where(cd > 1e-6, cl / cd, 0.0)

        axes[0, 0].plot(alpha, cl, "o-", color=colour, label=name, markersize=3)
        axes[0, 1].plot(alpha, cd, "o-", color=colour, label=name, markersize=3)
        axes[1, 0].plot(cd, cl, "o-", color=colour, label=name, markersize=3)
        axes[1, 1].plot(alpha, ld, "o-", color=colour, label=name, markersize=3)

    axes[0, 0].set_xlabel(r"$\alpha$ [deg]"); axes[0, 0].set_ylabel(r"$C_L$")
    axes[0, 0].set_title("Lift Coefficient"); axes[0, 0].grid(True, alpha=0.3); axes[0, 0].legend()

    axes[0, 1].set_xlabel(r"$\alpha$ [deg]"); axes[0, 1].set_ylabel(r"$C_D$")
    axes[0, 1].set_title("Drag Coefficient"); axes[0, 1].grid(True, alpha=0.3); axes[0, 1].legend()

    axes[1, 0].set_xlabel(r"$C_D$"); axes[1, 0].set_ylabel(r"$C_L$")
    axes[1, 0].set_title("Drag Polar"); axes[1, 0].grid(True, alpha=0.3); axes[1, 0].legend()

    axes[1, 1].set_xlabel(r"$\alpha$ [deg]"); axes[1, 1].set_ylabel(r"$C_L / C_D$")
    axes[1, 1].set_title("Lift-to-Drag Ratio"); axes[1, 1].grid(True, alpha=0.3); axes[1, 1].legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"  Figure saved to {save_path}")

    plt.show()


# ---------------------------------------------------------------------------
# Main: generate data and plot
# ---------------------------------------------------------------------------
if __name__ == "__main__":

    print("=" * 60)
    print("Generating polars: airfoil comparison at Re = 200,000")
    print("=" * 60)

    airfoils = ["NACA 0012", "NACA 2412", "NACA 4412", "NACA 23012"]
    polars_by_airfoil = {}

    for foil in airfoils:
        print(f"  Running {foil} ...")
        # Build a filesystem-safe filename from the airfoil command
        fname = foil.replace(" ", "_").lower()
        polar = run_xfoil_polar(
            airfoil_cmd=foil,
            reynolds=200_000,
            alpha_min=-6,
            alpha_max=14,
            alpha_step=0.5,
            polar_path=os.path.join("xfoil_demos", f"polar_{fname}_re200k.txt"),
            n_iter=200,
        )
        polars_by_airfoil[foil] = polar

    plot_airfoil_comparison(
        polars_by_airfoil,
        reynolds=200_000,
        save_path=os.path.join(RESULTS_DIR, "xfoil_demos", "polars_airfoil_comparison_re200k.png"),
    )

    print()
    print("Done. PNG figure saved to the 'results' folder.")
