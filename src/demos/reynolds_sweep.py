"""
Reynolds sweep for a single airfoil.

Runs XFOIL for NACA 2412 across a range of Reynolds numbers and plots a 2x2 grid:
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
def plot_polars_by_reynolds(polars_by_re, airfoil_name, save_path=None):
    """
    Plot a 2x2 grid of aerodynamic curves for a single airfoil across multiple Re.

    Parameters
    ----------
    polars_by_re : dict
        Dictionary mapping Reynolds number (int) to polar array (Nx7 NumPy array,
        as returned by run_xfoil_polar). Entries with None values are skipped.
    airfoil_name : str
        Used in the figure title and as the saved-file stem.
    save_path : str or None
        If provided, the figure is saved to this path. Otherwise just displayed.
    """

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Aerodynamic Polars — {airfoil_name}", fontsize=14, fontweight="bold")

    # Pick a colourmap with one distinct colour per Reynolds number
    re_values = sorted(re for re, p in polars_by_re.items() if p is not None)
    colours = plt.cm.viridis(np.linspace(0, 0.9, len(re_values)))

    for re, colour in zip(re_values, colours):
        polar = polars_by_re[re]
        alpha = polar[:, 0]
        cl = polar[:, 1]
        cd = polar[:, 2]
        # Avoid division-by-zero where cd is tiny
        ld = np.where(cd > 1e-6, cl / cd, 0.0)

        label = f"Re = {re:,}"

        # CL vs alpha
        axes[0, 0].plot(alpha, cl, "o-", color=colour, label=label, markersize=3)

        # CD vs alpha
        axes[0, 1].plot(alpha, cd, "o-", color=colour, label=label, markersize=3)

        # Drag polar (CL vs CD)
        axes[1, 0].plot(cd, cl, "o-", color=colour, label=label, markersize=3)

        # L/D vs alpha
        axes[1, 1].plot(alpha, ld, "o-", color=colour, label=label, markersize=3)

    # --- Formatting ---
    axes[0, 0].set_xlabel(r"$\alpha$ [deg]")
    axes[0, 0].set_ylabel(r"$C_L$")
    axes[0, 0].set_title("Lift Coefficient")
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=9)

    axes[0, 1].set_xlabel(r"$\alpha$ [deg]")
    axes[0, 1].set_ylabel(r"$C_D$")
    axes[0, 1].set_title("Drag Coefficient")
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=9)

    axes[1, 0].set_xlabel(r"$C_D$")
    axes[1, 0].set_ylabel(r"$C_L$")
    axes[1, 0].set_title("Drag Polar")
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=9)

    axes[1, 1].set_xlabel(r"$\alpha$ [deg]")
    axes[1, 1].set_ylabel(r"$C_L / C_D$")
    axes[1, 1].set_title("Lift-to-Drag Ratio")
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(fontsize=9)

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
    print("Generating polars: NACA 2412 across Reynolds sweep")
    print("=" * 60)

    reynolds_list = [50_000, 100_000, 200_000, 500_000, 1_000_000]
    polars_by_re = {}

    for re in reynolds_list:
        print(f"  Running Re = {re:,} ...")
        polar = run_xfoil_polar(
            airfoil_cmd="NACA 2412",
            reynolds=re,
            alpha_min=-4,
            alpha_max=14,
            alpha_step=0.5,
            polar_path=os.path.join("xfoil_demos", f"polar_naca2412_re{re}.txt"),
            n_iter=200,
        )
        polars_by_re[re] = polar

    plot_polars_by_reynolds(
        polars_by_re,
        airfoil_name="NACA 2412",
        save_path=os.path.join(RESULTS_DIR, "xfoil_demos", "polars_naca2412_reynolds_sweep.png"),
    )

    print()
    print("Done. PNG figure saved to the 'results' folder.")
