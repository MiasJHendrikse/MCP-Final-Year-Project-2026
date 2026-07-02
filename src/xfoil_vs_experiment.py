"""
XFOIL validation against published experimental data.

Runs XFOIL for NACA 0012 at Re = 6,000,000 and overlays the result on digitised
wind-tunnel data from Abbott & von Doenhoff, "Theory of Wing Sections". Produces a
2x2 grid (lift, drag, drag polar, lift error) and prints the mean/max lift error.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from xfoil_runner import run_xfoil_polar, RESULTS_DIR


# ---------------------------------------------------------------------------
# Experimental data loader
# ---------------------------------------------------------------------------
def _parse_two_column_dat(path):
    """
    Parse a whitespace-separated two-column .dat file into two NumPy arrays.

    Lines starting with "#" (comments) and the "variables=" header line are
    skipped, as are blank lines and trailing whitespace.

    Parameters
    ----------
    path : str
        Path to the .dat file.

    Returns
    -------
    tuple of numpy.ndarray
        (col0, col1) as 1D float arrays.
    """

    col0, col1 = [], []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("variables="):
                continue
            parts = line.split()
            col0.append(float(parts[0]))
            col1.append(float(parts[1]))

    return np.array(col0), np.array(col1)


def load_abbott_data(cl_path, cd_path):
    """
    Load digitised Abbott & von Doenhoff experimental data for NACA 0012.

    The two files use different formats and are NOT a single combined table:
      - cl_path : rows of (alpha_deg, cl)
      - cd_path : rows of (cl, cd)  — drag is a polar tabulated against lift
        coefficient, NOT against alpha, so its cl column is kept separate.

    Parameters
    ----------
    cl_path : str
        Path to the lift file (alpha vs cl).
    cd_path : str
        Path to the drag-polar file (cl vs cd).

    Returns
    -------
    dict
        Keys "alpha", "cl" (from the lift file) and "cl_for_cd", "cd" (from the
        drag file), each a 1D NumPy array.

    Raises
    ------
    FileNotFoundError
        If either file is missing.
    """

    if not os.path.exists(cl_path):
        raise FileNotFoundError(f"Experimental lift file not found: {cl_path}")
    if not os.path.exists(cd_path):
        raise FileNotFoundError(f"Experimental drag file not found: {cd_path}")

    alpha, cl = _parse_two_column_dat(cl_path)
    cl_for_cd, cd = _parse_two_column_dat(cd_path)

    return {"alpha": alpha, "cl": cl, "cl_for_cd": cl_for_cd, "cd": cd}


def load_ladson_data(path):
    """
    Load Ladson (NASA TM 4074) experimental data for NACA 0012.

    Unlike the Abbott files, this is a single whitespace-separated 3-column file
    (alpha, cl, cd) split into several roughness cases by "zone, t=..." separator
    lines. Only the first zone ("80 grit") is returned, since cd here is tabulated
    directly against alpha. Lines starting with "#" (comments) and the "variables="
    header line are skipped, as are the "zone" separators themselves.

    Parameters
    ----------
    path : str
        Path to the Ladson .dat file.

    Returns
    -------
    dict
        Keys "alpha", "cl", "cd", each a 1D NumPy array (first zone only).

    Raises
    ------
    FileNotFoundError
        If the file is missing.
    """

    if not os.path.exists(path):
        raise FileNotFoundError(f"Ladson experimental file not found: {path}")

    alpha, cl, cd = [], [], []
    zone_count = 0
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("variables="):
                continue
            if line.startswith("zone"):
                zone_count += 1
                # Stop once we leave the first zone ("80 grit").
                if zone_count > 1:
                    break
                continue
            parts = line.split()
            alpha.append(float(parts[0]))
            cl.append(float(parts[1]))
            cd.append(float(parts[2]))

    return {"alpha": np.array(alpha), "cl": np.array(cl), "cd": np.array(cd)}


# ---------------------------------------------------------------------------
# Plotting helper
# ---------------------------------------------------------------------------
def plot_xfoil_vs_experiment(xfoil_polar, exp_data, airfoil_name, reynolds,
                             ladson_data=None, save_path=None):
    """
    Overlay an XFOIL polar against digitised experimental data on a 2x2 grid.

    XFOIL results are drawn as solid lines; experimental points as discrete
    open markers (no connecting line). Abbott drag is only available against cl,
    so it appears on the drag-polar panel alone; the optional Ladson data has cd
    against alpha directly, so it is also shown on the Cd-vs-alpha panel.

    Parameters
    ----------
    xfoil_polar : numpy.ndarray
        Polar array (Nx7) as returned by run_xfoil_polar; columns
        [alpha, CL, CD, ...].
    exp_data : dict
        Abbott experimental data as returned by load_abbott_data, with keys
        "alpha", "cl", "cl_for_cd", "cd".
    airfoil_name : str
        Used in the figure title.
    reynolds : float
        Reynolds number of the datasets (used in the title).
    ladson_data : dict or None
        Optional Ladson experimental data as returned by load_ladson_data, with
        keys "alpha", "cl", "cd". If None, only Abbott is overlaid.
    save_path : str or None
        If provided, the figure is saved to this path. Otherwise just displayed.
    """

    alpha = xfoil_polar[:, 0]
    cl = xfoil_polar[:, 1]
    cd = xfoil_polar[:, 2]

    # Experiment marker styles: open circles (Abbott) and open squares (Ladson),
    # each a distinct colour, with no connecting line.
    exp_style = dict(marker="o", linestyle="none", markerfacecolor="none",
                     markeredgecolor="C1", markersize=6, label="Experiment (Abbott)")
    ladson_style = dict(marker="s", linestyle="none", markerfacecolor="none",
                        markeredgecolor="C2", markersize=6,
                        label="Experiment (Ladson, 80 grit)")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"XFOIL vs. Experiment — {airfoil_name} at Re = {reynolds:,}",
        fontsize=14, fontweight="bold",
    )

    # CL vs alpha
    axes[0, 0].plot(alpha, cl, "-", color="C0", label="XFOIL")
    axes[0, 0].plot(exp_data["alpha"], exp_data["cl"], **exp_style)
    if ladson_data is not None:
        axes[0, 0].plot(ladson_data["alpha"], ladson_data["cl"], **ladson_style)
    axes[0, 0].set_xlabel(r"$\alpha$ [deg]")
    axes[0, 0].set_ylabel(r"$C_L$")
    axes[0, 0].set_title("Lift Coefficient")
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=9)

    # CD vs alpha — XFOIL plus Ladson (which tabulates cd against alpha).
    # Abbott has no direct alpha-cd mapping, so it appears in the drag-polar panel only.
    axes[0, 1].plot(alpha, cd, "-", color="C0", label="XFOIL")
    if ladson_data is not None:
        axes[0, 1].plot(ladson_data["alpha"], ladson_data["cd"], **ladson_style)
    axes[0, 1].set_xlabel(r"$\alpha$ [deg]")
    axes[0, 1].set_ylabel(r"$C_D$")
    axes[0, 1].set_title("Drag Coefficient\n(Abbott $C_D$ shown in drag-polar panel only)")
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=9)

    # Drag polar (CL vs CD) — experiment plotted on its own (cd, cl) axes
    axes[1, 0].plot(cd, cl, "-", color="C0", label="XFOIL")
    axes[1, 0].plot(exp_data["cd"], exp_data["cl_for_cd"], **exp_style)
    if ladson_data is not None:
        axes[1, 0].plot(ladson_data["cd"], ladson_data["cl"], **ladson_style)
    axes[1, 0].set_xlabel(r"$C_D$")
    axes[1, 0].set_ylabel(r"$C_L$")
    axes[1, 0].set_title("Drag Polar")
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=9)

    # CL error: interpolate XFOIL CL onto experimental alpha, over the overlap only.
    # np.interp needs the XFOIL alpha ascending (ASEQ produces it ascending).
    exp_alpha = exp_data["alpha"]
    exp_cl = exp_data["cl"]
    mask = (exp_alpha >= alpha.min()) & (exp_alpha <= alpha.max())
    alpha_overlap = exp_alpha[mask]
    cl_interp = np.interp(alpha_overlap, alpha, cl)
    cl_err = np.abs(cl_interp - exp_cl[mask])

    axes[1, 1].plot(alpha_overlap, cl_err, "o-", color="C3", markersize=4,
                    label=r"$|C_{L,\mathrm{XFOIL}} - C_{L,\mathrm{Abbott}}|$")
    axes[1, 1].set_xlabel(r"$\alpha$ [deg]")
    axes[1, 1].set_ylabel(r"$|\Delta C_L|$")
    axes[1, 1].set_title("Lift Coefficient Error")
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(fontsize=9)

    if cl_err.size:
        print(f"  Mean |CL error| : {cl_err.mean():.4f}")
        print(f"  Max  |CL error| : {cl_err.max():.4f} at alpha = {alpha_overlap[np.argmax(cl_err)]:+.2f} deg")
    else:
        print("  [WARN] No overlap between experimental and XFOIL alpha ranges — no CL error computed.")

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
    print("Validation: NACA 0012 vs. Abbott data at Re = 6,000,000")
    print("=" * 60)

    # Alpha extended to +20 deg to capture the stall region, where the
    # XFOIL-vs-experiment divergence is most informative. n_iter raised to 300 to
    # help near-stall convergence; XFOIL simply omits any alpha that fails to
    # converge, so non-converged points are skipped rather than causing a crash.
    print("  Running NACA 0012 at Re = 6,000,000 ...")
    val_polar = run_xfoil_polar(
        airfoil_cmd="NACA 0012",
        reynolds=6_000_000,
        alpha_min=-6,
        alpha_max=20,
        alpha_step=0.5,
        polar_path="polar_naca0012_re6e6.txt",
        n_iter=300,
    )

    if val_polar is None:
        print("  [WARN] XFOIL produced no polar for NACA 0012 — skipping validation.")
    else:
        print(f"  Converged points: {len(val_polar)} "
              f"(alpha {val_polar[:, 0].min():+.1f} to {val_polar[:, 0].max():+.1f} deg)")

        # Digitised experimental data lives in the repo's top-level data/ folder.
        exp_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
        cl_path = os.path.join(exp_dir, "0012_abbottdata_cl.dat")
        cd_path = os.path.join(exp_dir, "0012_abbottdata_cd.dat")
        ladson_path = os.path.join(exp_dir, "CLCD_Ladson_expdata.dat")

        try:
            exp_data = load_abbott_data(cl_path, cd_path)
        except FileNotFoundError as err:
            print(f"  [WARN] {err}")
            print("  Skipping XFOIL-vs-experiment comparison.")
        else:
            # Ladson is optional: a missing file just drops that series, leaving
            # the Abbott comparison intact.
            try:
                ladson_data = load_ladson_data(ladson_path)
            except FileNotFoundError as err:
                print(f"  [WARN] {err}")
                print("  Continuing with Abbott data only.")
                ladson_data = None

            plot_xfoil_vs_experiment(
                val_polar,
                exp_data,
                airfoil_name="NACA 0012",
                reynolds=6_000_000,
                ladson_data=ladson_data,
                save_path=os.path.join(RESULTS_DIR, "naca0012_xfoil_vs_abbott_re6e6.png"),
            )

    print()
    print("Done. PNG figure saved to the 'results' folder.")
