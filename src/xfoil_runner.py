"""
XFOIL automation script.

Runs XFOIL as a subprocess, sweeps angle of attack at a given Reynolds number,
and returns the polar as a NumPy array.

Author: MJ Hendrikse
Project: DSP810S — Inverse Design of Small Wind Turbine Blades
"""

import os
import subprocess
import numpy as np


# XFOIL is installed outside the repo; resolve it absolutely so the wrapper
# works no matter what the caller's working directory is.
_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_XFOIL = r"C:\Users\miash\Documents\XFOIL6.99\xfoil.exe"

# All generated output (polars, figures) is collected in the repo's top-level
# results/ folder, one level up from this script.
RESULTS_DIR = os.path.abspath(os.path.join(_HERE, "..", "results"))


def run_xfoil_polar(
    airfoil_cmd,
    reynolds,
    alpha_min,
    alpha_max,
    alpha_step,
    polar_path="polar.txt",
    n_iter=100,
    mach=0.0,
    timeout=120,
    xfoil_executable=None,
):
    """
    Run XFOIL for a single airfoil at a single Reynolds number across an alpha sweep.

    Parameters
    ----------
    airfoil_cmd : str
        XFOIL command to load the airfoil. Either:
          - "NACA 2412"       (built-in NACA 4-digit)
          - "NACA 23012"      (built-in NACA 5-digit)
          - "LOAD myfoil.dat" (custom coordinate file)
    reynolds : float
        Reynolds number for viscous analysis.
    alpha_min, alpha_max, alpha_step : float
        Alpha sweep bounds and step, in degrees.
    polar_path : str, optional
        Where to write the polar file. A bare filename or relative path is placed
        inside the top-level "results" folder; an absolute path is used as given.
        Will be deleted first if it already exists.
    n_iter : int, optional
        Maximum viscous-solver iterations per alpha (default 100).
    mach : float, optional
        Mach number (default 0 for incompressible).
    timeout : float, optional
        Maximum wall-clock seconds before XFOIL is killed (default 120).
    xfoil_executable : str, optional
        Path to the XFOIL executable. Defaults to the XFOIL 6.99 install in
        Documents. Pass a full path here to use an XFOIL installed elsewhere.

    Returns
    -------
    numpy.ndarray or None
        Polar data with columns [alpha, CL, CD, CDp, CM, Top_Xtr, Bot_Xtr].
        Returns None if XFOIL timed out or failed to produce a polar file.
    """

    if xfoil_executable is None:
        xfoil_executable = _DEFAULT_XFOIL

    # A bare filename (or any relative path) is collected into RESULTS_DIR so all
    # output lands in one place. An absolute path is honoured as given.
    if not os.path.isabs(polar_path):
        polar_path = os.path.join(RESULTS_DIR, polar_path)
    os.makedirs(os.path.dirname(polar_path), exist_ok=True)

    # XFOIL refuses to overwrite an existing polar file — it would open a
    # confirmation prompt that we are not answering, and the process would hang.
    if os.path.exists(polar_path):
        os.remove(polar_path)

    # XFOIL 6.99 mangles file paths containing spaces or long absolute paths
    # (e.g. "...\XFOIL Results\polar.txt" gets truncated). So we run XFOIL with
    # its working directory set to the output folder and hand it only the bare
    # filename, which it always parses correctly.
    polar_dir = os.path.dirname(polar_path)
    polar_name = os.path.basename(polar_path)

    # The blank lines below are intentional:
    #   - After "G" : exits the PLOP submenu
    #   - After polar_name : skips the optional dump-filename prompt
    commands = f"""PLOP
G

{airfoil_cmd}
PANE
OPER
VISC {reynolds}
MACH {mach}
ITER {n_iter}
PACC
{polar_name}

ASEQ {alpha_min} {alpha_max} {alpha_step}
PACC

QUIT
"""

    try:
        subprocess.run(
            [xfoil_executable],
            input=commands,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=polar_dir,
        )
    except subprocess.TimeoutExpired:
        print(f"  [WARN] XFOIL timed out after {timeout}s for {airfoil_cmd} at Re={reynolds}")
        return None
    except FileNotFoundError:
        raise RuntimeError(
            f"Could not find XFOIL executable '{xfoil_executable}'. "
            "Either add XFOIL to your PATH or pass the full path via xfoil_executable=..."
        )

    if not os.path.exists(polar_path):
        print(f"  [WARN] XFOIL produced no polar file for {airfoil_cmd} at Re={reynolds}")
        return None

    # XFOIL polar files have a 12-line header before the numerical table.
    try:
        data = np.loadtxt(polar_path, skiprows=12)
    except (ValueError, StopIteration):
        print(f"  [WARN] Polar file empty or malformed for {airfoil_cmd} at Re={reynolds}")
        return None

    # If only one alpha converged, np.loadtxt returns a 1D array — promote to 2D
    # so the column-indexing convention is consistent for the caller.
    if data.ndim == 1:
        data = data.reshape(1, -1)

    return data


def summarise_polar(polar, label=""):
    """Print a quick readable summary of a polar array."""
    if polar is None or len(polar) == 0:
        print(f"{label}: no converged points")
        return

    alpha = polar[:, 0]
    cl = polar[:, 1]
    cd = polar[:, 2]
    ld = cl / cd  # Lift-to-drag ratio

    i_best = int(np.argmax(ld))

    print(f"{label}")
    print(f"  Converged points  : {len(polar)}")
    print(f"  Alpha range       : {alpha.min():+.2f} to {alpha.max():+.2f} deg")
    print(f"  CL range          : {cl.min():+.3f} to {cl.max():+.3f}")
    print(f"  Best L/D          : {ld[i_best]:.2f} at alpha = {alpha[i_best]:+.2f} deg")
    print(f"                      (CL = {cl[i_best]:.3f}, CD = {cd[i_best]:.4f})")


# Airfoils offered at the interactive prompt. Each value is the exact command
# string XFOIL uses to generate that airfoil, so the menu selection is passed
# straight through to XFOIL untouched.
AIRFOIL_CHOICES = {
    "1": "NACA 2412",  # cambered 4-digit
    "2": "NACA 0012",  # symmetric
    "3": "NACA 4412",  # higher camber
}


def prompt_airfoil_choice():
    """Prompt the user to pick one of three NACA airfoils.

    Returns the XFOIL airfoil command (e.g. "NACA 2412") for the chosen option.
    """
    print("Select a NACA airfoil:")
    for key, cmd in AIRFOIL_CHOICES.items():
        print(f"  {key}) {cmd}")

    while True:
        choice = input("Enter choice [1-3]: ").strip()
        if choice in AIRFOIL_CHOICES:
            airfoil_cmd = AIRFOIL_CHOICES[choice]
            print(f"Selected: {airfoil_cmd}")
            return airfoil_cmd
        print("  Invalid choice — please enter 1, 2, or 3.")


# ----------------------------------------------------------------------------
# Run when this file is executed directly
# ----------------------------------------------------------------------------
if __name__ == "__main__":

    airfoil_cmd = prompt_airfoil_choice()

    reynolds = 200_000

    # Build a filesystem-safe filename from the selection, e.g.
    # "NACA 2412" -> "polar_naca2412_re200k.txt"
    fname = f"polar_{airfoil_cmd.replace(' ', '').lower()}_re200k.txt"

    print()
    print("=" * 60)
    print(f"Running {airfoil_cmd} at Re = {reynolds:,}")
    print("=" * 60)

    polar = run_xfoil_polar(
        airfoil_cmd=airfoil_cmd,
        reynolds=reynolds,
        alpha_min=-4,
        alpha_max=14,
        alpha_step=0.5,
        polar_path=fname,
        n_iter=200,
    )
    summarise_polar(polar, label=f"{airfoil_cmd} @ Re={reynolds:,}")

    print()
    print(f"Done. Polar saved to: {os.path.join('results', fname)}")