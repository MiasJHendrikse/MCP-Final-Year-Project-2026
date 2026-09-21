"""
Task 3 report figures: the C1 interpolant's value and derivative curves, and
the staircase figure that compares its analytic derivative with the bilinear
lookup it replaced.

Per the work order: "the value plot looks fine with almost any scheme -- the
derivative plot is where a bad interpolant reveals itself."

Everything here reads the committed cache CSVs under `data/polars/` -- no
XFOIL, no BEM run. The default cache is the SG6043 at Re = 200 000 (the
design Reynolds number of the committed rotor); the S809 staircase at
Re = 600 000, which the polar work was first measured on, is kept for the
record under its own name. Both caches can be selected:

    python verification/polar_interpolant/generate_plots.py
    python verification/polar_interpolant/generate_plots.py --airfoil s809 --re 600000

Writes into this directory (through `plotting.figstyle.save`):

    cl_vs_alpha.png             Cl against alpha, one curve per Reynolds number
    dcl_dalpha_vs_alpha.png     the analytic derivative of the C1 interpolant
    staircase_before_after.png  bilinear central-FD derivative vs the C1
                                interpolant's analytic derivative, one zoomed
                                alpha range (the cache and Re selected above)
    staircase_s809_re600k.png   the S809 Re 600 000 staircase, kept for the
                                record (only when the selection above is a
                                different cache or Reynolds number)

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
_SRC = os.path.join(_REPO_ROOT, "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from plotting import figstyle  # noqa: E402
from polars.cache import PolarGrid  # noqa: E402
from polars.interpolant import PolarInterpolant  # noqa: E402
from xfoil.polar_lookup import PolarLookup  # noqa: E402

#: The two committed polar caches this directory's figures are drawn from.
CACHE_DIRS = {
    "s809": os.path.join(_REPO_ROOT, "data", "polars", "s809"),
    "sg6043": os.path.join(_REPO_ROOT, "data", "polars", "sg6043"),
}

DEFAULT_AIRFOIL = "sg6043"
DEFAULT_RE = 200_000

#: A spread of Reynolds numbers across each cache, including both ends.
RE_SAMPLES = {
    "s809": (100_000, 300_000, 600_000, 900_000, 1_300_000),
    "sg6043": (40_000, 100_000, 200_000, 400_000, 800_000),
}

#: The alpha window of the staircase figure: attached flow, where the
#: bilinear lookup's central-FD derivative is a staircase of four or five
#: steps and the C1 interpolant's analytic derivative is one smooth curve.
ZOOM_DEG = (4.0, 5.6)

#: The measurement the staircase figure was first made on, kept under its
#: own name for the record.
S809_RECORD = {"airfoil": "s809", "re": 600_000,
               "filename": "staircase_s809_re600k.png"}


def derivative_label():
    """dCl/dalpha, per degree: the interpolant's alpha axis is in degrees."""

    return figstyle.label(r"$\partial C_l/\partial\alpha$", None, r"deg$^{-1}$")


def build_interpolant(cache_dir):
    grid = PolarGrid(cache_dir)
    return grid, PolarInterpolant(grid)


def plot_cl_vs_alpha(interpolant, re_samples, path):
    """Cl against alpha, one curve per cached Reynolds number."""

    figstyle.apply()
    alphas = np.linspace(-5.0, 20.0, 2000)  # attached flow through deep stall
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    for re in re_samples:
        ax.plot(alphas, [interpolant.cl(a, re) for a in alphas],
                label=rf"$Re$ = {re:,.0f}")
    ax.set_xlabel(figstyle.LABELS["angle_of_attack"])
    ax.set_ylabel(figstyle.LABELS["lift_coefficient"])
    ax.legend()
    figstyle.save(fig, path)
    plt.close(fig)
    print(f"wrote {path}")


def plot_dcl_dalpha(interpolant, re_samples, path):
    """The analytic derivative, fine alpha resolution: no staircase."""

    figstyle.apply()
    alphas = np.linspace(-5.0, 20.0, 2000)
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    for re in re_samples:
        ax.plot(alphas, [interpolant.dcl_dalpha(a, re) for a in alphas],
                label=rf"$Re$ = {re:,.0f}")
    ax.set_xlabel(figstyle.LABELS["angle_of_attack"])
    ax.set_ylabel(derivative_label())
    ax.legend()
    figstyle.save(fig, path)
    plt.close(fig)
    print(f"wrote {path}")


def plot_staircase(interpolant, cache_dir, re, path):
    """The bilinear lookup's central-FD derivative against the C1 one."""

    figstyle.apply()
    lookup = PolarLookup(cache_dir)
    alphas = np.linspace(ZOOM_DEG[0], ZOOM_DEG[1], 400)
    step = 1e-6
    bilinear = [
        (lookup(a + step, re)[0] - lookup(a - step, re)[0]) / (2.0 * step)
        for a in alphas
    ]
    analytic = [interpolant.dcl_dalpha(a, re) for a in alphas]

    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    ax.plot(alphas, bilinear, drawstyle="steps-mid",
            label="bilinear interpolant (central FD)")
    ax.plot(alphas, analytic, label=r"C$^1$ interpolant (analytic)")
    ax.set_xlabel(figstyle.LABELS["angle_of_attack"])
    ax.set_ylabel(derivative_label())
    ax.legend()
    figstyle.save(fig, path)
    plt.close(fig)
    print(f"wrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--airfoil", choices=sorted(CACHE_DIRS),
                        default=DEFAULT_AIRFOIL,
                        help="polar cache to draw from (default: sg6043)")
    parser.add_argument("--re", type=int, default=DEFAULT_RE,
                        help="Reynolds number of the staircase figure")
    args = parser.parse_args(argv)

    cache_dir = CACHE_DIRS[args.airfoil]
    print(f"cache: {os.path.relpath(cache_dir, _REPO_ROOT)}; "
          f"staircase at Re = {args.re:,}")

    _, interpolant = build_interpolant(cache_dir)
    re_samples = RE_SAMPLES[args.airfoil]
    plot_cl_vs_alpha(interpolant, re_samples,
                     os.path.join(_HERE, "cl_vs_alpha.png"))
    plot_dcl_dalpha(interpolant, re_samples,
                    os.path.join(_HERE, "dcl_dalpha_vs_alpha.png"))
    plot_staircase(interpolant, cache_dir, args.re,
                   os.path.join(_HERE, "staircase_before_after.png"))

    if (args.airfoil, args.re) != (S809_RECORD["airfoil"], S809_RECORD["re"]):
        record_dir = CACHE_DIRS[S809_RECORD["airfoil"]]
        _, record_interpolant = build_interpolant(record_dir)
        plot_staircase(record_interpolant, record_dir, S809_RECORD["re"],
                       os.path.join(_HERE, S809_RECORD["filename"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
