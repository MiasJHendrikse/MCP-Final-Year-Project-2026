"""
Task 3 report figures: Cl and dCl/dalpha vs alpha, several Reynolds numbers.

Per the work order: "the value plot looks fine with almost any scheme -- the
derivative plot is where a bad interpolant reveals itself." Run from
`src/` (matches the rest of the project's scripts):

    python ../verification/polar_interpolant/generate_plots.py

Writes two PNGs into this directory: `cl_vs_alpha.png` (value) and
`dcl_dalpha_vs_alpha.png` (analytic derivative, fine alpha resolution). The
second is the actual test -- no staircase, no sawtooth, one smooth curve per
Reynolds number.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from polars.cache import PolarGrid
from polars.interpolant import PolarInterpolant

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(HERE, "..", "..", "data", "polars", "s809")


def main():
    grid = PolarGrid(CACHE_DIR)
    interp = PolarInterpolant(grid)

    # A spread of Reynolds numbers across the cache, including its endpoints.
    re_samples = [100_000, 300_000, 600_000, 900_000, 1_300_000]
    alphas = np.linspace(-5.0, 20.0, 2000)  # attached flow through deep stall

    fig, ax = plt.subplots(figsize=(8, 5))
    for re in re_samples:
        cl = [interp.cl(a, re) for a in alphas]
        ax.plot(alphas, cl, label=f"Re={re:,.0f}")
    ax.set_xlabel("alpha (deg)")
    ax.set_ylabel("Cl")
    ax.set_title("S809: Cl vs alpha (PolarInterpolant)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "cl_vs_alpha.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    for re in re_samples:
        dcl = [interp.dcl_dalpha(a, re) for a in alphas]
        ax.plot(alphas, dcl, label=f"Re={re:,.0f}")
    ax.set_xlabel("alpha (deg)")
    ax.set_ylabel("dCl/dalpha (1/deg)")
    ax.set_title("S809: dCl/dalpha vs alpha (analytic, PolarInterpolant)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "dcl_dalpha_vs_alpha.png"), dpi=150)
    plt.close(fig)

    # Zoomed comparison against the audit's specific measurement: dCl/dalpha
    # across alpha in [4.0, 5.6] deg at Re=600k, old bilinear lookup vs new.
    from xfoil.polar_lookup import PolarLookup

    old = PolarLookup(CACHE_DIR)
    alphas_zoom = np.linspace(4.0, 5.6, 400)
    h = 1e-6
    old_deriv = [
        (old(a + h, 600_000)[0] - old(a - h, 600_000)[0]) / (2 * h)
        for a in alphas_zoom
    ]
    new_deriv = [interp.dcl_dalpha(a, 600_000) for a in alphas_zoom]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(alphas_zoom, old_deriv, label="old bilinear lookup (central FD)",
            drawstyle="steps-mid")
    ax.plot(alphas_zoom, new_deriv, label="new C1 interpolant (analytic)")
    ax.set_xlabel("alpha (deg)")
    ax.set_ylabel("dCl/dalpha (1/deg)")
    ax.set_title("S809 @ Re=600k: the audit's staircase measurement, before/after")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "staircase_before_after.png"), dpi=150)
    plt.close(fig)

    print("Wrote cl_vs_alpha.png, dcl_dalpha_vs_alpha.png, "
          "staircase_before_after.png to", HERE)


if __name__ == "__main__":
    main()
