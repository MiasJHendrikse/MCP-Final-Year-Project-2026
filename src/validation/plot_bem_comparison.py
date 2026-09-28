"""
Combine compare_pybemt.py's and compare_ccblade.py's results into
the two-panel three-way Cp-lambda / Ct-lambda plot and deviation tables for
docs/validation/bem-cross-validation.md.

Third of three scripts in the cross-check pipeline (see
compare_pybemt.py's module docstring) -- this one needs only numpy and
matplotlib, both already used elsewhere in this repo, so it runs with the
plain repo python (no external virtualenv). It runs no solver: its only
inputs are the two committed result files written by the other two scripts,
and its only output is the figure (and the markdown tables on stdout).

    cd <repo root>
    python src/validation/plot_bem_comparison.py

(or, as the doc originally stated, `python -m validation.plot_bem_comparison`
from `src/` -- both work.)

Reads docs/validation/pybemt_case/results.json and
docs/validation/ccblade_case/results.json (written by the other two
scripts -- run those first) and writes
docs/validation/bem_cross_validation_comparison.png plus prints the
markdown tables used in the doc above. The figure is written through
`plotting.figstyle.save` (one style for every committed plot).

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, ".."))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from validation.compare_pybemt import DOCS_VALIDATION_DIR  # noqa: E402


def load_results():
    with open(os.path.join(DOCS_VALIDATION_DIR, "pybemt_case", "results.json")) as f:
        pybemt_rows = json.load(f)
    with open(os.path.join(DOCS_VALIDATION_DIR, "ccblade_case", "results.json")) as f:
        ccblade_rows = json.load(f)

    by_v_inf = {}
    for row in pybemt_rows:
        by_v_inf[row["v_inf"]] = {
            "v_inf": row["v_inf"], "tsr": row["tsr"],
            "Cp_ours": row["Cp_ours"], "Ct_ours": row["Ct_ours"],
            "Cp_pybemt": row["Cp_pybemt"], "Ct_pybemt": row["Ct_pybemt"],
        }
    for row in ccblade_rows:
        entry = by_v_inf.setdefault(row["v_inf"], {"v_inf": row["v_inf"], "tsr": row["tsr"]})
        entry.setdefault("Cp_ours", row["Cp_ours"])
        entry.setdefault("Ct_ours", row["Ct_ours"])
        entry["Cp_ccblade"] = row["Cp_ccblade"]
        entry["Ct_ccblade"] = row["Ct_ccblade"]

    rows = sorted(by_v_inf.values(), key=lambda r: -r["tsr"])
    for r in rows:
        assert "Cp_pybemt" in r and "Cp_ccblade" in r, (
            f"v_inf={r['v_inf']} missing a result -- did both compare_pybemt.py "
            "and compare_ccblade.py run over the same wind-speed sweep?"
        )
    return rows


def make_plot(rows, out_path):
    """The committed two-panel figure: power and thrust coefficient vs TSR."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from plotting import figstyle

    figstyle.apply()

    tsr = [r["tsr"] for r in rows]
    fig, (ax_cp, ax_ct) = plt.subplots(1, 2, figsize=figstyle.DOUBLE)

    for ax, key, title, ylabel in (
        (ax_cp, "Cp", "Power coefficient", figstyle.LABELS["power_coefficient"]),
        (ax_ct, "Ct", "Thrust coefficient", figstyle.LABELS["thrust_coefficient"]),
    ):
        ax.plot(tsr, [r[f"{key}_ours"] for r in rows], "o-", label="present solver")
        ax.plot(tsr, [r[f"{key}_ccblade"] for r in rows], "^-.", label="CCBlade")
        ax.plot(tsr, [r[f"{key}_pybemt"] for r in rows], "s--", label="pyBEMT")
        ax.set_xlabel(figstyle.LABELS["tip_speed_ratio"])
        ax.set_ylabel(ylabel)
        figstyle.title(ax, title)
        ax.legend()

    figstyle.save(fig, out_path)
    plt.close(fig)
    print(f"Wrote plot: {out_path}")


def make_full_table(rows):
    lines = [
        "| Wind speed (m/s) | TSR | Cp (ours) | Cp (pyBEMT) | Cp (CCBlade) | Ct (ours) | Ct (pyBEMT) | Ct (CCBlade) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['v_inf']:.1f} | {r['tsr']:.3f} | {r['Cp_ours']:.4f} | {r['Cp_pybemt']:.4f} | "
            f"{r['Cp_ccblade']:.4f} | {r['Ct_ours']:.4f} | {r['Ct_pybemt']:.4f} | {r['Ct_ccblade']:.4f} |"
        )
    return "\n".join(lines)


def make_deviation_table(rows):
    peak = max(rows, key=lambda r: r["Cp_ours"])
    off_design = [r for r in rows if r["v_inf"] in (10.0, 25.0)]
    highlighted = [peak] + off_design

    lines = [
        "| Wind speed (m/s) | TSR | Cp (ours) | Cp dev. vs pyBEMT | Cp dev. vs CCBlade | "
        "Ct (ours) | Ct dev. vs pyBEMT | Ct dev. vs CCBlade |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in highlighted:
        cp_dev_pybemt = r["Cp_ours"] - r["Cp_pybemt"]
        cp_pct_pybemt = 100.0 * cp_dev_pybemt / r["Cp_pybemt"]
        cp_dev_ccblade = r["Cp_ours"] - r["Cp_ccblade"]
        cp_pct_ccblade = 100.0 * cp_dev_ccblade / r["Cp_ccblade"]
        ct_dev_pybemt = r["Ct_ours"] - r["Ct_pybemt"]
        ct_pct_pybemt = 100.0 * ct_dev_pybemt / r["Ct_pybemt"]
        ct_dev_ccblade = r["Ct_ours"] - r["Ct_ccblade"]
        ct_pct_ccblade = 100.0 * ct_dev_ccblade / r["Ct_ccblade"]
        tag = " (peak $C_p$ in this sweep)" if r is peak else ""
        lines.append(
            f"| {r['v_inf']:.1f}{tag} | {r['tsr']:.3f} | {r['Cp_ours']:.4f} | "
            f"{cp_dev_pybemt:+.4f} ({cp_pct_pybemt:+.1f}%) | {cp_dev_ccblade:+.4f} ({cp_pct_ccblade:+.1f}%) | "
            f"{r['Ct_ours']:.4f} | {ct_dev_pybemt:+.4f} ({ct_pct_pybemt:+.1f}%) | "
            f"{ct_dev_ccblade:+.4f} ({ct_pct_ccblade:+.1f}%) |"
        )
    return "\n".join(lines)


def main():
    rows = load_results()

    plot_path = os.path.join(DOCS_VALIDATION_DIR, "bem_cross_validation_comparison.png")
    make_plot(rows, plot_path)

    print("\n--- Full sweep table (markdown) ---\n")
    print(make_full_table(rows))
    print("\n--- Deviation table: peak Cp + off-design points (markdown) ---\n")
    print(make_deviation_table(rows))

    mean_abs_pct_pybemt = sum(abs(100 * (r["Cp_ours"] - r["Cp_pybemt"]) / r["Cp_pybemt"]) for r in rows) / len(rows)
    mean_abs_pct_ccblade = sum(abs(100 * (r["Cp_ours"] - r["Cp_ccblade"]) / r["Cp_ccblade"]) for r in rows) / len(rows)
    print(f"\nMean |Cp deviation| across the sweep: vs pyBEMT = {mean_abs_pct_pybemt:.2f}%, "
          f"vs CCBlade = {mean_abs_pct_ccblade:.2f}%")


if __name__ == "__main__":
    main()
