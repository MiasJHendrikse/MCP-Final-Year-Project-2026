"""
Parser for the UIUC LSAT Volume 3 SG6043 clean-model files.

data/sg6043_uiuc_lsat/SG6043_lift_clean.txt and SG6043_drag_clean.txt are
committed byte-for-byte from the UIUC distribution (see that directory's
README). Their block format repeats "Average Reynolds #:" / "Number of
angles of attack:" / a header line / the alpha rows, once per run, with no
delimiter between runs other than a "Tabulated from ..." provenance line —
the outer "::::::::::::::" markers appear once per *file*, not once per run,
so a parser has to scan for the repeating markers rather than the file
header.

The drag file's 9 runs are not 9 Reynolds numbers: several are short
continuation runs at a nominal Re a fuller run already covers (see the
directory README). This module buckets every run to its nearest nominal
Re in {100k, 150k, 200k, 300k, 400k, 500k} and merges same-bucket runs,
because comparing against a two-point continuation run in isolation would
misread it as a curve.

Author: MJ Hendrikse
Project: MCP820S — Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
UIUC_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "data", "sg6043_uiuc_lsat"))

NOMINAL_RE = [100_000, 150_000, 200_000, 300_000, 400_000, 500_000]


def _nearest_nominal_re(re_value):
    return min(NOMINAL_RE, key=lambda r: abs(r - re_value))


def _parse_runs(path):
    """Parse every "Average Reynolds # -> data table" run in a UIUC file.

    Returns a list of (measured_re, rows) tuples, rows a 2D array whose
    columns are exactly as printed (alpha, Cl[, Cm] or alpha, Cl, Cd, ...).
    """
    with open(path) as f:
        lines = f.readlines()

    n = len(lines)
    runs = []
    i = 0
    while i < n:
        if lines[i].strip().startswith("Average Reynolds"):
            re_value = float(lines[i + 1].strip())
            j = i + 2
            while j < n and not lines[j].strip().startswith("Number of angles"):
                j += 1
            n_alpha = int(lines[j + 1].strip())
            j += 2
            while j < n and not lines[j].strip().lower().startswith("alpha"):
                j += 1
            j += 1  # step past the "alpha / Cl / ..." header line
            rows = np.array([
                [float(x) for x in lines[k].split()] for k in range(j, j + n_alpha)
            ])
            runs.append((re_value, rows))
            i = j + n_alpha
        else:
            i += 1
    return runs


def _bucket_and_merge(runs, cols):
    """Group runs by nearest nominal Re and merge, deduplicating shared alphas.

    cols selects which columns to keep from each run's row (e.g. [0, 1] for
    alpha, Cl). Duplicate alphas within a bucket (rounded to 0.01 deg, as can
    happen where a continuation run overlaps the run it extends) are averaged.
    Returns {nominal_re: array sorted by alpha}.
    """
    buckets = {re_nom: [] for re_nom in NOMINAL_RE}
    for measured_re, rows in runs:
        buckets[_nearest_nominal_re(measured_re)].append(rows[:, cols])

    merged = {}
    for re_nom, parts in buckets.items():
        if not parts:
            continue
        stacked = np.vstack(parts)
        alpha_key = np.round(stacked[:, 0], 2)
        unique_alphas = np.unique(alpha_key)
        out = np.array([
            [a] + list(stacked[alpha_key == a, 1:].mean(axis=0))
            for a in unique_alphas
        ])
        merged[re_nom] = out[np.argsort(out[:, 0])]
    return merged


def load_lift_clean():
    """Clean-model lift curves, columns [alpha, cl, cm], one array per nominal Re.

    One measured run per nominal Re already (no continuation runs in this
    file), so the bucketing is an identity op here — kept for a single
    consistent code path with load_drag_clean.
    """
    runs = _parse_runs(os.path.join(UIUC_DIR, "SG6043_lift_clean.txt"))
    return _bucket_and_merge(runs, cols=[0, 1, 2])


def load_drag_clean():
    """Clean-model drag polars, columns [alpha, cl, cd], one array per nominal Re.

    Drops the spanwise Cd replicate columns the file also carries — only the
    first (representative) Cd column is used, matching FORMAT03.TXT.
    """
    runs = _parse_runs(os.path.join(UIUC_DIR, "SG6043_drag_clean.txt"))
    return _bucket_and_merge(runs, cols=[0, 1, 2])


if __name__ == "__main__":
    lift = load_lift_clean()
    drag = load_drag_clean()
    print("Lift (clean):")
    for re_nom, arr in sorted(lift.items()):
        print(f"  Re={re_nom:>7,}  n={len(arr):3d}  "
              f"alpha=[{arr[:,0].min():+.2f}, {arr[:,0].max():+.2f}]  "
              f"cl_max={arr[:,1].max():.3f}")
    print("Drag (clean):")
    for re_nom, arr in sorted(drag.items()):
        print(f"  Re={re_nom:>7,}  n={len(arr):3d}  "
              f"alpha=[{arr[:,0].min():+.2f}, {arr[:,0].max():+.2f}]  "
              f"cd_min={arr[:,2].min():.4f}")
