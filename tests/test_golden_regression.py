"""
Golden-file regression: the safety net under every change to the solver.

This test recomputes the snapshot defined in `golden_reference.build_golden()`
and asserts it still matches the committed files in `tests/golden/`. It runs
after every remediation task. A task is not complete while it is failing,
unless the change is intended -- in which case the golden files are
regenerated in their own commit (see `tests/generate_golden.py`).

What it covers: the full Cp(lambda) and Ct(lambda) arrays for the NREL Phase
VI rotor across the cross-validated lambda range; per-station a, a', phi,
alpha, Re, F, Cl, Cd and local loads dT/dr, dQ/dr at four operating points
spanning attached flow through deep stall; and our solver's side of all three
cross-tool comparisons. See `golden_reference` for why each piece is included.


Tolerance, and why it is this number
=====================================

    RTOL = 1e-10, with a per-field ATOL of 1e-12 * max|golden| for that field.

Justification, from the bottom up.

*The solver is deterministic today.* The audit measured `solve_rotor` repeated
in-process as bitwise identical, and `test_snapshot_is_bitwise_reproducible`
below re-asserts that directly. So on unmodified code any tolerance at all
would pass, and the tolerance is doing a different job: distinguishing
"nothing changed" from "something changed" across a six-task refactor.

*The floor is set by the root-finder.* `station.solve_station` calls
`brentq(..., xtol=1e-12, rtol=1e-12)`. Those are tolerances on the root
location phi, so phi is only pinned to a relative 1e-12. Cp, Ct and the
spanwise quantities have O(1) sensitivity to phi, so a change that is pure
floating-point reassociation upstream -- reordering a sum, threading a
constant through a config object instead of a module global, moving an
identical value across a function boundary -- can still move the recorded
numbers by order 1e-12 without any behaviour having changed. A tolerance at
or below 1e-12 would therefore produce false failures on pure rewiring
changes that should be numerically inert.

*The ceiling is set by what has to be caught.* The smallest change that
matters physically is far larger. A 0.01 % shift in Cp is 1e-4 relative --
eight orders of magnitude above this tolerance. Every real change the
remaining tasks make (a C1 interpolant replacing bilinear, removal of the
alpha/Re clamps, Ning's residual form replacing the multiplied form) moves
these numbers by 1e-3 or more.

RTOL = 1e-10 sits two orders above the root-finder's own position tolerance
and eight below the smallest physically meaningful change. That gap is wide
enough that the choice is not delicate: anything from 1e-11 to 1e-8 would
behave the same way. 1e-10 is the round number in the middle of it.

*Why the per-field ATOL is scale-relative.* Several recorded quantities pass
through or sit near zero -- a' at the outermost stations, alpha at the near-
zero-twist stations, dQ/dr where the tangential load changes sign. A pure
relative tolerance is meaningless there and a fixed absolute one is either
too tight for dT/dr (order 1e3 N/m) or far too loose for a' (order 1e-3). So
each field gets ATOL = 1e-12 * max|golden value in that field|, which is the
same 1e-12 floor expressed against the field's own scale.

The QBlade-pinned block is the one place a tighter check is available and
used: `test_qblade_pinned_matches_committed_csv` compares against the
committed data/qblade/*.csv bitwise, because those files reproduce exactly
from today's tree.


A note on the cross-tool deviation figures
===========================================

`cross_tool_summary.json` records mean |Cp| deviations of 13.32 % (CCBlade)
and 15.79 % (pyBEMT) against the stored external values -- not the 0.51 % and
1.83 % published in docs/validation/bem-cross-validation.md. This is not a
regression. The stored results.json files predate both the Ncrit=5 S809 cache
rebuild and the extension of the Reynolds bucket list past 500k, so their
external columns come from polar tables the repo no longer contains. The
figures are pinned here as *tracking* values so a refactor cannot move them
unnoticed; they are not asserted to be validation results. See
tests/golden/README.md.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import csv
import math
import os

import pytest

import golden_reference as gr

#: See the tolerance discussion in the module docstring.
RTOL = 1e-10
ATOL_SCALE = 1e-12


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

def _flatten(value, prefix=""):
    """
    Flatten a nested dict/list of numbers into {dotted.path: number}.

    Only leaves reachable from the *golden* structure are produced, so a task
    that adds a new field to a station dict (such as `converged`,
    `residual`, `iterations`) does not fail this test for that reason alone.
    Strings and booleans are skipped -- nothing in the snapshot uses them as a
    numeric quantity.
    """

    out = {}
    if isinstance(value, dict):
        for key, sub in value.items():
            out.update(_flatten(sub, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(value, list):
        for i, sub in enumerate(value):
            out.update(_flatten(sub, f"{prefix}[{i}]"))
    elif isinstance(value, bool) or value is None or isinstance(value, str):
        pass
    elif isinstance(value, (int, float)):
        out[prefix] = float(value)
    return out


def _field_of(path):
    """
    The field name a flattened path belongs to, for the per-field ATOL scale.

    "phase_vi_spanwise.points[2].a_prime[7]" -> "points[].a_prime"
    -- i.e. list indices dropped, so every element of one recorded array
    shares one scale.
    """

    parts = []
    for chunk in path.split("."):
        parts.append(chunk[: chunk.index("[")] + "[]" if "[" in chunk else chunk)
    return ".".join(parts)


def _assert_matches(actual, golden, label):
    """Compare a recomputed structure against its golden counterpart."""

    flat_golden = _flatten(golden)
    flat_actual = _flatten(actual)

    missing = sorted(set(flat_golden) - set(flat_actual))
    assert not missing, (
        f"{label}: {len(missing)} value(s) present in the golden file but not "
        f"recomputed, first few: {missing[:5]}. The snapshot's shape changed; "
        f"if that is intended, regenerate the golden files in their own commit."
    )

    # Per-field absolute floor, expressed against that field's own scale.
    scale = {}
    for path, value in flat_golden.items():
        field = _field_of(path)
        scale[field] = max(scale.get(field, 0.0), abs(value))

    failures = []
    for path, expected in flat_golden.items():
        got = flat_actual[path]
        atol = ATOL_SCALE * scale[_field_of(path)]
        if not math.isclose(got, expected, rel_tol=RTOL, abs_tol=atol):
            delta = got - expected
            rel = abs(delta) / abs(expected) if expected else float("inf")
            failures.append(f"  {path}: golden {expected!r} -> got {got!r} "
                            f"(abs {delta:+.3e}, rel {rel:.3e})")

    assert not failures, (
        f"{label}: {len(failures)} of {len(flat_golden)} recorded values moved "
        f"beyond rtol={RTOL:g} / per-field atol={ATOL_SCALE:g}*scale.\n"
        + "\n".join(failures[:25])
        + (f"\n  ... and {len(failures) - 25} more" if len(failures) > 25 else "")
        + "\n\nIf this change is intended, regenerate the golden files with "
          "`python tests/generate_golden.py` in a SEPARATE commit from the "
          "code change, stating the reason. If it is not intended, this is "
          "the bug the current task introduced."
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def snapshot():
    """
    Recompute the whole snapshot once for the module.

    Fast now that one operating point costs about 15 ms (it used to take
    about three minutes at ~9 s per point).
    """

    return gr.build_golden()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_cp_lambda_matches_golden(snapshot):
    """Rotor-level Cp(lambda) and Ct(lambda) across the cross-validated range."""

    _assert_matches(snapshot["phase_vi_cp_lambda"],
                    gr.load_golden("phase_vi_cp_lambda"),
                    "Phase VI Cp(lambda)")


def test_spanwise_matches_golden(snapshot):
    """
    Per-station state and local loads at four points spanning the envelope.

    This is the check that catches what a mean Cp hides: a refactor can shift
    the spanwise a, a' and phi distributions materially while barely moving
    the integrated coefficients.
    """

    _assert_matches(snapshot["phase_vi_spanwise"],
                    gr.load_golden("phase_vi_spanwise"),
                    "Phase VI spanwise")


def test_cross_tool_matches_golden(snapshot):
    """
    Our side of the three cross-tool comparisons.

    Includes the two deviation figures against the stored external values --
    pinned as tracking numbers, not asserted as validation results. See the
    module docstring.
    """

    _assert_matches(snapshot["cross_tool_summary"],
                    gr.load_golden("cross_tool_summary"),
                    "Cross-tool summary")


def test_snapshot_is_bitwise_reproducible():
    """
    The rotor-level answer is bitwise identical when recomputed in-process.

    This held before the solver was restructured and it must stay true: the
    AEP determinism requirement (`AEP(d)` bitwise reproducible) rests
    on the solver underneath it being bitwise reproducible first. Kept cheap
    -- two operating points, not the whole snapshot.
    """

    from bem.rotor import phase_vi_geometry, solve_rotor

    geometry = phase_vi_geometry()
    for tsr in (3.0, 6.5):
        first = solve_rotor(geometry, tsr=tsr, v_inf=7.0,
                            air_density=gr.AIR_DENSITY,
                            kinematic_viscosity=gr.KINEMATIC_VISCOSITY)
        second = solve_rotor(geometry, tsr=tsr, v_inf=7.0,
                             air_density=gr.AIR_DENSITY,
                             kinematic_viscosity=gr.KINEMATIC_VISCOSITY)
        assert first["Cp"] == second["Cp"], f"Cp not bitwise stable at tsr={tsr}"
        assert first["Ct"] == second["Ct"], f"Ct not bitwise stable at tsr={tsr}"
        for a, b in zip(first["stations"], second["stations"]):
            for field in ("phi", "a", "a_prime", "Cl", "Cd"):
                assert a[field] == b[field], (
                    f"station r={a['r']} field {field} not bitwise stable at tsr={tsr}"
                )


def test_qblade_pinned_matches_committed_csv(snapshot):
    """
    The QBlade-side per-station results still reproduce the committed CSVs.

    Unlike the CCBlade and pyBEMT reference files, `data/qblade/
    bem_solver_phase_vi_result_Re*.csv` postdate the Ncrit=5 cache rebuild and
    reproduce bitwise from today's tree -- so this one is asserted exactly,
    not to a tolerance. It is a second, independent anchor on the same solver
    path, held in `data/` rather than in `tests/golden/`.
    """

    deg = 180.0 / math.pi
    for entry in snapshot["cross_tool_summary"]["qblade_pinned"]:
        reynolds = entry["fixed_reynolds"]
        path = os.path.join(gr.QBLADE_DIR,
                            f"bem_solver_phase_vi_result_Re{reynolds}.csv")
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f))

        assert len(rows) == len(entry["r"]), (
            f"station count changed for the Re={reynolds:,} QBlade pair"
        )
        for i, row in enumerate(rows):
            for csv_key, value in (
                ("r_m", entry["r"][i]),
                ("phi_deg", entry["phi"][i] * deg),
                ("alpha_deg", entry["alpha"][i] * deg),
                ("a", entry["a"][i]),
                ("a_prime", entry["a_prime"][i]),
                ("Cl", entry["Cl"][i]),
                ("Cd", entry["Cd"][i]),
                ("F", entry["F"][i]),
            ):
                assert float(row[csv_key]) == value, (
                    f"Re={reynolds:,} station {i} ({csv_key}): committed CSV has "
                    f"{row[csv_key]}, solver now gives {value!r}. If intended, "
                    f"re-run `python -m validation.compare_qblade --reynolds "
                    f"{reynolds}` from src/ in its own commit."
                )
