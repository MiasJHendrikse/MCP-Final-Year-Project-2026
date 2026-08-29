"""
Task 2 acceptance: the polar caches are gap-free, span the full circle, keep
provenance, and the promoted Viterna module changed nothing.

Each test here corresponds to one of the work order's "done when" clauses for
Task 2. The one clause not covered by a test is `validate_polars.py`'s 5 checks
still giving 4/5 on S809 -- that script prints a report and is run by hand
(`python -m validation.validate_polars s809`); its result is recorded in
config/polars_s809.yaml.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import glob
import os

import numpy as np
import pytest

from config import load_polar_cache
from polars import cache_format, envelope, viterna
from validation import check_stitch_continuity, export_qblade
from xfoil.polar_lookup import PolarLookup

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
POLAR_DIR = os.path.join(REPO_ROOT, "data", "polars")
QBLADE_DIR = os.path.join(REPO_ROOT, "data", "qblade")

#: Caches committed to the repo, as cache-name -> whether Task 2 extended it.
#: NACA 4412 is a retained secondary reference dataset (see
#: build_polar_cache.AIRFOILS); it is deliberately left as built -- extending
#: it would be a rebuild of validated data nothing currently reads.
CACHES = {"s809": True, "naca4412": False}


def _paths(cache):
    return sorted(glob.glob(os.path.join(POLAR_DIR, cache, "*_Re*.csv")))


# ---------------------------------------------------------------------------
# "Neither cache contains a NaN after rectangularisation"
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cache", sorted(CACHES))
def test_no_nan_after_rectangularisation(cache):
    """
    A NaN here is what Task 3 cannot fit a spline over.

    PolarLookup reindexes every Reynolds curve onto the union alpha axis and
    leaves NaN wherever a curve does not reach. That is the rectangularisation
    the acceptance criterion names, so it is checked through the real thing
    rather than by re-deriving it.
    """

    lookup = PolarLookup(os.path.join(POLAR_DIR, cache))
    for label, interp in (("cl", lookup._cl_interp),
                          ("cd", lookup._cd_interp),
                          ("cm", lookup._cm_interp)):
        grid = np.asarray(interp.values)
        n_nan = int(np.count_nonzero(np.isnan(grid)))
        assert n_nan == 0, f"{cache}: {n_nan} NaN in the rectangularised {label} grid"


# ---------------------------------------------------------------------------
# "Both span -180..180 deg"
# ---------------------------------------------------------------------------

def test_s809_spans_full_circle():
    step = load_polar_cache("s809").alpha_step_deg
    expected = np.round(np.arange(-180.0, 180.0 + step / 2, step), 6)

    for path in _paths("s809"):
        table, _source = cache_format.read_polar_csv(path)
        assert np.array_equal(np.round(table[:, 0], 6), expected), (
            f"{os.path.basename(path)} does not cover -180..180 on a "
            f"{step} deg grid"
        )


# ---------------------------------------------------------------------------
# The holes are closed, and what filled them is still visible
# ---------------------------------------------------------------------------

def test_s809_measured_band_is_gap_free():
    """The 8 documented interior holes are gone, at every Reynolds number."""

    config = load_polar_cache("s809")
    n = int(round((config.alpha_max_deg - config.alpha_min_deg)
                  / config.alpha_step_deg)) + 1
    expected = np.round(config.alpha_min_deg
                        + config.alpha_step_deg * np.arange(n), 6)

    for path in _paths("s809"):
        table, _source = cache_format.load_xfoil_band(path)
        missing = np.setdiff1d(expected, np.round(table[:, 0], 6))
        assert not len(missing), (
            f"{os.path.basename(path)} still missing alpha {missing.tolist()}"
        )

    assert config.known_gaps == (), (
        "config/polars_s809.yaml still lists outstanding gaps"
    )


def test_filled_points_stay_distinguishable():
    """
    A filled point must not be mistakable for a converged one.

    The cache carries provenance per row, so this is a property of the data,
    not of the README. Rows that are neither converged XFOIL output nor a
    retry are counted against the record in results/ and the cache README.
    """

    filled = {}
    for path in _paths("s809"):
        table, source = cache_format.load_xfoil_band(path)
        flagged = source != cache_format.SOURCE_XFOIL
        if flagged.any():
            filled[os.path.basename(path)] = [
                (float(a), cache_format.SOURCE_LABELS[int(s)])
                for a, s in zip(table[flagged, 0], source[flagged])
            ]

    # 8 holes: 1 recovered by the XFOIL retry, 7 filled from a local fit.
    total = sum(len(v) for v in filled.values())
    assert total == 8, f"expected 8 flagged rows, found {total}: {filled}"
    labels = [label for rows in filled.values() for _a, label in rows]
    assert labels.count("gap_retry") == 1
    assert labels.count("local_fit") == 7


# ---------------------------------------------------------------------------
# "export_qblade.py produces byte-identical .plr output through the promoted
#  module -- proof the promotion changed no behaviour"
# ---------------------------------------------------------------------------

def _data_lines(path):
    """The .plr's rows, minus the three lines that carry a timestamp."""

    with open(path) as f:
        lines = f.read().splitlines()
    return [line for i, line in enumerate(lines) if i not in (2, 3, 4)]


def test_export_qblade_reproduces_committed_plr(tmp_path):
    """
    Re = 500k has no filled rows, so its `.plr` must come back unchanged.

    Any difference here is the promotion having changed behaviour, which is
    the thing this test exists to rule out.
    """

    export_qblade.export(str(tmp_path), reynolds=500_000, plr_name="S809_Re500k.plr")
    generated = _data_lines(str(tmp_path / "S809_Re500k.plr"))
    committed = _data_lines(os.path.join(QBLADE_DIR, "S809_Re500k.plr"))
    assert generated == committed


def test_export_qblade_re100k_differs_only_where_a_gap_was_filled(tmp_path):
    """
    Re = 100k had one hole, at alpha = 13.0 deg.

    The committed `.plr` bridged it by linear interpolation; the cache now
    carries a documented local fit there instead. Exactly one row may move,
    and it must be that one -- anything else would be the promotion leaking.
    """

    export_qblade.export(str(tmp_path), reynolds=100_000, plr_name="S809.plr")
    generated = _data_lines(str(tmp_path / "S809.plr"))
    committed = _data_lines(os.path.join(QBLADE_DIR, "S809.plr"))

    assert len(generated) == len(committed)
    differing = [g for g, c in zip(generated, committed) if g != c]
    assert len(differing) == 1, f"{len(differing)} rows moved, expected 1"
    assert differing[0].split()[0] == "13.000000"


# ---------------------------------------------------------------------------
# The three reviewed constants
# ---------------------------------------------------------------------------

def test_reversed_amplitude_is_derived_and_reproduces_the_fitted_value():
    """0.7 * 1.8 / 2 = 0.63 exactly -- the fitted constant, now derived."""

    assert viterna.reversed_cl_amplitude() == 0.63


def test_cd_max_finite_blade_is_below_the_fitted_value_for_this_rotor():
    """
    The design rotor's derived CD_MAX, and the S809 default it replaces there.

    Recorded as a test so the two cannot silently converge: the S809 cache and
    the committed QBlade files must keep the fitted 1.8, and a design-rotor
    cache must not.
    """

    aspect_ratio = envelope.blade_aspect_ratio()
    derived = viterna.cd_max_finite_blade(aspect_ratio)
    assert 1.2 < derived < 1.5
    assert viterna.CD_MAX_QBLADE_MATCHED == 1.8
    assert derived < viterna.CD_MAX_QBLADE_MATCHED


# ---------------------------------------------------------------------------
# "Stitch value/slope discontinuities quantified per Re row and recorded"
# ---------------------------------------------------------------------------

def test_stitch_value_continuity_is_exact():
    """
    The extrapolation is anchored at the last converged alpha, so it must
    reproduce that point exactly. The slope jump is the interesting quantity
    and is not asserted -- it is measured and recorded (see the cache README);
    a threshold here would be a tuning knob, which is not what Task 2 asks for.
    """

    rows = check_stitch_continuity.check_cache("s809")
    assert len(rows) == len(load_polar_cache("s809").reynolds)
    for row in rows:
        assert row["value_jump_cl_upper"] == 0.0
        assert row["value_jump_cd_upper"] == 0.0


# ---------------------------------------------------------------------------
# "SG6043 Re bounds demonstrably cover the computed operating envelope with
#  margin; envelope recorded in config"
# ---------------------------------------------------------------------------

def test_sg6043_bounds_cover_the_computed_envelope_with_margin():
    """
    Recomputed from config, not read back from the file that records it.

    Task 4 makes an out-of-range Reynolds lookup raise, so this is the check
    that stands between a design-rotor sweep and a hard crash in the low
    wind-speed bins.
    """

    bounds = load_polar_cache("sg6043").reynolds
    computed = envelope.reynolds_envelope()
    re_min = computed["overall_min"]["reynolds"]
    re_max = computed["overall_max"]["reynolds"]

    assert bounds[0] < re_min, f"cache floor {bounds[0]:,} is above the envelope"
    assert bounds[-1] > re_max, f"cache ceiling {bounds[-1]:,} is below the envelope"
    assert bounds[0] <= 0.9 * re_min, "less than 10 % margin below the envelope"
    assert bounds[-1] >= 1.1 * re_max, "less than 10 % margin above the envelope"
