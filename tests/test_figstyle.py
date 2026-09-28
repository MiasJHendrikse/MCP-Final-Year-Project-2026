"""
The figure-style gate: `src/plotting/figstyle.py`.

`save()` is the mechanical guarantee behind the report's figure rules -- at
most two panels, side by side, no project history inside an image, and no
text or legend on the data -- so the refusal paths are tested one by one:
three axes, a super-title, each forbidden string, the two forbidden whole
words, and a legend or text sitting on a line or on another text.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from plotting import figstyle  # noqa: E402


def _axes(n):
    fig, axes = plt.subplots(1, n, figsize=figstyle.DOUBLE)
    return fig, list(axes)


def test_label_builds_the_vocabulary_strings():
    assert figstyle.label("Radius", "r", "m") == "Radius $r$ [m]"
    assert figstyle.label("Chord", "c", "mm") == "Chord $c$ [mm]"
    assert figstyle.label("Twist", r"\theta", "°") == "Twist $\\theta$ [°]"
    assert (figstyle.label("Tip-speed ratio", r"\lambda", "--")
            == "Tip-speed ratio $\\lambda$ [--]")
    assert (figstyle.label("Wind speed", "V", r"m s$^{-1}$")
            == "Wind speed $V$ [m s$^{-1}$]")
    assert (figstyle.label("Lift coefficient", "C_l", "--")
            == "Lift coefficient $C_l$ [--]")
    assert figstyle.label(r"$1/c(r)^3$", None, r"m$^{-3}$") == "$1/c(r)^3$ [m$^{-3}$]"
    assert (figstyle.label(r"$\|g_j(h) - g_j(h^*_j)\|$", None,
                           r"MWh yr$^{-1}$")
            == r"$\|g_j(h) - g_j(h^*_j)\|$ [MWh yr$^{-1}$]")


def test_labels_dict_carries_every_vocabulary_string():
    assert figstyle.LABELS["radius"] == "Radius $r$ [m]"
    assert figstyle.LABELS["residual"] == "Relative residual $\\|R\\|/R_0$ [--]"
    assert figstyle.LABELS["material_saved"] == "Shell material saved [%]"
    assert (figstyle.LABELS["shell_material"]
            == "Shell material relative to $\\mathbf{x}_0$ [%]")
    assert figstyle.LABELS["flexibility_weight"] == "$1/c(r)^3$ [m$^{-3}$]"


def test_blades_style_names_the_thing_not_its_history():
    assert set(figstyle.BLADES) == {"x0", "x_c", "x_m", "x_star"}
    for entry in figstyle.BLADES.values():
        figstyle.check(_titled(entry["label"]))
    assert figstyle.BLADES["x_m"]["label"] == (
        r"$\mathbf{x}_m$ minimum-material")


def _titled(text):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label=text)
    ax.legend()
    return fig


@pytest.mark.parametrize("text", [
    "resolved 2026-09-19",
    "grounded bounds",
    "Phase 4 optimum",
    "the audit measurement",
    "case A3",
    "case A4",
    "case B1",
    "case B2",
    "case B3",
    "case B5",
    "provisional range",
    "placeholder",
    "scratch run",
    "old bilinear lookup",
    "new C1 interpolant",
    "before/after",
    "fitted Schmitz",
    "chord_max_m = 0.30 m",
    "chord_min_m",
    "twist_min",
    "twist_max",
    "manufacturing.min_chord_m",
    "under the configured bounds",
    "a re-run",
    "journal",
    "that decision",
])
def test_save_refuses_forbidden_text(text):
    fig = _titled(text)
    with pytest.raises(figstyle.FigureStyleError):
        figstyle.save(fig, os.path.join(os.environ.get("TEMP", "/tmp"),
                                        "_figstyle_should_not_exist.png"))
    plt.close(fig)


def test_save_refuses_more_than_two_axes():
    fig, _axes_list = _axes(3)
    with pytest.raises(figstyle.FigureStyleError):
        figstyle.save(fig, os.path.join(os.environ.get("TEMP", "/tmp"),
                                        "_figstyle_3axes.png"))
    plt.close(fig)


def test_twin_axes_count_as_one_panel(tmp_path):
    fig, (ax, _) = plt.subplots(1, 2)
    twin = ax.twinx()
    twin.plot([0, 1], [1, 2])
    ax.plot([0, 1], [0, 1])
    path = tmp_path / "twin.png"
    assert figstyle.save(fig, str(path)) == str(path)
    assert path.exists()
    plt.close(fig)


def test_save_refuses_a_super_title():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    fig.suptitle("a super-title")
    with pytest.raises(figstyle.FigureStyleError):
        figstyle.save(fig, os.path.join(os.environ.get("TEMP", "/tmp"),
                                        "_figstyle_suptitle.png"))
    plt.close(fig)


def test_save_writes_a_clean_two_panel_figure(tmp_path):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=figstyle.DOUBLE)
    ax.plot([0, 1], [0, 1])
    bx.plot([0, 1], [1, 0])
    figstyle.title(ax, "Chord control points")
    figstyle.title(bx, "Twist control points")
    ax.set_xlabel(figstyle.LABELS["radius"])
    path = tmp_path / "clean.png"
    figstyle.save(fig, str(path))
    assert path.exists() and path.stat().st_size > 0
    plt.close(fig)


def test_title_refuses_long_parenthesised_or_precise_text():
    _fig, ax = plt.subplots()
    figstyle.title(ax, "Chord control points")  # three words: fine
    figstyle.title(ax, "SG6043 at 60 mm chord")  # a name's digits are not a number
    with pytest.raises(figstyle.FigureStyleError):
        figstyle.title(ax, "One two three four five six seven eight")
    with pytest.raises(figstyle.FigureStyleError):
        figstyle.title(ax, "Chord (mm)")
    with pytest.raises(figstyle.FigureStyleError):
        figstyle.title(ax, "Error at 12.3456")
    plt.close(_fig)


def test_sci_leaves_log_axes_alone_and_formats_linear_ones():
    import matplotlib.ticker
    _fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_xscale("log")
    figstyle.sci(ax)
    assert isinstance(ax.yaxis.get_major_formatter(),
                      matplotlib.ticker.LogFormatterSciNotation)
    assert isinstance(ax.xaxis.get_major_formatter(),
                      matplotlib.ticker.LogFormatterSciNotation)

    _fig2, ax2 = plt.subplots()
    ax2.set_yscale("linear")
    figstyle.sci(ax2)
    formatter = ax2.yaxis.get_major_formatter()
    assert isinstance(formatter, matplotlib.ticker.ScalarFormatter)
    assert formatter.get_useMathText()
    plt.close("all")


def test_apply_sets_the_stated_parameters():
    figstyle.apply()
    assert matplotlib.rcParams["font.family"] == ["serif"]
    assert matplotlib.rcParams["savefig.dpi"] == 300
    assert matplotlib.rcParams["savefig.bbox"] == "tight"
    # A closed box with inward ticks and framed legends.
    assert matplotlib.rcParams["legend.frameon"] is True
    assert matplotlib.rcParams["axes.spines.top"] is True
    assert matplotlib.rcParams["axes.spines.right"] is True
    assert matplotlib.rcParams["xtick.direction"] == "in"
    assert matplotlib.rcParams["ytick.direction"] == "in"
    assert (matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
            == figstyle.PALETTE)
    assert figstyle.SINGLE == (3.35, 2.6)
    assert figstyle.DOUBLE == (6.3, 2.6)


def test_palette_leads_with_the_blade_colours():
    assert figstyle.PALETTE[0] == figstyle.BLADES["x_m"]["color"]
    assert figstyle.PALETTE[1] == figstyle.BLADES["x_c"]["color"]
    assert figstyle.PALETTE[2] == figstyle.BLADES["x_star"]["color"]


def _diagonal(**legend):
    """A line corner to corner, with a legend placed as asked."""
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    ax.plot([0, 1], [0, 1], label="diagonal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(**legend)
    return fig, ax


def test_overlaps_finds_a_legend_on_a_line():
    fig, _ax = _diagonal(loc="center")
    found = figstyle.overlaps(fig)
    assert found and "legend sits on line 'diagonal'" in found[0]
    with pytest.raises(figstyle.FigureStyleError):
        figstyle.save(fig, os.path.join(os.environ.get("TEMP", "/tmp"),
                                        "_figstyle_overlap.png"))
    plt.close(fig)


def test_overlaps_passes_a_legend_in_a_free_corner(tmp_path):
    fig, _ax = _diagonal(loc="upper left")
    assert figstyle.overlaps(fig) == []
    figstyle.save(fig, str(tmp_path / "clear.png"))
    plt.close(fig)


def test_overlaps_passes_a_legend_below_the_axes():
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    ax.plot([0, 1], [0, 1], label="diagonal")
    figstyle.legend_below(ax, ncol=1)
    assert figstyle.overlaps(fig) == []
    plt.close(fig)


def test_overlaps_finds_text_on_a_line_and_text_outside_the_axes():
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    ax.plot([0, 1], [0, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, 0.5, "on the line", ha="center", va="center")
    ax.text(0.98, 0.05, "past the edge")
    found = " ".join(figstyle.overlaps(fig))
    assert "'on the line' sits on line" in found
    assert "'past the edge' runs outside its axes" in found
    plt.close(fig)


def test_overlaps_finds_two_texts_on_each_other():
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.1, 0.8, "first label")
    ax.text(0.12, 0.8, "second label")
    assert any("overlaps" in item for item in figstyle.overlaps(fig))
    plt.close(fig)


def test_overlap_ok_gid_exempts_an_artist():
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)
    ax.plot([0, 1], [0, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, 0.5, "curve label", ha="center", va="center",
            gid=figstyle.OVERLAP_OK)
    assert figstyle.overlaps(fig) == []
    plt.close(fig)


def test_save_can_skip_the_overlap_check_for_a_schematic(tmp_path):
    fig, _ax = _diagonal(loc="center")
    figstyle.save(fig, str(tmp_path / "schematic.png"), check_overlap=False)
    assert (tmp_path / "schematic.png").exists()
    plt.close(fig)


def test_twin_leaves_one_set_of_ticks_on_the_right_spine():
    fig, ax = plt.subplots()
    other = figstyle.twin(ax)
    assert not ax.yaxis._major_tick_kw.get("tick2On", True)
    assert other.get_position().bounds == ax.get_position().bounds
    plt.close(fig)
