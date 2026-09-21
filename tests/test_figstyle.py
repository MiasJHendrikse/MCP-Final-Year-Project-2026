"""
The figure-style gate: `src/plotting/figstyle.py`.

`save()` is the mechanical guarantee behind the report's two figure rules --
at most two panels, side by side, and no project history inside an image --
so the refusal paths are tested one by one: three axes, a super-title, each
forbidden string, and the two forbidden whole words.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
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
    assert matplotlib.rcParams["legend.frameon"] is False
    assert matplotlib.rcParams["axes.spines.top"] is False
    assert matplotlib.rcParams["axes.spines.right"] is False
    assert figstyle.SINGLE == (3.35, 2.6)
    assert figstyle.DOUBLE == (6.3, 2.6)
