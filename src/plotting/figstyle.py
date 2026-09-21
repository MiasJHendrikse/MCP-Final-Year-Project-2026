"""
One figure style for every committed plot (the figure-style module).

Every committed figure of this project is written through `save()`. It is the
mechanical guarantee behind two rules the report states: **a figure has at
most two panels, side by side**, and **no figure carries project history**
(no dates, no phase numbers, no "old"/"new", no bounds labels). A figure that
breaks either rule cannot be written -- the call raises instead.

What the module provides
------------------------

* `apply()`      -- sets `matplotlib.rcParams` once: serif text, 9 pt base,
                    no top/right spines, a light grid, 300 dpi PNG output.
* `SINGLE`, `DOUBLE` -- the two figure widths (inches). The A4 text block is
                    160 mm wide; a two-panel figure fills it, a single panel
                    takes a little over half.
* `BLADES`       -- the one style and legend wording for the three committed
                    blades, so they look the same in every chapter.
* `label()`      -- builds every axis label: sentence case, symbol in
                    mathtext, unit in square brackets.
* `LABELS`       -- the same labels under short names, built through
                    `label()`.
* `title()`      -- a short panel title, used only where a two-panel figure
                    needs its panels distinguished.
* `sci()`        -- math-text tick formatting, so exponents read 10^-6 and
                    not 1e-06.
* `save()`       -- asserts the panel count, the absence of a super-title and
                    the forbidden-string list, then writes the PNG.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import os
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text  # noqa: E402
import matplotlib.ticker  # noqa: E402

# ---------------------------------------------------------------------------
# size and style
# ---------------------------------------------------------------------------

#: A single panel: 85 mm, just over half the 160 mm A4 text block.
SINGLE = (3.35, 2.6)
#: Two panels side by side, filling the text block.
DOUBLE = (6.3, 2.6)

DPI = 300

SERIF = ["Palatino Linotype", "TeX Gyre Pagella", "DejaVu Serif"]


def apply():
    """Set the project's `matplotlib.rcParams`. Idempotent."""

    matplotlib.rcParams.update({
        "font.family": "serif",
        "font.serif": list(SERIF),
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
        "lines.linewidth": 1.4,
        "lines.markersize": 4,
        "figure.figsize": SINGLE,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


# ---------------------------------------------------------------------------
# the three committed blades, and the one permitted fourth
# ---------------------------------------------------------------------------

#: Style and legend wording for the blades. Every figure that draws them uses
#: these entries, so `x0`, `x_c` and `x_m` read the same in Chapters 6, 8, 9.
BLADES = {
    "x0": dict(color="0.35", ls="-", label=r"$\mathbf{x}_0$ Schmitz reference"),
    "x_c": dict(color="#b3452a", ls="--", label=r"$\mathbf{x}_c$ energy optimum"),
    "x_m": dict(color="#1f3b8b", ls="-", lw=2.0, label=r"$\mathbf{x}_m$ minimum-material"),
    "x_star": dict(color="#2a7f62", ls=":", label=r"$\mathbf{x}^*$ unconstrained optimum"),
}


# ---------------------------------------------------------------------------
# axis labels
# ---------------------------------------------------------------------------

def label(quantity, symbol=None, unit=None):
    """
    Build an axis label: sentence case, symbol in mathtext, unit in brackets.
    `label("Radius", "r", "m")` returns `"Radius $r$ [m]"`;
    `label("Tip-speed ratio", r"\\lambda", "--")` returns
    `"Tip-speed ratio $\\lambda$ [--]"`.

    `[--]` is the unit of a dimensionless quantity. A label whose whole text
    is a formula (the gradient-deviation axis, for instance) passes
    `symbol=None` and the formula as `quantity`.
    """

    text = quantity if symbol is None else f"{quantity} ${symbol}$"
    if unit is not None:
        text = f"{text} [{unit}]"
    return text


#: The label vocabulary, exact strings, built through `label()`.
LABELS = {
    "radius": label("Radius", "r", "m"),
    "chord": label("Chord", "c", "mm"),
    "twist": label("Twist", r"\theta", "°"),
    "tip_speed_ratio": label("Tip-speed ratio", r"\lambda", "--"),
    "wind_speed": label("Wind speed", "V", r"m s$^{-1}$"),
    "rotor_speed": label("Rotor speed", r"\Omega", "rpm"),
    "angle_of_attack": label("Angle of attack", r"\alpha", "°"),
    "lift_coefficient": label("Lift coefficient", "C_l", "--"),
    "drag_coefficient": label("Drag coefficient", "C_d", "--"),
    "power_coefficient": label("Power coefficient", "C_P", "--"),
    "thrust_coefficient": label("Thrust coefficient", "C_T", "--"),
    "fd_step": label("Step", "h", "--"),
    "gradient_deviation": label(r"$\|g_j(h) - g_j(h^*_j)\|$", None,
                                r"MWh yr$^{-1}$"),
    "residual": label("Relative residual", r"\|R\|/R_0", "--"),
    "evaluations": label("Residual evaluations", None, "--"),
    "design_variables": label("Design variables", "n", "--"),
    "wall_time": label("Wall time", None, "s"),
    "energy_given_up": label("Energy given up", None, "%"),
    "material_saved": label("Shell material saved", None, "%"),
    "spanwise_moment": label("Flapwise moment", "M(r)", "N m"),
    "flexibility_weight": label("$1/c(r)^3$", None, r"m$^{-3}$"),
    "probability_density": label("Probability density", None, r"s m$^{-1}$"),
    "bin_mass": label("Bin mass", "m_b", "--"),
    "ratio_to_reference": label("Ratio to", r"\mathbf{x}_0", "--"),
    "shell_material": label(r"Shell material relative to $\mathbf{x}_0$", None, "%"),
    "objective": label("Objective", "J", "--"),
    "wind_speed_bin": label("Wind speed bin", None, r"m s$^{-1}$"),
}


# ---------------------------------------------------------------------------
# panel titles
# ---------------------------------------------------------------------------

def title(ax, text):
    """
    Set a short panel title. Sentence case, at most seven words, no
    parentheses, no number with more than three significant figures. Used
    only where a two-panel figure needs its panels distinguished; a
    single-panel figure takes its title from the LaTeX caption instead.
    """

    words = text.split()
    if len(words) > 7:
        raise FigureStyleError(
            f"panel title has {len(words)} words (limit 7): {text!r}")
    if "(" in text or ")" in text:
        raise FigureStyleError(
            f"panel title carries parentheses: {text!r}")
    for number in re.findall(r"\d+(?:\.\d+)?", text):
        digits = number.replace(".", "").lstrip("0")
        if len(digits) > 3:
            raise FigureStyleError(
                f"panel title carries a number with more than three "
                f"significant figures: {number!r} in {text!r}")
    ax.set_title(text)
    return ax


# ---------------------------------------------------------------------------
# tick formatting
# ---------------------------------------------------------------------------

def sci(ax, which="both"):
    """
    Math-text tick formatting for the given axes, so powers of ten read
    `10^{-6}` rather than `1e-06`. Returns the axes.
    """

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((0, 0))
    if which in ("both", "x"):
        ax.xaxis.set_major_formatter(formatter)
    if which in ("both", "y"):
        ax.yaxis.set_major_formatter(formatter)
    return ax


# ---------------------------------------------------------------------------
# refusal rules
# ---------------------------------------------------------------------------

class FigureStyleError(Exception):
    """A figure breaks the style rules and was not written."""


#: Substrings that may not appear in any text of any committed figure.
#: Matched case-insensitively.
FORBIDDEN = (
    "2026-", "resolved", "grounded", "phase", "audit",
    "a3", "a4", "b1", "b2", "b3", "b5",
    "provisional", "placeholder", "scratch",
    "before/after", "fitted schmitz",
    "(energy optimum, phase",
    "chord_max_m", "chord_min_m", "twist_min", "twist_max",
    "manufacturing.min_chord_m", "configured bounds",
    "re-run", "journal", "decision",
)

#: Whole words that may not appear in any text of any committed figure. A
#: legend names the thing, not its history: "bilinear interpolant" /
#: "C^1 interpolant", never "old"/"new".
FORBIDDEN_WORDS = ("old", "new")

_WORD_RE = re.compile(r"(?<![A-Za-z0-9])(" + "|".join(FORBIDDEN_WORDS)
                      + r")(?![A-Za-z0-9])", re.IGNORECASE)


def _panel_count(fig):
    """
    The number of panels: distinct axes positions. Twin axes share one
    position rectangle and count once.
    """

    positions = []
    for ax in fig.axes:
        box = tuple(round(v, 3) for v in ax.get_position().bounds)
        if box not in positions:
            positions.append(box)
    return len(positions)


def _texts(fig):
    for text in fig.findobj(matplotlib.text.Text):
        value = text.get_text()
        if value:
            yield value


def check(fig):
    """
    Apply the refusal rules to `fig`: at most two panels, no super-title, no
    forbidden string in any text object. Raises `FigureStyleError`.
    """

    panels = _panel_count(fig)
    if panels > 2:
        raise FigureStyleError(
            f"figure has {panels} panels; a figure has at most two axes, "
            "side by side. Split it into two figure files.")

    suptitle = getattr(fig, "_suptitle", None)
    if suptitle is not None and suptitle.get_text():
        raise FigureStyleError(
            f"figure carries a super-title: {suptitle.get_text()!r}. "
            "A super-title is not part of the style.")

    for value in _texts(fig):
        lowered = value.lower()
        for needle in FORBIDDEN:
            if needle in lowered:
                raise FigureStyleError(
                    f"figure text carries a forbidden string {needle!r}: "
                    f"{value!r}")
        match = _WORD_RE.search(value)
        if match:
            raise FigureStyleError(
                f"figure text carries the forbidden word {match.group(0)!r}: "
                f"{value!r}")


def save(fig, path):
    """
    The only way a committed figure is written. Applies the refusal rules
    (`check`), writes the PNG at `DPI`, and returns `path`.
    """

    check(fig)
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fig.savefig(path, dpi=DPI)
    return path
