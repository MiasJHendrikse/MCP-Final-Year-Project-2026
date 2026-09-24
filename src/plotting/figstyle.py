"""
One figure style for every committed plot (the figure-style module).

Every committed figure of this project is written through `save()`. It is the
mechanical guarantee behind two rules the report states: **a figure has at
most two panels, side by side**, and **no figure carries project history**
(no dates, no phase numbers, no "old"/"new", no bounds labels). A figure that
breaks either rule cannot be written -- the call raises instead.

What the module provides
------------------------

* `apply()`      -- sets `matplotlib.rcParams` once: serif text (Palatino, the
                    report's body face), 9 pt base, a closed axes box with
                    inward ticks, a light grid, framed legends, one colour
                    cycle (`PALETTE`), 300 dpi PNG output.
* `PALETTE`      -- the colour cycle. Its first three entries are the blade
                    colours of `BLADES`, so a series drawn from the default
                    cycle and a blade drawn from `BLADES` never disagree.
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
* `twin()`       -- a right-hand y-axis sharing x, with the host's right-side
                    ticks switched off so the two sets do not double up.
* `merged_legend()` -- one legend carrying a host axes' and its twin's entries.
* `legend_below()` -- a legend under the axes, for a panel with no free corner.
* `overlaps()`   -- every place a text or legend sits on plotted data, or on
                    another text.
* `save()`       -- asserts the panel count, the absence of a super-title,
                    the forbidden-string list and (by default) that no text or
                    legend sits on the data, then writes the PNG.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import os
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.collections  # noqa: E402
import matplotlib.lines  # noqa: E402
import matplotlib.patches  # noqa: E402
import matplotlib.path  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text  # noqa: E402
import matplotlib.ticker  # noqa: E402
import matplotlib.transforms  # noqa: E402
import numpy as np  # noqa: E402
from cycler import cycler  # noqa: E402

# ---------------------------------------------------------------------------
# size and style
# ---------------------------------------------------------------------------

#: A single panel: 85 mm, just over half the 160 mm A4 text block.
SINGLE = (3.35, 2.6)
#: Two panels side by side, filling the text block.
DOUBLE = (6.3, 2.6)

DPI = 300

SERIF = ["Palatino Linotype", "TeX Gyre Pagella", "DejaVu Serif"]


#: The colour cycle: dark, print-safe hues that stay distinct for the common
#: colour-vision deficiencies. The first three are the `BLADES` colours (navy,
#: rust, teal); the rest follow in contrast order.
PALETTE = ["#1f3b8b", "#b3452a", "#2a7f62", "#c7881c", "#6a4c93",
           "#4a90c8", "#555555"]

#: A lighter partner of the first two palette colours, for a paired series
#: that belongs to the same thing (shell and solid material, for instance).
LIGHT = {"#1f3b8b": "#8fa4d8", "#b3452a": "#e3a894"}


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
        # A closed box with inward ticks: the engineering-report convention.
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.linewidth": 0.7,
        "axes.edgecolor": "0.15",
        "axes.prop_cycle": cycler(color=PALETTE),
        "axes.axisbelow": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.pad": 4,
        "ytick.major.pad": 4,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        # A framed, near-opaque legend: where it must sit near data, the data
        # passes behind the box, never through the words.
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "0.7",
        "legend.fancybox": False,
        "legend.borderpad": 0.4,
        "legend.handlelength": 1.8,
        "legend.labelspacing": 0.3,
        "lines.linewidth": 1.4,
        "lines.markersize": 4,
        "figure.figsize": SINGLE,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def twin(ax):
    """
    A right-hand y-axis sharing `ax`'s x-axis. The host's right-side ticks are
    switched off so the twin's are the only ones on that spine; the twin draws
    no grid of its own.
    """

    ax.tick_params(axis="y", which="both", right=False)
    other = ax.twinx()
    other.tick_params(axis="x", which="both", top=False)
    other.grid(False)
    return other


def merged_legend(ax, other, **kwargs):
    """One legend on `ax` carrying the entries of `ax` and its twin `other`."""

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = other.get_legend_handles_labels()
    return ax.legend(h1 + h2, l1 + l2, **kwargs)


def legend_below(ax, handles=None, labels=None, ncol=2, offset=0.22):
    """
    A legend centred under `ax`, clear of the tick labels and the x label, for
    a panel whose data leaves no free corner. `offset` is the gap below the
    axes as a fraction of the axes height.
    """

    kwargs = dict(loc="upper center", bbox_to_anchor=(0.5, -offset),
                  ncol=ncol, frameon=False)
    if handles is None:
        return ax.legend(**kwargs)
    return ax.legend(handles, labels, **kwargs)


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
    # A digit run that starts a word: "60" and "9.5" are numbers, the digits in
    # a name like "SG6043" are not.
    for number in re.findall(r"(?<![A-Za-z0-9])\d+(?:\.\d+)?", text):
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

    A logarithmic axis already carries `LogFormatterSciNotation`, which
    renders the exponent as math text; forcing `ScalarFormatter` on it
    mislabels the ticks, so a log axis is left alone.
    """

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((0, 0))
    if which in ("both", "x") and ax.xaxis.get_scale() != "log":
        ax.xaxis.set_major_formatter(formatter)
    if which in ("both", "y") and ax.yaxis.get_scale() != "log":
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


# ---------------------------------------------------------------------------
# text on data
# ---------------------------------------------------------------------------

#: An artist given this gid is exempt from the overlap check: a curve label
#: placed on its own curve by design, for instance.
OVERLAP_OK = "overlap-ok"

_PAD = 1.0  # display pixels shaved off every text box before it is tested


def _shrink(box):
    return matplotlib.transforms.Bbox.from_extents(
        box.x0 + _PAD, box.y0 + _PAD, box.x1 - _PAD, box.y1 - _PAD)


def _polyline_enters(points, box):
    """Whether the polyline `points` (N x 2, display units) enters `box`."""

    x0, y0, x1, y1 = box.x0, box.y0, box.x1, box.y1
    finite = np.all(np.isfinite(points), axis=1)
    inside = (finite & (points[:, 0] > x0) & (points[:, 0] < x1)
              & (points[:, 1] > y0) & (points[:, 1] < y1))
    if inside.any():
        return True
    ok = finite[:-1] & finite[1:]
    if not ok.any():
        return False
    a, b = points[:-1][ok], points[1:][ok]
    d = b - a
    # Liang-Barsky clipping of every segment against the box at once.
    t0 = np.zeros(len(a))
    t1 = np.ones(len(a))
    hit = np.ones(len(a), dtype=bool)
    for p, q in ((-d[:, 0], a[:, 0] - x0), (d[:, 0], x1 - a[:, 0]),
                 (-d[:, 1], a[:, 1] - y0), (d[:, 1], y1 - a[:, 1])):
        parallel = p == 0
        hit &= ~(parallel & (q < 0))
        r = np.divide(q, p, out=np.zeros_like(q), where=~parallel)
        t0 = np.where(~parallel & (p < 0), np.maximum(t0, r), t0)
        t1 = np.where(~parallel & (p > 0), np.minimum(t1, r), t1)
    return bool(np.any(hit & (t0 <= t1)))


def _line_points(line):
    x, y = line.get_xdata(orig=False), line.get_ydata(orig=False)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    style = line.get_drawstyle()
    if style not in (None, "default"):
        x, y = matplotlib.lines.STEP_LOOKUP_MAP[style](x, y)
    return line.get_transform().transform(np.column_stack([x, y]))


def _data_hit(ax, box):
    """A description of the first piece of `ax`'s data inside `box`, or None."""

    for line in ax.lines:
        if not line.get_visible() or line.get_gid() == OVERLAP_OK:
            continue
        pts = _line_points(line)
        if not len(pts):
            continue
        if line.get_linestyle() in ("None", "none", "", " "):
            # Markers only: test the points, not the joins between them.
            p = pts[np.all(np.isfinite(pts), axis=1)]
            if np.any((p[:, 0] > box.x0) & (p[:, 0] < box.x1)
                      & (p[:, 1] > box.y0) & (p[:, 1] < box.y1)):
                return f"markers {line.get_label()!r}"
            continue
        if _polyline_enters(pts, box):
            return f"line {line.get_label()!r}"
    for patch in ax.patches:
        if (not patch.get_visible() or patch.get_gid() == OVERLAP_OK
                or not isinstance(patch, matplotlib.patches.Rectangle)):
            continue
        extent = patch.get_window_extent()
        if extent.width > 0 and extent.height > 0 and extent.overlaps(box):
            return f"bar {patch.get_label()!r}"
    for coll in ax.collections:
        if (not coll.get_visible() or coll.get_gid() == OVERLAP_OK
                or not isinstance(coll, matplotlib.collections.PathCollection)):
            continue
        offsets = coll.get_offset_transform().transform(coll.get_offsets())
        if len(offsets) and np.any(
                (offsets[:, 0] > box.x0) & (offsets[:, 0] < box.x1)
                & (offsets[:, 1] > box.y0) & (offsets[:, 1] < box.y1)):
            return f"points {coll.get_label()!r}"
    return None


def _text_boxes(fig, renderer):
    """(description, display box, axes) for every text and legend in the axes."""

    for ax in fig.axes:
        if getattr(ax, "name", "") == "3d":
            continue
        legend = ax.get_legend()
        if (legend is not None and legend.get_visible()
                and legend.get_gid() != OVERLAP_OK):
            yield "legend", legend.get_window_extent(renderer), ax
        for text in ax.texts:
            if (text.get_visible() and text.get_text().strip()
                    and text.get_gid() != OVERLAP_OK):
                yield (f"text {text.get_text()!r}",
                       text.get_window_extent(renderer), ax)


def overlaps(fig):
    """
    Every place a text or legend on `fig` sits on plotted data (a line, a bar,
    a marker) of its own panel -- the host axes and any twin -- or on another
    text. Returns a list of strings; empty means clean.

    A legend placed outside the axes passes: the data is clipped at the axes
    edge and cannot reach it.
    """

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = list(_text_boxes(fig, renderer))
    found = []
    for what, box, ax in boxes:
        box = _shrink(box)
        if box.width <= 0 or box.height <= 0:
            continue
        # A text drawn in the axes must stay inside them, clear of the spines
        # and tick labels; a legend placed under the axes is outside by design.
        frame = ax.get_window_extent(renderer)
        if what != "legend" and not (
                frame.x0 <= box.x0 and box.x1 <= frame.x1
                and frame.y0 <= box.y0 and box.y1 <= frame.y1):
            found.append(f"{what} runs outside its axes")
            continue
        panel = tuple(round(v, 3) for v in ax.get_position().bounds)
        for other in fig.axes:
            if (getattr(other, "name", "") == "3d" or panel != tuple(
                    round(v, 3) for v in other.get_position().bounds)):
                continue
            clip = matplotlib.transforms.Bbox.intersection(
                box, other.get_window_extent(renderer))
            if clip is None:
                continue
            hit = _data_hit(other, clip)
            if hit:
                found.append(f"{what} sits on {hit}")
                break
    for i, (what_a, a, _) in enumerate(boxes):
        for what_b, b, _ in boxes[i + 1:]:
            if _shrink(a).overlaps(b):
                found.append(f"{what_a} overlaps {what_b}")
    return found


def save(fig, path, check_overlap=True):
    """
    The only way a committed figure is written. Applies the refusal rules
    (`check`) and, unless `check_overlap` is False, refuses a figure whose
    text or legend sits on its data (`overlaps`). Writes the PNG at `DPI` and
    returns `path`.

    `check_overlap=False` is for a schematic, whose labels sit beside drawn
    lines by design; its layout is checked by eye.
    """

    check(fig)
    if check_overlap:
        found = overlaps(fig)
        if found:
            raise FigureStyleError(
                f"{os.path.basename(path)}: text on data -- "
                + "; ".join(found))
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    fig.savefig(path, dpi=DPI)
    return path
