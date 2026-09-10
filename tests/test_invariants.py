"""
Structural invariants of the source tree.

Work order Task 7 lists these; the first one is pulled forward to Task 4
because Task 4 is what makes it true. `bem/` stopped importing anything under
`xfoil/` when `S809Polar` was deleted, and an invariant that is merely true
today is worth two lines to make permanent.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import ast
import os

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(_HERE, "..", "src"))

#: Packages that must never reach the XFOIL side. `objective/` and `design/`
#: do not exist yet (Phase 1.5/1.6); they are named now so the invariant
#: applies from their first commit rather than being remembered later.
SOLVE_PATH_PACKAGES = ["bem", "polars", "objective", "design"]


def _imported_modules(path):
    """Every module name imported by a source file, per its AST."""

    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=path)

    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names.add(node.module)
    return names


def _python_files(package):
    directory = os.path.join(SRC, package)
    if not os.path.isdir(directory):
        return []
    return [
        os.path.join(root, name)
        for root, _dirs, files in os.walk(directory)
        for name in files
        if name.endswith(".py")
    ]


@pytest.mark.parametrize("package", SOLVE_PATH_PACKAGES)
def test_solve_path_never_imports_xfoil(package):
    """
    Plan section 3.4: XFOIL never runs in the loop.

    Checked over the AST rather than by importing and inspecting
    `sys.modules`, so it catches an import that is present in the source but
    guarded, lazy, or inside a function body -- which is exactly the shape a
    reintroduction would take.

    Note this is stricter than "never calls `xfoil_runner`": nothing in the
    solve path may import *anything* from `xfoil/`, including the retired
    `polar_lookup`. Splitting cache-consuming code (`polars/`) from
    cache-generating code (`xfoil/`) is what makes the invariant checkable at
    all, and a partial ban would leave the boundary ambiguous.
    """

    offenders = []
    for path in _python_files(package):
        for module in _imported_modules(path):
            if module == "xfoil" or module.startswith("xfoil."):
                offenders.append(f"{os.path.relpath(path, SRC)} imports {module}")

    assert not offenders, (
        f"src/{package}/ must not import anything under src/xfoil/:\n  "
        + "\n  ".join(offenders)
    )


def test_the_retired_polar_lookup_globals_are_gone():
    """
    Task 4 retired `ACTIVE_AIRFOIL` / `set_active_airfoil()` / `get_polar()`
    and the `_lookup_cache` dict.

    A single global "current airfoil" made the answer depend on call order
    once two caches were live at the same time. The `PolarLookup` class itself
    is deliberately retained as the pre-remediation baseline for the committed
    staircase figure -- it is the module-level state that had to go, not the
    bilinear interpolation it sat on.
    """

    import xfoil.polar_lookup as polar_lookup

    for name in ("ACTIVE_AIRFOIL", "set_active_airfoil", "get_polar", "_lookup_cache"):
        assert not hasattr(polar_lookup, name), (
            f"xfoil.polar_lookup.{name} is back; the airfoil is a property of "
            "the blade (bem.rotor.RotorGeometry.polar_cache), not module state"
        )

    assert hasattr(polar_lookup, "PolarLookup"), (
        "PolarLookup is retained on purpose -- "
        "verification/polar_interpolant/generate_plots.py draws the committed "
        "before/after figure from it"
    )
