"""
The polar layer: everything that *consumes* a polar cache, with no XFOIL
dependency.

`src/xfoil/` generates caches (it shells out to XFOIL); this package reads
them, extrapolates them and reasons about the Reynolds envelope they have to
cover. Keeping the two apart is what makes plan section 3.4's "XFOIL never runs
in the solve loop" a testable invariant rather than an intention -- nothing
importable from here can start an XFOIL process.

Task 2 adds:
  `viterna`       -- post-stall extrapolation to +/-180 deg, promoted out of
                     validation/export_qblade.py.
  `cache_format`  -- the cache CSV schema, including the provenance column
                     that keeps a filled or extrapolated point distinguishable
                     from a converged one.
  `envelope`      -- the design rotor's per-station Reynolds envelope, which
                     is what sets the SG6043 cache bounds.

Task 3 adds:
  `cache`         -- CSV load and ragged->rectangular reindexing into a
                     `PolarGrid`, lifted from `xfoil.polar_lookup`.
  `interpolant`   -- `PolarInterpolant`, the C1 analytically differentiable
                     (alpha, Reynolds) -> (Cl, Cd, Cm) surface built on top
                     of a `PolarGrid`.

Task 4 adds:
  `polar`         -- `CachedPolar`, the BEM solver's per-station view of that
                     surface at one fixed Reynolds number. Replaces
                     `bem.airfoil.S809Polar`; clamps nothing, substitutes
                     nothing, raises out of range.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

from polars.cache import PolarGrid, PolarGridError
from polars.cache_format import (
    HEADER,
    SOURCE_GAP_RETRY,
    SOURCE_LABELS,
    SOURCE_LOCAL_FIT,
    SOURCE_VITERNA,
    SOURCE_XFOIL,
    extend_to_full_range,
    load_xfoil_band,
    read_polar_csv,
    write_polar_csv,
)
from polars.interpolant import PolarDomainError, PolarInterpolant
from polars.polar import (
    CachedPolar,
    CachedPolarFactory,
    interpolant_for,
    polar_factory_for,
)
from polars.viterna import (
    CD_MAX,
    CD_MAX_QBLADE_MATCHED,
    NEGATIVE_CL_SCALE,
    REVERSED_CD_FLOOR,
    build_full_range_polar,
    cd_max_finite_blade,
    reversed_cl_amplitude,
    viterna_coefficients,
)

# `viterna` (the single-point evaluator) is deliberately NOT re-exported: it
# would shadow the `polars.viterna` submodule on this package, so
# `from polars import viterna` would hand a caller a function where it asked
# for a module. Import it from `polars.viterna` directly.

__all__ = [
    "CD_MAX",
    "CD_MAX_QBLADE_MATCHED",
    "CachedPolar",
    "CachedPolarFactory",
    "HEADER",
    "NEGATIVE_CL_SCALE",
    "PolarDomainError",
    "PolarGrid",
    "PolarGridError",
    "PolarInterpolant",
    "REVERSED_CD_FLOOR",
    "SOURCE_GAP_RETRY",
    "SOURCE_LABELS",
    "SOURCE_LOCAL_FIT",
    "SOURCE_VITERNA",
    "SOURCE_XFOIL",
    "build_full_range_polar",
    "cd_max_finite_blade",
    "extend_to_full_range",
    "interpolant_for",
    "load_xfoil_band",
    "polar_factory_for",
    "read_polar_csv",
    "reversed_cl_amplitude",
    "viterna_coefficients",
    "write_polar_csv",
]
