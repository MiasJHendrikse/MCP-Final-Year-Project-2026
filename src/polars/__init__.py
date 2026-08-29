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

Task 3 adds `cache.py` and `interpolant.py` here.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

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
    "HEADER",
    "NEGATIVE_CL_SCALE",
    "REVERSED_CD_FLOOR",
    "SOURCE_GAP_RETRY",
    "SOURCE_LABELS",
    "SOURCE_LOCAL_FIT",
    "SOURCE_VITERNA",
    "SOURCE_XFOIL",
    "build_full_range_polar",
    "cd_max_finite_blade",
    "extend_to_full_range",
    "load_xfoil_band",
    "read_polar_csv",
    "reversed_cl_amplitude",
    "viterna_coefficients",
    "write_polar_csv",
]
