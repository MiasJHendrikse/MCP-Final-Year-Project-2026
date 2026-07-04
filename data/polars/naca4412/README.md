NACA 4412 was the original primary airfoil target for the polar cache. The primary
target has since moved to S809 (`../s809/`) — see PROJECT_PLAN.md and the
2026-07-04 journal entry for the rationale (S809 is the actual NREL Phase VI rotor
airfoil and has real low-Re tunnel data from Delft/OSU/CSU, unlike this NACA
section). This cache is retained as a validated secondary reference dataset; the
pipeline (`src/xfoil/build_polar_cache.py`) and lookup interface
(`src/xfoil/polar_lookup.py`) are airfoil-agnostic, so nothing here needed to
change to make room for S809.
