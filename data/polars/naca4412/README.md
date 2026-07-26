NACA 4412 was the original primary airfoil target for the polar cache. The primary
target has since moved to S809 (`../s809/`) — see PROJECT_PLAN.md and the
2026-07-04 journal entry for the rationale (S809 is the actual NREL Phase VI rotor
airfoil and has real low-Re tunnel data from Delft/OSU/CSU, unlike this NACA
section). This cache is retained as a validated secondary reference dataset; the
pipeline (`src/xfoil/build_polar_cache.py`) and lookup interface
(`src/xfoil/polar_lookup.py`) are airfoil-agnostic, so nothing here needed to
change to make room for S809.

**This cache is at `Ncrit = 9`** (XFOIL's default, clean wind tunnel), whereas
the primary S809 cache was rebuilt at `Ncrit = 5` on 2026-07-26 — the 21%-thick
S809 could not be converged onto a physical solution branch at low Re under the
clean-tunnel assumption, while NACA 4412 at 12% thickness has no such trouble.
That difference is deliberate and this set is staying as it is: it is a
reference dataset, not what the solver reads. `validate_polars.py naca4412`
scores 3/5 on marginal Re-trend non-monotonicity (Cl_max dips 0.8%, Cl(0)
spreads 0.051 against a 0.05 threshold) — recorded rather than tuned away.
See the 2026-07-26 journal entry.
