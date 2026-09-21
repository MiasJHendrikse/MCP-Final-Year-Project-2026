# The polar interpolant: the value, the derivative, and the staircase it removed

The figures behind Task 3's claim: the C¹ interpolant
(`polars.interpolant.PolarInterpolant`) replaces the bilinear lookup's
piecewise-constant alpha-derivative — the staircase a gradient-based
optimiser sees as noise — with one smooth, analytically differentiable
curve. The derivative plot is the actual test; the value plot would look
fine under almost any scheme.

## What was run

    python verification/polar_interpolant/generate_plots.py
    python verification/polar_interpolant/generate_plots.py --airfoil s809 --re 600000

No XFOIL and no BEM run: every figure reads the committed cache CSVs under
`data/polars/` (the default cache is the SG6043; `--airfoil` selects the
S809, `--re` the staircase Reynolds number). Figures are written through
`src/plotting/figstyle.py::save`.

## The figures

| file | shows |
|---|---|
| `cl_vs_alpha.png` | `Cl` against alpha, one curve per cached Reynolds number (SG6043: 40 k, 100 k, 200 k, 400 k, 800 k) |
| `dcl_dalpha_vs_alpha.png` | the analytic `dCl/dalpha`, per degree, at the same Reynolds numbers — no staircase, no sawtooth |
| `staircase_before_after.png` | the bilinear lookup's central finite-difference derivative against the C¹ interpolant's analytic derivative, on one axis, for the **SG6043 cache at Re = 200 000** — the design Reynolds number of the committed rotor |
| `staircase_s809_re600k.png` | the same measurement on the **S809 cache at Re = 600 000**, the cache and Reynolds number this check was first made on; kept for the record under its own name |

The staircase window is alpha in [4.0, 5.6] deg, attached flow through the
first half of the lift curve.

## What would make this wrong

- **A re-built cache.** Every curve is a read of the committed CSVs; a cache
  re-run at a different `N_crit`, envelope or gap policy changes every
  figure, and the JSON provenance of that build is the record, not this
  README.
- **Reading the low-Re SG6043 curves as interpolant artefacts.** At
  Re = 40 k–100 k the cached polar itself carries scatter and stall-region
  jumps, so the derivative curves grow wiggles between data points there;
  that is the data, reproduced faithfully (no staircase, no sawtooth), not a
  defect the C¹ construction introduced.
- **The plot is not a fit metric.** It compares shapes on one axis; the
  cache-vs-measurement RMSEs live in `results/sg6043_uiuc_validation/`.
