# The Schmitz reference blade x₀

Every result in the project is compared against one reference: the classical
Schmitz blade for this rotor, fitted to the same spline the optimiser uses.
`x0.json` is where every optimisation starts, and `baseline_reference.json`
holds the numbers everything else is compared against.

```
python verification/baseline/generate_baseline.py
```

This writes `x0.json`, `baseline_reference.json`, `baseline_geometry.png` and
`baseline_operating_line.png`. (`baseline.png`, an earlier three-panel figure,
is no longer produced; it's kept as the output of the run that made it.)

Re-running this should be rare and deliberate, as with the golden files in
`tests/golden/`. If these numbers change, something in the solver or the polar
layer has changed, and the commit message needs to say what and why.

## How it's built

1. The analytic **Schmitz** distribution (with wake rotation) at the SG6043
   best-L/D point.
2. Fitted by linear least squares to the **same spline the optimiser uses**
   (5 + 5 control points, the number chosen in the representation study).
   x₀ is the *fitted* design vector, not the analytic blade, because the
   optimiser can only start from a shape it can represent.
3. Checked against the design-variable bounds.
4. Evaluated with the same solver and settings as everything else.

| design point | |
|---|---|
| Reynolds number | 200,000 |
| design angle of attack | 5.36° |
| C_L | 1.2825 |
| max L/D | 98.5 |

## Fitting error

| | RMS |
|---|---|
| chord | **0.515 mm** (0.389 % of the 132 mm mean chord) |
| twist | **0.106°** |

That is well below any plausible manufacturing tolerance, so x₀ is a faithful
representation of the Schmitz blade. `verification/spline_fit_error/` shows
the fit changes AEP by only −0.015 %.

The blade runs from a **263 mm chord at the root to 69 mm at the tip**, with
twist from **23.07° to 0.57°**. A 263 mm root chord is wide for a 2.0 m blade,
which is typical of Schmitz designs, but it fits inside the 0.30 m chord bound
with 36 mm to spare. `x0.json`'s `feasibility` block records `checked: true`,
`clipped: false` and `violations: []`. If the blade had been outside the box,
`DesignBounds.clip_physical` would have clipped it and listed the violation.

## Reference numbers

| | |
|---|---|
| Cp at the design point (V = 11 m/s, λ = 5.711986642890533, 300 rpm) | **0.4516** |
| Ct at the design point | **0.7029** |
| aerodynamic power at the design point | **3656.9 W** |
| peak Cp on the fixed-λ curve | **0.4720 at λ = 6.5** |
| root bending moment | **177.38 N·m per blade** |
| peak thrust | **517.51 N at 11 m/s** |
| AEP | **10.2477 MWh/yr** |

The AEP is inside the 8–12 MWh/yr sanity band (`aep_in_sanity_band: true`;
see `docs/OUTSTANDING-INPUTS.md` §7).

**The peak Cp landing exactly at the design λ is the strongest check on the
construction.** The blade is designed for λ = 6.5, so if the Schmitz formulas,
the twist convention or the solver wiring were wrong, the peak would be
somewhere else. `tests/test_baseline.py` asserts it. The peak is on the
*fixed-λ* curve, which is what this check needs. In operation, the rotor is
capped at 300 rpm from 9.67 m/s, so at the 11 m/s design point it runs at
λ = 5.71 and Cp = 0.4516, deliberately below its design tip-speed ratio.

A Cp of 0.472, against the Betz limit of 0.593, is what a real blade with drag
achieves. The roughly 20 % shortfall is profile drag plus tip loss, not an
error.

The root bending moment is **per blade**, not per rotor, because the root
attachment carries one blade's load; a rotor total would overstate it by a
factor of B.

### Effect of the rotor-speed ceiling

| | no ceiling (λ = 6.5 everywhere) | 300 rpm ceiling |
|---|---|---|
| Cp at 11 m/s | 0.4720 | 0.4516 (λ = 5.712) |
| root bending moment | 206.7 N·m | 177.38 N·m |
| peak thrust | 597.2 N | 517.51 N |
| AEP | not computed | 10.2477 MWh/yr |

The drop is expected. Capping the rotor speed lowers λ above 9.67 m/s, which
lowers Cp and Ct at the same wind speed, and so both loads and the energy.
`tests/test_baseline.py` pins the current numbers.

The operating-line sweep runs from 3.0 to 11.0 m/s in 0.5 m/s steps (17
points, all converged), at λ = 6.5 up to 9.67 m/s and at 300 rpm above. Three
points sit on the ceiling (10.0, 10.5 and 11.0 m/s, at λ = 6.28, 5.98 and
5.71), which shows as the break in the λ curve in
`baseline_operating_line.png`.

## A bug this blade found

The first version of this baseline had a station fail to converge at
V = 3 m/s and Re = 70,995. The cause was a real value discontinuity in the
polar interpolant at angle-of-attack knots, for any Reynolds number between
table rows. That is fixed (see `tests/golden/README.md`), and the numbers
above are from after the fix. Every operating point and every λ in the sweep
now converges, which `test_every_operating_point_converges` asserts.

## Where the reference blade's Cp actually peaks

`cp_peak.json`, from `run_cp_peak.py`. The Schmitz formula ignores tip loss
and accounts for drag and Reynolds number only through one design point on
the polar, so x₀ is the exact optimum of a simplified problem, not of the full
BEM model. This measures the difference: the tip-speed ratio at which x₀'s Cp
peaks along fixed-λ lines, in the site atmosphere.

| V [m/s] | peak λ | peak Cp | Cp at λ = 6.5 |
|---|---|---|---|
| 5 | 6.25 | 0.4507 | 0.4493 |
| 7 | 6.41 | 0.4636 | 0.4634 |
| 9 | 6.46 | 0.4689 | 0.4689 |
| 11 | 6.49 | 0.4720 | 0.4720 |

At rated wind speed the peak is within 0.01 of the design ratio. It drifts
down to 6.25 at 5 m/s as the station Reynolds numbers fall away from the
200,000 design point, but the Cp given up there is only 0.3 %.

    python verification/baseline/run_cp_peak.py     # ~20 s
