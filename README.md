# Adjoint-based optimisation of a small wind turbine blade

This is my final-year Mechanical Engineering project (MCP820S) at the Namibia
University of Science and Technology, supervised by Prof. Hannes van der Walt.

**Title:** *Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine
Blade for Low-Wind-Speed Conditions Using a Discrete-Adjoint BEM Framework*

I wrote a blade-element momentum (BEM) solver from scratch in Python, derived
and implemented its discrete adjoint by hand, and used the resulting exact
gradients to design a lighter blade for a 2 m, three-bladed small wind turbine.
The site is a low-wind, high-altitude ridge in the Khomas Hochland west of
Windhoek (mean wind speed 6.5 m/s at the 20 m hub, air density about 21 %
below sea level).

![The three blades: the Schmitz reference, the energy optimum and the minimum-material blade](verification/mass_optimisation/blades_rendered_iso.png)

## What the project found

I started out trying to maximise annual energy yield (AEP). The optimiser
converged cleanly, but it only found about +0.15 % over the classical Schmitz
blade. Working out why took some time. Under a fixed tip-speed-ratio operating
law with a rotor-speed ceiling, the Schmitz blade is already within 0.15 % of
the best energy any blade of this rotor can deliver. There was almost nothing
left to gain on energy.

So I turned the problem around. Annual energy became a constraint, and the
optimiser minimises blade material instead:

```
minimise    shell material   ∫ c dr
subject to  AEP                 ≥ AEP of the Schmitz reference
            root bending moment ≤ reference     (KS aggregate over the load cases)
            root stress proxy   ≤ reference
            tip deflection      ≤ reference
            monotone chord and twist, 60 mm minimum tip chord,
            Reynolds numbers inside the validated polar range
```

At exactly the reference energy, the optimised blade uses **3.4 % less shell
material** (about 1.62 kg instead of 1.68 kg per blade for a 2 mm glass/epoxy
skin). Loads, stress and deflection all stay at or below the reference. If
0.5 % of the energy can be given up, the saving grows to **5.3 %**. Ten
different starting points all converge to the same design.

![Material saved against energy given up, and the chord distributions along the front](verification/mass_optimisation/pareto_front.png)

## How it works

```
XFOIL polars → smooth (C¹) polar interpolant → BEM solver (Ning 2014 form)
  → AEP, loads and material functionals → discrete adjoint → SLSQP
```

- **Aerodynamic data.** SG6043 airfoil polars from XFOIL over the Reynolds range
  the blade sees, extended to ±180° with Viterna. I checked them against UIUC
  wind-tunnel measurements and ran an n_crit sensitivity study.
- **BEM solver.** A single-residual formulation with Prandtl tip/hub loss and a
  Glauert/Buhl high-induction correction. It is guaranteed to converge because
  it brackets the root instead of iterating.
- **Discrete adjoint.** Every partial derivative is derived by hand
  ([`docs/adjoint_derivation.md`](docs/adjoint_derivation.md)) and checked
  against the complex step. Adjoints are implemented for energy, root bending
  moment and tip deflection.
- **Optimisation.** A B-spline chord and twist parameterisation (5 + 5 control
  points), solved with SciPy's SLSQP.

## Verification

I tried to make every number in the project checkable.

- **Solver cross-validation.** On the NREL Phase VI rotor, my Cp–λ curve agrees
  with CCBlade to 0.51 % on average. It agrees with pyBEMT to 1.8 % and with
  QBlade to within a few percent. These are comparisons against other BEM
  codes, not against measured data. The CCBlade and pyBEMT runs used an
  earlier version of the S809 polar table and haven't been repeated since.
  ([`docs/validation/bem-cross-validation.md`](docs/validation/bem-cross-validation.md))
- **Gradients.** The adjoint gradients agree with the complex step to about
  1e-13 and with finite differences down to the finite-difference noise floor.
  One variable is off by a small amount. I traced it to round-off in the
  objective itself and documented it rather than loosening the tolerance.
- **Cost.** Each adjoint gradient costs about 1.1 objective evaluations, however
  many design variables there are. Finite differences cost 2n.
- **Smoothness.** 3,000 objective evaluations show the energy function is
  smooth enough to differentiate. It is C¹, with one curvature break at the Buhl
  threshold.

Each folder under [`verification/`](verification/) has a README explaining
what was run, what came out, how to reproduce it and what would prove it
wrong. [`verification/README.md`](verification/README.md) gives the order to
re-run them in.

## Repository layout

```
src/
  bem/           BEM solver: station residual, corrections, rotor, power curves
  polars/        polar cache loading, C¹ interpolant, Viterna extension
  xfoil/         XFOIL wrapper and polar cache generation
  design/        B-spline parameterisation, bounds, Schmitz reference blade
  objective/     wind resource (Weibull), AEP, loads, material proxies
  adjoint/       discrete adjoint systems (energy, root moment, deflection)
  gradients/     finite-difference and adjoint gradients, the optimisation problem
  validation/    cross-checks against CCBlade, pyBEMT, QBlade and UIUC data
  config/        typed loader for config/
config/          site, rotor and polar-cache configuration (YAML)
data/            airfoil coordinates, polar caches, reference data
docs/            design basis, assumptions, adjoint derivation, validation write-up
tests/           pytest suite, including a golden-file regression
verification/    scripts, results (JSON) and figures behind every reported number
```

Two airfoils have two different jobs. **SG6043** is the design airfoil. **S809**
is used only to validate the solver on the NREL Phase VI rotor, since other
people have computed that rotor too.

A few rules I kept throughout the code:

- **No hidden defaults.** Air density and viscosity are required arguments.
  Config values that haven't been decided yet raise an error instead of falling
  back to a guess.
- **No silent clamping.** Asking the polar data for an angle or Reynolds number
  outside its range raises an error. Clamping to the edge of the table would
  quietly change Cp by as much as 60 %.
- **`src/bem` never imports XFOIL.** A test checks this, so the solver runs from
  the committed polar cache alone.

## Running it

Requires Python 3.12+.

```bash
python -m pip install -r requirements.txt
```

Run the tests from the repository root:

```bash
pytest
```

A few results are expected and documented rather than hidden. Five tests are
marked `xfail`, because the polar-cache sanity checks don't fit a thin,
high-camber, low-Reynolds airfoil (the reasoning is in
[`data/polars/sg6043/README.md`](data/polars/sg6043/README.md)). One test is
left failing on purpose: the gradient case above, kept visible so that any
change to it shows up.

To reproduce the main result:

```bash
python verification/mass_optimisation/run_mass_slsqp.py --delta 0     # ~20 s
python verification/mass_optimisation/run_multistart.py               # ~3 min
python verification/mass_optimisation/run_sweep.py                    # ~4 min
```

Using the solver directly (run from `src/`):

```python
from config import load_phase_vi_rotor
from bem.rotor import phase_vi_geometry, PHASE_VI_RATED_RPM
from bem.powercurve import power_curve

rotor = load_phase_vi_rotor()
air = dict(air_density=rotor.air_density,
           kinematic_viscosity=rotor.kinematic_viscosity)

curve = power_curve(phase_vi_geometry(), [5, 7, 10, 13, 15],
                    rpm=PHASE_VI_RATED_RPM, **air)
for p in curve["points"]:
    print(f"V={p['v_inf']:5.1f} m/s  Cp={p['Cp']:.4f}  P={p['power']:8.1f} W")
```

XFOIL, CCBlade, pyBEMT and QBlade are only needed to regenerate the polar
cache or re-run the cross-validation. The caches are committed, so none of
them is needed to run the solver or the optimisation.

## Further reading

- [`docs/DESIGN-BASIS.md`](docs/DESIGN-BASIS.md): the site, machine and model,
  why the objective is material rather than energy, and the results with a
  source for every number
- [`docs/OUTSTANDING-INPUTS.md`](docs/OUTSTANDING-INPUTS.md): each modelling
  assumption, the real data that would replace it, and what would need
  re-running; also the laminate behind the kilogram and stress figures
- [`docs/adjoint_derivation.md`](docs/adjoint_derivation.md): the full adjoint
  derivation

## Author

Mias J. Hendrikse, Mechanical Engineering, NUST (2026)
