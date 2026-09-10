# Objective smoothness gate

**Empty. Lands with plan step 1.8.** Created by work order Task 8 so the
location is fixed before the work starts, rather than being invented under
time pressure when the figures are needed.

## What goes here

Plan step 1.8, described in the brief as the single most valuable step in
Phase 1: sweep the objective `J(d)` along one design variable at a time at fine
resolution and look at it, before any gradient code exists. The failure it
exists to catch is a "staircase" — an objective that looks continuous at plot
scale but is piecewise constant or kinked underneath, which a gradient-based
optimiser reads as noise and a finite-difference check reads as a wrong
derivative.

Expected artefacts:

- `J` against each design variable at fine resolution, and the first
  difference of each — the derivative plot is where a bad objective reveals
  itself, exactly as it was for the polar interpolant in
  `../polar_interpolant/`.
- The same sweeps for any variable that leaves the polar table or the Reynolds
  envelope, since those now raise rather than clamp (Task 4) and the gate
  should show where the objective stops being defined.
- Per-station convergence status across the sweep, since `solve_rotor` reports
  non-convergence rather than raising (Task 5) and a point that failed must be
  visibly distinct from a point that converged to a bad value.

## Why this is now feasible

The audit measured the gate at roughly 3000 evaluations × 18 wind-speed bins.
At the 9.0 s per operating point it recorded, that is **~5.5 days of compute
per pass**, which is why the audit called Task 5's cost target the dependency
that made the whole step urgent. After Tasks 4 and 5 an operating point costs
about 15 ms, so the same sweep is **~13 minutes**.

That changes what this directory should contain: a sweep this cheap can be run
per design variable, at high resolution, and re-run whenever the objective
changes — so the artefacts here should be regenerable in one command rather
than being a one-off run that nobody dares repeat.
