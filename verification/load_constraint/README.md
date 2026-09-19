# Phase 4: the flapwise root-moment constraint (KS aggregate, adjoint Jacobian)

**New 2026-09-19 (Phase 4, Step 2d), `src/` commit `a220e1b`** — the commit the
JSON files stamp as `src_commit`. Run under the 300 rpm operating law
(`lambda(V) = min(6.5, Omega_max R / V)`, `Omega_max = 300 rpm`, `V_c = 9.67 m/s`)
and the project bounds in `config/rotor_design.yaml` (`chord_max_m = 0.30 m`,
`chord_min_m = 0.045 m`, `twist_min_deg = -2`, `twist_max_deg = 35`, local
solidity cap 0.5). The fixed generator rating is `3822.189755449124 W` (B2
provisional). The `--eps 0` run below is the **Phase 5 production optimum for
now**; the Pareto runs price the load reduction against it.

## What was run

    python verification/load_constraint/run_load_checks.py                    # ~1 min
    python verification/load_constraint/run_constrained_slsqp.py --eps 0       # ~30 s
    python verification/load_constraint/run_constrained_slsqp.py --eps 0.02
    python verification/load_constraint/run_constrained_slsqp.py --eps 0.05
    python verification/load_constraint/run_constrained_slsqp.py --eps 0.10
    python verification/load_constraint/run_pareto.py                          # ~4 min (8 runs)

All of them read `verification/baseline/x0.json`, `verification/fd_step_size/sweep.json`
(the global `h* = 3.162277660168379e-06`) and
`verification/adjoint_optimisation/result.json`'s `u*`; `load_constraint` runs
after `adjoint_optimisation` in `verification/README.md`'s order.

## Definition

**Load.** The flapwise root bending moment **per blade**, the same integrand
`design.baseline.root_bending_moment` computes from a solved rotor:

    m_{b,i} = 1/2 rho w_{b,i}^2 c_i Cn_{b,i} (r_i - r_hub)     [N m/m, per blade]
    M_b     = sum_i t_i m_{b,i}                                 [N m]     (trapezoid t_i)

with `w = V (1-a)/sin phi`, `Cn = Cl cos phi + Cd sin phi`, `r_hub = root_fraction R`.
No blade count (per blade) and no `Omega_b` (the power integrand has both; this
does not). The model is out-of-plane moment, not rotated into the section flap
axis; no centrifugal relief, no gravity, no dynamic amplification, no in-plane
component.

**Load operating set `L`** — the B3-independent points, fixed at construction
from `x0` and never recomputed per design (a set that moved with `d` would be a
non-smooth constraint):

    L = { (V_b, lambda_b) : bin midpoint b with P_b(x0) <= P_rated }
        union { (11.0, lambda(11.0) = 5.711986642890533) }
      = 3.5 .. 9.5 m/s at lambda = 6.5 ;  10.5 at lambda = 5.984 ;  11.0 at 5.712
      -> nine points, 9 x 25 = 225 states

**KS aggregate**, on moments normalised by `M_ref = M(x0)` at the design
condition (`177.3755406092969 N m`, so `rho` is dimensionless and the constraint
is O(1) for SLSQP):

    Mhat_b     = M_b / M_ref
    KS_rho(Mhat) = max_b Mhat_b + ln( sum_b exp(rho (Mhat_b - max)) ) / rho
    dKS/dMhat_b  = softmax weights, sum to 1

`rho = 100` by default. **The constraint**, as SciPy sees it, for reduction
fraction `eps`:

    g_eps(u) = (1 - eps) KS_rho(x0) - KS_rho(u) >= 0

written in the code as the regrouped `(KS0 - KS(u)) - eps KS0`, algebraically
identical and exact to the bit at `u0`. The limit is `KS_rho(x0)`, **not**
`Mhat_max(x0) = 1`: at 1 the baseline would be infeasible by the conservatism at
`eps = 0`. The Jacobian is `dg/du = -(dKS/dd) * span`, the discrete adjoint of
`adjoint/loads.py`.

**The problem** every constrained run solves:

    minimise   fun(u) = J(u)/|J(u0)|,   J = -AEP
    over       u in [0, 1]^10
    s.t.       polar-cache Reynolds envelope (50 rows) and solidity cap (25 rows)
               g_eps(u) >= 0

## Tiers 1-3 at `x0` and at the unconstrained optimum `u*`

`run_load_checks.py` (`checks.json`). Tier 1 is a complex step (`h = 1e-30`) of
the code's own value path, reported as the worst mixed error
`|estimate - partial| / max(1, |partial|)`; Tier 2 is the forward-mode tangent
against the adjoint direction and the adjoint assembly identity; Tier 3 is the
constraint Jacobian against central FD at the global `h*` and `h*/sqrt(10)`,
`h* sqrt(10)`, with `eps_j` re-measured as the three-step local jitter and the
measured round-off floor (below).

| point | Tier 1 `m` (worst) | Tier 1 KS (worst) | Tier 2 tangent | Tier 2 assembly | Tier 3 worst `|diff|/eps` | Taylor min ratio |
|---|---|---|---|---|---|---|
| `x0` | 3.19e-14 (`dm_dphi`) | 8.33e-16 (`dKS_dx`) | 8.9e-16 | 0.0 | **0.396** (`twist_0`) | 99.0 |
| `u*` | 3.62e-14 (`dm_dphi`) | 1.17e-15 (`dKS_dx`) | 8.8e-16 | 0.0 | **22.64** (`chord_1`) — see below | 94.8 |

Tier 3 at `x0`, all ten variables, at `h* = 3.162e-6` (`|diff|` and the floor in
units of `KS` per unit `u`; the floor is Tier 4's `delta_g/h*` from nine samples
of `g` along `u + t e_j`, `t = -4e-12 .. 4e-12`):

| variable | `|adj - fd|` | `eps_j` | floor | `|diff|/eps_j` | `|diff|/floor` |
|---|---|---|---|---|---|
| chord_0 | 9.42e-12 | 3.50e-11 | 2.26e-11 | 0.27 | 0.42 |
| chord_1 | 3.20e-14 | 1.16e-10 | 2.78e-11 | 0.00 | 0.00 |
| chord_2 | 2.88e-11 | 8.45e-11 | 2.78e-11 | 0.34 | 1.04 |
| chord_3 | 6.51e-12 | 7.21e-11 | 3.44e-11 | 0.09 | 0.19 |
| chord_4 | 4.49e-11 | 1.14e-10 | 2.78e-11 | 0.39 | 1.62 |
| twist_0 | 2.21e-11 | 5.58e-11 | 3.87e-11 | 0.40 | 0.57 |
| twist_1 | 1.65e-11 | 9.46e-11 | 2.19e-11 | 0.17 | 0.75 |
| twist_2 | 1.05e-11 | 1.11e-10 | 4.02e-11 | 0.09 | 0.26 |
| twist_3 | 1.22e-10 | 9.17e-10 | 2.52e-11 | 0.13 | 4.83 |
| twist_4 | 3.66e-11 | 4.31e-10 | 3.27e-11 | 0.08 | 1.12 |

At `u*` the disagreement is 0.03 … 4.8 x the measured floor for every variable;
the one variable above `3 eps_j` is the subject of the next section.

**One variable at `u*` fails its acceptance scale, and it is reported as found.**
`chord_1` at `u*`:

| quantity | value |
|---|---|
| `|adj - fd|` | 7.74e-12 |
| `eps_j` (the acceptance scale) | **3.42e-13** (the smallest of the ten, ~2 decades below the median) |
| `|diff|/eps_j` | **22.64** (limit 3) |
| measured round-off floor `delta_g/h*` | 4.41e-11 |
| `|diff|` / that floor | **0.18** |
| Taylor remainder, `u*`, min over 3 draws | 94.8 (limit 30) |

This is the same fragility `verification/gradient_verification/README.md`
records for the objective at `x0`: a three-step flatness estimate `eps_j` that
lands on a degenerate neighbour pair and comes out below the true round-off
floor. The adjoint is not the suspect — the disagreement is *below* the measured
floor, the absolute disagreement is among the smallest of the ten, and the
whole-chain Taylor test passes at 94.8. No tolerance and no estimator was
changed; `checks.json` records `passes: false`, `fragile_only: true`,
`resolved: true` for this point. It is in the hand-off's BLOCKING list for MJ.
It does **not** touch the constrained runs: they start at `x0`, where the same
Jacobian passes with worst ratio 0.396.

## The constrained optima

`result_eps{E}.json`, `E = 0, 0.02, 0.05, 0.1`; SLSQP `ftol = 1e-8`,
`maxiter = 200`, from `x0`, exit status 0 at every `eps`. `x0` has
`M = 177.376 N m` and `517.51 N` thrust at the design condition; the
unconstrained `u*` has `177.919 N m` (+0.306 %) and `526.66 N` (+1.77 %).

| `eps` | AEP [MWh/yr] | vs `x0` | cost vs `u*` | `KS0 -> KS*` | KS reduction | rated M `x_c` [N m] | M vs `x0` | thrust `x_c` [N] | thrust vs `x0` | active | moment slack | KKT residual | `nit` | moment solves | wall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 10.262719 | **+0.1465 %** | 0.00017 % | 1.000436 -> 1.000436 | 0.000 % | 177.411 | +0.020 % | 525.49 | +1.54 % | `chord_0` upper | 5.1e-8 | 9.6e-6 | 34 | 34 | 28 s |
| 0.02 | 10.260624 | +0.1260 % | 0.0206 % | -> 0.980427 | 2.000 % | 173.833 | -1.997 % | 516.13 | -0.27 % | `chord_0` upper, `chord_4` lower | 2.2e-8 | 4.3e-6 | 39 | 39 | 33 s |
| 0.05 | 10.249582 | +0.0183 % | 0.1282 % | -> 0.950414 | 5.000 % | 168.441 | -5.037 % | 501.99 | -3.00 % | same | 1.2e-8 | 6.0e-6 | 37 | 38 | 30 s |
| 0.10 | 10.203370 | -0.4327 % | **0.5785 %** | -> 0.900392 | 10.000 % | 159.480 | -10.089 % | 478.52 | -7.54 % | same | 2.0e-9 | 6.0e-5 | 29 | 31 | 24 s |

Every run leaves the constraint **active** (`KS` reduction equals `eps KS0` to
five digits) and `chord_0` on its upper bound. Cost relative to the
unconstrained optimum is 0.00017 / 0.0206 / 0.1282 / 0.5785 % — inside the
stated bands (`eps = 0`: +0.10 … +0.14666 %; `eps = 0.10`: 0.3 … 3 %), and far
below the 5 % defect gate. The `eps = 0` result is the statement that the
unconstrained optimum's extra +0.285 % KS bought essentially no AEP: capping the
load at the baseline's own level costs 0.0002 % of the +0.1467 % gain.

`eps = 0` active set: `chord_0` at its upper bound; no envelope or solidity row
active (tightest envelope row `floor r = 1.9660 m`, 16.6 mm slack; `sigma_max`
0.409 against the 0.5 cap). KKT check (SLSQP exposes no multipliers): the
component of `grad fun` in the span of the active normals has norm 2.46e-3 and
the residual after removing it is 9.6e-6 (scaled), on the moment row and the
`chord_0` bound.

## The Pareto front

`pareto.json` / `pareto.png`. Each `eps` solved **cold** from `x0` and **warm**
from the previous `eps`'s optimum; the pair must agree to the multi-start
study's spread `0.0166` in `u`.

| `eps` | cold AEP | warm AEP | `||u_cold - u_warm||_inf` | agrees | kept |
|---|---|---|---|---|---|
| 0 | 10.262719 | 10.262719 | 0 | yes | cold |
| 0.02 | 10.260624 | 10.260624 | 2.43e-04 | yes | cold |
| 0.05 | 10.249582 | 10.249582 | 1.40e-03 | yes | cold |
| 0.10 | 10.203370 | 10.203371 | 2.58e-03 | yes | warm (higher by 1e-6 MWh/yr) |

All four agree, so no optimum is ambiguous. The cost relative to `u*` is
monotone in `eps` (0.00017 -> 0.0206 -> 0.1282 -> 0.5785 %). The Pareto figure
plots the AEP change against `x0` versus rated-point moment reduction, with `x0`
at the origin and `u*` at **-0.306 % reduction / -0.147 % change** — the wrong
side of the origin, paying a *smaller* moment to gain AEP; the constrained curve
runs to +10.1 % reduction at +0.43 % AEP change.

## What the KS actually did

On the normalised moments the 11 m/s rated point carries softmax weight
**0.957** at `x0` (and the 10.5 m/s point 0.043; the other seven sum to <1e-37),
so at `rho = 100` the KS *is* the rated-point moment to within the
conservatism. The aggregate is kept because it is the principled single smooth
scalar with one adjoint right-hand side, not because the eight lower points
shape the answer. At `u*` the rated weight rises to 0.978.

| `rho` | `KS(x0)` | conservatism `KS - max` | `KS(u*)` | conservatism |
|---|---|---|---|---|
| 30 | 1.0115041 | 0.0115041 | 1.0126510 | 0.0095869 |
| 100 | 1.0004357 | 0.0004357 | 1.0032871 | 0.0002230 |
| 300 | 1.0000003 | 2.944e-7 | 1.0030641 | 3.824e-8 |

Per-point moments and weights over `L` (N m):

| V [m/s] | lambda | rpm | `M(x0)` | `Mhat(x0)` | `w(x0)` | `M(u*)` | `Mhat(u*)` | `w(u*)` |
|---|---|---|---|---|---|---|---|---|
| 3.5 | 6.500 | 108.6 | 19.090 | 0.1076 | ~0 | 18.735 | 0.1056 | ~0 |
| 4.5 | 6.500 | 139.7 | 33.069 | 0.1864 | ~0 | 32.598 | 0.1838 | ~0 |
| 5.5 | 6.500 | 170.7 | 50.463 | 0.2845 | ~0 | 49.869 | 0.2812 | ~0 |
| 6.5 | 6.500 | 201.7 | 71.278 | 0.4018 | ~0 | 70.402 | 0.3969 | ~0 |
| 7.5 | 6.500 | 232.8 | 95.455 | 0.5381 | ~0 | 94.264 | 0.5314 | ~0 |
| 8.5 | 6.500 | 263.8 | 122.990 | 0.6934 | ~0 | 121.466 | 0.6848 | ~0 |
| 9.5 | 6.500 | 294.8 | 153.900 | 0.8677 | ~0 | 152.038 | 0.8572 | ~0 |
| 10.5 | 5.984 | 300.0 | 171.857 | 0.9689 | 0.0426 | 171.193 | 0.9651 | 0.0221 |
| 11.0 | 5.712 | 300.0 | 177.376 | 1.0000 | 0.9574 | 177.919 | 1.0031 | 0.9779 |

## What the optimiser did to the blade

The load comes off the outer blade. Relative to `x0`, every constrained optimum
pins `chord_0` at the 0.30 m upper bound (as the unconstrained one does) and
then, as `eps` rises: `chord_1` and `chord_2` fall (0.182 -> 0.165 and
0.097 -> 0.109), `chord_4` drops to and then sits at its 0.045 m lower bound,
`chord_3` falls (0.079 -> 0.063 at `eps = 0.10`), and the tip twist rises from
-1.94 deg toward -0.60 deg while `twist_0` rises slightly. Control points
(physical, from `result_eps{E}.json::x_c`):

| | c0 | c1 | c2 | c3 | c4 [m] | th0 | th1 | th2 | th3 | th4 [deg] |
|---|---|---|---|---|---|---|---|---|---|---|
| `x0` | 0.2758 | 0.1816 | 0.0971 | 0.0786 | 0.0670 | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 |
| `u*` | 0.3000 | 0.1796 | 0.1304 | 0.0769 | 0.0491 | 24.35 | 10.42 | 5.60 | 2.10 | -1.85 |
| `eps = 0` | 0.3000 | 0.1817 | 0.1262 | 0.0784 | 0.0477 | 24.27 | 10.63 | 5.22 | 2.23 | -1.94 |
| `eps = 0.02` | 0.3000 | 0.1809 | 0.1158 | 0.0782 | 0.0450 | 24.38 | 10.61 | 4.51 | 1.93 | -1.87 |
| `eps = 0.05` | 0.3000 | 0.1772 | 0.1096 | 0.0730 | 0.0450 | 24.60 | 10.33 | 4.34 | 1.00 | -1.19 |
| `eps = 0.10` | 0.3000 | 0.1651 | 0.1094 | 0.0626 | 0.0450 | 25.01 | 9.52 | 4.68 | -0.02 | -0.60 |

`pareto.png`'s second and third panels show the same in chord and twist.

## Cut-out post-check (reported, never constrained)

At 20 m/s, `lambda = 62.831853071796 / 20 = 3.14159`, 300 rpm:

| | Ct | thrust [N] | moment [N m] | alpha [deg] |
|---|---|---|---|---|
| `x0` | 0.2593 | 631.1 | 216.46 | 10.9 … 30.1 |
| `u*` | 0.2814 | 684.9 | 228.54 | 9.6 … 29.3 |

**B3-dependent: the model holds `P = P_rated` with no mechanism, so the state
here is not the machine's; alpha reaches ~30 deg on the Viterna extrapolation.**
It is computed through the forward path only, never in a constraint and never in
a tier. The cut-out thrust is *larger* than at the design condition because the
rotor is held at 300 rpm into 20 m/s; a real above-rated mechanism (B3) would
change it.

## What would make this wrong

- **The load model.** Out-of-plane moment only (not rotated into the section
  flap axis), no centrifugal relief, no gravity, no dynamic amplification, no
  in-plane component. It is the same integrand the baseline quotes, so
  `x0` vs `x_c` comparisons are consistent, but the absolute number is not a
  structural load.
- **The fixed load set `L`.** If a design pushed a currently-uncapped bin over
  the rating, or the rated point ceased to dominate the KS, the constraint would
  be measuring a different thing. At `x0` the parent's `limited` mask is all
  `False` on `L` (`tests/test_loads.py`), and the 0.957 rated weight makes the
  aggregate the rated point; both are asserted where they can be.
- **Above rated (B3).** The constraint is at <= rated by construction; the
  cut-out case is a labelled post-check. If B3 lands, this artefact must be
  re-run.
- **The `u*` Tier 3 fragility.** The one failing variable is an acceptance-scale
  artefact, reported above; it does not touch `x0` or any constrained run.
- **The objective's known degeneracy.** A +0.147 % unconstrained gain is the
  subject of `docs/AEP_GAIN_AUDIT.md`; this artefact does not re-open it.

## Reproduce

    python verification/load_constraint/run_load_checks.py
    python verification/load_constraint/run_constrained_slsqp.py --eps 0
    python verification/load_constraint/run_pareto.py

## Files

- `run_load_checks.py` — Tiers 1-3, the moment tables, the cut-out post-check.
- `checks.json` — everything above, plus the per-variable Tier 3 tables, the
  measured floors, the per-point moments and weights, timings, the bounds/law
  labels, `M_ref`, `KS0`, the load set, and the `src_commit`.
- `moments_x0_xstar.png` — per-point moment and softmax weight at `x0` and `u*`.
- `run_constrained_slsqp.py` — the eps-constrained SLSQP (recorder, KKT check).
- `result_eps{E}.json`, `iterates_eps{E}.json`, `constrained_blade_eps{E}.png`
  for `E = 0, 0.02, 0.05, 0.1`.
- `run_pareto.py` — the cold/warm eps-sweep.
- `pareto.json`, `pareto.png` — the front and the blade shapes.
