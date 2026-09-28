# The root-moment constraint

This folder adds a constraint on the flapwise root bending moment to the
energy optimisation, checks its adjoint gradient, and prices it: how much
energy each per cent of moment reduction costs.

The `eps = 0` run gives **x_c, the energy optimum with the root moment capped
at the Schmitz blade's value**: +0.1465 % AEP for +6.2 % shell material. It
is the energy reference that the minimum-material blade in
`verification/mass_optimisation/` is compared against.

All runs use the 300 rpm operating law (`λ(V) = min(6.5, Ω_max R / V)`,
ceiling reached at 9.67 m/s), the design bounds in
`config/rotor_design.yaml` (chord 0.045–0.30 m, twist −2° to 35°, local
solidity at most 0.5) and the fixed generator rating of 3822.189755449124 W.

## Reproducing it

    python verification/load_constraint/run_load_checks.py                    # ~1 min
    python verification/load_constraint/run_constrained_slsqp.py --eps 0       # ~30 s
    python verification/load_constraint/run_constrained_slsqp.py --eps 0.02
    python verification/load_constraint/run_constrained_slsqp.py --eps 0.05
    python verification/load_constraint/run_constrained_slsqp.py --eps 0.10
    python verification/load_constraint/run_pareto.py                          # ~4 min (8 runs)

These read `verification/baseline/x0.json`,
`verification/fd_step_size/sweep.json` (the global
`h* = 3.162277660168379e-06`) and the u\* in
`verification/adjoint_optimisation/result.json`, so this folder runs after
`adjoint_optimisation`.

## Definition

**The load** is the flapwise root bending moment per blade, the same
integrand `design.baseline.root_bending_moment` uses on a solved rotor:

    m_{b,i} = 1/2 rho w_{b,i}^2 c_i Cn_{b,i} (r_i - r_hub)     [N m/m, per blade]
    M_b     = sum_i t_i m_{b,i}                                 [N m]     (trapezoid weights t_i)

with `w = V (1-a)/sin phi`, `Cn = Cl cos phi + Cd sin phi` and
`r_hub = root_fraction R`. Unlike the power integrand, it has no blade count
(it's per blade) and no `Omega_b`. It is the out-of-plane moment, not rotated
into the section's flap axis, and it ignores centrifugal relief, gravity,
dynamic amplification and in-plane loads.

**The load cases `L`** are the operating points at or below rated. They are
fixed once from x₀ and never recomputed for a new design, because a set that
changed with the design would make the constraint non-smooth:

    L = { (V_b, lambda_b) : bin midpoint b with P_b(x0) <= P_rated }
        union { (11.0, lambda(11.0) = 5.711986642890533) }
      = 3.5 .. 9.5 m/s at lambda = 6.5 ;  10.5 at lambda = 5.984 ;  11.0 at 5.712
      -> nine points, 9 x 25 = 225 states

**The KS aggregate** combines the nine moments into one smooth value. The
moments are normalised by `M_ref = M(x0)` at the design condition
(177.3755406092969 N m), so `rho` is dimensionless and the constraint is of
order one for SLSQP:

    Mhat_b       = M_b / M_ref
    KS_rho(Mhat) = max_b Mhat_b + ln( sum_b exp(rho (Mhat_b - max)) ) / rho
    dKS/dMhat_b  = softmax weights, summing to 1

with `rho = 100`. For a reduction fraction `eps`, the constraint SciPy sees is

    g_eps(u) = (1 - eps) KS_rho(x0) - KS_rho(u) >= 0

The code writes it as `(KS0 - KS(u)) - eps KS0`, which is algebraically the
same and exact at `u0`. The limit is `KS_rho(x0)`, **not**
`Mhat_max(x0) = 1`; with a limit of 1, the reference blade itself would be
infeasible at `eps = 0`, because KS slightly overestimates the maximum. The
Jacobian is `dg/du = -(dKS/dd) * span`, from the discrete adjoint in
`adjoint/loads.py`.

**The optimisation problem:**

    minimise   fun(u) = J(u)/|J(u0)|,   J = -AEP
    over       u in [0, 1]^10
    s.t.       polar-table Reynolds envelope (50 rows) and solidity limit (25 rows)
               g_eps(u) >= 0

## Gradient checks at x₀ and at the energy optimum u\*

`run_load_checks.py` writes `checks.json`:

- **Tier 1** is a complex step (`h = 1e-30`) through the code's own value
  path, reported as the worst mixed error `|estimate - partial| / max(1, |partial|)`.
- **Tier 2** compares the forward-mode tangent with the adjoint, and checks
  the adjoint assembly identity.
- **Tier 3** compares the constraint Jacobian with central finite differences
  at `h*`, `h*/sqrt(10)` and `h* sqrt(10)`, with `eps_j` re-measured as the
  local three-step jitter, alongside the measured round-off floor.

| point | Tier 1 `m` (worst) | Tier 1 KS (worst) | Tier 2 tangent | Tier 2 assembly | Tier 3 worst `|diff|/eps` | Taylor minimum ratio |
|---|---|---|---|---|---|---|
| x₀ | 3.19e-14 (`dm_dphi`) | 8.33e-16 (`dKS_dx`) | 8.9e-16 | 0.0 | **0.396** (`twist_0`) | 99.0 |
| u\* | 3.62e-14 (`dm_dphi`) | 1.17e-15 (`dKS_dx`) | 8.8e-16 | 0.0 | **22.64** (`chord_1`), see below | 94.8 |

(The Taylor column is the minimum over three random directions;
`tests/test_loads.py` uses one direction at x₀ and gets 100.4.)

Tier 3 at x₀, for all ten variables at `h* = 3.162e-6`. `|diff|` and the
floor are in units of KS per unit `u`. The floor is `delta_g/h*`, from nine
samples of `g` along `u + t e_j` with `t = -4e-12 .. 4e-12`:

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

At u\*, the disagreement is 0.03 to 4.8 times the measured floor for every
variable. One variable exceeds 3 eps_j:

| `chord_1` at u\* | value |
|---|---|
| `|adj - fd|` | 7.74e-12 |
| `eps_j` (the acceptance scale) | **3.42e-13** (the smallest of the ten, about two decades below the median) |
| `|diff|/eps_j` | **22.64** (limit 3) |
| measured round-off floor `delta_g/h*` | 4.41e-11 |
| `|diff|` / floor | **0.18** |
| Taylor remainder ratio at u\* (minimum of 3) | 94.8 (limit 30) |

This is the same weakness `verification/gradient_verification/README.md`
describes for the objective at x₀. The three-step flatness estimate `eps_j`
lands on a degenerate pair of neighbouring steps and comes out below the true
round-off floor. The adjoint isn't at fault: the disagreement is *below* the
measured floor, it's among the smallest of the ten in absolute terms, and the
whole-chain Taylor test passes at 94.8. I didn't change any tolerance or
estimator. `checks.json` records `passes: false` and `floor_limited: true` for
this point and lists it under `tier3_failures`, so the failure stays visible.
It doesn't affect the constrained runs, which start at x₀, where the Jacobian
passes with a worst ratio of 0.396.

## The constrained optima

`result_eps{E}.json` for `E = 0, 0.02, 0.05, 0.1`. SLSQP (`ftol = 1e-8`,
`maxiter = 200`) from x₀ exits normally at every `eps`. At the design
condition, x₀ has `M = 177.376 N m` and 517.51 N thrust; the unconstrained
u\* has 177.919 N m (+0.306 %) and 526.66 N (+1.77 %).

| `eps` | AEP [MWh/yr] | vs x₀ | cost vs u\* | `KS0 -> KS*` | KS reduction | rated M of x_c [N m] | M vs x₀ | thrust of x_c [N] | thrust vs x₀ | active bounds | moment slack | KKT residual | `nit` | moment solves | time |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 10.262719 | **+0.1465 %** | 0.00017 % | 1.000436 -> 1.000436 | 0.000 % | 177.411 | +0.020 % | 525.49 | +1.54 % | `chord_0` upper | 5.1e-8 | 9.6e-6 | 34 | 34 | 28 s |
| 0.02 | 10.260624 | +0.1260 % | 0.0206 % | -> 0.980427 | 2.000 % | 173.833 | −1.997 % | 516.13 | −0.27 % | `chord_0` upper, `chord_4` lower | 2.2e-8 | 4.3e-6 | 39 | 39 | 33 s |
| 0.05 | 10.249582 | +0.0183 % | 0.1282 % | -> 0.950414 | 5.000 % | 168.441 | −5.037 % | 501.99 | −3.00 % | same | 1.2e-8 | 6.0e-6 | 37 | 38 | 30 s |
| 0.10 | 10.203370 | −0.4327 % | **0.5785 %** | -> 0.900392 | 10.000 % | 159.480 | −10.089 % | 478.52 | −7.54 % | same | 2.0e-9 | 6.0e-5 | 29 | 31 | 24 s |

In every run the constraint is **active** (the KS reduction equals `eps KS0`
to five digits) and `chord_0` sits on its upper bound. The cost relative to
the unconstrained optimum is 0.00017 / 0.0206 / 0.1282 / 0.5785 %, inside the
ranges I set out beforehand (+0.10 to +0.14666 % gain at `eps = 0`, 0.3–3 %
cost at `eps = 0.10`) and far below the 5 % level that would suggest a bug.
The `eps = 0` result shows that the unconstrained optimum's extra +0.285 % KS
bought almost no energy: capping the load at the reference level costs
0.0002 % of the +0.1467 % gain.

At `eps = 0`, only `chord_0`'s upper bound is active. No envelope or solidity
row is active: the tightest envelope row (`floor r = 1.9660 m`) has 16.6 mm
slack, and the largest local solidity is 0.409 against the 0.5 limit.

**KKT check and the shadow price.** SLSQP doesn't report multipliers, so
`kkt_report` estimates them as the least-squares coefficients of the
objective gradient on the active constraint and bound normals (each written
so that a minimiser's multiplier is ≥ 0). At `eps = 0` the projection has norm
2.46e-3 and the residual after it is 9.6e-6 (scaled). The multipliers are
**1.67e-03** for the moment and 4.2e-05 for `chord_0`'s upper bound, both
positive. `all_multipliers_nonnegative` is true at every `eps`, so all four
optima are genuine KKT points. The moment multiplier gives the shadow price
of the cap, `|J0| lambda` in MWh/yr per unit KS:

| `eps` | `lambda_moment` (scaled) | shadow price [MWh/yr per unit KS] | AEP per 1 % KS at the margin |
|---|---|---|---|
| 0 | 1.67e-03 | 0.0171 | 0.0017 % |
| 0.02 | 1.99e-02 | 0.2035 | 0.0199 % |
| 0.05 | 5.39e-02 | 0.5524 | 0.0539 % |
| 0.10 | 1.29e-01 | 1.3181 | 0.1286 % |

This turns "the extra moment bought almost no energy" into a measurement: at
`eps = 0` the marginal price is 0.0017 % AEP per 1 % KS, and it grows
convexly. Integrating the multiplier from `eps = 0` to 0.02 (trapezoid rule)
predicts a cost of 0.0215 %, against the measured 0.0206 %.

**Why the rated moment rises 0.020 % at `eps = 0`.** The cap is on `KS_rho`,
not on the rated-point moment. `KS0` exceeds `Mhat_max(x0) = 1` by 0.00044,
contributed by the 10.5 m/s point, and the optimiser uses some of that
margin. At x_c the 10.5 m/s moment is lower relative to the rated one, the
overestimate shrinks, and the rated moment can sit 0.020 % above x₀'s while
`KS = KS0` to within 5e-8. For the same reason, the steps of `eps` are in KS;
the rated-point moment reductions they give are 1.997 / 5.037 / 10.089 %
(the "M vs x₀" column).

## The energy–moment trade-off

`pareto.json`, `pareto.png`. Each `eps` is solved from x₀ ("cold") and from
the previous optimum ("warm"), and the pair must agree to within the
multi-start spread of 0.0166 in `u`.

| `eps` | cold AEP | warm AEP | `||u_cold - u_warm||_inf` | agree | kept |
|---|---|---|---|---|---|
| 0 | 10.262719 | 10.262719 | 0 | yes | cold |
| 0.02 | 10.260624 | 10.260624 | 2.43e-04 | yes | cold |
| 0.05 | 10.249582 | 10.249582 | 1.40e-03 | yes | cold |
| 0.10 | 10.203370 | 10.203371 | 2.58e-03 | yes | warm (higher by 1e-6 MWh/yr) |

All four agree, so none of the optima is ambiguous, and the cost rises
steadily with `eps`. The figure plots AEP cost relative to x₀ (a gain is
negative) against the reduction in rated-point moment, with x₀ at the origin.
The unconstrained u\* sits at −0.306 % reduction and −0.147 % cost, on the
wrong side of the origin: it pays a *larger* moment for its extra energy. The
constrained front runs from `eps = 0` at (−0.02 %, −0.1465 %), next to u\*, up
to `eps = 0.10` at (+10.09 %, +0.4327 %).

## What the KS aggregate actually does

At x₀, the 11 m/s rated point carries a softmax weight of **0.957**, the
10.5 m/s point 0.043, and the other seven together less than 1e-37. So at
`rho = 100` the KS value is essentially the rated-point moment. I kept the
aggregate because it gives one smooth value with a single adjoint right-hand
side, not because the lower points shape the result. At u\* the rated weight
rises to 0.978.

| `rho` | `KS(x0)` | overestimate `KS - max` | `KS(u*)` | overestimate |
|---|---|---|---|---|
| 30 | 1.0115041 | 0.0115041 | 1.0126510 | 0.0095869 |
| 100 | 1.0004357 | 0.0004357 | 1.0032871 | 0.0002230 |
| 300 | 1.0000003 | 2.944e-7 | 1.0030641 | 3.824e-8 |

Moments and weights at each load case (N m):

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

`moments_x0_xstar.png` plots the same.

## What happens to the blade

The load comes off the outer blade. Every constrained optimum puts `chord_0`
on its 0.30 m upper bound, as the unconstrained one does. As `eps` rises,
`chord_1` falls (0.182 → 0.165), `chord_2` moves from 0.097 to 0.109, `chord_4`
drops to its 0.045 m lower bound and stays there, `chord_3` falls (0.079 →
0.063 at `eps = 0.10`), and the tip twist rises from −1.94° toward −0.60°
while `twist_0` rises slightly. Control points (physical units, from
`result_eps{E}.json::x_c`):

| | c0 | c1 | c2 | c3 | c4 [m] | th0 | th1 | th2 | th3 | th4 [deg] |
|---|---|---|---|---|---|---|---|---|---|---|
| x₀ | 0.2758 | 0.1816 | 0.0971 | 0.0786 | 0.0670 | 24.53 | 9.73 | 4.04 | 1.00 | 0.59 |
| u\* | 0.3000 | 0.1796 | 0.1304 | 0.0769 | 0.0491 | 24.35 | 10.42 | 5.60 | 2.10 | −1.85 |
| `eps = 0` | 0.3000 | 0.1817 | 0.1262 | 0.0784 | 0.0477 | 24.27 | 10.63 | 5.22 | 2.23 | −1.94 |
| `eps = 0.02` | 0.3000 | 0.1809 | 0.1158 | 0.0782 | 0.0450 | 24.38 | 10.61 | 4.51 | 1.93 | −1.87 |
| `eps = 0.05` | 0.3000 | 0.1772 | 0.1096 | 0.0730 | 0.0450 | 24.60 | 10.33 | 4.34 | 1.00 | −1.19 |
| `eps = 0.10` | 0.3000 | 0.1651 | 0.1094 | 0.0626 | 0.0450 | 25.01 | 9.52 | 4.68 | −0.02 | −0.60 |

The second and third panels of `pareto.png` show the chord and twist.

## Cut-out check (reported, never constrained)

At 20 m/s and 300 rpm (`lambda = 62.831853071796 / 20 = 3.14159`):

| | Ct | thrust [N] | moment [N m] | alpha [deg] |
|---|---|---|---|---|
| x₀ | 0.2593 | 631.1 | 216.46 | 10.9 … 30.1 |
| u\* | 0.2814 | 684.9 | 228.54 | 9.6 … 29.3 |

**The model holds rated power above 11 m/s without any limiting mechanism, so
this isn't the real machine's state; the angle of attack reaches about 30°,
in the Viterna extrapolation.** It is computed through the forward path only,
never in a constraint or a gradient check. The cut-out thrust is *larger*
than at the design condition because the rotor is held at 300 rpm into
20 m/s. A real above-rated mechanism would change it.

## What would make this wrong

- **The load model.** It is the out-of-plane moment only, not rotated into
  the flap axis, with no centrifugal relief, gravity, dynamic amplification or
  in-plane component. It's the same integrand the baseline uses, so
  comparisons between x₀ and x_c are consistent, but the absolute value isn't
  a structural design load.
- **The fixed load cases.** If a design pushed a bin that is currently below
  rated over the rating, or the rated point stopped dominating the KS, the
  constraint would be measuring something different. At x₀ no point in the
  load set is power-limited (`tests/test_loads.py`), and the 0.957 rated
  weight makes the aggregate the rated point; both are checked where they
  can be.
- **Above rated.** The constraint only covers operation at or below rated by
  construction, and the cut-out case is a separate check. Modelling a real
  above-rated mechanism would mean re-running this folder.
- **The Tier 3 weakness at u\*.** The one failing variable comes from the
  acceptance scale, as explained above, and doesn't affect x₀ or any
  constrained run.
- **The `eps = 0` cap isn't a cap on the rated moment.** It caps `KS_rho`,
  and the +0.020 % rise in rated moment is the aggregate's overestimate being
  used up. Capping the rated moment itself would need `rho -> infinity` (not
  smooth), or the limit set at `Mhat_max(x0)` with the overestimate added
  back explicitly.
- **The flat energy objective.** Why the unconstrained gain is only +0.147 %
  is explained in `docs/DESIGN-BASIS.md` §2.

## Files

- `run_load_checks.py`: Tiers 1–3, the moment tables and the cut-out check.
- `checks.json`: everything above, plus the per-variable Tier 3 tables, the
  measured floors, the moments and weights per point, timings, `M_ref`, `KS0`,
  the load set and the source commit.
- `moments_x0_xstar.png`: moment and softmax weight per point at x₀ and u\*.
- `run_constrained_slsqp.py`: the constrained SLSQP run, with the iteration
  recorder, the KKT check and the shadow price.
- `result_eps{E}.json`, `iterates_eps{E}.json`, `constrained_blade_eps{E}.png`
  for `E = 0, 0.02, 0.05, 0.1`.
- `run_pareto.py`: the cold and warm `eps` sweep.
- `pareto.json`, `pareto.png`: the trade-off and the blade shapes.

The figures here predate the shared figure style and are kept as the output
of the runs that produced them.
