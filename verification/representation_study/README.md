# How many control points does the blade need?

The blade's chord and twist are B-splines, and the number of control points is
a design choice. Rather than just pick one, I compared 6, 8, 10, 12 and 16
control points in total.

**Conclusion: 10 in total, 5 for chord and 5 for twist.** The reasons are
below, including why the count that fits best isn't the one I chose.

```
python verification/representation_study/run_study.py
```

This writes `representation_study.json` and `representation_study.png`.

The study only compares geometry: how well each number of control points can
represent the analytic Schmitz chord and twist. The bounds and operating law
don't enter it, and the Schmitz blade it uses hasn't changed, so it didn't
need re-running when those changed.

## Method

The target is the **analytic Schmitz blade** for this rotor (R = 2.0 m, B = 3,
λ = 6.5) at the SG6043 best-L/D point (L/D = 98.5 at α = 5.36°, Cl = 1.2825,
Re = 200,000). That is the right target because it's the shape the optimiser
starts from. What matters is whether the spline can represent *that* blade,
not some generic taper chosen to be easy to fit.

The blade runs from a 263 mm chord at the root to 69 mm at the tip (mean
132 mm), with twist from 23.07° to 0.57°.

The parameterisation is linear in its control points, so fitting is a single
linear least-squares solve: exact, with no starting guess and no local minima.
The error reported is therefore a property of the **basis**, not of how well an
optimiser happened to fit it.

## Results

| total | chord + twist | degree | chord RMS | (% of mean) | twist RMS | cond(N) | sawtooth reproduced |
|---|---|---|---|---|---|---|---|
| 6 | 3 + 3 | 2 | 4.112 mm | 3.111 % | 0.8626° | 3.17 | 9.8 % |
| 8 | 4 + 4 | 3 | 0.642 mm | 0.486 % | 0.2926° | 5.98 | 9.8 % |
| **10** | **5 + 5** | **3** | **0.515 mm** | **0.389 %** | **0.1061°** | **5.59** | **14.5 %** |
| 12 | 6 + 6 | 3 | 0.427 mm | 0.323 % | 0.0301° | 5.43 | 15.6 % |
| 16 | 8 + 8 | 3 | 0.122 mm | 0.093 % | 0.0041° | 5.60 | 20.7 % |

"Sawtooth reproduced" is the fraction of a pure strip-to-strip oscillation the
basis can represent. At 0 %, the optimiser can't see that kind of unphysical
zig-zag at all; 100 % is what you'd get with one design variable per station.
Keeping this low is one of the reasons for using a spline in the first place.

## Why 10

**6 is ruled out.** Three control points per distribution can't carry a
cubic, so the degree drops to 2, and it shows: 3.1 % chord error, 0.86° twist
error, and a visible miss at both root and tip in the left panel of the
figure. It can't represent the shape it would have to start from.

**Beyond 8, better fitting stops meaning anything.** At 10 control points the
chord error is **0.515 mm RMS on a 132 mm mean chord**. No composite blade of
this size is made to anything like that tolerance, so the extra accuracy of 12
(0.427 mm) or 16 (0.122 mm) buys nothing physical. The error curve in the
middle panel is steep from 6 to 8 and flat from 8 to 12, which shows the basis
has already captured the shape.

**Conditioning doesn't help decide.** cond(N) is between 5.4 and 6.0 for
every count from 8 up, because the basis functions are well separated at all
these resolutions on 25 strips. The lower value at 6 is a side effect of the
lower degree, not an advantage. I expected conditioning to be the deciding
factor, but it isn't, so the choice rests on the other criteria.

**16 fits best, and I still didn't choose it**, for two reasons. It
reproduces **20.7 %** of a strip-to-strip oscillation against 14.5 % at 10,
handing the optimiser noticeably more of the zig-zag freedom the spline is
there to remove. And it costs 60 % more chain-rule work through the
parameterisation in every adjoint evaluation, for a fitting improvement
already far below manufacturing tolerance.

**10 is at the knee.** Twist error at 10 (0.106°) is an order of magnitude
better than at 8 (0.293°) for one more control point per distribution, while
chord error barely changes. Twist is the harder of the two to represent here,
and 10 is where its error becomes negligible.

## What this study doesn't settle

- **Whether the optimisation result depends on the count.** That needs full
  optimisations at other counts, which I haven't done. If an optimum with 12
  control points turned out to be materially better than with 10, this choice
  should be revisited; the study itself is cheap to re-run.
  (`verification/cost_scaling/` shows that the gradient cost doesn't depend on
  the count, so cost isn't a reason to keep it low.)
- **How the count interacts with the bounds.** The bounds were set after this
  study (`docs/OUTSTANDING-INPUTS.md` §2). A tight limit on the root chord
  might favour more control points inboard, or an uneven split between chord
  and twist rather than the even one assumed here.
- **An uneven chord/twist split.** Every row splits the total evenly. Since
  twist error dominates at low counts, 4 chord + 6 twist might beat 5 + 5 at
  the same cost. I haven't explored that.
