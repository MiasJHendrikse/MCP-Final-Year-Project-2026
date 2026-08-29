# SG6043 wind-tunnel data — UIUC Low-Speed Airfoil Tests, Volume 3

The measured reference the SG6043 polar cache has to be validated against
(PROJECT_PLAN.md §1.2, §3.1). This is the experimental half of the pair: the
coordinates it goes with are `data/airfoils/sg6043.dat`, and the cache built
from them is specified in `config/polars_sg6043.yaml`.

Same role for SG6043 that `data/naca0012_validated/` plays for NACA 0012 —
an independent measurement to check XFOIL against, rather than another thing
XFOIL produced.

## Source and provenance

Downloaded 2026-08-29 from the UIUC Applied Aerodynamics Group's LSATs data
distribution:

- `https://m-selig.ae.illinois.edu/pd/pub/lsat/volume03/DRAG03.TXT`
  — sha256 `5e493c472975e5a84b26889b856a1c7bd37103ba33a43ec7b3281c4f831b094d`
- `https://m-selig.ae.illinois.edu/pd/pub/lsat/volume03/LIFT03.TXT`
  — sha256 `088d7aa0c902b3c173b035fa4c8b2c9271ae00a4350004eb70c78edeed6ef4f5`

Those two files hold every airfoil in Volume 3. What is committed here are the
**SG6043 blocks copied out byte-for-byte**, each keeping its original
`::::::::::::::` header so it can be located in the parent file, and each
keeping the run provenance lines (`Tabulated from data in file pg02397.dat …`)
the distribution attaches. Nothing was re-tabulated, re-interpolated,
re-formatted or rounded — checking a file here against the parent download is
a plain byte comparison.

One caveat on that comparison: the distribution uses LF line endings, and this
repository has no `.gitattributes`, so a checkout on Windows may hand you CRLF.
Compare with `git show HEAD:<path>` (which gives you what is stored), or
normalise line endings first, before concluding that a byte differs.

## Citation

> Lyon, C. A., Broeren, A. P., Giguère, P., Gopalarathnam, A., and Selig,
> M. S., *Summary of Low-Speed Airfoil Data, Volume 3*, SoarTech Publications,
> Virginia Beach, VA, March 1998. ISBN 0-9646747-3-4.
> Available from SoarTech Publications, c/o Herk Stokely, 1504 N. Horseshoe
> Circle, Virginia Beach, VA 23451.

**This data was produced under the UIUC Low-Speed Airfoil Tests program.**
That statement, and the citation above with the address, are required by the
distribution terms whenever the data is used or passed on — see `README03.TXT`
(copyright notice), `GPL.TXT` (licence) and `MANIFEST.TXT` (programme
manifesto), all included here verbatim for that reason. The report must carry
the citation too.

## Files

| File | Contents |
|---|---|
| `SG6043_drag_clean.txt` | **The primary reference.** Clean model, α / Cl / Cd, 9 runs. |
| `SG6043_lift_clean.txt` | Clean model, α / Cl at fine α resolution, 6 runs to ±20°. |
| `SG6043_drag_trip_typeC.txt` | Zigzag trip type C, h = 0.023″ at 2 % u.s. / 5 % l.s. |
| `SG6043_drag_trip_Re200k.txt` | Re = 200k: clean vs plain trips at x/c = 5…45 %. |
| `SG6043_drag_trip_Re300k.txt` | Re = 300k: same comparison. |
| `SG6043_lift_trip_typeC.txt` | Lift curves for the type C trip configuration. |
| `FORMAT03.TXT` | The distribution's own description of the file format. |
| `README03.TXT`, `GPL.TXT`, `MANIFEST.TXT` | Copyright notice, licence, manifesto. |

The trip files are here because boundary-layer tripping and `n_crit` are two
descriptions of the same physics — an `n_crit` choice is a statement about
freestream turbulence and transition location, and the trip data is measured
evidence about what moving transition does to this section. They are context
for the Step 2 study, not the baseline it validates against.

## What the clean data actually covers — read this before using it

**Re = 100k, 150k, 200k, 300k, 400k and 500k.** Not 60k.

`PROJECT_PLAN.md` §3.1 previously said the section was tested "at Re = 60k,
100k, 200k and 300k", and §1.4's justification for R = 2.0 m rested on
bracketing those points. Both halves were wrong in the same direction: **there
is no 60k run**, and there are two more test points at the top (400k, 500k)
than the plan credited. Volumes 1, 2 and 4 were checked — SG6043 appears only
in Volume 3, so the six runs above are the complete set. The plan's
enumerations were corrected in the same commit that added this directory; the
report still needs to follow.

The R = 2.0 m choice comes out *better* justified than argued — the envelope's
upper half is measured, not extrapolated — but the low-Re anchor the plan
leaned on does not exist.

The consequence for `config/polars_sg6043.yaml`: that cache is specified down
to **40k**, because the design rotor's computed envelope reaches 49k at the
root at cut-in. Nothing below 100k has experimental support. Those rows are
XFOIL alone, in the regime where XFOIL is least trustworthy, and the `n_crit`
study cannot validate them — it can only validate 100k–500k and assume the
calibration carries downward. That is a limitation to name in the report, not
one to discover later.

Summary of the clean drag runs, derived from the committed files:

| Nominal Re | α range | Cl_max | Cd_min | max L/D |
|---|---|---|---|---|
| 100k | −4.3 … +11.3° | 1.432 | 0.0218 | 59.4 |
| 150k | −4.0 … +12.4° | 1.538 | 0.0173 | 74.2 |
| 200k | −5.0 … +11.4° | 1.542 | 0.0143 | 86.6 |
| 300k | −6.2 … +13.3° | 1.607 | 0.0106 | 105.3 |
| 400k | −6.3 … +12.4° | 1.614 | 0.0088 | 118.0 |
| 500k | −6.3 … +12.3° | 1.627 | 0.0080 | 125.1 |

Cl_max rises and Cd_min falls monotonically with Re, which is the trend
`validate_polars.py` check 5 will look for in the cache built to match this.

Note that the nine drag runs are **not** nine Reynolds numbers: several are
short continuation runs at a nominal Re already listed (one is a single α
point), filling in a range the main run did not reach. Group by nominal Re
before comparing against anything, or a two-point run will look like a curve.

## How this gets used

1. **The `n_crit` sensitivity study** (plan §1.2, and the outstanding
   `build.ncrit` TODO in `config/polars_sg6043.yaml`). Sweep candidate `n_crit`
   values in XFOIL at these six Reynolds numbers, compare against the clean
   data, and *select* one — the way S809's `Ncrit = 5` was selected against
   published behaviour rather than assumed.
2. **Naming the low-Re uncertainty** (plan §3.3, risk R4). The disagreement
   that remains at the chosen `n_crit` is the measurement of how much to trust
   the cache, and it belongs in the limitations chapter with a number attached.
