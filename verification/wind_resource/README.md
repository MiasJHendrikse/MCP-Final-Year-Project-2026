# Wind resource at the 20 m hub height

This is where the values in `config/site.yaml:wind_resource` come from. The
most important thing about them isn't in the config file:

> They are **extrapolated**, not measured at hub height. The source data is at
> 50 m, and the rotor's hub is at 20 m.

To regenerate:

```
python verification/wind_resource/run_extrapolation.py
```

This writes `wind_resource_20m.json`. `tests/test_height_extrapolation.py`
checks it against its own recorded inputs, so a hand edit would be caught.

---

## The source data

A **GASP point-data** panel for the site, plus four Google Earth views pinning
the location. The panel is saved here as `gasp-point-data-50m.png` and the
regional location view as `site-location-regional.png`. I committed them on
purpose: the panel is a screenshot of a tool session that can't be re-created
from a URL, and an extrapolation whose input can't be re-read can't be checked.

| Shown | Value |
|---|---|
| Height | 50 m (100 m and 150 m were also offered; 20 m wasn't) |
| Sector | Total (all directions) |
| Omni A / k | **8.8 m/s / 1.87** |
| Mean wind speed U | 7.8 m/s |
| Coordinates | S22.428, W16.558 |
| Location pin | 22°25′40.7″S 16°33′28.2″E |

**The values are self-consistent.** `A · Γ(1 + 1/k) = 7.8131 m/s`, against the
displayed U of 7.8, a difference of 0.013 m/s that is just display rounding. A
misread digit in either parameter would move it well outside that, so this
confirms the panel was read correctly.

## Three things to note about the source

**1. It is GASP, not the Global Wind Atlas.** The project plan originally
named globalwindatlas.info, but the tool used was GASP. That's probably fine,
since it is a mesoscale wind atlas of the same kind, and the panel's other
outputs (IEC class I, a sensible 1.0 MW reference power curve, a plausible
wind rose) are coherent. But the report should name GASP.

**2. `location.gwa_area` stays unresolved, deliberately.** This is a single
**point** extraction, so there is no area selection to record. Filling the
field with the point's coordinates would claim an extraction that never
happened. `tests/test_config.py` checks that the field still behaves as
unresolved.

**3. The longitude sign contradicts itself.** GASP labels the point `W16.558`,
but the satellite views read `16°33′28.2″E`. West 16.558° is open ocean several
hundred kilometres off the Namibian coast, where the panel's terrain-based
statistics couldn't have been computed. **East is recorded**, and I treated
the W as a display error in the tool. Nothing downstream uses the longitude.

---

## The extrapolation

**Method: Justus & Mikhail (1976)**, implemented in
`src/objective/height_extrapolation.py`:

```
n    = (0.37 − 0.0881 ln A_ref) / (1 − 0.0881 ln(z_ref/10))
A(z) = A_ref · (z/z_ref)^n
k(z) = k_ref · (1 − 0.0881 ln(z_ref/10)) / (1 − 0.0881 ln(z/10))
```

I chose it over the two obvious alternatives for two reasons:

- **It adjusts `k` as well as `A`.** A plain power law on `A` alone assumes the
  distribution keeps its shape with height, but it doesn't: the distribution
  broadens closer to the ground, and `k` falls from 1.87 to 1.709 over these
  30 m. Keeping `k` at its 50 m value would bias the AEP in a way nothing
  downstream could detect. Here the effect is **−0.91 %** (10.2477 →
  10.1542 MWh/yr with `k` held at 1.87; `wind_resource_20m.json →
  held_k_penalty`, for the reference blade under the current operating law).
  That's smaller than I'd have guessed, but it's still worth getting right. An
  earlier version of the model gave −0.90 % (10.270 → 10.178).
- **It doesn't need a surface roughness.** A log-law central case would need a
  roughness length `z₀` for this terrain, and there's no survey to take one
  from. So `z₀` only appears below as a *range*, in the cross-check, where a
  range is the honest thing to use.

| | 50 m (source) | **20 m (hub, used)** |
|---|---|---|
| `k` | 1.87 | **1.709226** |
| `A` = `c` [m/s] | 8.8 | **7.273759** |
| `V̄` [m/s] | 7.813 | **6.487612** |

The shear exponent is `n = 0.207880`.

## Cross-check with the logarithmic profile

This uses an independent physical basis (surface-layer similarity rather than
an empirical correlation) over plausible roughness lengths for open thornbush
savanna in the Khomas Hochland:

| `z₀` [m] | `A(20 m)` [m/s] |
|---|---|
| 0.05 | 7.633 |
| 0.10 | 7.503 |
| 0.20 | 7.343 |
| 0.30 | 7.224 |
| 0.50 | 7.049 |

The range is **7.049–7.633 m/s**, and the central value of **7.274 m/s is
inside it**, at about `z₀ ≈ 0.25 m`, which is reasonable for this terrain.

This is a check, not a tolerance, and the two methods don't have to agree
exactly. But if the central value fell outside the range for every plausible
roughness, one of them would be wrong, and I'd want to know that before
quoting any AEP.

## Compared with my prior expectation, it's different on both counts

Before getting the data I expected roughly `k ≈ 1.8–2.4` and `c ≈ 6–7 m/s`,
intended only as a sanity check on the extraction:

| | expected | result | |
|---|---|---|---|
| `k` | 1.8–2.4 | **1.709** | 0.09 **below** |
| `c` | 6–7 m/s | **7.274** | 0.27 m/s **above** |

The site is windier, with a broader distribution, than I guessed. I didn't
widen the expected range to fit the result, or use it to adjust the result;
`tests/test_height_extrapolation.py` pins the difference so it stays visible.
Both directions make physical sense for this site: a ridge at 1,800 m is
windy, and 20 m is low enough in the surface layer for shear to broaden the
distribution noticeably.

## What this doesn't establish

**The source of the formula.** The Justus & Mikhail equations are
transcribed from general knowledge, not from a copy of the 1976 paper (the
same situation as for Buhl and Ning). The tests check internal consistency,
physical behaviour and agreement with the log law, which would catch a garbled
transcription. They can't catch a faithful transcription of the wrong
correlation. This is item 8 in `docs/OUTSTANDING-INPUTS.md`.

**How representative one point is.** GASP, like the Global Wind Atlas, uses a
mesoscale grid. The Khomas Hochland is complex terrain, so a single point
carries real uncertainty. This folder gives the central case. How much the
final blade depends on it is checked by evaluating the blades at four
resource corners (`verification/mass_optimisation/cross_evaluation.json`).

**The AEP sanity band.** The reference blade gives 10.27 MWh/yr, while the
first sanity band I set was 4–6 MWh/yr. The resource isn't the main reason for
that mismatch: even at the calmest end of my prior range the rotor gives 5.9,
and at the middle 7.7. The band itself was inconsistent with the rotor size, so
I revised it to **8–12 MWh/yr**, based on the resource uncertainty; the
derivation is in `config/rotor_design.yaml` and `docs/OUTSTANDING-INPUTS.md`
§7.
