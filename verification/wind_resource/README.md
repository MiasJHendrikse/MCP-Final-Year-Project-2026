# Wind resource at 20 m hub height

**This is the source of `config/site.yaml:wind_resource`.** Read it before
quoting `k`, `c` or `V̄` anywhere, because the headline fact about these
numbers is not in the config file:

> They are **extrapolated**, not extracted. The extraction is at 50 m. The
> rotor's hub is at 20 m.

Regenerate with:

```
python verification/wind_resource/run_extrapolation.py
```

which writes `wind_resource_20m.json`. `tests/test_height_extrapolation.py`
checks the artefact against its own recorded inputs, so a hand-edit is caught.

---

## What was supplied

MJ, 2026-09-13 — a **GASP point-data** panel, plus four Google Earth views
pinning the location.

**The panel is committed here as `gasp-point-data-50m.png`, and the regional
location view as `site-location-regional.png`.** They arrived in `misc/`,
which is gitignored as scratch — so the sole source of the project's headline
wind resource would not have been versioned. Copied in deliberately: an
extrapolation whose input cannot be re-read is not auditable, and the input
here is a screenshot of a tool session that cannot be re-run from a URL.

| Displayed | Value |
|---|---|
| Height | 50 m (panel selector; 100 m and 150 m also offered, 20 m not) |
| Sector | Total (omni-directional) |
| Omni A / k | **8.8 m/s / 1.87** |
| Avg. wind speed U | 7.8 m/s |
| Coordinates | S22.428, W16.558 |
| Location pin | 22°25′40.7″S 16°33′28.2″E |

**The input is self-consistent.** `A · Γ(1 + 1/k) = 7.8131 m/s` against a
displayed `U = 7.8`, a residual of +0.013 m/s — display rounding. That identity
is the evidence the panel was read correctly; a misread digit in either
parameter would move it well outside rounding.

## Three things about the source, recorded rather than smoothed over

**1. It is GASP, not the Global Wind Atlas.** Plan §1.3, `config/site.yaml` and
`docs/OUTSTANDING-INPUTS.md` all name globalwindatlas.info. The tool actually
used is different. The substitution may well be fine — it is a mesoscale wind
atlas of the same family, and the panel's other outputs (IEC class I, a
sensible 1.0 MW reference power curve, a plausible sector rose) are coherent —
but it is not the thing the plan specified, and the report should say GASP.

**2. `location.gwa_area` is still `TODO`, and deliberately.** This is a single
**point** extraction. There is no area selection to record. Filling that field
with the point's coordinates would assert an extraction that was never
performed, which is what ground rule 3 forbids. The field stays unresolved and
`tests/test_config.py` guards that it still behaves like a sentinel.

**3. The longitude sign disagrees with itself.** GASP labels the point
`W16.558`; the pinned satellite views read `16°33′28.2″E`. West 16.558 is open
ocean several hundred kilometres off the Namibian coast, where the panel's own
terrain-derived statistics could not have been computed. **East is recorded.**
Treated as a display artefact of that tool. Worth one confirming glance from
MJ; nothing downstream reads longitude today.

---

## The extrapolation

**Method: Justus & Mikhail (1976)**, implemented in
`src/objective/height_extrapolation.py`.

```
n    = (0.37 − 0.0881 ln A_ref) / (1 − 0.0881 ln(z_ref/10))
A(z) = A_ref · (z/z_ref)^n
k(z) = k_ref · (1 − 0.0881 ln(z_ref/10)) / (1 − 0.0881 ln(z/10))
```

Chosen over the two obvious alternatives for two specific reasons:

- **It moves `k` as well as `A`.** A bare power law on `A` alone implicitly
  asserts the distribution keeps its shape with height. It does not — the
  distribution broadens downward, `k` falling 1.87 → 1.709 over these 30 m.
  Holding `k` at a value read 30 m above hub height biases the AEP integral in
  a direction nothing downstream could detect. Measured here it is **−0.90 %**
  (10.270 → 10.178 MWh/yr with `k` pinned at 1.87) — smaller than one might
  guess, and worth removing regardless, since the point is that the shape is
  computed rather than assumed.
- **It needs no surface roughness.** A log-law central case would require a
  `z₀` for this terrain, and no roughness survey exists. Inventing one is
  exactly what ground rule 3 forbids. `z₀` appears below only as a *range*, in
  the cross-check, where a range is honest.

| | 50 m (extraction) | **20 m (hub, used)** |
|---|---|---|
| `k` | 1.87 | **1.709226** |
| `A` = `c` [m/s] | 8.8 | **7.273759** |
| `V̄` [m/s] | 7.813 | **6.487612** |

Shear exponent `n = 0.207880`.

## Cross-check: the logarithmic profile

An independent physical basis — surface-layer similarity rather than an
empirical correlation — across plausible roughness lengths for Khomas Hochland
open thornbush savanna:

| `z₀` [m] | `A(20 m)` [m/s] |
|---|---|
| 0.05 | 7.633 |
| 0.10 | 7.503 |
| 0.20 | 7.343 |
| 0.30 | 7.224 |
| 0.50 | 7.049 |

Band **7.049 – 7.633 m/s**; the central case **7.274 m/s sits inside it**,
around `z₀ ≈ 0.25 m`, which is a reasonable value for this terrain.

This is a check, not a tolerance. The two methods are not required to agree
exactly. A central case falling outside the band across the whole plausible
roughness range would mean one of them is wrong, and that is worth knowing
before any AEP figure is quoted.

## Against the plan's prior expectation — it disagrees, on both

Plan §1.3 records `k ≈ 1.8–2.4, c ≈ 6–7 m/s` explicitly "for sanity-checking
an extraction only".

| | prior | result | |
|---|---|---|---|
| `k` | 1.8 – 2.4 | **1.709** | 0.09 **below** |
| `c` | 6 – 7 m/s | **7.274** | 0.27 m/s **above** |

The site is windier, and its distribution broader, than the plan guessed. Per
ground rule 5 the band is **neither widened to admit the result nor used to
adjust it**. `tests/test_height_extrapolation.py` pins the disagreement so it
stays visible.

Both directions are physically coherent for this site: a 1 800 m highland ridge
is windy, and 20 m is low enough into the surface layer for shear to broaden
the distribution noticeably.

## What this does not establish

**Provenance.** The Justus & Mikhail formulae are transcribed from general
knowledge, not from a retrieved copy of the 1976 paper — the same position the
repo is in for Buhl and Ning. The tests establish internal consistency,
physical behaviour and log-law agreement, which would catch a *mangled*
transcription. They cannot catch a *faithful* transcription of the wrong
correlation. Open as item 8 in `docs/OUTSTANDING-INPUTS.md`.

**Terrain representativeness.** GASP, like GWA, runs on a mesoscale grid. The
Khomas Hochland is complex terrain and a single point carries real uncertainty
— plan §1.3 already commits to "a central case plus a sensitivity band" with
the optimisation re-run across it. This artefact is the central case. The band
is not yet chosen.

**That the resource explains the AEP sanity band's failure.** It did not, or
not mostly. The baseline returns 10.27 MWh/yr against what was a 4–6 MWh/yr
band — but at the *calmest* corner of the plan's own prior resource range this
rotor already returns 5.9, and at the middle 7.7. The inconsistency was between
plan §1.3 and plan §1.4 and predated this data. The band was revised to
**8–12 MWh/yr** on 2026-09-13 on MJ's decision; the derivation is in
`config/rotor_design.yaml` and the reasoning in the journal entry.
