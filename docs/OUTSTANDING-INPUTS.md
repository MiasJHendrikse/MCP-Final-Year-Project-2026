# Outstanding external inputs

**Things only MJ can supply.** Everything on this list blocks a specific,
named exit criterion, and none of it can be filled in by inference, convention
or a sensible default — that is ground rule 3 of the work order:

> **No invented values.** The site `TODO`s in the plan (GWA coordinates,
> Weibull `k` and `c`, material allowable stress) stay `TODO` until real data
> arrives.

Each entry says exactly what is needed, where it goes, and what unblocks when
it lands. Keep this file current: an item that gets filled should be deleted
from here in the same commit that fills it.

> **When one of these arrives, read the resumption checklist first:**
> `docs/journal/Session Notes/2026-09-10.md`, section *"When the data lands"*,
> and then the 2026-09-13 entry, which records where that checklist was
> wrong when it was first used in anger — two of the three tests it named
> were not the tests that actually went red, and one of them could not have.
> This file says **what** is missing; that section says **what to do** with it
> — which config line, which tests will deliberately go red and what to replace
> them with, what evidence to regenerate, and how to tell whether the result is
> actually right rather than merely green.

Last reviewed: **2026-09-13**.

---

## 1. Wind resource — Weibull `k` and `c` ✅ **RESOLVED 2026-09-13**

Kept in place rather than deleted, because what landed is not quite what this
entry asked for and the differences are load-bearing. Nothing inherits the 🔴 — the one
follow-on it raised (section 7) was decided the same day.

**What landed:** a **GASP point-data** extraction supplied by MJ
(`misc/Screenshot 2026-09-13 114138.png`), omni-directional, at **50 m**:
`A = 8.8 m/s`, `k = 1.87`, `U = 7.8 m/s`.

**What was done with it:** extrapolated 50 m → 20 m hub height by the Justus &
Mikhail (1976) correlation, per the resumption checklist's instruction that a
height mismatch "must be done and recorded, not fudged". Recorded in
`verification/wind_resource/` — **read that README before quoting these
numbers.**

`config/site.yaml` now carries `weibull_k = 1.709226`,
`weibull_c_ms = 7.273759`, `mean_wind_speed_ms = 6.487612`, plus
`latitude_deg` and `longitude_deg`.

**Three things this did NOT resolve:**

- **`location.gwa_area` is still `TODO`.** Not an oversight. The extraction is
  a single *point*, not a Global Wind Atlas *area*, so there is no area
  selection to record; filling it would assert an extraction never performed.
- **The source is GASP, not globalwindatlas.info,** which is what plan §1.3
  names. Probably fine — same family of mesoscale atlas — but the report should
  say GASP, and MJ may want to confirm the substitution was intended.
- **The longitude sign.** GASP labelled the point `W16.558`; the pinned
  satellite views read `16°33′28.2″E`. East is recorded, west being open ocean.
  Worth one confirming glance; nothing downstream reads longitude today.

**And one thing it exposed, which is now the more interesting problem:** the
baseline AEP is **10.27 MWh/yr** against plan §1.4's 4–6 MWh/yr sanity band.
The band has **not** been widened (ground rule 5). The diagnosis is in the
2026-09-13 journal entry; the short version is that the discrepancy is mostly
*not* caused by this data — at the calmest corner of plan §1.3's own prior
resource range the rotor already returns 5.9 MWh/yr — and is an inconsistency
between plan §1.3, plan §1.4 and the fact that the AEP model produces
aerodynamic shaft energy with no drivetrain efficiency applied.

MJ decided the same day: widen the band. Done, to **8–12 MWh/yr**, with the
derivation recorded rather than the numbers merely stretched to fit. See
section 7, which is closed.

So section 1 leaves nothing blocking behind it. **The 🟠 bounds entry
(section 2) is now the most blocking item on the list.**

---

## 2. Design-variable bounds 🟠 partially blocking

**Needed for:** plan step 1.6's bounds bullet, step 1.7's feasibility check,
and step 1.8's "sweep `d_i` across its **feasible** range".

**Where it goes:** `config/rotor_design.yaml`, under `parameterisation.bounds`

| field | needed |
|---|---|
| `chord_min_m` | manufacturability / root-attachment limit |
| `chord_max_m` | manufacturability / solidity limit |
| `twist_min_deg` | physically realistic range |
| `twist_max_deg` | physically realistic range |

**Source:** plan section 7.1's manufacturability study. The plan commits only
to the fact that bounds exist and are applied to the **scaled** variables; it
selects no numbers, and no study has.

**Context that may inform them:** the Schmitz baseline runs from **263 mm chord
at the root to 69 mm at the tip**, with twist **23.07° → 0.57°** (see
`verification/representation_study/`). A root chord of 263 mm on a 2.0 m blade
is wide — characteristic of Schmitz — so whether the baseline is feasible
against a real manufacturability envelope is a genuinely open question, not a
formality.

**What unblocks:** `DesignBounds.from_config()`, which currently raises; step
1.7's "clip and record if violated"; and the sweep ranges in step 1.8.

**Interim position — updated 2026-09-13, MJ's decision:**
`tests/test_parameterisation.py::PROVISIONAL_BOUNDS` now carries a *grounded*
set for three of the four values — `chord_min_m = 0.045` (SG6043 10 % t/c at
32.1 % chord → ~4.5 mm laminate minimum), `twist −2° … 35°` (Schmitz baseline
0.57° / 23.07° with margin) — and a *provisional* `chord_max_m = 0.45`, which
depends on the hub radius and root-attachment concept, neither yet decided.
`x0` is feasible against this set (control points 0.067–0.276 m,
0.59°–24.5°; zero violations). Bounds-dependent work (plan steps 1.6–1.8,
Phase 2 SLSQP bounds) proceeds by constructing `DesignBounds(**PROVISIONAL_BOUNDS)`
explicitly. **Nothing goes into `config/` and `from_config()` keeps raising**
until `chord_max_m` has a basis — see `PROJECT_DIRECTION_v2.md` §16.

---

## 3. Buhl (2005), NREL/TP-500-36834 🟡 provenance

**Needed for:** plan step 1.3's checkbox — verify the Glauert/Buhl constants
against the original paper, *not* from secondary sources or from prior notes in
this repo.

**Where it goes:** `docs/references/` (create it), then a short commit
recording the equation number and the stated `a_c`.

**Why it is not already done:** the report could not be retrieved on
2026-09-10. `nrel.gov` and `docs.nrel.gov` do not resolve from the tooling;
`web.archive.org` is unreachable; UNT's Digital Library copy sits behind a
scripted interstitial; and OSTI's full-text link redirects to **`nlr.gov`** —
note `nlr`, not `nrel` — which is a look-alike domain and must not be used as a
citation source.

**What is established meanwhile:** internal consistency, as committed tests.
`tests/test_corrections.py` asserts C⁰ and C¹ against momentum theory at
`a = 0.4` and the F-independent anchor `Ct(1) = 2` — the three conditions that
pin all three coefficients uniquely. That would catch a transcription error but
**not** constants faithfully transcribed from the wrong source, which is
precisely what the checkbox is for. See the PROVENANCE block in
`src/bem/corrections.py`.

---

## 4. Ning (2014), *Wind Energy* 17(9), 1327–1345 🟡 provenance

**Needed for:** the claim that the solver implements "the Ning formulation".

**Where it goes:** `docs/references/`, then a commit reconciling
`bem.station.momentum_region_bracket` against Ning's own region definitions.

**Why it is not already done:** paywalled Wiley journal article, no reachable
copy.

**What is established meanwhile:** the region classification was **derived**
from the residual's structure and verified numerically over 247 stations
(`verification/phase_vi/`). The mathematics is checked to machine precision.
What is *not* established is that it matches Ning's own definitions — and
Ning's propeller-brake region on `φ ∈ (−π/4, 0)` is deliberately not
implemented, with a station that would need it reported rather than solved by
an untested branch. See the PROVENANCE block in `src/bem/station.py`.

---

## 5. Root cut-out fraction 🟢 confirm only

**Needed for:** nothing is blocked; this is a confirmation.

**Where it lives:** `config/rotor_design.yaml`, `geometry.root_fraction: 0.15`.

A conventional root cut-out for a rotor this size, and a **modelling choice**
rather than a measured or derived quantity. It is already load-bearing: it
sized the SG6043 cache's Reynolds bounds during work order Task 2, and it sets
the inboard end of the parameterised span. Promoted from a duplicated private
constant in `polars/envelope.py` into config during step 1.6 so there is one
definition rather than two.

If the real hub/root attachment geometry differs materially, the SG6043 cache
bounds and the baseline should both be re-checked.

---

## 6. NREL Phase VI experimental performance data ⚫ blocked, external

**Status:** marked `[!]` in the plan and believed unobtainable — no numeric
measured Cp–λ dataset from any primary source.

**Consequence, stated plainly:** every Phase VI number this project produces is
a prediction cross-checked against other predictions, never against
measurement. The infrastructure to plot measurements as a fourth series exists
(`src/validation/plot_bem_comparison.py`) should a dataset ever surface.

---

## 7. The AEP sanity band ✅ **RESOLVED 2026-09-13 — MJ's decision**

Raised and closed the same day. Kept in place rather than deleted because it
records a deliberate revision of a stated expectation, which this repo's ground
rule 5 otherwise forbids, and the basis for that needs to stay findable.

**The problem:** when the wind resource landed, plan §1.4's 4–6 MWh/yr AEP
sanity band could be evaluated for the first time. The baseline returns
**10.27 MWh/yr**.

**MJ's decision, as given:** widen the band — if the blade already makes
substantially more power than the original band, the original band was simply
too narrow.

**Acted on.** `config/rotor_design.yaml` now carries **8–12 MWh/yr**, with the
full derivation written out beside the numbers. Summary:

- **Why 4–6 was wrong, on two counts, both predating the data.** It compared
  an *electrical* capacity-factor estimate (~3.1 kW at `Cp = 0.42`,
  `η = 0.90`) against a model that produces *aerodynamic* shaft energy at the
  solver's actual `Cp = 0.472` with no drivetrain efficiency in the chain —
  3.82 kW rated, not 3.1 kW. And it was never consistent with plan §1.3's own
  resource prior: at the calmest corner of that prior the rotor already returns
  5.9 MWh/yr, at the middle 7.7.
- **How 8–12 was derived** — from the resource uncertainty this project has
  recorded, not by fitting to 10.27. ±10 % on `c` and `k` spans 8.21–12.01; the
  log-law roughness band on the height extrapolation spans 9.73–11.10 and sits
  inside. Rounded outward. Capacity factor 0.24–0.36 on the aerodynamic rating.
- **It is still a check.** Verified to fail on unlimited power (14.6),
  sea-level density (13.0), unnormalised bin masses (12.9) and any factor of
  two. `test_the_sanity_band_is_still_narrow_enough_to_catch_a_bug` asserts
  exactly this, so a future widening cannot quietly turn the band into a note.

**One sub-question is deliberately left open**, because it is a report
decision rather than a blocker: whether the Phase 1 exit criterion should be
stated in **aerodynamic** or **delivered electrical** energy. The band as
written is aerodynamic, matching what the model computes and labelling itself
as such. Stating it electrically means an `η` entering the AEP chain and a
config field to hold it; `10.27 × 0.90 = 9.24 MWh/yr`, which is inside the new
band either way. Nothing is blocked on this.

---

## 8. Justus & Mikhail (1976) 🟡 provenance

**Raised 2026-09-13.** The third entry of the same shape as Buhl and Ning.

**Needed for:** the claim that `src/objective/height_extrapolation.py`
implements the Justus & Mikhail correlation, and therefore for the 20 m `k`
and `c` in `config/site.yaml`.

**Where it goes:** `docs/references/`, then a commit confirming the two
constants (`0.37` and `0.0881`), the 10 m anchor height, and the form of the
`k(z)` relation against the paper.

**Why it is not already done:** the formulae were transcribed from general
knowledge of the method; no primary copy was retrieved.

**What is established meanwhile:** `tests/test_height_extrapolation.py` — the
identity at `z = z_ref`, monotonicity in `z`, the *sign* of both height
relations (the `k` relation being upside down is the most plausible
transcription error and would pass any check that looked only at the scale),
and agreement with an independent log-law calculation across the plausible
roughness range for this terrain. That would catch a mangled transcription. It
would not catch a faithful transcription of the wrong correlation, which is
what this item is for.

**Weight it carries:** more than Buhl's or Ning's. Those two verify
implementations whose behaviour is separately pinned by 247 stations and a
golden regression. This one sets two numbers that multiply straight through
every AEP figure in the report.
