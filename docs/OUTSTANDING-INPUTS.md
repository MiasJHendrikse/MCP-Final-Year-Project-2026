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
entry asked for and the differences are load-bearing. The next entry down
inherits the 🔴 blocking status.

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

**This needs MJ's decision** and is the successor to this entry. See section 7.

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

**Interim position:** `DesignBounds(...)` accepts explicit values, so studies
can state their own provisional range locally without it entering `config/` and
silently constraining every later result. `tests/test_parameterisation.py` has
a `PROVISIONAL_BOUNDS` constant used for exactly this, deliberately kept out of
the configuration.

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

## 7. The AEP sanity band 🔴 blocking — **MJ's decision**

**Raised 2026-09-13**, when the wind resource landed and the band could be
evaluated for the first time. This is the successor to section 1 and is now the
item on the critical path.

**The discrepancy:** baseline AEP is **10.27 MWh/yr**. Plan §1.4's band is
**4–6 MWh/yr**. The band is unchanged and `evaluate_baseline` reports
`aep_in_sanity_band: false`; `tests/test_baseline.py` pins the failure so it
cannot be quietly tuned into agreement.

**It is mostly not the new data's fault.** Sweeping the plan's own prior
resource range from §1.3:

| resource | `V̄` [m/s] | AEP [MWh/yr] |
|---|---|---|
| plan prior, calmest (`k=2.4, c=6.0`) | 5.32 | 5.88 |
| plan prior, mid (`k=2.1, c=6.5`) | 5.76 | 7.74 |
| plan prior, windiest (`k=1.8, c=7.0`) | 6.23 | 9.53 |
| **site, 20 m (actual)** | **6.49** | **10.27** |

Only the *calmest corner* of the plan's own prior expectation lands inside the
band at all. Plan §1.3 and plan §1.4 were never consistent with each other for
this rotor, and the 2026-09-13 extraction merely made that visible.

**Two identified contributors, neither of which may be fixed unilaterally:**

1. **Aerodynamic vs electrical energy.** The band comes from a 0.15–0.22
   capacity factor on the ~3.1 kW rating in plan §1.4, which was computed with
   `Cp = 0.42` and `η = 0.90`. The AEP model applies **no drivetrain
   efficiency** — `src/objective/power.py` says so explicitly and always has —
   and the solver's actual design-point `Cp` is **0.472**, not 0.42. Rated
   aerodynamic power is therefore 3.82 kW, not 3.1 kW.
2. **The rotor is deliberately oversized for its rating.** 2.0 m radius with
   cut-in at 3.0 m/s is a low-wind machine by design (plan §1.5 fixes `R` on
   polar-validation grounds, explicitly *not* to hit a power target). A
   capacity factor typical of a generic small turbine is the wrong prior for
   it.

**What MJ has to decide** — this is a plan question, not a code question, and
nothing should be changed until it is answered:

- Is the exit criterion about **aerodynamic** or **delivered electrical**
  energy? If electrical, the AEP chain needs an `η` and the config needs a
  field for it; `10.27 × 0.90 = 9.24 MWh/yr`, still outside the band.
- Should the §1.4 band be **recomputed** from the rotor as actually specified
  and the resource as actually measured, rather than from the round numbers it
  was estimated with? A band derived from `Cp = 0.472`, `η` as decided, and
  `V̄ = 6.49 m/s` would be a real check; the present one is not.
- Or is the **rotor** wrong for the site — is a 3.1 kW rating simply under-set
  for a 2.0 m rotor at `V̄ = 6.49 m/s`?

**What must not happen:** the band widened to 4–11 so the test goes green.
That converts a genuine finding into a silent assumption, which is ground
rule 5's whole subject.

**Interim position:** the number is computed, reported, and flagged as
out-of-band everywhere it appears — `baseline_reference.json`, the generator's
stdout, and a test. Phase 2 is *not* blocked on the answer, because gradients
do not care about the band's absolute value; the report's Phase 1 exit
criterion is.

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
