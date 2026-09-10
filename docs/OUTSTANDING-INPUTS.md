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

Last reviewed: **2026-09-10**.

---

## 1. Wind resource — Weibull `k` and `c` 🔴 blocking

**Needed for:** plan step 1.5 (AEP model) and therefore step 1.8 (the objective
smoothness gate), because the plan defines the objective as

```
J(d) = −AEP(d) = −T · ∫ P(V; d) · f_Weibull(V; k, c) dV
```

Without `k` and `c` there is no `J`.

**Where it goes:** `config/site.yaml`

| field | current | needed |
|---|---|---|
| `wind.weibull_k` | `TODO: from GWA at 20 m hub height (plan 1.3)` | shape parameter, dimensionless |
| `wind.weibull_c_ms` | `TODO: from GWA at 20 m hub height (plan 1.3)` | scale parameter, m/s |
| `wind.mean_wind_speed_ms` | `TODO: derived, V̄ = c·Γ(1 + 1/k)` | derived — the loader can compute it, but record the GWA value too so they can be cross-checked |

**Source:** Global Wind Atlas extraction at the site, at **20 m hub height**
(`site.hub_height_m`). Also fills `site.gwa_area`, `latitude_deg` and
`longitude_deg`, which are separately `TODO`.

**What unblocks:** step 1.5's exit criterion "baseline geometry returns
4–6 MWh/yr" (the sanity band is already in `config/rotor_design.yaml`), and
step 1.8's `J` sweeps against the real objective rather than the surrogate.

**Interim position:** the AEP machinery is built as mechanism and raises
`UnresolvedConfigError` naming the missing field. Step 1.8 runs on a
**unit-weighted energy surrogate** — the same bins with equal weights instead
of Weibull weights. Because the Weibull weights are a fixed convex combination
independent of `d`, the surrogate exercises the identical numerical chain
(parameterisation → geometry → BEM → power) and reveals the same staircasing,
kinks and noise. Re-running against the real objective costs ~13 minutes once
the numbers land.

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
