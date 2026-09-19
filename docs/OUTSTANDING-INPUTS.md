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

Last reviewed: **2026-09-19 (evening)**. Phases 1–3 are complete. Since the
14th: the generator rating is **frozen in code** (audit recommendation 2,
section 9) at a provisional, baseline-derived value, so B2 asks for a number
to *replace* a placeholder. **On the 19th MJ decided B1 (maximum rotor
speed, 300 rpm, provisional with a stated basis) and the bounds (section 2,
`chord_max_m = 0.30 m` with a solidity cap)** on the evidence of the
isolated experiment now recorded in `verification/aep_optimisation_experiment/`.
Both are in `config/rotor_design.yaml`. Nothing on this list blocks Phase 4
or Phase 5 any more; B2 and B3 remain open as *replacements* and *scope
limits* (section 9), and section 10 is a confirmation.

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

So section 1 leaves nothing blocking behind it. The bounds entry (section
2) was the most blocking item until 2026-09-13's provisional set unblocked
Phases 2–3, and was resolved on 2026-09-19; section 9 carries the machine
facts, of which B1 was resolved provisionally the same day.

---

## 2. Design-variable bounds ✅ **RESOLVED 2026-09-19 — MJ's decision (O4)**

Kept in place rather than deleted, because the fourth value is an
engineering judgement with a recorded basis, not a measured hub/attachment
limit, and the basis has to stay findable.

**What landed:** `config/rotor_design.yaml`, `parameterisation.bounds`:

| field | value | basis |
|---|---|---|
| `chord_min_m` | 0.045 m | SG6043 10 % t/c at 32.1 % chord → ~4.5 mm laminate minimum (2026-09-13) |
| `chord_max_m` | **0.30 m** | local solidity at the root cut-out r = 0.3 m: 0.45 m gives σ = 0.72 (blades nearly touching; BEM's independent-annuli assumption and the Prandtl loss model no longer hold), 0.30 m gives σ = 0.48, the conventional edge of BEM validity; c/R = 0.15 is the top of the small-turbine commercial range (~0.08–0.15); the Schmitz root is 276 mm, so the baseline stays feasible and unclipped (audit §3.4) and the fairness argument is unchanged; measured neither generous nor punitive — +0.147 % at 300 rpm with either 0.30 or 0.45 m, and 2.33 → 2.28 % at 263 rpm |
| `twist_min_deg` | −2° | Schmitz tip twist 0.57° with margin (2026-09-13). **Active at the optimum** from every start |
| `twist_max_deg` | 35° | Schmitz root twist 23.07° with margin (2026-09-13) |

plus `constraints.max_local_solidity: 0.5` — the radius-aware form of the
same limit, `σ_i = B c_i / (2π r_i) ≤ 0.5` at every station, wired through
`ScaledProblem.solidity_constraint(cap)`. Honesty note, recorded in the YAML
too: with the 0.30 m box on the control points and the spline's convex-hull
property, the cap cannot be reached (at the first station, r = 0.334 m,
σ = 0.5 needs 350 mm), so it will report *inactive* in every run; the box is
what bites at the root.

**What changed in code:** `DesignBounds.from_config()` now loads instead of
raising; the provisional set in `tests/test_parameterisation.py::PROVISIONAL_BOUNDS`
(0.45 m placeholder, 2026-09-13) is retired and every script and test reads
the config. Three tests that were written to go red when this landed did,
and were replaced (`test_design_bounds_are_resolved_with_their_recorded_values`,
`test_bounds_from_config_carry_the_decided_values`,
`test_x0_is_feasible_against_the_configured_bounds`). Because the box
changed, the `u` scaling changed with it, so every gradient and optimisation
artefact was re-run — see the 2026-09-19 journal entry.

**What would replace it:** a hub radius and root-attachment concept that
fixes a different upper chord. Then the value is replaced, not reconciled,
and the artefacts re-run.

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

---

## 9. The machine — three facts that decide the Phase 5 objective 🟡 B1 resolved provisionally, B2 provisional, B3 open (scope limit)

**Raised 2026-09-13** by `docs/AEP_GAIN_AUDIT.md` §5, after the FD-driven
and adjoint-driven optimisations both stopped at **+0.217 %** over Schmitz
against the plan's 2–6 % expectation. The audit's finding: under plan §6.2's
strategy — fixed λ = 6.5 at every wind speed, ideal power hold above rated —
AEP is proportional to a single number, `Cp(λ = 6.5)`, up to Reynolds
effects, and the polar-consistent Schmitz `x0` is the analytic maximiser of
exactly that. The optimiser is right; the problem as posed has almost nothing
in it to optimise. A larger honest number requires the machine model to
contain a feature the real machine has. Which features it has is hardware,
not a modelling choice.

**Status 2026-09-19 (evening).** Nothing here blocks Phase 4 or Phase 5.

### B1 — maximum rotor speed ✅ **resolved provisionally 2026-09-19: 300 rpm**

`config/rotor_design.yaml`, `operating.max_rotor_speed_rpm: 300` (62.8 m/s
tip speed at R = 2 m; the ceiling bites at V_c = 9.67 m/s). The operating
law in `src/objective/power.py` is `λ(V) = min(6.5, Ω_max R / V)`. Basis,
recorded in the YAML:

- 341 rpm (= no ceiling below rated) is 71.4 m/s tip speed at rated — above
  the ~60–65 m/s noise-conscious upper limit generally applied to
  residential-scale small turbines.
- The closest commercial analogue in size and rating (Skystream 3.7, 2.4 kW,
  3.72 m rotor) runs to ~330 rpm ≈ 64 m/s tip speed; off-the-shelf
  direct-drive 3 kW PMGs cluster at 250–300 rpm rated.
- 263 and 239 rpm (55 / 50 m/s) were explicitly *not* chosen: they are the
  cases that inflate the gain and would amount to designing for an
  under-sped generator. 300 rpm gives **+0.147 %** (from +0.121 %), i.e. the
  decision does not flatter the optimiser.
- Evidence: `verification/aep_optimisation_experiment/` (2026-09-19).

**Provisional** in the sense that it is a stated engineering basis, not a
generator datasheet. A specified generator's rated rpm replaces it, and the
artefacts re-run. Section 10 is the confirmation item.

### B2 — generator nameplate rating 🟡 provisional (unchanged since the 19th, morning)

The cap no longer floats: `power_per_bin` holds power above rated at
`operating.rated_power_w`, and the adjoint carries the capped bins as a
constant (`adjoint.system.BEMSystem.J_capped`). The value is **provisional**
— `P_aero(11 m/s; x0) = 3822.189755449124 W`, the Schmitz baseline's own
aerodynamic power at the rated wind speed at full precision — and
`tests/test_baseline.py::test_configured_rating_is_the_baselines_aerodynamic_rated_power`
pins it to that basis. Note that under the 300 rpm ceiling the baseline runs
λ = 5.71 at 11 m/s, not 6.5, so `P_aero(11 m/s; x0)` on the *schedule* is no
longer exactly the configured rating; the test pins the basis as stated
(λ = 6.5), and the config comment says so. What B2 still needs is the
nameplate, which replaces the number and retires that test.

### B3 — above-rated limiting mechanism 🟠 open; a scope limit, not a blocker

Ideal hold (current model): a capped bin contributes a constant to the
objective and nothing to the gradient. That is defensible for the *energy*
axis. It is **not** a model of the *loads* above rated: at (V > 11 m/s,
Ω_max) the BEM state produces P_aero > P_rated with no mechanism shedding
the difference, so the thrust and root moment there are those of a machine
that is not the one modelled — and at 20 m/s / 300 rpm, λ = 3.1 and the
inboard stations are at α ≈ 30° on the Viterna extrapolation. Consequence
for Phase 4: the load constraint is applied at **rated wind speed at Ω_max
(11 m/s, 300 rpm, λ = 5.71)**, the highest B3-independent point, inside the
polar cache's α band; the cut-out moment is reported as a labelled
post-check only. B3 would be needed to constrain above rated, and would
change the adjoint structure (a trim equation). Not required for the
project's stated scope.

**Phase 4 delivered 2026-09-19.** The root-moment KS constraint at ≤ rated is
implemented and verified (`verification/load_constraint/`,
`docs/adjoint_derivation.md` §10): the `ε = 0` run is the Phase 5 production
optimum (+0.1465 % AEP over `x0` at the baseline's own moment cap), and the
2/5/10 % Pareto costs 0.0206/0.1282/0.5785 % AEP. Loads **above** rated remain
B3's, kept as a labelled cut-out post-check.

| # | fact | status | measured sensitivity |
|---|---|---|---|
| **B1** | maximum rotor speed | **300 rpm, provisional with basis** | none: +0.121 %; 300 rpm: **+0.147 %**; 286 rpm (60 m/s): +0.34 %; 263 rpm: +2.33 %; 239 rpm: +7.84 % |
| **B2** | generator nameplate rating | provisional, `P_aero(11 m/s; x0)` at λ = 6.5 | floating: +0.217 %; fixed: +0.121 % |
| **B3** | above-rated limiting mechanism | open; the ≤ rated load constraint is **done** (Phase 4, 2026-09-19), above-rated still B3 | not measured |

**What is NOT acceptable as a way of resolving this** (audit §3.5): choosing
a ceiling because it makes the number bigger; swapping the baseline for a
linear-taper blade (report it as context only); keeping the floating rating;
building Schmitz at a datasheet α; lowering `n_crit`. The 300 rpm decision
was taken with those exclusions in front of it, which is why 263/239 rpm are
recorded as rejected.

---

## 10. Kestrel e400 rated rotor speed 🟢 confirm only

**Raised 2026-09-19** with the B1 decision. Nothing is blocked; this is a
citation to strengthen the basis of `max_rotor_speed_rpm: 300`.

**Needed for:** the report's justification of the rotor-speed ceiling. The
Kestrel e400 (3 kW, 4.0 m rotor diameter, manufactured in Gqeberha) is the
nearest regional comparable in size and rating.

**What to pull:** the datasheet rated rotor speed (rpm) and, if given, the
rated tip speed or noise rating. If it supports ~300 rpm, cite it beside the
Skystream 3.7 figure in the config comment and the methodology chapter. If it
does not, record what it says and whether the basis should move; the
artefacts re-run only if the config value changes.
