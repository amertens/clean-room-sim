# Paper Polishing Report: *cleanTMLE — Staged, Outcome-Blind Targeted Learning with Plasmode Data-Quality Stress Testing*

Constructive advisor pass over the manuscript (`reports/manuscript_outcome_blind_dq.qmd`) and the two package vignettes (`cleanTMLE/vignettes/`). Every issue is paired with a concrete fix. Line numbers are to the `.qmd` sources.

Two items raised in the earlier peer review are already being addressed in this work cycle and are therefore only cross-referenced here, not re-argued: the `run_clean_tmle()` enforcement gap (fixed via the new `run_clean_tmle_preoutcome()`/`run_clean_tmle_primary()` split) and the Scenario C detection-vs-fragility calibration (recalibration in `sandbox/scenarioC_recalibration/`).

## Overall Assessment

**Recommendation: Major Revision (mechanical, not conceptual).** The methodology and its unusually honest framing are sound; the FIORD boundary is drawn cleanly and the "what passing does and does not establish" discipline is exemplary. What holds the paper back is a dense cluster of internal-consistency drift from many piecemeal revisions (capability claims and numbers that were true at different package versions) plus a few notation-precision points. Reconciling the drift and finalizing two provisional result sections would clear the highest-value issues.

## Top 3 Strengths

1. **Honest epistemics.** "A STOP is evidence of fragility but a GO is not evidence of robustness" (`line 1308-1333`) and the personnel-vs-simulation blinding tradeoff (`line 794-883`) are more careful than most applied TL-RWE papers.
2. **A falsifiable contribution.** The paper commits to a testable claim (the `min_max_rmse` rule can lock a different candidate than `min_rmse`) and demonstrates it twice (`@sec-divergence`, case study).
3. **Clean positioning against FIORD.** The division of labour (software implementation + DQ extension + audit record) is stated repeatedly and never claims the outcome-blind framework as original (`line 364-372`, `line 2932-2942`).

---

## Critical Issues (fix before submission)

1. **MAR/MNAR support is stated both ways.** "now degrades covariates under MAR and MNAR" (`line 332`, `line 1122`, `line 1641`, `line 3234`) vs "The covariate-missingness scenario assumes MCAR. MAR and MNAR ... are not yet supported" (`line 3114-3116`), with MCAR-only pseudocode at `line 1469` and a MCAR-only assumption list at `line 691`. *(Independently verified: the code implements `.degrade_missingness_mar/_mnar`.)* **Fix:** state the true state once — the stress test *degrades* W under MCAR/MAR/MNAR but the only *handling* exercised is median imputation; rewrite `line 3114-3116` accordingly and fix the `line 1469` pseudocode and `line 2897` list.

2. **"Four threats" vs "five threat families."** Four in the abstract (`line 128`), `@tbl-threat-matrix` (`line 1110`, four rows), and `line 3024`; five at `line 1641`, `line 2074`, `line 2897`, `line 3312`. The fifth (near-positivity) is called a "first-class threat" (`line 2102`) yet is absent from the one table that enumerates what each scenario probes. **Fix:** add a near-positivity row to `@tbl-threat-matrix` (probes positivity; operationalised by amplifying the covariate→treatment association ×2/×3/×4) and make the abstract say five; sweep the wording globally.

3. **Function-reference vignette will not knit (relocated functions).** `cleanTMLE-functions.qmd` calls `attach_estimand()` (L177), `declare_sensitivity_plan()` (L183), `summarize_stage_path()` (L551), `build_stage_manifest()` (L555) as **live** chunks, but these moved to the companion `cleanroomGov` package in the estimation/governance split and are not exported by `cleanTMLE` (`export()=0`), which does not `library(cleanroomGov)` or list it in `Suggests`. *(Independently verified.)* **Fix:** either add `cleanroomGov` to `Suggests`, `library(cleanroomGov)` in setup, and guard the chunks with `eval = requireNamespace("cleanroomGov", quietly = TRUE)`; or replace the calls with surviving `cleanTMLE` primitives. Either way remove the four rows from the Function Index (L74, L75, L112, L117) and the two event-process rows (L139-140) that also relocated.

## Major Issues

4. **The six "identifiability assumptions" mislabel two conditions** (`line 666-698`). "Correctly specified nuisance models" is an *estimation* condition (TMLE's point is that identification does not require it), and outcome-measurement / covariate-missingness are measurement conditions, not identification assumptions. As written the text implies the g-formula fails to identify Ψ if a SuperLearner library is misspecified. **Fix:** present three identification assumptions (exchangeability, positivity, consistency) that identify Ψ, then list the measurement/missingness/estimation conditions the DQ test probes separately. Keep the elegant "probed by" mapping.

5. **Main-simulation configuration described incompatibly.** Design section: 200 replicates (`line 1666`), 19 severity cells (`line 1653`), three PS candidates at `t ∈ {0.001, 0.025, 0.20}` (`line 1656`). Cost section: 30 replicates/cell, nine degraded cells, two candidates (`line 3168`). Limitations: a two-candidate `{0.01, 0.05}` grid (`line 3154`). **Fix:** recompute the cost paragraph against the real 19-cell × 3-candidate × 200-rep design (or label it a reduced illustrative config), and settle on one truncation grid (`{0.001, 0.025, 0.20}` matches the divergence study).

6. **Replicate counts contradict within results.** Case-study DQ reported as 200 (`line 2685`, `line 3392`) and as 50 with "confirmation in progress" (`line 2712`, `line 2763`, `line 2698`). Variance sub-study caption says "100 MC replicates" (`line 1885`) while TMLE/TMLE_CF rows are 60 (`line 1853`). **Fix:** pick the final counts; if 200, delete the "in progress" hedges; footnote the 60/100 split or re-run to a common count.

7. **Overlap-regime letters collide with scenario letters.** Main scenarios A/B/C = good/marginal/**unmeasured confounding** (`line 1622`); variance regimes A/B/C = good/marginal/**very-good overlap** (`line 54-59`, `line 1802`). "Scenario C" means two different things. **Fix:** rename the variance regimes (e.g. overlap tiers VG/G/M).

8. **Bias threshold 0.01 vs 0.02.** DQ defaults and case study use 0.01 (`line 1231`, `line 2684`); the simulation figures and interpretation use 0.02 (`line 1746`, `line 1988`, `line 2138`). **Fix:** state once that the simulation locks a looser 0.02 than the 0.01 package default and why; otherwise the "touching the threshold at 0.010" story (`line 1248`) reads against the 0.02 line.

9. **Case-study N bookkeeping.** Three "complete-case" estimators (`line 2642-2647`) report N = 1,277 (Crude) but N = 1,693 (IPTW, full-cohort TMLE) at `line 2822`; the matched cohort is reported as both 615 (pairs) and 1,230 (individuals). **Fix:** report design N and analyzed N consistently across rows; annotate "615 matched pairs (1,230 individuals)."

10. **Two provisional result sections.** Gate-OC sensitivity/specificity appear only via inline code with "reps set higher for a final run" (`line 2492`); estimator validation says "Provisional results, to be confirmed" (`line 3351`). The gate's operating characteristics are the headline deliverable of `@sec-gate-oc`. **Fix:** run the final reps and state literal numbers with Monte Carlo SE, or mark both subsections "preliminary" up front. *(The 200-rep Scenario C confirmation in this cycle supplies part of this.)*

11. **Plasmode fidelity never quantified for the case study.** Synthetic-outcome fidelity is named as *the* central residual risk (`line 228`, `line 831`, `line 3068`) but no real-vs-synthetic covariate SMD or distributional distance is reported for rescueCo. **Fix:** report one fidelity diagnostic on the case study, or state explicitly that none was computed and why.

12. **Abstract omits Scenario D while over-weighting the variance supplement.** The abstract spends a dense paragraph on the supplementary variance comparison (`line 134-158`) but never mentions Scenario D — the paper's strongest result, where CV-TMLE recovers a misspecified surface that leaves GLM-based estimators on the wrong side of zero (`line 1775-1795`). **Fix:** add one abstract sentence on Scenario D; trim the variance paragraph to two sentences.

13. **Missing seminal citations.** MCAR/MAR/MNAR taxonomy (Rubin 1976 / Little & Rubin) is never cited though the concept is load-bearing from `line 318`; AIPW / double robustness (Robins–Rotnitzky–Zhao 1994, Bang & Robins 2005) is uncited despite `estimate_aipwrisk()` and repeated "doubly robust" claims (`line 420`). **Fix:** add both at first use.

---

## Section-by-Section Feedback

**Abstract (`line 103-177`)** — Add Scenario D; trim the variance paragraph; ensure the 0.007/threshold framing matches the recalibrated Scenario C.

**Introduction (`line 179-431`)** — Strong and well-cited. Consolidate the ~3 restatements of the outcome-blind definition (keep `line 250-260`) and ~2 contribution statements; soften "central contribution" (`line 240`) toward the narrower `line 344` wording. One sentence on FIORD's publication/availability status would preempt the "building on an unpublished 2026 framework" objection.

**Roadmap (`line 433-484`)** — Clean; seven-step enumeration and Step-5/6/7 mapping are accurate.

**Walkthrough (`line 486-598`)** — Uses `fuel_wood` as a surviving negative control (`line 529`, "both near null") but it is later dropped by the near-zero-variance filter (`line 2574`); use a survivor such as `chronic_hypertension`. Also update the walkthrough to show the new enforced two-pass entry point (`run_clean_tmle_preoutcome()` → `authorize_outcome_analysis()` → `run_clean_tmle_primary()`) now that `run_clean_tmle()` is the unguarded convenience path.

**Estimand (`line 600-707`)** — Fix the six-assumption framing (Issue 4). Optionally formalize the IPCW estimator: `line 318` claims unbiasedness "under conditional independent censoring" but no censoring-weight formula appears. Ψ is used for both causal and statistical estimands (`line 610`/`614`); add half a sentence noting they coincide under `@sec-id`.

**Framework (`line 708-909`)** — A highlight. Only needs acronym expansion.

**DQ (`line 911-1181`)** — Add near-positivity to `@tbl-threat-matrix`; reconcile MCAR-only pseudocode with MAR/MNAR claims; fix the 0.01-vs-0.02 default.

**Implementation (`line 1335-1600`)** — The candid bug-disclosure ethos is a strength. Reconcile the two-candidate `{0.01,0.05}` grid with the design section.

**Simulation (`line 1601-2358`)** — Surface Scenario D in the abstract and `@tbl-workflow-contrast`. Fix the regime-letter collision, the 60/100 and 30/200 replicate discrepancies, and the truncation-grid mismatch.

**Gate-OC (`line 2414-2493`)** — Report final sensitivity/specificity as literal numbers; add the promised gate figure or fix the dangling "gate-OC figures" reference (`line 3380`).

**Case study (`line 2495-2926`)** — Fix N bookkeeping and the 50/200 replicate count; use a surviving NC; label the two Love-plot figures.

**Discussion (`line 2928-3305`)** — Reconcile the MCAR-only limitation with the MAR/MNAR capability claims; recompute or re-scope the computational-cost model count (`line 3168`); add the MCAR/MAR/MNAR and double-robustness citations.

## Equation / Notation Review

| Location | Issue | Suggestion |
|---|---|---|
| `line 666-698` | Six "identifiability assumptions" include estimation + measurement conditions | Split into 3 identification assumptions + separate probed conditions |
| `line 610` / `614` | Ψ denotes both causal and statistical estimand | One clause noting they coincide under `@sec-id` |
| `line 318-321` | IPCW unbiasedness asserted, estimator never written | Add a one-line censoring-weight definition or pointer |

## Figure & Table Review (selected)

| Item | Current | Suggested improvement |
|---|---|---|
| `@fig-mc-boxplots` (`1934`) | RD distributions by estimator×scenario | Caption should highlight the Scenario-D panel (crude/GLM on wrong side of zero) |
| `@fig-coverage-bar` (`1953`) | Coverage bars | Add Monte Carlo SE error bars (already computed) |
| `@fig-se-calibration` (`1973`) | SE/SD ratio, line at 1.0 | Shade the (0.8, 1.2) calibration band |
| `@fig-dq-degradation` (`2113`) | Five-threat gradients | Colorblind-safe palette; direct line labels instead of legend |
| `@fig-divergence` (`2289`) | Truncation-candidate gradients | Annotate the crossing point where `min_max_rmse` overtakes `min_rmse` |
| `@tbl-threat-matrix` (`1108`) | Four rows | Add near-positivity row (Issue 2) |
| `@tbl-workflow-contrast` (`2018`) | A/B/C only | Add a Scenario-D row or explain its omission |
| `@tbl-boot-variance` (`1879`) | Mixed 60/100 reps | Footnote the per-row replicate counts |
| `tbl-spifd2-map` (`1020`) | Caption hardcodes "Table 1:" | Add `#| label:` and drop the literal "Table 1:" |
| Love-plot chunks (`2734`, `2738`) | `fig.cap` but no `#fig-` label | Add labels so text can cross-reference |

## Missing Citations

| Location | Claim | Suggested citation |
|---|---|---|
| `line 691` / `2075` | MCAR/MAR/MNAR taxonomy | Rubin (1976); Little & Rubin |
| `line 420` | "doubly robust" / AIPW | Robins–Rotnitzky–Zhao (1994); Bang & Robins (2005) |
| `line 678` / `1060` | practical positivity, truncation-changes-estimand | Petersen et al. (2012) |
| `line 1547` | "Kish formula" | Kish (1965) |
| `references.bib:175` | `@shi2020selecting` in bib, uncited | Cite in the NCO discussion (`line 3033`) or drop |
| `line 141`, `1711` | free-text "(Abadie and Imbens 2008)", "(Gruber 2023 ...)" | Convert to `@abadie2008bootstrap`, `@gruber2023future` |

## Vignette-specific findings

- **Critical:** function-reference relocated-function calls (Issue 3 above).
- **Major:** experimental time-to-event helpers (`estimate_ipwrisk/gcomprisk/aipwrisk/ipwhr`) run as ordinary live examples in `cleanTMLE-functions.qmd` L462-534 despite carrying runtime "outside the tested v0.1 scope" warnings — add a `callout-warning` before the "Model Specification DSL" header (L457) and/or set those chunks `eval = FALSE`.
- **Major:** the tested-scope limitation is prominent in the staged vignette (callouts at L270-278) but nearly absent up front in the function reference — add a `callout-important` after the intro stating binary point-treatment / binary outcome / marginal RD scope.
- **Major:** the hand-written `pick_ci()` helper (staged L130-141) exists only because CIs live in `$estimate/$ci_lower` for `run_*_workflow()` objects but `$estimates$ATE$…` for the modular extractor — point readers at `summarize_cleanroom_results()` (already exported) rather than modelling reconciliation glue.
- **Minor:** the staged vignette introduces itself twice (plain wrap-up at L179, then a formal "Purpose/Conceptual workflow" restart at L194) — mark the `---` at L192 as "Part II"; the stage ladder appears three times (L224, L253, L1362) — show it once. RWE/SPIFD2/FIORD unexpanded; two parallel decision-log APIs shown without saying which is canonical; a redundant hand-built forest plot where `forest_plot()` is exported.
- **Highest-value addition:** a short "applying the workflow to messy data" subsection (inject a factor covariate, ~10% covariate NAs, thin the outcome to a rare event) so `sanitize_covariates()`, `checkpoint_cohort_adequacy()`, and the IPCW path are exercised on data that actually stresses them.

## Revision Checklist (critical → minor)

- [ ] **P1 (Critical):** Reconcile MCAR vs MAR/MNAR capability across `line 332`/`1122`/`1641`/`3234` and `line 691`/`1469`/`2897`/`3114`.
- [ ] **P2 (Critical):** Make "four vs five threat families" consistent; add near-positivity row to `@tbl-threat-matrix`.
- [ ] **P3 (Critical):** Fix `cleanTMLE-functions.qmd` relocated-function calls (guard for `cleanroomGov` or replace) and the Function Index rows.
- [ ] **P4 (Major):** Recast the six "identifiability assumptions" as 3 identification + separate probed conditions (`line 666`).
- [ ] **P5 (Major):** Reconcile main-simulation config (candidates, truncation grid, cells, reps) between `line 1653-1666` and `line 3154-3168`.
- [ ] **P6 (Major):** Settle case-study (50/200) and variance-study (60/100) replicate counts; delete "in progress" hedges.
- [ ] **P7 (Major):** Rename variance overlap regimes to remove the A/B/C collision.
- [ ] **P8 (Major):** State the 0.01 (package default) vs 0.02 (simulation) threshold once, explicitly.
- [ ] **P9 (Major):** Fix case-study N bookkeeping (1,693 vs 1,277; 615 vs 1,230).
- [ ] **P10 (Major):** Finalize gate-OC and validation numbers (use the 200-rep Scenario C confirmation).
- [ ] **P11 (Major):** Report one plasmode-fidelity diagnostic for the case study, or state none was computed.
- [ ] **P12 (Major):** Add Scenario D to the abstract; trim the variance paragraph.
- [ ] **P13 (Major):** Add MCAR/MAR/MNAR and AIPW/double-robustness citations.
- [ ] **P14 (Major, vignette):** Caveat the experimental time-to-event helpers; add function-reference scope callout.
- [ ] **P15 (Minor):** acronym expansions; hardcoded "Table 1"; free-text citations → `@key`; `fuel_wood` → surviving NC in walkthrough; soften "lower edge" for 0.93 coverage; label the two Love-plot figures; fix the dangling "gate-OC figures" reference; STaRT-RWE capitalization.
- [ ] **P16 (Minor, vignette):** de-duplicate the staged-vignette intro and stage ladder; point to `summarize_cleanroom_results()`; add a messy-data worked example.
