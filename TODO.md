# cleanTMLE TODO

Long-tail development items, grouped. Updated 2026-05-11.

Items shipped or merged into the manuscript/presentation/vignettes are marked **[done]** and kept for traceability.

---

## A. Methodological extensions (cleanTMLE 0.2)

These require new R code and tests, and rerun no existing simulation.

### A.1 FIORD two-stage candidate selector
- Implement `select_tmle_candidate(..., rule = "fiord_two_stage")`.
- Stage 1: filter candidates to those whose Monte Carlo empirical-SE-based CI achieves nominal oracle coverage; among the survivors pick the lowest variance.
- Stage 2: choose the variance method (IF, TMLE-robust, bootstrap) that achieves nominal coverage on the locked point estimator.
- Documentation, example chunk, and a unit test against a small known-truth simulation.
- Cite: Nance et al. 2026 (FIORD). Manuscript and vignette already reference this as planned for 0.2.

### A.2 Nonparametric bootstrap variance
- Add a `variance_method = "bootstrap"` argument to `run_clean_tmle()`, `run_ipcw_tmle()`, `run_matched_tmle()`, `estimate_ipwrisk()`, and the plasmode loop.
- Defaults: `B = 200`, optional block-bootstrap over the matching draw for matched TMLE.
- Surface the bootstrap SE alongside the IF-based SE in `print()` methods.
- Becomes the recommended default for IPTW, matched TMLE, and any TMLE on a non-i.i.d. sample.

### A.3 Unmeasured-confounding severity gradient
- Extend `run_plasmode_dq_stress()` so the latent-U scenario accepts a grid of strengths, not a single fixed value.
- Default grid: `unmeasured_U_strength = c(0.1, 0.3, 0.5, 0.7)`, or a logged risk-ratio sweep.
- Output the per-strength bias / RMSE / coverage as a tidy degradation gradient that the existing `summarize_dq_degradation()` can plot.
- Documented in the manuscript Future-extensions list and in the functions vignette Planned-helpers list.

### A.4 DGP-spectrum argument
- Add `dgp_mode = c("hybrid", "structural_equation", "parametric_bootstrap")` to `run_plasmode_feasibility()` and `run_plasmode_dq_stress()`.
- "hybrid" = current default (empirical W + A, fitted outcome model).
- "structural_equation" = user-specified DGP via `simcausal`-style spec.
- "parametric_bootstrap" = resample from the fitted joint $(W, A, Y)$ — outcome-blind only with sample splitting; document the trade-off.
- Lock fingerprint should record the DGP mode used.

### A.5 Synthetic-data fidelity metrics
- Add `assess_dgp_fidelity()` that computes Hellinger or Wasserstein distance (or, as a simpler default, SMDs and KS distances) between real and synthetic covariates and treatment.
- Surface a fidelity warning in the plasmode print methods when the metrics exceed prespecified thresholds.

### A.6 Monte Carlo SE per estimator
- Add an explicit `mc_se` column to the plasmode summary tables.
- Update the simulation manuscript figures to display per-estimator MC SEs as error bars.

### A.7 Oracle-coverage diagnostic plot
- New `plot_oracle_coverage()` that draws empirical oracle coverage by candidate against the nominal level, with the FIORD-style screen marked.
- Becomes a standard plasmode visualisation in the dossier.

### A.8 Matched-cohort vs full-cohort plasmode comparison
- Run both estimators on the same plasmode replicates under a known truth.
- Report the gap as bias-vs-truth for each, and whether it is consistent with sampling noise + the efficiency cost of matching (vs effect modification or positivity strain).
- Will close out the rescueCo open question about the 0.026 vs 0.031 gap.

### A.9 Pre-protocol stress-test mode
- New `run_preprotocol_plasmode()`: takes user-specified covariate distributions, treatment-mechanism strength, outcome event rate, and the same DQ scenarios; produces a candidate-selection report without requiring real data.
- Use case: deciding which dataset to use before a DUA is signed.
- Documented in the manuscript and functions vignette as planned for 0.2.

### A.10 MAR / MNAR missingness in the DQ stress test
- Add an MAR mechanism (missingness depends on W and A) and an MNAR mechanism (missingness also depends on Y).
- Pair with an IPCW or G-computation handling step inside the candidate run, so the stress test actually exercises the missingness model rather than just the imputation.
- Currently the DQ covariate-missingness scenario is MCAR + median imputation, which the manuscript and reviewer both correctly noted is a weak test.

### A.11 Time-window confounders
- Add a built-in DQ scenario for time-window confounders (strike-window in the rescueCo motivation, calendar-time shocks more generally).
- Operationalised as a covariate available only outside a calendar window; the workflow shows whether candidate estimators degrade when the window-period covariate is dropped.

### A.12 Inter-facility transfer / restriction imbalance
- A built-in DQ scenario for differential restriction by treatment arm (motivated by the rescueCo inter-facility-transfer pattern).

### A.13 Longitudinal / iterated TMLE
- Wrap `ltmle` so the same lock + plasmode + DQ + gate machinery works for time-varying treatments with censoring.
- This is the FIORD case-study setting and the main next-step extension Carrie / Zhiwei flagged.

---

## B. Documentation, governance, and operational

### B.1 Repo move to a Berkeley / group GitHub org
- Travis's suggestion in the GHEP methods meeting.
- Migrate `amertens/cleanTMLE` to a group-owned org (e.g. `berkeley-rwe/cleanTMLE`).
- Update README, DESCRIPTION URLs, pkgdown site, manuscript resources section, presentation resources slide.
- Set up a redirect on the old repo.

### B.2 Spin out the outcome-blind simulation core
- Andrew's own action item in the meeting: package the plasmode + DQ stress test as a separate, lighter-weight package (working name `obsim` or similar).
- Keeps the clean-room / lock / gate machinery in `cleanTMLE`; lets analysts who only want the simulation loop adopt it without buying the full staged-analysis governance frame.
- Will require deciding the dependency direction (cleanTMLE depends on obsim, or vice versa).

### B.3 Calibration helpers for severity ranges
- A YAML reader that ingests validation-literature parameters (PPV, sensitivity, specificity, claims-accuracy ranges) and translates them into `data_quality_scenarios` lists with citations attached to the audit log.
- Already listed in the manuscript Future-extensions list (item 11).

### B.4 `clean_target_trial_table()` and `checkpoint_dashboard()` helpers
- Both currently documented as planned-for-0.2 in the functions vignette.
- Each can be built from documented primitives today; promoting them to exported helpers reduces user-side boilerplate.

---

## C. Manuscript / presentation / vignettes

### C.1 Re-render after each batch of A or B items
- Manuscript: HTML + DOCX.
- Presentation: PPTX.
- Vignettes: both HTML.
- Tutorial: HTML + DOCX.

### C.2 Confirm the small-font slide template reads cleanly at presentation magnification
- The font-size shrink (30% on `sz` attributes in the slide master) was applied; visually verify before any external presentation.

### C.3 Survey of related software for the manuscript
- The current side-by-side table compares TMLE-only packages, causalRisk, and clean-room methods papers.
- Consider adding rows for: `tmle3` (target-learning ecosystem), `lmtp`, `survtmle`, `WeightIt`, `MatchIt`. The point is operational — what each does and does not provide — not a method ranking.

### C.4 Tutorial follow-up
- The `cleanTMLE_for_applied_analysts.qmd` tutorial currently renders to HTML and DOCX. Consider adding PDF output once the cross-reference style stabilises.

---

## D. Housekeeping

### D.1 Remove leftover local Word files
- `reports/manuscript_outcome_blind_dq_old.docx` and `reports/~$nuscript_outcome_blind_dq.docx` are local-only and gitignored, but should be deleted once Word releases them.

### D.2 Re-verify DOCX render after Word closes the lock
- Documented in earlier commit messages; the workflow is now stable but worth noting.

---

## Already shipped (kept for traceability)

### Manuscript
- **[done]** Cite Nance et al. 2026 (FIORD) at the start of §6 and note the DGP spectrum.
- **[done]** Reframe `select_tmle_candidate(rule = "min_rmse")` in §6 as a single-step approximation of the FIORD two-stage selector; flag the full two-stage rule as planned for 0.2.
- **[done]** Adopt the "all the data except the outcome" phrasing in §6.4 (Why outcome-blind matters here).
- **[done]** Add a variance-estimator-selection paragraph in the simulation interpretation (§8.2.4) citing FIORD and Mertens-Zhang and recommending nonparametric bootstrap as the principled fix for IPTW / matched-cohort SE conservativeness.
- **[done]** Expand the manuscript Future-extensions list to 11 items including the FIORD two-stage selector, bootstrap variance, unmeasured-U gradient, pre-protocol mode, DGP spectrum, synthetic-data fidelity, matched-vs-full plasmode, longitudinal extension via `ltmle`.
- **[done]** Note the matched-vs-full plasmode diagnostic as planned in §9 to disambiguate the 0.026 vs 0.031 case-study gap.

### Presentation
- **[done]** Plain-English five-section body restructure (Why this matters, How the package works, Simulation results, Case study, How cleanTMLE improves on existing tools).
- **[done]** Appendix A–D (conceptual tables, case-study setup, diagnostic figures, helpers + outcome-access taxonomy).
- **[done]** Appendix E: 16-step function-by-function Rescue.Co walkthrough.
- **[done]** Matched vs full-cohort slide with the corrected framing (sampling noise + efficiency cost in the rescueCo case; same-estimand under standard assumptions).
- **[done]** Six-methods plain-language table, traditional-vs-OB-workflow contrast, method recommendations for colleagues, what-the-simulation-does-and-does-not-show slide.
- **[done]** Slide-master font-size shrink (30% on `sz` attributes).

### Vignettes / README
- **[done]** Pre-outcome study dossier section in README and staged-analysis vignette.
- **[done]** Outcome-access taxonomy in functions vignette.
- **[done]** Low-replicate warning callout.
- **[done]** "What external governance must still provide" section in README.
- **[done]** Planned-helpers section in functions vignette, expanded with FIORD-driven items (A.1, A.3, A.9).

### Package code
- **[done]** `run_ipcw_tmle()` Delta-path fallback now warns and exposes `method_used`.
- **[done]** Outcome-access checks on Stage 4 estimators; `authorize_outcome_analysis()` GO/FLAG/STOP + override audit semantics in roxygen.

### Tutorial
- **[done]** `reports/cleanTMLE_for_applied_analysts.qmd` — 12-section plain-English numbered outline, ~1-hour read for an applied analyst with no prior TMLE or clean-room background.
