# cleanTMLE TODO

Long-tail development items, grouped. Updated 2026-05-11.

Items shipped or merged into the manuscript/presentation/vignettes are marked **[done]** and kept for traceability.

---

## Z. Five-reviewer audit action items (2026-07-01)

From the independent five-reviewer audit (Muntner design, FDA regulatory, TL/R
code, outcome-blind simulation, applied analyst) plus the consensus sign-off.
Severity in brackets. Items checked off were fixed in the 2026-07-01 pass and
verified; the rest are staged by the joint sign-off's four categories. Existing
A-section items are cross-referenced where they overlap.

### (a) Statistical / code correctness
- [x] Expose model-based (TMLE plug-in) arm risks from `extract_tmle_estimate()`
  as `estimates$risk_treated` / `estimates$risk_control`; switch the
  staged-analysis vignette risk-report chunk off the `crude +/- ATE/2`
  reconstruction. [important] — verified: `risk_treated - risk_control == ATE`.
- [x] Fix the vignette E-value chunk to use the confounding-adjusted arm risks
  instead of the crude control risk (`crude$r0`). [important]
- [ ] Keep survival estimates strictly secondary and document that
  `.surv_tmle_fallback` (`tmle.R:587-694`) excludes rather than IPCW-weights
  early-censored subjects; scope a real discrete-time / IPCW survival TMLE.
  [important] (rel. A.13)

### (b) Clean-room / regulatory representation
- [x] Correct the Muntner citation year 2020 -> 2024 in
  `cleanTMLE-staged-analysis.qmd` (only surviving 2020 instance). [important]
- [ ] Wire `assert_outcome_authorized()` into the Stage 4 estimators so the gate
  enforces recorded authorization, not just outcome masking
  (`.check_outcome_access`, `cleanroom.R`). Gate only when
  `cleanroom_enabled = TRUE`; update the staged vignette and tests to authorize
  first. [critical]
- [ ] Soften the manuscript enforcement prose (`manuscript_outcome_blind_dq.qmd`
  ~304-305, "withholding ... until checkpoints recorded") to match what the code
  enforces. [important]
- [ ] Re-document the RescueCo inter-facility-transfer exclusion as a dated,
  pre-outcome protocol amendment with recorded rationale; acknowledge that
  injury-to-arrival transport time is an outcome-adjacent (proxy) quantity; log
  the CP1-STOP override that the primary analysis proceeded past. [critical]
  (rel. A.12)
- [ ] Regenerate `case_study_metadata.json` under 0.1.5 (records "0.1.1"); fix
  provenance by re-running, not by hand-editing. [minor]

### (c) Package scope / organization — consensus: TRIM now, SPLIT direction
- [ ] Execute A.17 dead-export cleanup: internalise `expit`, `logit`,
  `wrap_ps_fit` (and candidates); test or deprecate the rest. None of A.17 has
  landed. [important] (rel. A.17)
- [ ] Generate the missing `man/*.Rd` for `plot_dq_heatmap` (exported but
  undocumented + untested). [important]
- [ ] Add tests for the ~40 untested exports, or deprecate them. [important]
- [ ] Decide TRIM vs SPLIT into an estimation package + a companion governance
  package/`cleanroomGov` (move the note-taking metadata family to docs; keep the
  real enforcers in code). See the audit's Reviewer 3 merge/internalise/remove
  lists. (rel. B.2)

### (d) Usefulness / adoption
- [ ] Add a positivity-strained DGP (reuse the in-package `near_positivity`
  mechanism) so the truncation-varying candidate grid is identifiable; the
  current five DGPs omit positivity strain. [important]
- [ ] Re-run the plasmode study at an inferential budget (outer >= 200,
  n >= 1000, inner >= 200 if FIORD/coverage selection is retained) or drop FIORD
  to a min-RMSE-only rule; the pilot cannot demonstrate value and FIORD collapses
  to parametric at inner_reps = 10 (coverage MC SE ~0.16 >> the +/-0.02
  tolerance). [important]
- [ ] Report per-estimator / per-candidate Monte Carlo SEs. [important]
  (rel. A.6)
- [ ] Implement the A.15 low-replicate guard (`cleantmle_workflow_demo` class +
  forced FLAG); currently prose-only. [important] (rel. A.15)
- [ ] Wire the existing MAR/MNAR degraders into the default DQ preset and
  exercise a candidate-side missingness handler (default is MCAR + median
  imputation). [important] (rel. A.10)
- [x] Manuscript factual fixes: pilot budget 40 reps/n=600 -> 10 reps/n=500;
  q0_library version 0.2.1 -> 0.1.3 (0.2.x versions do not exist). [minor]
- [ ] Revisit the manuscript's "FIORD only with inner_reps >= 30"
  recommendation; the audit's MC-SE analysis puts the usable floor near 200 when
  coverage drives selection. [minor]

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

### A.14 Estimand-revision hook after STOP
- When `gate_all()` returns STOP, surface a structured list of estimand-revision options the analyst can record in the decision log: "revise target population", "revise treatment strategy", "revise follow-up window", "switch to a matched / overlap-restricted target", "change missingness strategy", "stop before outcome access".
- Each revision becomes a logged decision rather than a silent override.
- Cite: targeted-learning RWE roadmap; the manuscript's new "estimand revision after STOP" subsection in §6.

### A.15 Automatic FLAG on low-replicate examples
- When `run_plasmode_feasibility()` or `run_plasmode_dq_stress()` are called with `reps`/`n_reps` below a threshold (e.g. < 50), tag the result object with class `cleantmle_workflow_demo` and force `gate_all()` to return FLAG with the rationale "workflow-demo-only replicate count".
- Prevents low-replicate bundled examples from being mistaken for inferential simulations.

### A.16 Structured negative-control attrition in the audit log
- When a prespecified negative-control endpoint is dropped at the data-sanitisation step (e.g. NZV filter), record the drop as a structured `decision_log` entry rather than only printing a console message.
- Capture: NC name, filter that fired (NZV / collinearity / other), prevalence in the cohort, the SAP rationale for including the NC, and the impact ("narrower NC panel; residual-confounding probe partial").
- Closes the gap the rescueCo case study identified.

### A.17 Dead-export cleanup (code audit, 2026-05)

The function-usage audit identified 14 exported functions with no
caller in tests, vignettes, manuscript, or case study. Each needs
either an executable usage example or a deprecation.

Advertised in README and `_pkgdown.yml` but not exercised:
- `hr_data()`
- `update_censoring()`, `update_outcome()`, `update_treatment()`
- `identify_missing()`
- `inspect_ipw_weights()`
- `re_estimate()`

Name-mentioned in vignettes but not invoked:
- `refine_ps_after_nco()`
- `get_final_cohort()`
- `load_audit()`
- `fit_ps_parallel()`
- `create_analysis_lock_from_yaml()`

Only referenced in manuscript prose:
- `tipping_point_sensitivity()` (rescueCo case study uses the
  output value but does not call the function explicitly)

Internal-only candidates:
- `validate_superlearner_spec()` — only its `.Rd` references it
- `expit`, `logit`, `wrap_ps_fit` — used internally; mark
  `@keywords internal`

Action: for each "advertised but unused" function, either add an
executable example in `cleanTMLE-functions.qmd` or move to
`@keywords internal` and drop from `_pkgdown.yml`.

### A.18 Pedagogical pass on `cleanTMLE-staged-analysis.qmd`

Current vignette opens with governance vocabulary before any
clinical anchor and runs the DQ stress test with `eval = FALSE`.
Pedagogically less accessible than the comparator's vignettes.
Concrete actions:

- Add a real-dataset opening with one clinical question (port an
  ACTG / WIHS-style example or use a clean simulated analogue with
  a clinical narrative).
- Add a "concept primer" block at the head defining clever
  covariate, EIC, plasmode, ESS, SMD, MDD, IPCW in one place.
- Render the four-step TMLE output between code chunks rather
  than only at the end.
- Add a stacked-estimators risk-curve overlay (crude → IPW →
  g-comp → AIPW → TMLE) — the single most pedagogically effective
  figure in the comparator's vignettes.
- Demote the GO / FLAG / STOP vocabulary to where it is actually
  earned (after the reader has seen a Table 1).
- Eliminate forward references to stages that have not yet been
  introduced.

### A.19 Comparator feature-parity items (2026-05 parity audit)

Gaps identified by the function-by-function parity audit. Each is
a concrete API or dataset addition:

- `estimate_ipwcount()` and the cumulative-count family
  (`identify_count`, `update_count`, `plot.cumcount`).
- `compare_protocols()` analogue for > 2-arm cumulative-risk
  comparison.
- `subgroup()` as an exported model-spec verb (current workaround
  uses `subset_idx` arguments).
- `update_missing()`, `update_count()`, `update_label()`.
- SMR weight scheme (`wt_type = 1`) in `estimate_ipwrisk()`.
- Public reproducible datasets equivalent to `actg`, `wihs`,
  `leukemia` (license-permitting) or clean simulated analogues
  with realistic structure.
- Standalone `trim_ps()` exported helper.
- `hist.ipw` method on weight objects for one-liner muscle memory.

### A.20 Block-bootstrap variance for matched TMLE

The matched-cohort TMLE EIC currently treats the matched subset as
an iid sample. The reported SE is therefore not a paired-design SE.
Add a block-bootstrap variant that resamples matched pairs and
returns the bootstrap SE alongside the IF-based SE. See manuscript
§8.2.4 and the in-source NOTE in `run_matched_tmle()`.

### A.21 Plasmode candidate-grid expansion (workshop-driven, 2026-05)

Five PDFs from the March workshop materials (Phillips et al. 2023
*Practical considerations for specifying a super learner*; Gruber et al.
2023 *Developing a TL-based SAP* + supplement; Gruber et al. 2023
*Evaluating and improving RWE with TL*; Gruber & van der Laan 2009
*Gentle intro to TMLE*) identify seven axes that the plasmode loop
should expose so candidate selection probes the dimensions that
actually move finite-sample bias and CI coverage:

**A.21.1 Variance menu.** Add `variance_method = c("IF","cv_IF","robust","bootstrap_HAL","targeted_bootstrap")`
to `tmle_candidate()` and surface it as a candidate-grid axis.
IF (current default) is anti-conservative under near-positivity;
`robust` calls upstream `tmle::tmle(variance.method = "tmle")`
plug-in variance (Tran et al. 2018); `bootstrap_HAL` forces HAL on
both Q and g and bootstraps; `targeted_bootstrap` implements Coyle
& van der Laan 2018. SAP p.472.

**A.21.2 CV-TMLE vs sample-split vs no-CV.** Add `cv_scheme = c("none","cv_tmle","sample_split")`.
The current `n_folds/fold_vec` arguments collapse into
`cv_scheme = "cv_tmle"`. Phillips Step 3.

**A.21.3 SuperLearner library presets.** New helper
`build_sl_library(role = c("Q","g","Delta"), n_eff, p, preset = c("small_n","default","rich","very_rich"))`
returning role-aware learner lists keyed by effective sample size
and dimensionality. Adopt the SAP-supplement Checklists C–E spec
for the three roles (binomial family / NNLS / stratifyCV per role).
Phillips Steps 2-5.

**A.21.4 Discrete vs convex super learner.** Add `discrete_sl = FALSE/TRUE`;
when TRUE, build the SL with `method.NNLS` and include the eSL
itself in the dSL pool. Phillips Step 5.

**A.21.5 Screeners.** Add `screener = c("All","corP","corRank","glmnet","randomForest")`;
couple a screener with every Q, g learner inside CV. Phillips p.1283.

**A.21.6 Truncation rule names.** Allow `truncation_rule` to accept
either a numeric bound or a named rule: `"sqrt_n_ln_n"` (= `5/(sqrt(n)*ln(n))`,
the SAP-supplement default), `"fixed_001"`, `"fixed_025"`, `"fixed_05"`.
Document running the same TMLE under at least three rules and
reporting all three. SAP-supp §8.3.2.

**A.21.7 Estimator family.** Add `estimator = c("tmle","aipw","onestep","ctmle","drtmle")`.
AIPW and one-step are wrappers around existing nuisances; C-TMLE
delegates to `ctmle::ctmleDiscrete()`; drtmle delegates to
`drtmle::drtmle()`. C-TMLE is the high-priority addition because
of its identifiability-recovery behaviour under borderline ETA.
Gruber 2009 §3.2.

**A.21.8 Tmle-control pass-through.** Add `tmle_control = list(fluctuation, alpha, target.gwt, automate, min.retain, cv_Qinit)`
so upstream `tmle::tmle()` controls are surfaceable from the
candidate spec. SAP-supp B2.

### A.22 New plasmode helpers and post-processors (workshop-driven)

- `compute_n_eff(Y, family)` — implements `n_eff = min(n, 5*n_rare)`
  for binary outcomes. Drives default `cv_V` and library size.
- `compute_G_value(fit)` — the size of causal gap that would flip
  the conclusion: `G = min(|psi - 1.96*se - null|, |psi + 1.96*se - null|)`.
  RWE Gruber et al. 2023 p.5.
- `run_delta_sensitivity(fit, delta_grid)` — varies the causal-gap
  delta over a grid and returns the corresponding point estimate,
  CI, and p-value. SAP §8.3.4.
- `run_positivity_diagnostics(fit)` — PS C-statistic, propensity
  histograms by arm, % truncated, bounded-g quantiles, and a
  near-positivity-violation flag using
  `eps(n) = max(0.01, 5/(sqrt(n)*ln(n)))`. Already partly available;
  reorganise into a single object.
- `run_bootstrap_variance(fit, B, type = c("nonparam","targeted","HAL_targeted"))` --
  three bootstrap variance implementations on the modular path.

### A.23 Selection-criterion expansion in `select_tmle_candidate()`

Add `criterion = c("min_rmse","minmax_rmse","min_bias","ci_coverage","min_se","composite")`
so the selector can prioritise CI coverage rather than RMSE alone.
The workshop materials emphasise coverage under near-positivity.
This is also the second stage of the FIORD two-stage selector (A.1).

### A.24 Plasmode-loop Cartesian-product grid

Update `run_plasmode_feasibility()` to accept the new axes as named
lists (`variance_methods`, `estimators`, `cv_schemes`,
`truncation_rules`, `sl_presets`, `discrete_sl`, `screeners`) and
Cartesian-product them with a `max_candidates` cap so the grid
stays bounded. Emit a stable `tmle_candidate_id` from a hash of the
serialised args so `select_tmle_candidate()` can join results across
runs.

### A.25 Difficulty grouping for items A.21-A.24

- **Low (argument pass-through):** A.21.2, A.21.4, A.21.5, A.21.6,
  A.21.8, A.22 (compute_n_eff, compute_G_value), A.23.
- **Medium (new code path):** A.21.1 (IF/cv_IF/robust variance),
  A.21.3 (build_sl_library), A.21.7 (AIPW/one-step), A.22 (delta-
  sensitivity, positivity diagnostics consolidation), A.24
  (Cartesian-product grid).
- **High (new estimator family or wrapping a separate package):**
  A.21.1 (bootstrap_HAL, targeted_bootstrap), A.21.7 (C-TMLE via
  ctmle pkg, drtmle via drtmle pkg), A.22 (run_bootstrap_variance).

The low and medium items can be staged across two cleanTMLE 0.2
point releases; the C-TMLE and bootstrap variance work is the
natural 0.3 release.

### A.26 Shipped in 2026-05-11 (status update)

**Shipped (in `R/plasmode_extensions.R` and `cleanroom.R`):**

- `compute_n_eff()`, `recommend_cv_V()` — Phillips Steps 2-3.
- `build_sl_library(role, n_eff, p, preset, include_screeners)` —
  four presets keyed by `n_eff`.
- `resolve_truncation_rule()` — named rules `"sqrt_n_ln_n"`,
  `"fixed_001"`, `"fixed_025"`, `"fixed_05"`, plus numeric pass-through.
- `run_positivity_diagnostics(ps_fit)` — single tidy object with C-stat,
  per-arm summary, % truncated, bounded-g quantiles, near-violation flag.
- `compute_G_value()` — pre-QBA tipping-point readout.
- `run_delta_sensitivity()` — causal-gap delta sweep on top of an
  existing fit.
- `compute_aipw(g_fit, Q_fit)` — AIPW estimate from the same nuisances
  the modular TMLE uses.
- `tmle_candidate()` extended with `variance_method`, `cv_scheme`,
  `cv_V`, `estimator`, `discrete_sl`, `screener`, `tmle_control`,
  `match_spec`. Back-compatible (all new args default to current
  behaviour).
- `expand_tmle_candidate_grid()` extended with the new axes as named
  lists; Cartesian-products them with a `max_candidates` cap.
- `select_tmle_candidate()` extended with `rule = "ci_coverage"`,
  `"min_se"`, `"composite"`, `"fiord_two_stage"`. FIORD two-stage:
  Stage 1 screens candidates by oracle coverage within
  `fiord_coverage_tol` of nominal, Stage 2 picks the smallest mean SE.
- `fit_tmle_treatment_mechanism()` now resolves named truncation
  rules at fit time (when sample size is known) and auto-sets V from
  the candidate's `cv_scheme = "cv_tmle"` flag using the Phillips rule.

**Still pending (the work that requires deeper code paths):**

- A.21.1 `variance_method = "robust"` dispatch to
  `tmle::tmle(variance.method = "tmle")` inside the modular path.
- A.21.1 `variance_method ∈ {"bootstrap_HAL","targeted_bootstrap"}`
  end-to-end implementations.
- A.21.7 dispatch of `estimator = "aipw"` / `"onestep"` inside the
  plasmode candidate loop (currently the candidate spec records the
  choice but the fitter still runs TMLE). `compute_aipw()` is
  available as a separate post-processor.
- A.21.7 C-TMLE wrapper (`run_ctmle_candidate()`) and drtmle wrapper.
- A.21.2 `cv_scheme = "sample_split"` (single 50/50 holdout) path
  inside the modular TMLE.
- A.21.5 screener wrapping inside the SL library construction (the
  candidate spec records the choice but the SL fit does not yet wrap
  each base learner with the screener).
- A.21.8 forwarding of `tmle_control` arguments to `tmle::tmle()` on
  the delegation path.

These pending items are tracked here and in the "Experimental /
planned extensions" section of README.md.

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
