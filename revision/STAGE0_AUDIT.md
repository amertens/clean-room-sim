# Stage 0 audit: cleanTMLE revision

Audit date: 2026-09-14. Repository commit: `863efc3` (clean-room-sim, branch main, working tree clean at audit start); cleanTMLE submodule at `9f22d50`. Environment: R 4.4.2 (Windows 11), Quarto 1.6.42, cleanTMLE 0.2.0 reinstalled from the current source before any run. Everything below distinguishes what was run in this session from what was read statically. No manuscript prose, package code, or simulation script was modified.

Layout deviations from the assumed layout in the revision prompt: `revision/` did not exist and was created for these outputs; the vignettes live in `cleanTMLE/vignettes/`, not a top-level `vignettes/`; `run_gate_operating_characteristics.R` sits at the repository root, not in `reports/`; the four-scenario simulation outputs live in `results_new/`; the presentation source is `reports/cleanTMLE_presentation.qmd`. Two stray directories named `C2, C3, C4` and `PRIMARY, C1` sit at the repository root (apparently created by a comma-separated argument reaching a file-system call); they are repository hygiene, not results.

---

## 2.1 Repository and build state

### R CMD check --as-cran (run in this session)

`R CMD build` followed by `R CMD check --as-cran --no-manual` on the built tarball (`cleanTMLE_0.2.0.tar.gz`), in a scratch directory outside the repository. Result: **2 ERRORs, 6 WARNINGs, 4 NOTEs**.

ERRORs:

1. **Examples fail.** The example for `assert_outcome_authorized` calls `create_audit_log()`, which the 0.2.0 consolidation unexported; the example halts with `could not find function "create_audit_log"`. The superseded governance functions kept their Rd files and executable examples after their dependencies were unexported, so the first such example kills the whole examples run. Other superseded examples after it are untested.
2. **Vignette rebuilding fails** for both vignettes inside the check environment. The proximate error is `unable to start png() device` when writing `cleanTMLE-staged-analysis_files/figure-html/...`. Both vignettes render cleanly from a clean session outside the check (see below), so this is at least partly a check-environment path issue, but as-cran is the standard CRAN environment and the failure stands as reported.

WARNINGs:

1. Non-ASCII characters in `R/plots.R`.
2. Undeclared dependencies: `rpart` used via `::` and `requireNamespace` without a Suggests entry; `earth` and `gam` likewise; and a `cleanTMLE:::.plasmode_fit_one_candidate` call to the package's own namespace (in the callr child function of `.fit_candidates_bounded`).
3. Rd cross-references: `support_thresholds.Rd` contains literal square brackets (`[0.043, 0.937]`, `[0.004, 0.940]`) that Rd parses as links.
4. **`estimate_effect` is undocumented** ("Undocumented code objects: 'estimate_effect'"). The roxygen block exists in `R/estimate_effect.R` but `man/` has no Rd for it: the documentation was not regenerated after the source changed.
5. Code/documentation mismatches: `build_sl_library` (argument `y` missing from docs), `declare_estimand_ladder` (argument `candidate` missing), `design_report` (argument `protocol` missing), `estimand_feasibility` (defaults differ). All consistent with a stale `man/` directory: running `devtools::document()` is the single fix for warnings 4 and 5.
6. Rd usage: `assess_support.Rd` does not document its `balance` argument.

NOTEs: CRAN incoming feasibility (new submission; Suggests `cleanroomGov` and `survtmle` are not in mainstream repositories); hidden file `inst/decision_logs/.gitkeep`; future file timestamps; possible problems in R code (`estimate_ipwhr` has no visible binding for `.weights`; `plot.support_assessment` calls `hist` without importing it from graphics).

The test suite inside the check passed (293 s).

### Tests and coverage (run in this session)

`devtools::test("cleanTMLE")`: **662 passed, 0 failed, 3 skipped, 11 warnings**, wall time 315 s. The warnings are the intentional experimental-scope and superseded-function messages. `covr::package_coverage()`: **total 68.3 %**, wall time 560 s. By file (ascending):

| File | Coverage % |
|---|---:|
| sl_and_sensitivity.R | 0.0 |
| utils.R | 0.0 |
| plasmode_extensions.R | 17.9 |
| integration.R | 23.7 |
| plots.R | 25.6 |
| estimate_effect.R | 49.5 |
| utils_weights.R | 56.8 |
| tmle_clean_room_wrapper.R | 58.2 |
| clean_tmle_staged.R | 61.4 |
| zzz.R | 66.7 |
| tmle.R | 68.3 |
| plasmode_dq.R | 69.4 |
| cleanroom.R | 70.8 |
| estimate_gcomprisk.R | 71.7 |
| sl_learners.R | 72.1 |
| identify.R | 76.6 |
| estimate_aipwrisk.R | 76.6 |
| estimand_estimators.R | 76.6 |
| design_report.R | 78.9 |
| specify_models.R | 79.0 |
| variance_methods.R | 79.9 |
| design_tools.R | 80.4 |
| event_process_helpers.R | 82.7 |
| simulate_support.R | 84.4 |
| stages.R | 85.2 |
| tables.R | 85.9 |
| assess_support.R | 85.9 |
| estimand_ladder.R | 86.4 |
| estimate_ipwhr.R | 89.7 |
| estimate_ipwrisk.R | 92.3 |
| sl_glmnet_bounded.R | 92.3 |
| data.R | 100.0 |
| superseded.R | 100.0 |

The front door `estimate_effect()` at 49.5 % is the notable gap: the single function every analysis is supposed to route through is half untested. There are no tests asserting that design-stage functions never touch the outcome column (the planned `test-blinding.R` does not exist yet).

### Exported surface

`getNamespaceExports("cleanTMLE")` returns exactly **54 exports**, matching the NAMESPACE and the manuscript claim (lines 431 and 3455) and the consolidation commit ("API consolidation: 149 exports to 54"). Classification:

* Core workflow verbs (16): `create_analysis_lock`, `create_contrast_locks`, `declare_estimand_ladder`, `fit_ps`, `assess_support`, `estimand_feasibility`, `simulate_support`, `run_plasmode_feasibility`, `run_plasmode_dq_stress`, `select_tmle_candidate`, `tmle_candidate`, `run_negative_control_ladder`, `check_process_indicators`, `design_report`, `run_estimand_ladder`, `estimate_effect`.
* Lock and masking utilities (5): `mask_outcome`, `unmask_outcome`, `save_lock`, `load_lock`, `estimate_design_precision`.
* Estimator layer and grammar (17): `run_clean_tmle`, `estimate_ipwrisk`, `estimate_gcomprisk`, `estimate_aipwrisk`, `estimate_ipwhr`, `estimate_surv_tmle`, `estimate_lmtp`, `specify_models`, `identify_outcome`, `identify_treatment`, `identify_censoring`, `identify_competing_risk`, `identify_subject`, `identify_interval`, `select_variance_method`, `run_delta_sensitivity`, `compute_evalue`.
* Learners and helpers (9): `SL.glmnet.bounded`, `SL.glmnet.ridge`, `SL.glmnet.enet`, `SL.glm.pca`, `SL.glm.pca5`, `SL.glm.pca10`, `build_sl_library`, `resolve_truncation_rule`, `sanitize_covariates`.
* Reporting and plotting (6): `make_table1`, `make_table2`, `attrition_table`, `love_plot`, `forest_plot`, `clean_weight_diagnostics`.
* Data (1): `sim_func1`.

No exported function emits a superseded message; the 27 superseded governance helpers are unexported and reachable only via `:::` (they emit once-per-session notes when called).

### Reconciliation against the functions vignette

`cleanTMLE-functions.qmd` index tables name **95 distinct functions**. Of these, **62 are not cleanTMLE exports**: the superseded governance layer (`gate_check`, `gate_all`, `authorize_outcome_analysis`, `assert_outcome_authorized`, `create_audit_log`, `record_stage`, `record_checkpoint`, `export_audit_trail`, `record_decision_log_entry`, `export_decision_log`, `checkpoint_*`, `refine_ps_after_nco`, `summarize_dq_degradation`, `summarize_event_support`, and others), internal workers (`fit_ps_glm`, `fit_ps_superlearner`, `run_negative_control_tmle`, `run_matched_tmle`, `fit_tmle_*`, `expand_tmle_candidate_grid`, ...), and six functions exported by the companion package cleanroomGov (`attach_estimand`, `declare_sensitivity_plan`, `build_stage_manifest`, `summarize_stage_path`, `clean_event_process_table`, `clean_check_event_processes`). Conversely, **21 of the 54 exports never appear in the vignette's index tables**, including the entire cumulative-risk grammar (`specify_models`, all seven `identify_*`, the four `estimate_*risk`/`estimate_ipwhr` verbs), `build_sl_library`, `select_variance_method`, `make_table1`, `make_table2`, `forest_plot`, and `SL.glmnet.bounded`. So the manuscript's "54 functions" is correct for the NAMESPACE; the function dictionary both over-covers (documents 62 non-exports without always distinguishing their status) and under-covers (omits 21 exports). Full lists: `scratchpad` file `vignette_index_reconciliation.txt` (regenerate with `revision`-adjacent audit scripts).

### Vignette builds

Both vignettes render from a clean R session with the installed package (Quarto CLI, fresh directory, no cache): `cleanTMLE-staged-analysis.qmd` in **5 min 14 s**, `cleanTMLE-functions.qmd` in **1 min 9 s**, both exit code 0. The as-cran vignette failure above is therefore specific to the check environment, but must still be fixed for a CRAN submission (WP3).

### cleanroomGov

Exists in the workspace at `cleanroomGov/` (version 0.1.0, installed). It exports **9 functions**: `attach_estimand`, `build_stage_manifest`, `clean_check_event_processes`, `clean_event_process_table`, `clean_missing_data_plan`, `clean_risk_report_table`, `clean_target_population`, `declare_sensitivity_plan`, `summarize_stage_path` (plus 4 print methods). That is below the fewer-than-ten threshold at which WP3 recommends folding it back into cleanTMLE.

---

## 2.2 Blinding audit of the code paths

Method: every read of the outcome column across `cleanTMLE/R/` was located by pattern search (`lock$outcome`, `data[[outcome]]`, and variants), each hit was read in context, and the six named items were additionally verified by executing the functions on `sim_func1()` data in this session (probe transcript reproduced below where it matters). Y denotes the primary outcome column of `lock$data`; NC denotes registered negative-control columns.

### Access table (exported functions)

| Function | reads Y | reads A | reads Y by A | reads NC | simulated Y only | Notes |
|---|:--:|:--:|:--:|:--:|:--:|---|
| `create_analysis_lock` | stores | stores | no | stores | no | `lock$data` is the full data frame including Y (in memory) |
| `create_contrast_locks` | stores | stores | no | stores | no | one full-data lock per contrast |
| `mask_outcome` | writes | no | no | no | no | sets column to NA on a copy; column retained; see below |
| `unmask_outcome` | writes | no | no | no | no | copies Y back from `original_lock` |
| `estimate_design_precision` | **yes** | yes | **yes** | no | no | events and crude rates per arm; refuses masked locks; see below |
| `fit_ps` (glm / superlearner / external) | no | yes | no | no | no | A and W only |
| `assess_support` | no | yes | no | no | no | g, A, W only |
| `estimand_feasibility` | no | yes | no | no | no | g, A, W; `unsupported` profile is aggregate SMD rows |
| `simulate_support` (defaults) | marginal mean | yes | no | optional | yes | reads `mean(Y)` for the base rate anchor; `q0_source = "heldout"`/`"primary_outcome"` fit real Y, the latter gated by `allow_outcome_q0` |
| `run_plasmode_feasibility` | **yes (Q0)** | yes | no | no | candidates: yes | fits Q0(W) on real Y; refuses masked locks; returns the lock (with Y) on the result; see below |
| `run_plasmode_dq_stress` | **yes (Q0)** | yes | no | no | candidates: yes | same pattern; returns the lock (with Y) on the result |
| `select_tmle_candidate` | via input | via input | no | no | yes | operates on plasmode metrics, but the input object carries the lock |
| `tmle_candidate` | no | no | no | no | no | specification constructor |
| `run_negative_control_ladder` | no | yes | no | **yes** | no | default `method = "tmle"`; silent unadjusted fallback (see below) |
| `check_process_indicators` | no | yes | no | no | no | indicators, A, conditioning W |
| `design_report` | no | via support | no | via nc_ladder | no | no Y; carries per-patient g and A vectors (see below) |
| `declare_estimand_ladder` | no | no | no | no | no | metadata only |
| `run_estimand_ladder` | **yes** | yes | yes | no | no | Stage 4 by design; note `anyNA(Y)` is read before the outcome guard fires |
| `estimate_effect` | **yes** | yes | yes | no | no | Stage 4 by design; guard via `.check_outcome_access` |
| `run_clean_tmle`, `estimate_*risk`, `estimate_ipwhr`, `estimate_surv_tmle`, `estimate_lmtp` | **yes** | yes | yes | no | no | estimation layer, post-outcome by design |
| `estimate_design_precision`-adjacent unexported `summarize_event_support` | **yes** | yes | **yes** | no | no | see below |
| `attrition_table`, `make_table1`, `love_plot`, `clean_weight_diagnostics`, `sanitize_covariates`, learners, `build_sl_library` | no | some | no | no | no | design-side utilities |
| `compute_evalue`, `run_delta_sensitivity`, `select_variance_method` | no (inputs are estimates) | no | no | no | no | post-outcome sensitivity layer |
| `save_lock` / `load_lock` | serialise | serialise | no | serialise | no | the serialised lock contains whatever the lock contains |

One further internal fact worth recording: `fit_tmle_treatment_mechanism()` (the Stage 4 g-fit) reads the outcome column to choose the cross-validation fold count (`recommend_cv_V(compute_n_eff(data[[lock$outcome]], ...))`, cleanroom.R:2370). It is a Stage 4 function, so no blinding rule is violated, but a treatment-mechanism fitter reading Y is worth knowing during WP3.

### 2.2.1 `estimate_design_precision()` and `summarize_event_support()`

Confirmed: both tabulate events **by arm**. `estimate_design_precision()` computes `events_treated <- sum(Y[A == 1L])` and `events_control` (stages.R:963-964), returns `events_per_arm` and `crude_rates` (crude event rate per arm), and its print method prints "Events (treated)" and "Events (control)" rows (stages.R:1024). `summarize_event_support()` (unexported; attached to the precision object as `$event_support`, stages.R:1004) returns a data frame with one row per arm containing `events` and `event_rate`. Live probe on `sim_func1(500)`:

```
Events (treated)   126      Events (control)   125
dp$crude_rates:  treated 0.4772727   control 0.5296610
summarize_event_support: Treated 264 126 0.4773 / Control 236 125 0.5297
```

Two crude arm rates are the crude treatment-outcome association. This contradicts, verbatim: the staged vignette, "`estimate_design_precision()` reads only marginal event counts (never the treatment-outcome association)" (cleanTMLE-staged-analysis.qmd:594-595); the manuscript's stage table, which labels Stage 1b "Events only" (manuscript:860) and the dossier table's "Marginal Y only" (manuscript:1048); the demo code comment "(marginal event counts only)" (manuscript:1715); and the functions vignette, which labels `summarize_event_support()` "Marginal Y only" in the same row that describes its purpose as "Event counts by arm" (cleanTMLE-functions.qmd:181-182), a self-contradiction. Both functions refuse masked locks, so in the documented workflow they run on the unmasked lock at Stage 1b. The roxygen of `estimate_design_precision` itself says "uses outcome counts for power/precision summaries only, not for estimation", which is accurate about intent and silent about the per-arm split.

### 2.2.2 The plasmode Q0 and what is returned

`run_plasmode_feasibility()` (cleanroom.R:1346-1379) and `run_plasmode_dq_stress()` (plasmode_dq.R:701-717) both fit `Q0(W) = E[Y|W]` on the real outcome (covariate-only logistic GLM by default; a SuperLearner library when `q0_library` is supplied). Neither stores the fitted Q0 object or `p_base` predictions on the returned object. However, **both return the entire lock, including `lock$data` with the real outcome column, as `result$lock`** (cleanroom.R:1590; plasmode_dq.R:1091). Verified live: `identical(plas$lock$data$event_24, dat$event_24)` is `TRUE`. So the question of whether the analyst could recover `E[Y | g(W)]` by propensity stratum is moot at the object level: the returned object contains the raw outcome wholesale. Both functions also refuse masked locks outright (error message: "If the outcome is masked, call unmask_outcome() before ..."), so the plasmode stage runs unblinded at the software level; the run_simulation.R pipeline masks only after the plasmode (run_simulation.R:551-553 and 636-641 confirm this ordering, with a comment acknowledging it).

On leakage under separation: with the default covariate-only GLM Q0, the analyst-visible metrics tables (bias, RMSE, coverage per candidate) do not expose `E[Y|W]` by propensity stratum. The leak path is the returned lock and, in principle, a user calling `predict` on a Q0 they refit themselves; the latter is unavoidable given the design, the former is a fixable object-hygiene defect (WP1 item 2, WP3 item 3).

The manuscript claims a configurable alternative exists: "a stricter reading of outcome blinding would require E[Y|W] to come from external pilot data, and the package supports both modes, recording which mode was used" (manuscript:1010-1013) and "Depending on configured mode, the plasmode generator may use ... or it may take E[Y|W] from external pilot data" (manuscript:1637-1642). **No such mode exists** for `run_plasmode_feasibility()` or `run_plasmode_dq_stress()`: the only knob is `q0_library`, which changes the learner, not the data source, and nothing records a mode. (`simulate_support()` does have `q0_source` with `synthetic`, `negative_control`, `auxiliary`, `heldout`, `primary_outcome` values, the last gated by `allow_outcome_q0 = TRUE`; none is an external-pilot mode.) The manuscript's own Future extensions item 5 lists the DGP-spectrum modes as future work, contradicting the two passages above. Recorded as FINDINGS F4.

### 2.2.3 `mask_outcome()` / `unmask_outcome()`

`mask_outcome()` sets `lock$data[[lock$outcome]] <- NA` on the returned copy (stages.R:1755). The column is **not removed**; it is NA-filled and retained, and `lock$.outcome_masked` is set. Because R copies on modify, the caller's original lock (and any other copy) still holds the real outcome in memory, and `unmask_outcome()` requires exactly that (`original_lock` supplies the column back). `lock$data` is an in-memory copy of the full dataset including the outcome at every stage unless the analyst masked at creation and physically kept the outcome elsewhere. The rescueCo multi-arm pipeline does the stronger thing outside the package: outcome columns live in `multiarm_outcomes.rds` and the locks carry masked placeholders, which is what the manuscript describes for the case study. The package itself provides no `outcome_store` separation (that is WP3 item 3).

### 2.2.4 `run_negative_control_ladder()`

Default is `method = c("tmle", "unadjusted")`, so **"tmle"** (design_tools.R:134-141), with `ps_method = "glm"` refit per rung; the vignettes' explicit `method = "tmle"` calls match the default. The manuscript's case study reports the **unadjusted** ladder as primary ("fits all five registered controls, unadjusted, on each nested cohort", manuscript:3153-3154) with the TMLE variant in the supplement, i.e. the documented primary is the non-default. Two engineering findings: (a) when the per-rung TMLE fit fails, the code **silently falls back** to the unadjusted two-proportion estimate (design_tools.R:191-198) and the output row still says `status = "estimated"` with no marker of which method produced it; (b) there are no prespecified decision criteria on the lock (no null band, no per-domain minimum, no consistency rule); the only rule is `p < alpha` per rung. The superseded `checkpoint_residual_bias` did implement an equivalence/TOST rule with a null band, but nothing in the current path uses it.

### 2.2.5 `refine_ps_after_nco()`

Unexported and superseded (integration.R:1218-1318; the superseded note points to `run_negative_control_ladder()`). What it does: runs the registered NCs with the current PS, creates a **new lock** with the expanded covariate set via `create_analysis_lock()` (which computes a fresh `lock_hash`, so re-fingerprinting does happen implicitly), refits a GLM PS, re-runs the NCs, and returns a before/after comparison. It **re-fits the propensity model after seeing negative-control results** by construction. A deviation record is written **only if** the caller passes an `audit` object (default `NULL`: nothing is recorded anywhere); no design-log entry is written on the lock; `design_report()` has no notion of a "protocol amendment" flag. WP1 item 4 applies as written.

### 2.2.6 `design_report()`

Neither the returned object nor the printed output contains individual-level *rows*: `feasibility$unsupported` (from `who_is_unsupported()`) is aggregate (population x variable: n, means, SMDs), and the printed report is summary-only. However, the returned object embeds the `support_assessment`, which carries **per-patient vectors** `support$g` (propensity scores) and `support$A` (treatment indicators), length n (verified live: length 500 on a 500-row lock; assess_support.R:347-348, kept for `plot()`). These are individual-level data on two variables, without covariates or outcome. Whether that violates the Muntner summary-statistics-only requirement is a judgment for WP1 item 6; the audit records the fact. The report does not embed the lock or the raw data frame.

---

## 2.3 Gate-logic audit

### What actually produces GO / FLAG / STOP

Three generations of verdict machinery coexist:

1. **Superseded checkpoint layer** (used by `run_simulation.R`): `checkpoint_cohort_adequacy`, `checkpoint_balance`, `checkpoint_residual_bias` each yield GO/FLAG/STOP; `gate_all()`/`authorize_outcome_analysis()` combine them with STOP > FLAG > GO precedence, plus a safety net that any recorded STOP checkpoint blocks authorisation. `gate_check()` evaluates only the **baseline** plasmode row of the selected candidate (bias, coverage, RMSE, SE/SD window). `gate_dq()` converts DQ degraded rows into a checkpoint (STOP if any degraded row for the candidate breaches `max_abs_bias = 0.02`, `min_coverage = 0.85`, or `rmse_ratio > 1.5`).
2. **Current design-stage path** (0.2.0): `assess_support()` PASS/FLAG/SEVERE/FAIL, `estimand_feasibility()` per-estimand verdicts, `simulate_support()` pass rule (`|bias| <= tolerance` and `coverage >= floor` per cell), combined narratively by `design_report()` into a recommendation string. No numeric DQ-based gate exists on this path.
3. The DQ metrics themselves, which are tables, not decisions.

### The Scenario A GO versus Scenario C STOP question

The recorded artifacts in `results_new/` (read in this session; `done_*.rds`, `audit_*.csv`, `decision_log_*.csv`) show that **the executed pipeline returned GO and authorised outcome access in all four scenarios, including Scenario C**:

```
unmeasured_conf: cp1 GO, cp2 GO, cp3 GO (nc_outcome est -0.016, p 0.44, not flagged),
Composite Gate GO, "Outcome unblinding authorised by gate. All required checkpoints passed."
```

The reasons are structural, not accidental:

* In `run_simulation.R` the pre-outcome gate is `gate_all(cp1, cp2, cp3)` (line 664-665): cohort adequacy, balance, and the single negative control. **The DQ stress results never enter the gate**; `gate_dq()` is never called anywhere in the script; the DQ run is saved and logged as "executed" but produces no checkpoint. `gate_check()` on the baseline plasmode (line 593) is printed but also never recorded as a checkpoint.
* The negative control cannot detect Scenario C's confounder even in principle: in `dgp_scenarios.R:64-66` `nc_outcome` depends only on age, sex, and biomarker, with no U term and no treatment term, so its weighted association with treatment is null in truth in every scenario. In the realised reference dataset it was not flagged (p = 0.44).
* With a single registered NC, `checkpoint_residual_bias` can only return GO or STOP (one flag exceeds m/2), never FLAG.

Against this, the manuscript asserts: "the prespecified pre-outcome decision correctly returned STOP" (abstract, line 176-181); "the staged decision rule returns STOP and the locked primary analysis does not run" (Section 1.5, line 508-510); "The pre-outcome decision returns STOP for Scenario C" (Interpretation, line 2121); "**STOP** under the prespecified unmeasured-U DQ scenario; the locked analysis does not run as planned ... The workflow blocks publication of the misleading null" (@tbl-workflow-contrast, line 2396); "the workflow's prespecified rule refuses to authorise the Stage 4 analysis once the unmeasured-U row of the DQ stress test exceeds the locked bias threshold" (line 2400). The verdicts in these sentences are asserted in prose; no chunk in the manuscript computes them.

**Answer to the (a)/(b)/(c) question.** The Scenario C STOP is (b) and (c) combined:

* (b) The textual justification leans on the realised-bias oracle: the quoted evidence is the Monte Carlo mean bias of the realised TMLE against the true effect (`scC_bias` = 0.037, coverage 0.42; setup chunk lines 105-119 read them from `simulation_results.rds`), quantities that require the true RD and 200 unblinded replicates and are unavailable pre-outcome in practice.
* (c) It is not the output of any executed rule: the recorded pipeline authorised Scenario C.
* Decisively, **the same locked rule applied to the recorded DQ rows returns STOP for Scenario A as well.** The manuscript locks `max_abs_bias = 0.02, min_coverage = 0.90` (lines 2090-2097). Recorded DQ rows for the selected candidates: Scenario A (candidate `aggressive`) unmeasured_U bias 0.0103 / 0.0264 / 0.0417 / 0.0756 / 0.0992 at OR 2 / 3 / 4 / 6 / 8, coverage 0.90 / 0.60 / 0.30 / 0.00 / 0.00; Scenario C (candidate `robust`) 0.0135 / 0.0325 / 0.0534 / 0.0850 / 0.1033, coverage 0.833 / 0.433 / 0.267 / 0.00 / 0.00. Every cell from OR 3.0 upward breaches the bias bound in **both** scenarios (the prompt's four quoted numbers are exactly these cells). There is no restriction of the grid under which a locked bias rule separates A (GO) from C (STOP): at OR 2.0 both pass on bias; at OR 3.0 and above both fail. The only cell-level difference is A's OR 2.0 coverage 0.90 versus C's 0.833, a two-replicates-in-thirty difference at 30 DQ replicates (binomial MC SE about 0.06), far inside noise. The DQ tables also compare different candidates across scenarios (A selected `aggressive`, C selected `robust`), so even the appearance of a comparison is confounded by the selection.

Consequently the abstract, Section 1.5, the Interpretation subsection, @tbl-workflow-contrast, and presentation slides "Built-in confounding stops the workflow" and "Same fits, different pipeline" overclaim, exactly as anticipated by the revision prompt. The exact sentences are catalogued in `revision/FINDINGS.md` (F1, F2). One nuance in the manuscript's favour: several passages already state the correct reading ("A STOP decision is a failure of prespecified operating-characteristic thresholds under the simulated DGP, not detected real-data bias", lines 2073-2085), and the DQ verdict as a tipping-point statement ("the workflow flags a prespecified unmeasured-confounding scenario under which the planned analysis would no longer meet the locked operating-characteristic thresholds", abstract). The defect is that the surrounding sentences then claim an executed asymmetric STOP that the recorded pipeline did not produce and that no locked rule reproduces. The honest framing (WP2 item 2) is already half-present and needs to displace the detection framing everywhere, with the additional disclosure that under the locked grid the tipping point is crossed in Scenario A as well.

### Section 11 (gate operating characteristics)

`run_gate_operating_characteristics.R` injects a what-if confounder of strength s into the plasmode screen at exactly the strength s of the true confounder in the oracle arm (`screen_bias(d, s, ...)`, lines 89-97; per-cell STOP at `|bias| > 0.02`, lines 148-150). Sensitivity is therefore P(screen mean bias crosses the threshold | realised mean bias crosses it) **at matched, known threat strength**. On the linear surface this is precisely the statement that a plasmode with an injected OR-s confounder produces (approximately) the same bias as real data with an OR-s confounder: a fidelity property of the generator, not a detection property of a screen that must find an unknown threat. The manuscript already concedes the operative caveats (specificity rests on two GO cells; the nonlinear s = 1 cell STOPs from nuisance misspecification alone with realised bias about 0.07, lines 2861-2887), but retains the sensitivity/specificity vocabulary. Reframing per WP2 item 3: present Section 11 as plasmode fidelity (does predicted bias match realised bias for a matched threat, and how does the match degrade under nuisance misspecification), and either drop sensitivity/specificity or keep it only with the explicit statement that the decision problem is defined over the named threat grid, not over the data. Also note a provenance detail: the script writes `results/gate_oc.rds` (config line 40) while the manuscript reads `results_new/gate_oc.rds`; the artifact (mtime 2026-05-26) was evidently moved by hand.

---

## 2.4 Manuscript numeric audit

### Machine-written versus hand-typed numbers

The manuscript is better instrumented than the prompt feared in some places and worse in others. The four-scenario Monte Carlo tables, the DQ degradation tables, the gate-OC table, the variance table (Table 5 analogue, `tbl-boot-variance`), and all case-study design tables are computed in chunks from named result files at render time, and the setup chunk hard-fails when `bootstrap_variance.csv` is missing rather than substituting placeholders. The hand-typed numbers concentrate in: Section 12.6 (the Stage 2b candidate table and the worst-case RMSE ratios in prose), the simulation-design prose (Scenario C parameters), Scenario D interpretive prose (bias/coverage values restated from the tables), the case-study narrative (effect sizes restated from `ladder_estimates.csv`), and the presentation prose.

### SOURCE_MAP and UNMATCHED files

`revision/SOURCE_MAP.csv` (429 rows) and `revision/UNMATCHED_NUMBERS.csv` were machine-generated in this session by an extraction-and-matching script (scratchpad `audit_numbers.R`; to be superseded by `revision/build_source_map.R` in WP0). Method: literal numeric tokens were extracted from the manuscript and presentation sources outside code chunks (inline `r` expressions are machine-computed by construction and were excluded); tokens were classified (years, citations, section references, versions, small counts, declared parameters, candidates); candidates and declared parameters were matched against 155,515 numeric values pooled from `results_new/`, `rescueCo/results/` (including `multiarm/`), `sandbox/*/`, and `plasmode_selection_paper/results/`, with rounding tolerance and x100 variants. Interpretation caveat, stated plainly: a **match asserts only that a numerically identical value exists in some result file**. Each row carries a `match_strength` grade; of the manuscript's 348 matched tokens, 96 are strong (3 or more decimals, or magnitude at least 100), 8 moderate, and 244 weak (low-precision values that collide with the pool by arithmetic necessity). `UNMATCHED_NUMBERS.csv` is empty, which reflects the permissiveness of value-level matching for low-precision tokens, not verified provenance; the numbers whose provenance targeted checking could **not** establish are itemised in FINDINGS F6 and belong in any honest reading of the unmatched set. WP0's `build_source_map.R` must therefore match at claim level (this table cell comes from this field of this file), not value level.

### Word count

Parsed from the Quarto source with code chunks, the abstract, headings, pipe tables, table captions, and image lines removed, and inline expressions counted as one token each: **22,207 body words** (23,661 including the abstract and headings). The 33,800 figure circulating previously evidently counted rendered front and back matter. Either way the body is roughly 2.5 times the 9,000-word WP2 target.

### Internal inconsistencies (manuscript, vignettes, slides, artifacts)

1. **Scenario C design description versus executed configuration.** The manuscript specifies U prevalence 0.20 and OR 2.0 on both A and Y (lines 1973-1975); the script and every recorded artifact use prevalence 0.40 and OR 3.0 (`run_simulation.R:146-153`; `done_unmeasured_conf.rds`). The July recalibration (see `sandbox/scenarioC_recalibration/`, `results_new/_regen_scenarioC.log`) was never propagated to the design prose. The E-value comparison chunk also prints "True RD (U present, OR = 2.0)" (line 2764) for the OR 3.0 scenario.
2. **Scenario C STOP claims versus recorded GO** (Section 2.3 above; FINDINGS F1).
3. **Same locked DQ rule fails Scenario A** and the manuscript presents A as the clean GO benchmark without noting its own OR 3.0 and 4.0 rows breach the same threshold (FINDINGS F2). The E-value section even prints Scenario A's unmeasured-U bias ladder (0.0103 to 0.0992) without remarking that the locked rule would flag it.
4. **Slide "Marginal overlap strains the variance" versus Table 5.** The slide says matched TMLE intervals are too wide under marginal overlap, "roughly doubling the true spread", and that both IPTW and matched TMLE are "fixable with a bootstrap variance". `results_new/bootstrap_variance.csv`: marginal-overlap Match_TMLE SE/SD (IF) 1.041, coverage 0.97 (approximately calibrated, not doubled); marginal IPTW SE/SD 1.229 with bootstrap 1.220 (the bootstrap does not fix it, as the manuscript itself states at lines 2192-2210); the anti-conservative matched-TMLE regime is **good** overlap (SE/SD 0.931, coverage 0.89), where the re-matching bootstrap raises coverage to 0.95. Three distinct contradictions on one slide (FINDINGS F9). The "roughly doubling" phrase loosely matches the four-scenario run's marginal Match_TMLE SE/SD of 1.452, a different study from the one the manuscript's variance section reports, and 1.45 is not 2.
5. **Time to arrival.** Section 12.7 states time to arrival "is estimated here only through bounded-weight estimands" (lines 3339-3342). `ladder_estimates.csv` contains ATE rows for time_to_arrival_minutes on C1 (the reconciliation override), C3, and C4; only C2 lacks one. The slide "The disaster outcome, re-analysed" prints the full ladder including the C1 ATE override row and its caption text acknowledges it, so the slide matches the artifact and the manuscript sentence does not (FINDINGS F8).
6. **Section 12.6 candidate table is unverifiable and its artifacts were overwritten.** The Stage 2b table (glm_t01 +0.00001 / 0.0134 / 0.92 and three more candidates) matches no file in the repository: `plasmode_metrics.csv` holds a two-candidate 50-rep run with different values (glm_t01 bias 0.00016, RMSE 0.01306, coverage 0.94), and `plasmode_summary.csv` and both `*_dq_stress.rds` files hold a four-candidate **3-replicate smoke run** (mtimes 2026-07-29). The worst-case RMSE ratios in prose (1.82 / 1.73 / 1.53 / 1.27) trace **only** to `_dq_degradation_50rep_backup.csv` (1.818 / 1.731 / 1.525 / 1.268; mtime 2026-06-07), while the reproducibility table points to `plasmode_dq_degradation.csv`, which now contains 3-rep values (worst ratios 1.14 to 1.22). Additionally the prose attributes every candidate's worst case to near-positivity at slope x4; the backup shows ensemble_t05's worst at slope x3. And the manuscript states Stage 2b "feasibility coverage of 0.94 ... across 50 replicates" (line 3387) three paragraphs after the table showing 0.92 for all four candidates; 0.94 matches `plasmode_metrics.csv` and 0.92 matches nothing on disk (FINDINGS F6).
7. **Case-study package-version claim versus released metadata.** The manuscript: "All stages were run on cleanTMLE 0.2.0, with the package version and the lock fingerprints recorded in the released case-study metadata" (lines 2913-2916). The only released metadata file, `rescueCo/results/manuscript_artifacts/case_study_metadata.json`, records `"cleanTMLE_version": "0.1.5"` (the earlier single-contrast run, rendered 2026-07-29); the multi-arm tracked artifacts (`lock_summary.csv` and the other CSVs) carry neither a version nor a fingerprint. Also, Section 12.6's candidate selection and DQ stress come from that earlier 0.1.x-era run by the manuscript's own account (FINDINGS F7).
8. **Simulation artifacts predate the package the paper documents.** `results_new/` summaries: Scenarios A, B, D completed 2026-06-08/09; Scenario C on 2026-07-11. The 0.2.0 consolidation happened in September; `run_simulation.R` in the working tree (modified 2026-09-13, now calling superseded functions via `:::`) is not the exact script that produced the July artifacts. The four-scenario results were produced by cleanTMLE 0.1.x.
9. **"54 functions"** is correct against the NAMESPACE; the functions vignette's index names 95 (Section 2.1). No manuscript correction needed on the count itself; the vignette needs status labelling.
10. **Kent (2026).** `kent2026fiord` is an unpublished conference presentation ("Conference presentation, 4th Annual FIORD ... Day 1 Lunch session", references.bib:56-62), cited four times, including for quantitative results of the evolocumab example (roughly -4 percentage points at four years, lines 3595-3604). The peer-review-track source exists: Jin R, Hurwitz K, Kent S, et al., Circulation 2025;152(Suppl_3): Abstract 4363106 (verified resolvable at ahajournals.org). Whether the staged-decision details (negative-control flag, review-team conditional proceed) appear in the abstract rather than only in the talk needs verification before the citation is swapped: [[VERIFY]] (FINDINGS F10).
11. **Version 0.2.1 reference.** "cleanTMLE 0.2.1 adds a q0_library argument" (line 3686) while DESCRIPTION and the rest of the paper say 0.2.0 (the argument exists in 0.2.0).
12. **Claimed but absent plasmode DGP modes** (Section 2.2.2; FINDINGS F4).
13. **Stage table and dossier table blinding labels** for Stage 1b (Section 2.2.1; FINDINGS F3).
14. **Threshold provenance prose.** Section "What the DQ supplement adds" claims the decision-threshold values are recorded "as part of the lock fingerprint" (line 1582); in code the base `lock_hash` never covers thresholds; the superseded `attach_decision_thresholds()` records a separate `thresholds_hash`, and the current path records thresholds only inside the objects that used them. It also names an "ESS floor at 30% of n" default that exists nowhere in the code (the superseded balance checkpoint uses 50%).
15. **`gate_oc.rds` script/output path mismatch** (Section 2.3).
16. **`run_simulation_full.R`** is a stale variant (unbounded SL.glmnet, no environment overrides, 425 diff lines against `run_simulation.R`); nothing marks which is canonical.

### Terminology and style counts (run in this session)

| Document | em-dashes | contractions | "synthetic" | "rather than" | banned vocabulary |
|---|---:|---:|---:|---:|---|
| manuscript | 0 | 0 | 64 | 62 | 1 ("underscore-prefixed scripts", a filename description) |
| presentation | 0 | 0 | 4 | 2 | 0 |
| staged vignette | 0 | 0 | 8 | 7 | 0 |
| functions vignette | 0 | 0 | 7 | 2 | 0 |
| workflow_contrast_c | 0 | 0 | 0 | 0 | 0 |
| cleanTMLE/R (code + roxygen) | - | - | 38 | - | - |

Em-dashes and contractions are already at zero everywhere, and the banned-vocabulary list is effectively absent. The two systematic style items are "synthetic" (121 occurrences across manuscript, vignettes, slides, and roxygen, plus DESCRIPTION) and the corrective "X rather than Y" construction at 62 manuscript occurrences (the "not X, but Y" comma form appears once by strict pattern; reversed forms such as "a failure of ..., not ..." add several more). The "supplement, not substitute" caveat appears in 7 distinct manuscript passages (lines 257, 301, 391, 480, 961, 1303, 3517), close to the 8 locations the prompt lists; WP2 item 9 consolidation applies.

---

## 2.5 Case-study provenance

**Which version and scripts produced Section 12.** The multi-arm design and estimation artifacts (`support_by_contrast.csv`, `estimand_feasibility.csv`, `nc_ladder.csv`, `ladder_estimates.csv`, support maps, design comparisons) were produced 2026-09-13/14 by `rescueCo/scripts/10-14_multiarm_*.R` under the 0.2.0-era package (design stage 14:43-18:05, estimation 18:21-22:50 for the 10 contrast-outcome pairs on 4 workers, reconciliation and render overnight; `design_run.log`, `estimation_run1.log`, `stage14_run.log`). The final step's attempt to render the companion case-study report failed (`quarto ... had status 1`, stage14_run.log). The candidate-selection and DQ subsection (12.6) comes from the earlier single-contrast run under cleanTMLE 0.1.5 (per `case_study_metadata.json`), as the manuscript states in general terms; Section 2.4 item 6 records that its numbers are no longer reproducible from the repository. Reconciliation table verdicts (`build_reconciliation_table.R`, 275 rows): agree 17, benign_difference 10, missing_in_cleanroom 219, missing_in_main 25, not_comparable 4.

**DQ stress replicates and rerun cost.** Confirmed at 50 replicates for the reported numbers, surviving only in `_dq_degradation_50rep_backup.csv`. Timeline from `rescueCo/logs/pipeline.log`: the 50-rep run on 2026-06-04 took 27 minutes of feasibility plus **85 minutes for the five-threat DQ sweep** (10:12 to 11:36). A 200-replicate run (`dq_reps=200, feas_reps=50`) was launched 2026-06-07 07:29, at which moment the 50-rep backup was written; the log records the sweep starting 07:49 and **no completion line**, so it appears to have died. A 3-replicate smoke run on 2026-07-29 then overwrote `plasmode_dq_stress.rds` and `plasmode_dq_degradation.csv`. Estimated wall time for the 200-replicate rerun that WP2 item 7 contemplates: DQ sweep about 4 x 85 = 340 minutes, plus the 50-rep feasibility (about 27 minutes), so **roughly six hours** on this machine. `rescueco_dq_full.R` already defaults to `dq_reps = 200`, so the rerun is one command; it requires the restricted registry data on disk (present on this machine, per the Stage 1/2 artifacts it reuses).

**Ethics placeholders.** Still unfilled, verbatim at manuscript lines 2928-2931: "The study was approved by [approving ethics committee and approval number to be completed], and participation was covered by [informed consent / an approved waiver of consent to be completed]." Nothing in the repository supplies these values (searched the rescueCo config, READMEs, and docs).

**Personnel framing.** The record itself establishes that the case-study analyst is the author of the parallel unblinded analysis: the manuscript reconciles against "the parallel unblinded analysis" throughout (lines 198, 3329, 3345, 3465), and `design_log_notes.md` records the design-stage decisions as made by "A. Mertens / automated pipeline assistant". The multi-arm run is therefore a re-execution and reconciliation demonstration under data-level masking, by an analyst who had previously seen the outcomes of the same registry; it is not personnel blinding, and the software cannot make it so. Sentences that could be read as implying otherwise and need the WP2 item 7 disclosure next to them: "The outcome columns are physically absent from the design stage: the locks carry masked placeholders, and the outcome store is joined only at estimation" (2991-2994); "Every recommendation was fixed, and the design log exported, before the outcome store was joined" (3198-3199); "All stages were run on cleanTMLE 0.2.0 ..." (2913-2914); the walkthrough claim "Everything before run_estimand_ladder() runs with no outcome access, and the outcome columns are physically absent from the design stage in the case study" (622-624); and the slide "Every verdict, declaration, and switch was fixed and logged before the outcome store was joined" (presentation:427). Each is true as a statement about the software record and silent about the analyst's prior outcome knowledge; one honest sentence in 12.1 fixes all of them. It is worth adding on the positive side: the design-log notes were recorded before Stage 13 ran in the September rebuild, and the manuscript already reports the earlier run's silently overridden STOP candidly (lines 3360-3383).

---

## 2.6 Three-reviewer preliminary read

### Reviewer A: scientific R programmer

The package is in better shape than most research software (662 passing tests, a real consolidation behind it, informative error messages, a thoughtful callr-based timeout guard in the DQ loop), but it is not submittable to CRAN and its central guarantees are asserted rather than tested. Five problems in priority order. (1) The blinding contract is untested and partly false: nothing asserts that design-stage functions ignore the outcome column, and the two violations this audit confirmed (`estimate_design_precision`/`summarize_event_support` returning per-arm crude rates, stages.R:963-1111; `run_plasmode_*` returning the unmasked lock on their results, cleanroom.R:1590, plasmode_dq.R:1091) would have been caught by the planned `test-blinding.R` on day one. The lock also never fingerprints the data values, only dimensions and column names (`.compute_lock_hash`, cleanroom.R:186-195), so "the same lock" can be validated against different data of the same shape, which undercuts the reconciliation story. (2) `R CMD check --as-cran`: 2 errors, 6 warnings, 4 notes; most warnings are one stale-`man/` regeneration away, but broken examples in superseded Rd files (assert_outcome_authorized calling the unexported `create_audit_log`) show the consolidation was not finished, and `estimate_effect`, the advertised front door, has no Rd and 49.5 % coverage. (3) `set.seed()` is called inside at least six package functions (`run_plasmode_feasibility` cleanroom.R:1405, `run_plasmode_dq_stress` plasmode_dq.R:916, `simulate_support` simulate_support.R:372/427/661, `estimand_feasibility` estimand_ladder.R:127, `run_estimand_ladder` estimand_ladder.R:403/462, `estimate_effect` estimate_effect.R:122/173), clobbering the caller's RNG state; `withr::with_seed` or explicit restore is the fix, and `estimand_feasibility`'s hard-coded `set.seed(1)` is a bug hiding as reproducibility. (4) Parallelism is ad hoc: `fit_ps_parallel` manages its own cluster, the DQ loop is serial with a callr session, and the analysis scripts build PSOCK clusters by hand; one future/furrr path would replace all three. (5) The dead weight is still on board: 62 unexported functions still documented in the vignette index, `run_simulation_full.R` as an unmarked stale twin, `tmle_clean_room_wrapper.R` at 58 % coverage, and `sl_and_sensitivity.R`/`utils.R` at 0 %. On the cleanroomGov split: nine exports, one consumer, and the manuscript never mentions it as a separate package; fold it back.

### Reviewer B: applied pharmacoepidemiologist

The package takes Muntner et al. seriously and gets several things right that most software ignores (graded verdicts phrased as review-team recommendations, the design log written before estimation, negative controls declared with domains at lock time, inestimable controls reported rather than dropped). Five problems. (1) Check Point 3 as executed is not Muntner Section 11: the primary negative-control reading in the case study is the unadjusted ladder, while Muntner specifies negative controls analysed in the final adjusted population with the primary-analysis estimator; the package's own default (`method = "tmle"`) does this, so the manuscript's primary/supplement assignment is backwards, and the silent TMLE-to-unadjusted fallback inside the ladder (design_tools.R:191-198) means a reader cannot even tell which rows are adjusted. (2) No prespecified negative-control decision criteria exist on the lock: the rule is per-rung p < 0.05, exactly the unstandardised practice Tan 2024 and Li 2026 document; the equivalence-band machinery already written into the superseded `checkpoint_residual_bias` shows the authors know better, and `nc_criteria` (WP1 item 5) is the fix. (3) Role separation exists only as prose: the software records no roles, and the case study's analyst previously ran the unblinded analysis of the same registry, so the demonstration is re-execution, not blinding; Section 12.1 must say so and a minimal `roles` field (WP1 item 7) should record the personnel structure the framework assumes. (4) Deviations are not loggable in a reviewable form: `refine_ps_after_nco()` re-fits the propensity after seeing NC results and records nothing unless handed an optional audit object, overrides in `run_estimand_ladder` are logged as free text, and nothing exports a Muntner Table S1-shaped decision log; the earlier run's silently overridden STOP, candidly reported in the manuscript, is the case in point. (5) The verdict thresholds (1/25/60 % outside the band; weights 30/150/1000; DQ 0.02/0.85/1.5) are presented as defaults with, at best, a two-fit anecdotal calibration recorded in a roxygen comment (assess_support.R:39-46); the manuscript should either give their basis or say plainly they are conventions, and the estimand ladder needs the E9(R1) paragraph (WP2 item 6): the software currently executes the switch itself at the declared trigger, whereas the framework's own governance story requires the masked review team to approve it.

### Reviewer C: causal inference and targeted learning methodologist

The estimability layer is genuinely useful and the generate-treatment correction is exactly right, but the paper's headline decision claim does not survive contact with its own artifacts, and several statistical presentation choices need repair. Five problems. (1) The Scenario C STOP is not produced by any locked, outcome-blind rule: the executed gate authorised all four scenarios, the asserted rule (bias 0.02 on the DQ rows) also fails Scenario A, and the quoted evidence is realised bias against the true effect, an oracle. The tipping-point reframing (WP2 item 2) is the only defensible presentation, and it is a perfectly publishable one. (2) "TMLE" in Scenarios A to C is TMLE with GLM-only nuisances (candidates constructed with g_library "SL.glm" and q_library defaulting to g_library, cleanroom.R:1002; run_simulation.R:202-209), with two-fold cross-fitting for TMLE_CF; a targeted-learning package should label the row "TMLE (GLM nuisances)" everywhere (the support-map labels already do this), and either justify V = 2 against the customary 5 to 10 or rerun. (3) The IPTW variance discussion is internally inconsistent with the literature it cites: treating the estimated propensity as known makes the IF variance conservative (Robins/Rotnitzky, Lunceford and Davidian 2004), the sandwich correction is standard and implementable in an afternoon, and the manuscript both reports the over-coverage as a finding and defers the standard fix to future work while the slide claims a bootstrap fixes it, which the manuscript's own table refutes. (4) Scoring truncated candidates against the fixed marginal-RD truth mixes nuisance quality with estimand shift; the manuscript concedes this in two places (lines 2717-2732) but the selection rules still rank on it; the grid should separate estimand-preserving from estimand-shifting knobs, or score truncated candidates against their own implied estimand, or route heavy truncation to the ATO explicitly. (5) The DQ stress test is a prospective, parametric quantitative bias analysis and should be named and cited as such (Lash et al.), which would also clarify its relation to the support map: the map's harshest confounding rung and the DQ unmeasured-U threat are the same perturbation direction applied through two devices, as the design-log note of 2026-09-13 already recognises. Additional items for WP2: Monte Carlo SEs exist in the four-scenario tables but coverage claims from 30-replicate DQ cells and 60 to 100-replicate variance rows need explicit exploratory labelling; the ATT complete-case and IPCW statements are correct as written; add Crump 2009 (present), Sturmer trimming, and Li-Morgan-Zaslavsky 2018 (present) alongside Petersen 2012 in the positivity paragraph, plus the missing AIPW/double-robustness and MCAR/MAR/MNAR citations flagged by TODO comments in the source (lines 540, 828).

---

## Audit deliverables written in this session

`revision/STAGE0_AUDIT.md` (this file), `revision/FINDINGS.md`, `revision/DECISIONS_PENDING.md`, `revision/SOURCE_MAP.csv`, `revision/UNMATCHED_NUMBERS.csv`. Audit scripts and logs remain in the session scratchpad and will be superseded by the WP0 tooling (`style_check.R`, `build_source_map.R`, `render_all.R`, `test-blinding.R`). Nothing in the repository outside `revision/` was modified, and nothing has been committed.
