# API consolidation review (2026-09-13)

The count at the moment of this review: 149 exported functions, 61 registered
S3 methods, about 55 internal helpers, 265 function definitions in R/. The
honest diagnosis is that the package accreted one function per idea for four
release cycles, and the September revision continued the pattern: it added
seventeen net new exports while soft-deprecating rather than deleting. Three
forces drove the growth. Organizational functions were written to record and
declare rather than to compute (the checkpoint, gate, audit, decision-log and
threshold families, 38 exports at peak, most of which the package's own case
study bypassed). Families were written as sibling functions where one
function with an argument would do (four ways to fit a propensity score,
seven estimator entry points, three lock constructors, two love plots, two
decision logs that never talked to each other). And helpers were exported by
default when internal would have served (the dt_* extractors, validators,
one-line wrappers). This review names every cut, merge and internalization,
and the target API.

## 1. Delete: organizational, not computational

These functions exist to record process, not to compute statistics. The
0.2.0 position is that recording belongs in the lock's design log (a plain
data frame the workflow functions append to), in cleanroomGov for formatted
reporting, or in the analyst's protocol, and that the one honest software
enforcement is outcome masking. Everything below leaves the public API. What
a repo driver still references survives as an internal (reachable by
`cleanTMLE:::` from the in-repo legacy scripts and run_simulation.R, which
get patched); what nothing references goes to attic/.

The checkpoint and gate family: checkpoint_cohort_adequacy,
checkpoint_balance, checkpoint_residual_bias, checkpoint_weights,
new_checkpoint, record_checkpoint, gate_all, gate_check, gate_dq. The
verdicts now live on the objects themselves (assess_support,
estimand_feasibility, simulate_support results). The audit family:
create_audit_log, record_stage, export_audit_trail, save_audit, load_audit.
The token layer: authorize_outcome_analysis, assert_outcome_authorized, and
the enforced two-pass entry (run_clean_tmle_preoutcome,
run_clean_tmle_primary, build_dossier); the two-pass path was the
enforcement showcase, and the enforcement is exactly what the case study
routed around. Both decision logs: init_decision_log, log_decision_entry,
save_decision_log, record_decision_log_entry, export_decision_log. The
threshold registry: decision_thresholds, attach_decision_thresholds,
dt_cohort, dt_balance, dt_plasmode, dt_dq, dt_nco. The stragglers:
print_locked_spec, create_analysis_lock_from_yaml, validate_superlearner_spec,
get_final_cohort, refine_ps_after_nco, run_residual_confounding_stage,
run_negative_control (IPTW version), sensitivity_truncation,
tipping_point_sensitivity, run_positivity_diagnostics (absorbed),
summarize_cleanroom_results, summarize_plasmode_results, fit_final_workflows,
fit_tmle_candidate_set, make_wt_summary_table, extreme_weights,
inspect_ipw_weights.

Cost of the cut: the three manuscripts describe and execute the old layer
(gate_all 18 call sites in the main draft alone), so they do not re-render
until their 0.2 rewrite, which is happening anyway; the rendered .docx and
.html stand in the meantime. run_simulation.R and the legacy rescueCo stage
scripts are patched with ::: references and marked legacy.

## 2. Merge: one function, arguments instead of siblings

- fit_ps(lock, method = c("superlearner", "glm", "external"), scores,
  truncate, cv_folds, cluster) absorbs fit_ps_superlearner, fit_ps_glm,
  fit_ps_parallel (its cluster argument) and wrap_ps_fit (method
  "external"). The four become internal workers.
- create_analysis_lock(..., negative_controls = NULL, mask = FALSE,
  enforce = FALSE) absorbs create_simple_lock (enforcement and masking are
  arguments, not constructors) and define_negative_control (registration at
  creation; the ladder also accepts controls directly).
- assess_support(ps_fit, thresholds = list(...), balance = TRUE) absorbs
  support_thresholds (a plain named list with documented defaults),
  compute_ps_diagnostics (the SMD and ESS balance table becomes $balance)
  and run_positivity_diagnostics. One object now answers the whole
  design-diagnostics question, and love_plot() reads it.
- estimand_feasibility(ps_fit, profile_vars = NULL) absorbs
  who_is_unsupported (the removed-versus-kept profile becomes $unsupported
  when profile_vars is given).
- simulate_support(lock, ps_fit, surface = list(...), design =
  c("generate_treatment", "sample_treatment", "parametric_bootstrap"))
  absorbs support_surfaces (a plain list) and check_locked_estimator
  (Petersen's parametric bootstrap is the third design, labelled optimistic
  in the output rather than in a separate function).
- estimate_effect(lock, ps_fit, estimand = c("ATE", "ATT", "ATO",
  "trimmed_ATE", "matched_ATT"), estimator = c("tmle", "iptw", "match",
  "crude"), missing = c("auto", "ipcw", "complete_case")) is the one
  estimation front door. It absorbs run_att_tmle, estimate_ato,
  run_trimmed_tmle, run_ipcw_tmle, run_iptw_workflow, run_match_workflow,
  run_crude_workflow, run_matched_tmle and estimate_tmle_risk_point, all of
  which become internal workers, and attaches the implausibility guard to
  every result. run_estimand_ladder() and run_clean_tmle() dispatch through
  it. The modular quartet (fit_tmle_treatment_mechanism,
  fit_tmle_outcome_mechanism, run_tmle_targeting_step,
  extract_tmle_estimate) becomes internal; estimate_effect(return_steps =
  TRUE) exposes the pieces for teaching.
- declare_estimand_ladder(..., candidate = NULL) absorbs
  lock_primary_tmle_spec and get_primary_tmle_spec: the selected candidate
  rides on the ladder declaration.
- run_plasmode_feasibility(candidates = NULL) builds the default grid
  itself; expand_tmle_candidate_grid and validate_tmle_candidates become
  internal. run_plasmode_dq_stress(scenarios = "regulatory_standard")
  accepts the preset name; default_dq_scenarios becomes internal;
  summarize_dq_degradation becomes summary() on the result;
  assess_dgp_fidelity becomes internal and its table rides on the result.
- run_negative_control_ladder() absorbs run_negative_control_tmle (internal
  worker).
- love_plot(x, matched = NULL) absorbs love_plot_threeway and
  compute_matched_smds; plot(tmle_fit, type = c("summary", "ic", "clever"))
  absorbs ic_histogram and clever_covariate_plot.
- build_sl_library(role, y = NULL, ...) computes the effective sample size
  and the fold recommendation itself; compute_n_eff and recommend_cv_V
  become internal. run_delta_sensitivity() reports the tipping value that
  compute_G_value computed; estimate_design_precision() carries
  summarize_event_support's table as $event_support. design_report(...,
  protocol = NULL) carries the TARGET table as $emulation; emulation_table
  becomes internal. select_variance_method() keeps bootstrap_rd_variance as
  its internal engine.
- expit and logit become internal (stats::plogis and stats::qlogis exist).

## 3. The target public API: 55 exports, a 15-function core

Core workflow (15): create_analysis_lock, mask_outcome, unmask_outcome,
save_lock, load_lock, declare_estimand_ladder, fit_ps, assess_support,
estimand_feasibility, simulate_support, estimate_effect,
run_estimand_ladder, run_clean_tmle, design_report, create_contrast_locks.

Design checks (2): run_negative_control_ladder, check_process_indicators.

Outcome-blind candidate selection (4): tmle_candidate,
run_plasmode_feasibility, run_plasmode_dq_stress, select_tmle_candidate,
plus select_variance_method (5).

Nuisance and utilities (4): build_sl_library, resolve_truncation_rule,
sanitize_covariates, sim_func1.

Learners (6, exported by technical necessity: tmle::tmle resolves learner
names on the search path): SL.glmnet.bounded, SL.glm.pca, SL.glm.pca5,
SL.glm.pca10, SL.glmnet.ridge, SL.glmnet.enet.

Sensitivity and precision (3): compute_evalue, run_delta_sensitivity,
estimate_design_precision.

Plots and tables (7): love_plot, forest_plot, clean_weight_diagnostics,
make_table1, make_table2, attrition_table, plus the S3 plot and print
methods, which are not counted as exports.

Cumulative-risk grammar, the causalRisk alignment (13): specify_models,
identify_outcome, identify_treatment, identify_censoring,
identify_competing_risk, identify_subject, identify_interval,
estimate_ipwrisk, estimate_gcomprisk, estimate_aipwrisk, estimate_ipwhr,
estimate_surv_tmle, estimate_lmtp. (identify_missing becomes internal.) A
future release can spin this layer into a sister package, which would take
the core package to about 42 exports.

Net: 149 exports fall to 55; the user-facing clean-room surface is the
15-function core; everything cut keeps working internally where an in-repo
driver still needs it, and dies in attic/ where nothing does.

## 4. Sequencing and risk

Implemented in this order, with devtools::test() (source-loaded, safe beside
the running background chain) after each family: the merges (fit_ps,
assess_support, estimand_feasibility, simulate_support, estimate_effect,
lock constructors), then the internalizations, then the deletions, then the
caller updates (package tests, rescueCo scripts 10 to 14, the two vignettes,
run_simulation.R and legacy stage scripts via :::), then NAMESPACE
regeneration and the reinstall once the background chain finishes. The
manuscripts are rewritten against this API, not the old one.
