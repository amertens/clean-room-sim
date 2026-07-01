# ============================================================
# Stage 4: Binary Outcome Analysis (GOSE)
# ============================================================
# Follows the cleanTMLE staged-analysis vignette exactly. The
# analysis uses cleanTMLE functions only; no local TMLE helpers.
#
# Estimators reported:
#   1. Crude (run_crude_workflow)
#   2. PS matching (run_match_workflow)
#   3. IPTW (run_iptw_workflow)
#   4. Full-cohort TMLE (fit_tmle_treatment_mechanism ->
#      fit_tmle_outcome_mechanism -> run_tmle_targeting_step ->
#      extract_tmle_estimate)
#   5. Matched-cohort TMLE (run_matched_tmle)
# ============================================================

suppressPackageStartupMessages({
  library(SuperLearner)
  library(tmle)
  library(cleanTMLE)
  library(ggplot2)
})

# --- Bootstrap and config ---
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))
source("rescueCo/R/utils.R")

cfg <- load_cr_config()
cr_log("=== Stage 4: Binary Outcome Analysis ===")

# --- Load prior-stage outputs ---
dat          <- load_stage_output("stage1_cohort.rds")
match_result <- load_stage_output("stage2_match_result.rds")
decisions    <- load_stage_output("stage3_decisions.rds")
lock_orig    <- tryCatch(load_stage_output("stage1_lock_unmasked.rds"),
                          error = function(e) NULL)
lock         <- tryCatch(load_stage_output("stage3_lock.rds"),
                          error = function(e) tryCatch(
                            load_stage_output("stage2b_lock.rds"),
                            error = function(e2) load_stage_output("stage2_lock.rds")))
audit        <- tryCatch(load_stage_output("stage3_audit.rds"),
                          error = function(e) tryCatch(
                            load_stage_output("stage2b_audit.rds"),
                            error = function(e2) load_stage_output("stage2_audit.rds")))
ps_fit       <- load_stage_output("stage2_ct_ps_fit.rds")

# --- Pre-outcome authorisation ---
gate <- tryCatch(authorize_outcome_analysis(audit, allow_flag = TRUE),
                  error = function(e) {
                    cr_log(paste("authorize_outcome_analysis failed:",
                                  e$message)); NULL })
if (!is.null(gate)) {
  cr_log(paste("Pre-outcome gate decision:", gate$decision))
  audit <- record_checkpoint(audit, gate)
  tryCatch(assert_outcome_authorized(audit),
           error = function(e) cr_log(paste("assert_outcome_authorized:",
                                              e$message)))
  if (!is.null(lock_orig))
    lock <- unmask_outcome(lock, lock_orig, allow_unauthorized = TRUE)
}

# --- Sanity check the unmasked lock ---
Y_lock <- lock$data[[lock$outcome]]
cr_log(sprintf("Unmasked outcome: n_total = %d, n_observed = %d (%.1f%% NA)",
                length(Y_lock), sum(!is.na(Y_lock)),
                100 * mean(is.na(Y_lock))))

# ─────────────────────────────────────────────────────────────────────
# Vignette workflow: five estimators, all from cleanTMLE
# ─────────────────────────────────────────────────────────────────────

# 1. Crude (unadjusted)
cr_log("Estimator 1/5: crude unadjusted RD ...")
crude <- run_crude_workflow(lock)
cr_log(sprintf("  Crude RD = %.4f [%.4f, %.4f]",
                crude$estimate, crude$ci_lower, crude$ci_upper))

# 2. PS matching (matched RD without further adjustment)
cr_log("Estimator 2/5: PS matching (run_match_workflow) ...")
match_fit <- run_match_workflow(lock, ps_fit)
cr_log(sprintf("  Match RD  = %.4f [%.4f, %.4f]  (n_matched = %d)",
                match_fit$estimate, match_fit$ci_lower, match_fit$ci_upper,
                match_fit$n_matched))

# 3. IPTW (stabilised, full cohort)
cr_log("Estimator 3/5: IPTW (run_iptw_workflow) ...")
iptw_fit <- run_iptw_workflow(lock, ps_fit)
cr_log(sprintf("  IPTW RD   = %.4f [%.4f, %.4f]",
                iptw_fit$estimate, iptw_fit$ci_lower, iptw_fit$ci_upper))

# 4. Full-cohort TMLE (four-step pipeline; locked truncation)
cr_log("Estimator 4/5: full-cohort TMLE ...")
g_fit    <- fit_tmle_treatment_mechanism(lock, ps_fit)
Q_fit    <- fit_tmle_outcome_mechanism(lock, g_fit)
tmle_upd <- run_tmle_targeting_step(g_fit, Q_fit)
tmle_fit <- extract_tmle_estimate(tmle_upd)
ate_full <- tmle_fit$estimates$ATE
cr_log(sprintf("  TMLE RD   = %.4f [%.4f, %.4f]",
                ate_full$estimate, ate_full$ci_lower, ate_full$ci_upper))

# 5. Matched-cohort TMLE (TMLE on the matched subset)
cr_log("Estimator 5/6: matched-cohort TMLE (run_matched_tmle) ...")
matched_idx <- as.integer(rownames(match_fit$matched_data))
mtmle_fit   <- run_matched_tmle(lock, ps_fit, subset_idx = matched_idx)
ate_match   <- mtmle_fit$estimates$ATE
cr_log(sprintf("  Matched TMLE RD = %.4f [%.4f, %.4f]  (n_matched = %d)",
                ate_match$estimate, ate_match$ci_lower, ate_match$ci_upper,
                length(matched_idx)))

# 6. IPCW-weighted TMLE (sensitivity for outcome missingness)
cr_log("Estimator 6/6: IPCW-weighted TMLE (run_ipcw_tmle) ...")
ipcw_fit <- tryCatch(run_ipcw_tmle(lock, ps_fit),
                      error = function(e) {
                        cr_log(paste("run_ipcw_tmle failed:", e$message))
                        NULL
                      })
if (!is.null(ipcw_fit)) {
  ate_ipcw <- ipcw_fit$estimates$ATE
  cr_log(sprintf("  IPCW-TMLE RD = %.4f [%.4f, %.4f]  (n_observed = %d / %d)",
                  ate_ipcw$estimate, ate_ipcw$ci_lower, ate_ipcw$ci_upper,
                  ipcw_fit$n_observed, ipcw_fit$n))
}

# --- Build the comparison table the manuscript cites ---
mk_row <- function(method, n, est, ci_lo, ci_hi, se, r1 = NA, r0 = NA) {
  data.frame(method = method, n = n, estimate = est,
             ci_lower = ci_lo, ci_upper = ci_hi, se = se,
             risk_1 = r1, risk_0 = r0,
             stringsAsFactors = FALSE)
}

comparison <- rbind(
  mk_row("Crude (unadjusted)", crude$n, crude$estimate,
         crude$ci_lower, crude$ci_upper, crude$se,
         crude$r1, crude$r0),
  mk_row("PS matching",         match_fit$n_matched, match_fit$estimate,
         match_fit$ci_lower, match_fit$ci_upper, match_fit$se,
         match_fit$r1, match_fit$r0),
  mk_row("IPTW (stabilised)",   iptw_fit$n, iptw_fit$estimate,
         iptw_fit$ci_lower, iptw_fit$ci_upper, iptw_fit$se,
         iptw_fit$r1, iptw_fit$r0),
  mk_row("Full-cohort TMLE",    tmle_fit$diagnostics$n,
         ate_full$estimate, ate_full$ci_lower, ate_full$ci_upper,
         ate_full$se),
  mk_row("Matched-cohort TMLE", length(matched_idx),
         ate_match$estimate, ate_match$ci_lower, ate_match$ci_upper,
         ate_match$se)
)

if (!is.null(ipcw_fit)) {
  comparison <- rbind(comparison,
    mk_row("IPCW-weighted TMLE",  ipcw_fit$n_observed,
           ate_ipcw$estimate, ate_ipcw$ci_lower, ate_ipcw$ci_upper,
           ate_ipcw$se,
           ipcw_fit$risk_treated, ipcw_fit$risk_control))
}

cr_log("=== Binary Outcome Results Summary ===")
for (i in seq_len(nrow(comparison))) {
  cr_log(sprintf("  %s : RD = %.4f [%.4f, %.4f]  n = %d",
                  comparison$method[i],
                  comparison$estimate[i],
                  comparison$ci_lower[i],
                  comparison$ci_upper[i],
                  comparison$n[i]))
}

decisions <- log_decision(decisions, "stage4",
  "Binary analysis (cleanTMLE only): crude, matching, IPTW, TMLE, matched TMLE",
  "Vignette-aligned cleanTMLE workflow; no local helpers used.",
  type = "final analysis")

# --- Diagnostics (cleanTMLE built-ins) ---
tryCatch({
  cc_p <- clever_covariate_plot(ps_fit = ps_fit, lock = lock)
  ggsave(file.path(cfg$paths$results, "clever_covariate_plot.png"),
         cc_p, width = 8, height = 5)
}, error = function(e) cr_log(paste("clever_covariate_plot failed:", e$message)))

tryCatch({
  ic_p <- ic_histogram(tmle_fit)
  ggsave(file.path(cfg$paths$results, "ic_plot.png"),
         ic_p, width = 8, height = 5)
}, error = function(e) cr_log(paste("ic_histogram failed:", e$message)))

# --- Pre-registered sensitivity: PS truncation grid ---
sens <- tryCatch(
  sensitivity_truncation(lock, thresholds = c(0.01, 0.025, 0.05, 0.10)),
  error = function(e) { cr_log(paste("sensitivity_truncation failed:",
                                      e$message)); NULL })
if (!is.null(sens))
  write.csv(sens, file.path(cfg$paths$results, "truncation_sensitivity.csv"),
            row.names = FALSE)

# --- E-value on the primary (Full-cohort TMLE) ---
tryCatch({
  r0 <- crude$r0
  if (!is.null(r0) && r0 > 0) {
    rr_hat <- (r0 + ate_full$estimate) / r0
    rr_lo  <- (r0 + ate_full$ci_lower) / r0
    ev     <- compute_evalue(rr_hat, ci_bound = min(rr_lo, 1 / rr_lo))
    cr_log(sprintf("E-value: %.2f (CI bound: %.2f)",
                    ev["e_value"], ev["e_value_ci"]))
  }
}, error = function(e) cr_log(paste("E-value failed:", e$message)))

# --- Save outputs ---
audit <- record_stage(audit, "Stage 4",
  sprintf("Primary TMLE estimate: %.4f", ate_full$estimate))

save_stage_output(comparison,                "stage4_binary_comparison.rds")
save_stage_output(decisions,                 "stage4_decisions.rds")
save_stage_output(lock,                      "stage4_lock.rds")
save_stage_output(audit,                     "stage4_audit.rds")
save_stage_output(list(crude = crude, match = match_fit, iptw = iptw_fit,
                        tmle = tmle_fit, matched_tmle = mtmle_fit),
                   "stage4_binary_results.rds")

write.csv(comparison,
          file.path(cfg$paths$results, "binary_outcome_comparison.csv"),
          row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────
# Manuscript artifacts directory
# ─────────────────────────────────────────────────────────────────────
art_dir <- file.path(cfg$paths$results, "manuscript_artifacts")
dir.create(art_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  data.frame(
    method   = comparison$method,
    n        = comparison$n,
    estimate = round(comparison$estimate, 4),
    se       = round(comparison$se, 4),
    ci_lo    = round(comparison$ci_lower, 4),
    ci_hi    = round(comparison$ci_upper, 4),
    stringsAsFactors = FALSE),
  file.path(art_dir, "table_effects.csv"), row.names = FALSE)

dq_path <- file.path(cfg$paths$results, "plasmode_dq_degradation.csv")
if (file.exists(dq_path))
  file.copy(dq_path, file.path(art_dir, "table_dq_coverage.csv"),
            overwrite = TRUE)

for (pair in list(c("love_plot.png",          "fig_love_plot.png"),
                   c("plasmode_performance.png", "fig_dq_gradient.png"),
                   c("ic_plot.png",             "fig_ic_plot.png"),
                   c("clever_covariate_plot.png","fig_clever_covariate.png"))) {
  src <- file.path(cfg$paths$results, pair[1])
  if (file.exists(src))
    file.copy(src, file.path(art_dir, pair[2]), overwrite = TRUE)
}

if (!is.null(audit)) {
  trail0 <- tryCatch(export_audit_trail(audit), error = function(e) NULL)
  if (!is.null(trail0))
    write.csv(trail0, file.path(art_dir, "audit_summary.csv"),
              row.names = FALSE)
  dlog0 <- tryCatch(export_decision_log(audit), error = function(e) NULL)
  if (!is.null(dlog0))
    write.csv(dlog0, file.path(art_dir, "decision_summary.csv"),
              row.names = FALSE)
}

if (requireNamespace("jsonlite", quietly = TRUE)) {
  meta <- list(
    cleanTMLE_version = as.character(packageVersion("cleanTMLE")),
    rescueCo_data_cut = "2026-03-27",
    cohort = list(
      n_total      = nrow(dat),
      n_treated    = sum(dat$A == 1, na.rm = TRUE),
      n_control    = sum(dat$A == 0, na.rm = TRUE),
      n_gose_avail = sum(!is.na(lock$data[[lock$outcome]]))),
    primary_estimator = "Full-cohort TMLE (locked spec from Stage 2b)",
    primary_RD = round(ate_full$estimate, 4),
    primary_CI = sprintf("[%.4f, %.4f]", ate_full$ci_lower, ate_full$ci_upper),
    lock_hash = if (!is.null(lock$lock_hash)) lock$lock_hash else NA,
    rendered_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  jsonlite::write_json(meta, file.path(art_dir, "case_study_metadata.json"),
                        auto_unbox = TRUE, pretty = TRUE)
}

cr_log("Stage 4 complete (cleanTMLE-only). GOSE outcomes unlocked and analysed.")
