# ============================================================
# Stage 2: Design Adequacy Checks
# ============================================================
# Propensity score estimation, matching, overlap diagnostics.
# NO outcome information used at this stage.
# ============================================================

library(SuperLearner)
library(MatchIt)
library(ggplot2)

# --- Source helpers ---
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/ps_matching.R")
source("rescueCo/R/plotting.R")

# --- Load config and Stage 1 outputs ---
cfg <- load_cr_config()
cr_log("=== Stage 2: Design Adequacy Checks ===")

dat      <- load_stage_output("stage1_cohort.rds")
W_matrix <- load_stage_output("stage1_W_matrix.rds")
decisions <- load_stage_output("stage1_decisions.rds")
lock     <- load_stage_output("stage1_lock.rds")
audit    <- load_stage_output("stage1_audit.rds")

A <- dat$A
n <- length(A)
cr_log(paste("Working with", n, "observations,", ncol(W_matrix), "covariates"))

# --- Set up parallel cluster (Windows PSOCK) ---
cl <- NULL
if (isTRUE(cfg$superlearner$parallel)) {
  n_cores <- cfg$superlearner$n_cores %||% 2L
  cr_log(paste("Setting up parallel cluster:", n_cores, "cores"))
  cl <- setup_sl_cluster(n_cores)
}

# --- Estimate propensity scores using SuperLearner ---
cr_log("Fitting propensity score model with SuperLearner...")

sl_lib <- cfg$superlearner$candidate_learners
cv_folds <- cfg$superlearner$cv_folds

ps_result <- estimate_propensity_score(
  A = A, W = W_matrix, sl_lib = sl_lib,
  cv_folds = cv_folds, seed = cfg$seed, cl = cl
)

# Clean up cluster
stop_sl_cluster(cl)

cr_log(paste("PS model fit. CV risks:", paste(round(ps_result$cv_risk, 4), collapse = ", ")))
cr_log(paste("SL coefficients:", paste(round(ps_result$coef, 3), collapse = ", ")))

decisions <- log_decision(decisions, "stage2",
                          paste("PS model: SuperLearner with", paste(sl_lib, collapse = ", ")),
                          "Pre-specified ensemble learner set",
                          type = "pre-specified")

# --- Truncate PS ---
ps_raw <- ps_result$ps_scores
ps_trunc <- truncate_ps(ps_raw,
                         lower = cfg$propensity_score$truncation_lower,
                         upper = cfg$propensity_score$truncation_upper)

n_truncated <- sum(ps_raw != ps_trunc, na.rm = TRUE)
cr_log(paste("PS truncated to [", cfg$propensity_score$truncation_lower, ",",
             cfg$propensity_score$truncation_upper, "]:",
             n_truncated, "values truncated"))

decisions <- log_decision(decisions, "stage2",
                          paste("PS truncation:", cfg$propensity_score$truncation_lower,
                                "-", cfg$propensity_score$truncation_upper),
                          "Standard practical positivity enforcement",
                          type = "pre-specified")

# --- PS overlap diagnostics ---
cr_log("Computing PS overlap diagnostics...")
ps_valid <- !is.na(ps_trunc) & !is.na(A)

ps_treated <- ps_trunc[ps_valid & A == 1]
ps_control <- ps_trunc[ps_valid & A == 0]

overlap_diagnostics <- list(
  ps_treated_mean   = mean(ps_treated),
  ps_treated_median = median(ps_treated),
  ps_treated_range  = range(ps_treated),
  ps_control_mean   = mean(ps_control),
  ps_control_median = median(ps_control),
  ps_control_range  = range(ps_control),
  ks_test           = ks.test(ps_treated, ps_control),
  overlap_region    = c(max(min(ps_treated), min(ps_control)),
                        min(max(ps_treated), max(ps_control)))
)

cr_log(paste("PS treated: mean =", round(overlap_diagnostics$ps_treated_mean, 3),
             ", range = [", round(overlap_diagnostics$ps_treated_range[1], 3), ",",
             round(overlap_diagnostics$ps_treated_range[2], 3), "]"))
cr_log(paste("PS control: mean =", round(overlap_diagnostics$ps_control_mean, 3),
             ", range = [", round(overlap_diagnostics$ps_control_range[1], 3), ",",
             round(overlap_diagnostics$ps_control_range[2], 3), "]"))

# --- PS Overlap Plot ---
results_dir <- cfg$paths$results
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

plot_ps_overlap(ps_trunc, A, file.path(results_dir, "ps_overlap_density.png"))

# --- IPTW weights and effective sample size ---
iptw_weights <- ifelse(A == 1, 1 / ps_trunc, 1 / (1 - ps_trunc))
iptw_weights[is.na(iptw_weights)] <- 0

ess_treated <- effective_sample_size(iptw_weights[A == 1 & ps_valid])
ess_control <- effective_sample_size(iptw_weights[A == 0 & ps_valid])

cr_log(paste("Effective sample size: treated =", round(ess_treated, 1),
             ", control =", round(ess_control, 1)))

plot_weight_distribution(iptw_weights[ps_valid], A[ps_valid],
                          file.path(results_dir, "weight_distribution.png"))

# --- IPTW weight diagnostics (Stage 2a) ---
cr_log("Computing IPTW weight diagnostics...")

wt_summary <- make_weight_summary_table(iptw_weights[ps_valid], A[ps_valid])
write.csv(wt_summary, file.path(results_dir, "weight_summary_table.csv"), row.names = FALSE)
cr_log("Weight summary by group:")
for (r in seq_len(nrow(wt_summary))) {
  cr_log(paste0("  ", wt_summary$group[r], ": median = ", round(wt_summary$median[r], 2),
                ", max = ", round(wt_summary$max[r], 2),
                ", SD = ", round(wt_summary$SD[r], 2)))
}

# Extreme weight detection
iptw_diag_cfg <- cfg$iptw_diagnostics
extreme_wt <- detect_extreme_weights(
  iptw_weights[ps_valid], A[ps_valid],
  pctl = iptw_diag_cfg$extreme_weight_percentile %||% 0.99,
  ratio_threshold = iptw_diag_cfg$max_min_ratio_threshold %||% 50
)

if (extreme_wt$flag) {
  cr_log(paste("WARNING: Extreme weights detected!",
               extreme_wt$n_extreme, "observations (",
               round(extreme_wt$pct_extreme, 1), "%) above threshold;",
               "max/min ratio =", round(extreme_wt$max_min_ratio, 1)))
  decisions <- log_decision(decisions, "stage2",
    paste("FLAG: Extreme IPTW weights detected (",
          extreme_wt$n_extreme, " obs, ratio =",
          round(extreme_wt$max_min_ratio, 1), ")"),
    "Consider trimming weights or adding covariates to improve overlap",
    type = "design-stage")
} else {
  cr_log("IPTW weight diagnostics: no extreme weight concerns")
}

# ESS check against minimum fraction
ess_min_frac <- iptw_diag_cfg$ess_min_fraction %||% 0.10
n_treated <- sum(A[ps_valid] == 1)
n_control <- sum(A[ps_valid] == 0)

if (ess_treated < ess_min_frac * n_treated) {
  cr_log(paste("WARNING: Low ESS for treated group:", round(ess_treated, 1),
               "< threshold", round(ess_min_frac * n_treated, 1)))
  decisions <- log_decision(decisions, "stage2",
    paste("FLAG: Low ESS treated =", round(ess_treated, 1)),
    "ESS below minimum fraction threshold",
    type = "design-stage")
}
if (ess_control < ess_min_frac * n_control) {
  cr_log(paste("WARNING: Low ESS for control group:", round(ess_control, 1),
               "< threshold", round(ess_min_frac * n_control, 1)))
  decisions <- log_decision(decisions, "stage2",
    paste("FLAG: Low ESS control =", round(ess_control, 1)),
    "ESS below minimum fraction threshold",
    type = "design-stage")
}

# --- Matching ---
cr_log("Performing 1:1 nearest-neighbor matching...")

match_result <- perform_ps_matching(
  A = A, ps = ps_trunc,
  caliper_sd_mult = cfg$propensity_score$caliper_sd_multiplier
)

cr_log(paste("Matching complete:",
             match_result$n_treated_matched, "treated matched of",
             match_result$n_treated_total, "total (",
             round(match_result$frac_treated_matched * 100, 1), "%)"))
cr_log(paste("Caliper:", round(match_result$caliper, 4)))

decisions <- log_decision(decisions, "stage2",
                          paste("1:1 NN matching, caliper =",
                                round(match_result$caliper, 4),
                                "(", cfg$propensity_score$caliper_sd_multiplier,
                                "x SD(logit PS))"),
                          "Standard caliper matching specification",
                          type = "pre-specified")

# --- Covariate balance ---
cr_log("Computing covariate balance (SMDs)...")

smd_pre  <- compute_all_smds(W_matrix, A)
smd_post <- compute_all_smds(W_matrix, A, match_result$matched_idx)

# Love plot data
love_data <- merge(smd_pre, smd_post, by = "covariate", suffixes = c("_pre", "_post"))

# Summary
n_imbalanced_pre  <- sum(abs(love_data$smd_pre) > 0.1)
n_imbalanced_post <- sum(abs(love_data$smd_post) > 0.1)
max_smd_pre  <- max(abs(love_data$smd_pre))
max_smd_post <- max(abs(love_data$smd_post))

cr_log(paste("Covariates with |SMD| > 0.1: pre =", n_imbalanced_pre,
             ", post =", n_imbalanced_post))
cr_log(paste("Max |SMD|: pre =", round(max_smd_pre, 3),
             ", post =", round(max_smd_post, 3)))

# Love plot (pre vs. post matching)
plot_love(love_data, file.path(results_dir, "love_plot.png"))

# IPTW-weighted SMDs for three-way comparison
smd_weighted <- compute_weighted_smds(W_matrix[ps_valid, , drop = FALSE],
                                       A[ps_valid],
                                       iptw_weights[ps_valid])

# Three-way love plot: unweighted vs. matched vs. IPTW
plot_weighted_love(smd_pre, smd_post, smd_weighted,
                   file.path(results_dir, "weighted_love_plot.png"))

n_imbalanced_weighted <- sum(abs(smd_weighted$smd) > 0.1)
max_smd_weighted <- max(abs(smd_weighted$smd))
cr_log(paste("IPTW-weighted balance: covariates with |SMD| > 0.1:", n_imbalanced_weighted,
             ", max |SMD| =", round(max_smd_weighted, 3)))

# Balance table (all three methods)
balance_table <- merge(love_data, smd_weighted, by = "covariate", suffixes = c("", "_weighted"))
names(balance_table)[names(balance_table) == "smd"] <- "smd_weighted"
balance_table <- balance_table[order(-abs(balance_table$smd_pre)), ]
write.csv(balance_table, file.path(results_dir, "balance_table.csv"), row.names = FALSE)

# --- cleanTMLE PS diagnostics (modernised: SuperLearner PS) ---
cr_log("Fitting cleanTMLE SuperLearner PS (Stage 2 modular API)...")

# unmask the lock (Stage 2 still doesn't access outcome; unmask_outcome
# only restores the column for downstream stages that legitimately need it)
lock_unmasked <- tryCatch(
  cleanTMLE::unmask_outcome(lock, load_stage_output("stage1_lock_unmasked.rds"),
                            allow_unauthorized = TRUE),
  error = function(e) { cr_log(paste("unmask_outcome failed:", e$message)); NULL }
)

ct_ps_fit <- NULL
ct_diag   <- NULL
if (!is.null(lock_unmasked)) {
  # 10-fold CV (Part 3 #3 of case-study prompt). fit_ps_parallel doesn't
  # take cv_folds; fit_ps_superlearner does. We use the sequential
  # SuperLearner version because it supports the cv_folds argument
  # explicitly. Wrap in tryCatch to fall back gracefully.
  cv_folds <- cfg$superlearner$cv_folds %||% 10L
  ct_ps_fit <- tryCatch(
    cleanTMLE::fit_ps_superlearner(lock_unmasked,
      truncate = cfg$propensity_score$truncation_lower,
      cv_folds = cv_folds),
    error = function(e) {
      cr_log(paste("fit_ps_superlearner failed:", e$message,
                    "— falling back to fit_ps_glm"))
      tryCatch(
        cleanTMLE::fit_ps_glm(lock_unmasked,
          truncate = cfg$propensity_score$truncation_lower),
        error = function(e2) NULL)
    })
  if (!is.null(ct_ps_fit))
    cr_log(paste("Stage 2 PS fit OK (cv_folds =", cv_folds, ")"))
}

if (!is.null(ct_ps_fit)) {
  ct_diag <- tryCatch(
    cleanTMLE::compute_ps_diagnostics(ct_ps_fit),
    error = function(e) { cr_log(paste("compute_ps_diagnostics failed:", e$message)); NULL }
  )

  ct_wt_summary <- tryCatch(cleanTMLE::make_wt_summary_table(ct_ps_fit),
    error = function(e) NULL)
  ct_extreme_wt <- tryCatch(cleanTMLE::extreme_weights(ct_ps_fit, k = 10),
    error = function(e) NULL)

  # Inspect IPW weights (built-in diagnostic plot/table)
  ct_ipw_inspect <- tryCatch(cleanTMLE::inspect_ipw_weights(ct_ps_fit),
    error = function(e) NULL)
  if (!is.null(ct_ipw_inspect)) saveRDS(ct_ipw_inspect,
    file.path(results_dir, "cleanTMLE_ipw_weight_inspection.rds"))

  # Three-way love plot (unmatched / IPTW / matched)
  matched_idx <- tryCatch(match_result$matched_idx, error = function(e) NULL)
  if (!is.null(matched_idx) && length(matched_idx) > 0) {
    m_smds <- tryCatch(
      cleanTMLE::compute_matched_smds(lock_data <- lock_unmasked$data,
        treatment = lock_unmasked$treatment,
        covariates = lock_unmasked$covariates,
        subset_idx = matched_idx),
      error = function(e) { cr_log(paste("compute_matched_smds failed:", e$message)); NULL })
    if (!is.null(m_smds)) {
      ct_love3 <- tryCatch(
        cleanTMLE::love_plot_threeway(ct_diag, matched_smds = m_smds),
        error = function(e) NULL)
      if (!is.null(ct_love3))
        ggsave(file.path(results_dir, "cleanTMLE_love_plot_threeway.png"),
               ct_love3, width = 10, height = 9)
    }
  }

  ct_love <- tryCatch(cleanTMLE::love_plot(ct_diag),
    error = function(e) NULL)
  if (!is.null(ct_love))
    ggsave(file.path(results_dir, "cleanTMLE_love_plot.png"),
           ct_love, width = 10, height = 8)

  # Baseline characteristics table 1 (cleanTMLE 0.1.1 accepts a lock directly,
  # per Issue #4 in the improvement prompt)
  ct_table1 <- tryCatch(cleanTMLE::make_table1(lock_unmasked),
    error = function(e) {
      cr_log(paste("make_table1(lock) failed:", e$message,
                    "— falling back to manual construction"))
      d <- lock_unmasked$data
      A_var <- lock_unmasked$treatment
      tab1 <- do.call(rbind, lapply(lock_unmasked$covariates, function(v) {
        x <- as.numeric(d[[v]])
        data.frame(variable = v,
                   mean_treated = mean(x[d[[A_var]] == 1], na.rm = TRUE),
                   mean_control = mean(x[d[[A_var]] == 0], na.rm = TRUE),
                   smd = NA_real_)
      }))
      tab1
    })
  if (!is.null(ct_table1))
    write.csv(ct_table1, file.path(results_dir, "table1_baseline.csv"),
              row.names = FALSE)
}

# Check Point 2: Balance
if (!is.null(ct_diag)) {
  cp2 <- tryCatch(
    cleanTMLE::checkpoint_balance(ct_diag,
      max_smd     = 0.10,
      min_ess_pct = 50,
      lock_hash   = lock$lock_hash),
    error = function(e) { cr_log(paste("checkpoint_balance failed:", e$message)); NULL }
  )
  if (!is.null(cp2)) {
    audit <- cleanTMLE::record_checkpoint(audit, cp2)
    cr_log(paste("Check Point 2 (balance):", cp2$decision))
    decisions <- log_decision(decisions, "stage2",
      paste("cleanTMLE Check Point 2:", cp2$decision),
      cp2$rationale %||% "Balance assessment via cleanTMLE",
      type = "design-stage")
  }
}

audit <- tryCatch(
  cleanTMLE::record_stage(audit, "Stage 2", "PS estimation and balance checks complete"),
  error = function(e) { cr_log(paste("record_stage failed:", e$message)); audit }
)

# --- Stage 2 summary ---
design_summary <- list(
  ps_model_cv_risk  = ps_result$cv_risk,
  ps_model_coef     = ps_result$coef,
  n_truncated       = n_truncated,
  overlap           = overlap_diagnostics,
  matching          = list(
    n_treated_matched  = match_result$n_treated_matched,
    n_control_matched  = match_result$n_control_matched,
    frac_matched       = match_result$frac_treated_matched,
    caliper            = match_result$caliper
  ),
  balance = list(
    n_imbalanced_pre      = n_imbalanced_pre,
    n_imbalanced_post     = n_imbalanced_post,
    n_imbalanced_weighted = n_imbalanced_weighted,
    max_smd_pre           = max_smd_pre,
    max_smd_post          = max_smd_post,
    max_smd_weighted      = max_smd_weighted
  ),
  ess_treated = ess_treated,
  ess_control = ess_control,
  iptw = list(
    weight_summary    = wt_summary,
    extreme_weights   = extreme_wt,
    smd_weighted      = smd_weighted
  )
)

# --- Save stage outputs ---
save_stage_output(ps_result, "stage2_ps_result.rds")
save_stage_output(ps_trunc, "stage2_ps_truncated.rds")
save_stage_output(match_result, "stage2_match_result.rds")
save_stage_output(love_data, "stage2_love_data.rds")
save_stage_output(design_summary, "stage2_design_summary.rds")
save_stage_output(decisions, "stage2_decisions.rds")
save_stage_output(lock, "stage2_lock.rds")
save_stage_output(audit, "stage2_audit.rds")
if (!is.null(ct_ps_fit)) save_stage_output(ct_ps_fit, "stage2_ct_ps_fit.rds")

cr_log("Stage 2 complete. NO outcomes were used.")
