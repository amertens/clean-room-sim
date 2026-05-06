# ============================================================
# Stage 2b: Negative Control Analysis
# ============================================================
# Tests TMLE on pre-treatment covariates that should NOT be
# affected by ambulance transport. Also runs plasmode simulation
# to compare TMLE vs. matching vs. g-computation under known ATE.
#
# NO real outcomes are used at this stage.
# ============================================================

library(SuperLearner)
library(tmle)

# --- Source helpers ---
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/negative_controls.R")
source("rescueCo/R/plotting.R")

# --- Load config and prior stage outputs ---
cfg <- load_cr_config()
cr_log("=== Stage 2b: Negative Control Analysis ===")

dat          <- load_stage_output("stage1_cohort.rds")
W_matrix     <- load_stage_output("stage1_W_matrix.rds")
ps_result    <- load_stage_output("stage2_ps_result.rds")
ps_trunc     <- load_stage_output("stage2_ps_truncated.rds")
match_result <- load_stage_output("stage2_match_result.rds")
decisions    <- load_stage_output("stage2_decisions.rds")
lock         <- load_stage_output("stage2_lock.rds")
audit        <- load_stage_output("stage2_audit.rds")
ct_ps_fit    <- load_stage_output("stage2_ct_ps_fit.rds")

A <- dat$A
matched_idx <- match_result$matched_idx
sl_lib <- cfg$superlearner$candidate_learners
results_dir <- cfg$paths$results

# --- Negative control outcome variables ---
nc_vars <- cfg$negative_controls$outcomes
cr_log(paste("Negative control outcomes:", paste(nc_vars, collapse = ", ")))

# Check which NC variables exist in W_matrix
nc_available <- nc_vars[nc_vars %in% colnames(W_matrix)]
nc_missing   <- nc_vars[!nc_vars %in% colnames(W_matrix)]

if (length(nc_missing) > 0) {
  cr_log(paste("WARNING: NC variables not found in W_matrix:", paste(nc_missing, collapse = ", ")))
}

if (length(nc_available) == 0) {
  cr_log("ERROR: No negative control variables found in covariate matrix")
  decisions <- log_decision(decisions, "stage2b",
    "SKIP: No NC variables available",
    "None of the configured NC variables exist in W",
    type = "design-stage")
  save_stage_output(decisions, "stage2b_decisions.rds")
  stop("No negative control variables available")
}

# --- cleanTMLE negative control analysis ---
# NOTE: cleanTMLE's run_negative_control() uses IPTW, not TMLE.
# This is a methodological difference but still serves the same diagnostic purpose.
cr_log("Running cleanTMLE negative control analysis...")

# Unmask lock for NC analysis (NCs don't use primary outcome)
lock_for_nc <- tryCatch(
  cleanTMLE::unmask_outcome(lock, load_stage_output("stage1_lock_unmasked.rds")),
  error = function(e) { cr_log(paste("unmask_outcome failed:", e$message)); NULL }
)

stage3_nc <- NULL
if (!is.null(lock_for_nc) && !is.null(ct_ps_fit)) {
  stage3_nc <- tryCatch(
    cleanTMLE::run_residual_confounding_stage(lock_for_nc, ct_ps_fit),
    error = function(e) {
      cr_log(paste("cleanTMLE run_residual_confounding_stage failed:", e$message))
      cr_log("Falling back to custom negative control analysis...")
      NULL
    }
  )
}

if (!is.null(stage3_nc)) {
  nc_results <- stage3_nc$summary_table
  cr_log("cleanTMLE NC results (IPTW):")
  print(nc_results)

  cp3 <- stage3_nc$checkpoint
  audit <- cleanTMLE::record_checkpoint(audit, cp3)
  cr_log(paste("Check Point 3 (residual bias):", cp3$decision))

  decisions <- log_decision(decisions, "stage2b",
    paste("cleanTMLE Check Point 3:", cp3$decision),
    cp3$rationale %||% "Negative control analysis (IPTW)",
    type = "design-stage")

  # Also run doubly-robust TMLE-based NC for each registered NC
  cr_log("Running doubly-robust TMLE negative-control analysis...")
  nc_tmle_rows <- list()
  for (nc_var in nc_available) {
    nc_tmle <- tryCatch(
      cleanTMLE::run_negative_control_tmle(lock_for_nc, nc_var, ct_ps_fit),
      error = function(e) { cr_log(paste("run_negative_control_tmle(",
                                           nc_var, ") failed:", e$message)); NULL })
    if (!is.null(nc_tmle)) {
      nc_tmle_rows[[length(nc_tmle_rows) + 1]] <- data.frame(
        nc_name = nc_var,
        estimator = nc_tmle$method %||% "TMLE",
        estimate = nc_tmle$estimate,
        se = nc_tmle$std_error %||% NA_real_,
        p_value = nc_tmle$p_value %||% NA_real_,
        flagged = isTRUE(nc_tmle$flagged),
        stringsAsFactors = FALSE)
    }
  }
  if (length(nc_tmle_rows) > 0) {
    nc_tmle_df <- do.call(rbind, nc_tmle_rows)
    write.csv(nc_tmle_df, file.path(results_dir,
              "negative_control_tmle_results.csv"), row.names = FALSE)
    cr_log("Doubly-robust NC TMLE results saved.")
    print(nc_tmle_df)
  }

  # Structured residual-bias checkpoint
  cp_resid <- tryCatch(
    cleanTMLE::checkpoint_residual_bias(stage3_nc,
      max_abs_bias = 0.02, min_p = 0.05,
      lock_hash = lock$lock_hash),
    error = function(e) NULL)
  if (!is.null(cp_resid)) {
    audit <- cleanTMLE::record_checkpoint(audit, cp_resid)
    cr_log(paste("Residual-bias checkpoint:", cp_resid$decision))
  }
} else {
  # Fallback to custom NC analysis (keep existing code as fallback)
  cr_log("Using custom negative control analysis as fallback...")
  cr_log(paste("Running TMLE on", length(nc_available), "negative control outcomes..."))

  nc_results <- run_all_negative_controls(
    W = as.matrix(W_matrix),
    A = A,
    nc_vars = nc_available,
    sl_lib = sl_lib
  )

  cr_log("Negative control results:")
  for (i in seq_len(nrow(nc_results))) {
    status <- ifelse(isTRUE(nc_results$pass[i]), "PASS", "FAIL")
    cr_log(paste0("  ", nc_results$nc_name[i], ": ATE = ",
                  round(nc_results$estimate[i], 4),
                  " [", round(nc_results$ci_lower[i], 4), ", ",
                  round(nc_results$ci_upper[i], 4), "] -> ", status))
  }

  decisions <- log_decision(decisions, "stage2b",
    paste("NC TMLE (custom fallback):", sum(nc_results$pass, na.rm = TRUE), "/",
          sum(!is.na(nc_results$pass)), "pass"),
    "Testing pre-treatment covariates for zero treatment effect",
    type = "design-stage")
}

# --- Plasmode simulation ---
cr_log("Running plasmode simulation...")

plasmode_n_sims <- cfg$negative_controls$plasmode_n_sims %||% 50
plasmode_ate    <- cfg$negative_controls$plasmode_true_ate %||% 0

# --- cleanTMLE plasmode feasibility ---
cr_log("Running cleanTMLE plasmode feasibility...")

candidates <- list(
  cleanTMLE::tmle_candidate("glm_t01", "GLM trunc=0.01",
    g_library = c("SL.glm"), truncation = 0.01),
  cleanTMLE::tmle_candidate("glm_t05", "GLM trunc=0.05",
    g_library = c("SL.glm"), truncation = 0.05)
)

ct_plasmode <- NULL
if (!is.null(lock_for_nc)) {
  ct_plasmode <- tryCatch(
    cleanTMLE::run_plasmode_feasibility(
      lock_for_nc,
      tmle_candidates = candidates,
      effect_sizes    = c(0.05),
      reps            = plasmode_n_sims,
      verbose         = FALSE
    ),
    error = function(e) {
      cr_log(paste("cleanTMLE plasmode failed:", e$message))
      NULL
    }
  )
}

if (!is.null(ct_plasmode)) {
  selected <- cleanTMLE::select_tmle_candidate(ct_plasmode, rule = "min_rmse")
  lock <- cleanTMLE::lock_primary_tmle_spec(lock, selected)
  cr_log(paste("Selected TMLE candidate:", selected$label))

  audit <- cleanTMLE::record_stage(audit, "Stage 2b",
    paste("Plasmode complete; selected:", selected$candidate_id))

  decisions <- log_decision(decisions, "stage2b",
    paste("TMLE spec selected:", selected$candidate_id),
    "Minimum RMSE rule from cleanTMLE plasmode",
    type = "design-stage")

  # Also report cleanTMLE plasmode metrics
  if (!is.null(ct_plasmode$metrics)) {
    cr_log("cleanTMLE plasmode metrics:")
    print(ct_plasmode$metrics)
  }
} else {
  # Fallback to custom plasmode
  cr_log("Using custom plasmode simulation as fallback...")
  plasmode_results <- run_plasmode_simulation(
    A = A,
    W = as.matrix(W_matrix),
    ps = ps_trunc,
    matched_idx = matched_idx,
    sl_lib = "SL.glm",    # Use GLM for speed in plasmode
    n_sims = plasmode_n_sims,
    true_ate = plasmode_ate,
    seed = cfg$seed
  )

  cr_log("Plasmode simulation metrics:")
  for (i in seq_len(nrow(plasmode_results$metrics))) {
    m <- plasmode_results$metrics[i, ]
    cr_log(paste0("  ", m$estimator, ": bias = ", round(m$bias, 4),
                  ", RMSE = ", round(m$rmse, 4),
                  ", coverage = ", round(m$coverage, 3)))
  }

  decisions <- log_decision(decisions, "stage2b",
    paste("Plasmode simulation (custom fallback):", plasmode_n_sims,
          "reps, true ATE =", plasmode_ate),
    "Comparing TMLE, matching, and g-computation under known effect",
    type = "design-stage")
}

# --- cleanTMLE gate check ---
if (!is.null(ct_plasmode)) {
  ct_gate <- tryCatch(
    cleanTMLE::gate_check(ct_plasmode$metrics, "plasmode",
      targets = list(max_abs_bias = 0.01, min_coverage = 0.90)),
    error = function(e) { cr_log(paste("gate_check failed:", e$message)); NULL }
  )
  if (!is.null(ct_gate)) {
    cr_log(paste("cleanTMLE gate check:", ct_gate$decision))
  }
}

# --- GO / FLAG / STOP decision ---
cr_log("Evaluating GO/FLAG/STOP gating criteria...")

gfs_criteria <- cfg$go_flag_stop
# Use custom plasmode_results if cleanTMLE plasmode was not available
if (!exists("plasmode_results")) {
  plasmode_results <- if (!is.null(ct_plasmode)) ct_plasmode else NULL
}
# Ensure nc_results has a 'pass' column (cleanTMLE uses 'flagged' = inverse)
if (!"pass" %in% names(nc_results) && "flagged" %in% names(nc_results)) {
  nc_results$pass <- !nc_results$flagged
}
gfs_result <- evaluate_go_flag_stop(nc_results, plasmode_results, gfs_criteria)

cr_log(paste("GATING DECISION:", gfs_result$decision))
cr_log(paste("  Rationale:", gfs_result$rationale))

decisions <- log_decision(decisions, "stage2b",
  paste("Stage 2b gating:", gfs_result$decision),
  gfs_result$rationale,
  type = "design-stage")

# --- Save outputs ---
save_stage_output(nc_results, "stage2b_nc_results.rds")
save_stage_output(plasmode_results, "stage2b_plasmode_results.rds")
save_stage_output(gfs_result, "stage2b_gating_decision.rds")
save_stage_output(decisions, "stage2b_decisions.rds")

# CSV exports
write.csv(nc_results, file.path(results_dir, "negative_control_results.csv"),
          row.names = FALSE)
write.csv(plasmode_results$metrics, file.path(results_dir, "plasmode_metrics.csv"),
          row.names = FALSE)

# Plasmode performance plot (reuse simulation plot function)
tryCatch(
  plot_simulation_performance(plasmode_results$metrics,
                              file.path(results_dir, "plasmode_performance.png")),
  error = function(e) cr_log(paste("Plasmode plot skipped:", e$message))
)

save_stage_output(lock, "stage2b_lock.rds")
save_stage_output(audit, "stage2b_audit.rds")

cr_log("Stage 2b complete. NO real outcomes were used.")
