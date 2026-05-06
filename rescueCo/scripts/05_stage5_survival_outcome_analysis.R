# ============================================================
# Stage 5: Survival Outcome Analysis
# ============================================================
# Unlock time-to-event analyses.
# Implements:
#   E. Cox proportional hazards regression (full cohort)
#   F. Cox regression in matched cohort
#   G. Survival-TMLE (discrete-time, full cohort)
#   H. Matched survival-TMLE
# ============================================================

library(survival)
library(SuperLearner)

# --- Source helpers ---
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/tmle_survival.R")
source("rescueCo/R/plotting.R")

# --- Load config and prior outputs ---
cfg <- load_cr_config()
cr_log("=== Stage 5: Survival Outcome Analysis ===")

dat          <- load_stage_output("stage1_cohort.rds")
W_matrix     <- load_stage_output("stage1_W_matrix.rds")
match_result <- load_stage_output("stage2_match_result.rds")
sel_libs     <- tryCatch(
  load_stage_output("stage3_selected_libraries.rds"),
  error = function(e) {
    list(Q_model      = cfg$superlearner$candidate_learners,
         g_model      = cfg$superlearner$candidate_learners,
         censor_model = cfg$superlearner$survival_learners %||%
                          c("SL.glm", "SL.mean"))
  })
decisions    <- load_stage_output("stage4_decisions.rds")

# --- Load cleanTMLE lock and audit ---
lock  <- tryCatch(load_stage_output("stage4_lock.rds"),
                  error = function(e) NULL)
audit <- tryCatch(load_stage_output("stage4_audit.rds"),
                  error = function(e) NULL)

A <- dat$A
time_var  <- dat$surv_time_days
event_var <- dat$surv_event
matched_idx <- match_result$matched_idx
target_times <- cfg$survival_outcome$time_horizons

cr_log(paste("Survival data:",
             "n_events =", sum(event_var == 1, na.rm = TRUE),
             ", n_censored =", sum(event_var == 0, na.rm = TRUE),
             ", n_missing =", sum(is.na(event_var))))
cr_log(paste("Target time horizons:", paste(target_times, collapse = ", "), "days"))

# --- Check survival data availability ---
if (all(is.na(time_var)) || all(is.na(event_var))) {
  cr_log("ERROR: Survival variables are all NA. Cannot proceed with Stage 5.")
  cr_log("Please define surv_time_days and surv_event in Stage 1.")
  save_stage_output(decisions, "stage5_decisions.rds")
  stop("Survival outcome not available")
}

decisions <- log_decision(decisions, "stage5",
                          "Survival analysis: time_to_event and death indicator",
                          paste("Time horizons:", paste(target_times, collapse = ", "), "days"),
                          type = "pre-specified")

# --- Progress tracker for 6 estimators ---
pb5 <- cr_progress(6, "Stage 5")

# --- E0. survtmle package (PRIMARY survival estimator) ---
cr_log("Running survtmle (primary survival estimator)...")

result_survtmle <- tryCatch({
  require(survtmle)

  # Prepare data: survtmle expects ftime (positive integer), ftype (0/1)
  cc <- complete.cases(time_var, event_var, A)
  ftime_cc  <- as.integer(pmax(time_var[cc], 1))
  ftype_cc  <- as.integer(event_var[cc])
  A_cc      <- A[cc]
  W_cc      <- as.data.frame(W_matrix[cc, , drop = FALSE])

  # Remove zero-variance covariates
  col_var <- apply(W_cc, 2, var, na.rm = TRUE)
  W_cc <- W_cc[, !is.na(col_var) & col_var > 1e-8, drop = FALSE]

  # Bin to target times (survtmle needs bounded ftime)
  max_t <- max(target_times)
  ftime_cc <- pmin(ftime_cc, max_t)

  # Run survtmle at each target time
  survtmle_fit <- survtmle::survtmle(
    ftime    = ftime_cc,
    ftype    = ftype_cc,
    trt      = A_cc,
    adjustVars = W_cc,
    t0       = max(target_times),
    SL.ftime = c("SL.glm", "SL.mean"),
    SL.ctime = c("SL.glm", "SL.mean"),
    SL.trt   = c("SL.glm", "SL.mean"),
    method   = "hazard",
    returnModels = TRUE
  )

  # Extract cumulative incidence at each target time using timepoints
  survtmle_tp <- survtmle::timepoints(survtmle_fit, times = target_times)

  # Build results table
  survtmle_rows <- list()
  for (i in seq_along(target_times)) {
    tt <- target_times[i]
    tp_i <- survtmle_tp[[i]]

    # tp_i$est has rows = treatment groups, cols = event types
    risk_0 <- tp_i$est[1, 1]  # control, cause 1
    risk_1 <- tp_i$est[2, 1]  # treated, cause 1
    rd <- risk_1 - risk_0
    rr <- risk_1 / max(risk_0, 1e-10)

    # SE from tp_i$var
    se_0 <- sqrt(tp_i$var[1, 1])
    se_1 <- sqrt(tp_i$var[2, 1])
    rd_se <- sqrt(se_0^2 + se_1^2)

    survtmle_rows[[i]] <- data.frame(
      method = "survtmle (primary)",
      time = tt,
      risk_1 = risk_1,
      risk_0 = risk_0,
      risk_difference = rd,
      risk_ratio = rr,
      rd_se = rd_se,
      rd_ci_lower = rd - 1.96 * rd_se,
      rd_ci_upper = rd + 1.96 * rd_se,
      n = sum(cc),
      stringsAsFactors = FALSE
    )
  }
  list(
    method = "survtmle",
    estimates = do.call(rbind, survtmle_rows),
    fit = survtmle_fit
  )
}, error = function(e) {
  cr_log(paste("survtmle FAILED:", e$message))
  NULL
})

if (!is.null(result_survtmle)) {
  cr_log("survtmle results (PRIMARY):")
  for (i in seq_len(nrow(result_survtmle$estimates))) {
    est <- result_survtmle$estimates[i, ]
    cr_log(paste("  t =", est$time, "days:",
                 "RD =", round(est$risk_difference, 4),
                 "[", round(est$rd_ci_lower, 4), ",",
                 round(est$rd_ci_upper, 4), "]",
                 "RR =", round(est$risk_ratio, 3)))
  }
}

pb5$tick()

# --- E. Cox PH regression (sensitivity) ---
cr_log("Running Cox PH regression (sensitivity)...")

result_cox_full <- tryCatch(
  run_cox_regression(time_var, event_var, A, W_matrix),
  error = function(e) {
    cr_log(paste("Cox full failed:", e$message))
    NULL
  }
)

if (!is.null(result_cox_full)) {
  cr_log(paste("Cox full HR:", round(result_cox_full$hr, 3),
               "[", round(result_cox_full$hr_ci[1], 3), ",",
               round(result_cox_full$hr_ci[2], 3), "]",
               "p =", round(result_cox_full$pvalue, 4)))

  # PH assumption
  if (!is.null(result_cox_full$ph_test)) {
    ph_global <- result_cox_full$ph_test$table
    if ("GLOBAL" %in% rownames(ph_global)) {
      ph_p <- ph_global["GLOBAL", "p"]
      cr_log(paste("PH global test p-value:", round(ph_p, 4)))
      if (ph_p < 0.05) {
        cr_log("WARNING: PH assumption may be violated (p < 0.05)")
        decisions <- log_decision(decisions, "stage5",
                                  "PH assumption violation detected in full Cox model",
                                  paste("Global PH test p =", round(ph_p, 4)),
                                  type = "design-stage")
      }
    }
  }
}

pb5$tick()

# --- F. Matched Cox regression ---
cr_log("Running Cox regression in matched cohort...")

result_cox_matched <- tryCatch(
  run_cox_regression(time_var, event_var, A, W_matrix, matched_idx),
  error = function(e) {
    cr_log(paste("Cox matched failed:", e$message))
    NULL
  }
)

if (!is.null(result_cox_matched)) {
  cr_log(paste("Cox matched HR:", round(result_cox_matched$hr, 3),
               "[", round(result_cox_matched$hr_ci[1], 3), ",",
               round(result_cox_matched$hr_ci[2], 3), "]",
               "p =", round(result_cox_matched$pvalue, 4)))
}

pb5$tick()

# --- G. Discrete-time Survival-TMLE (full cohort) — SENSITIVITY ---
# Note (Part 3 #1 of case-study prompt): the proxy survival times built in
# Stage 1 (hospital deaths=day 1; FU deaths=day 90; alive=180; LFTU=30
# censored) introduce informative censoring (lost ≠ MAR). Discrete-time
# TMLE on this scheme has historically given a 180-day RD that differs
# from `survtmle`'s by ~1.1 pp; the gap is artifact, not a real estimator
# disagreement. Therefore: SURVTMLE is primary survival; DISCRETE-TIME
# TMLE is reported as sensitivity ONLY.
cr_log("Running survival TMLE (full cohort) — sensitivity only, see header note...")

result_surv_tmle <- tryCatch(
  run_survival_tmle(
    time_var      = time_var,
    event_var     = event_var,
    A             = A,
    W             = W_matrix,
    target_times  = target_times,
    sl_lib        = sel_libs$Q_model,
    sl_lib_censor = sel_libs$censor_model,
    seed          = cfg$seed,
    cv_folds      = cfg$superlearner$cv_folds,
    n_boot        = cfg$simulation$n_sims
  ),
  error = function(e) {
    cr_log(paste("Survival TMLE failed:", e$message))
    NULL
  }
)

if (!is.null(result_surv_tmle)) {
  cr_log("Survival TMLE results:")
  for (i in seq_len(nrow(result_surv_tmle$estimates))) {
    est <- result_surv_tmle$estimates[i, ]
    cr_log(paste("  t =", est$time, "days:",
                 "RD =", round(est$risk_difference, 4),
                 "[", round(est$rd_ci_lower, 4), ",",
                 round(est$rd_ci_upper, 4), "]",
                 "RR =", round(est$risk_ratio, 3)))
  }
}

pb5$tick()

# --- H. Matched survival-TMLE ---
cr_log("Running matched survival TMLE...")

result_surv_tmle_matched <- tryCatch(
  run_matched_survival_tmle(
    time_var      = time_var,
    event_var     = event_var,
    A             = A,
    W             = W_matrix,
    matched_idx   = matched_idx,
    target_times  = target_times,
    sl_lib        = sel_libs$Q_model,
    sl_lib_censor = sel_libs$censor_model,
    seed          = cfg$seed,
    cv_folds      = cfg$superlearner$cv_folds,
    n_boot        = cfg$simulation$n_sims
  ),
  error = function(e) {
    cr_log(paste("Matched survival TMLE failed:", e$message))
    NULL
  }
)

if (!is.null(result_surv_tmle_matched)) {
  cr_log("Matched survival TMLE results:")
  for (i in seq_len(nrow(result_surv_tmle_matched$estimates))) {
    est <- result_surv_tmle_matched$estimates[i, ]
    cr_log(paste("  t =", est$time, "days:",
                 "RD =", round(est$risk_difference, 4),
                 "RR =", round(est$risk_ratio, 3)))
  }
}

pb5$tick()

decisions <- log_decision(decisions, "stage5",
                          "Survival: Cox full + Cox matched + TMLE full + TMLE matched",
                          "All four estimators run for comparison",
                          type = "final analysis")

# --- Combine survival results ---
surv_results <- list(
  survtmle_primary   = result_survtmle,
  cox_full           = result_cox_full,
  cox_matched        = result_cox_matched,
  surv_tmle_full     = result_surv_tmle,
  surv_tmle_matched  = result_surv_tmle_matched
)

# Build comparison table
surv_comparison <- summarize_survival_results(
  list(result_survtmle, result_cox_full, result_cox_matched,
       result_surv_tmle, result_surv_tmle_matched)
)

if (!is.null(surv_comparison)) {
  cr_log("=== Survival Analysis Comparison ===")
  print(surv_comparison)
}

# --- Survival curves plot (with CI bands where available) ---
plot_survival_curves(result_surv_tmle,
                      file.path(cfg$paths$results, "survival_curves_full.png"))
plot_survival_curves(result_surv_tmle_matched,
                      file.path(cfg$paths$results, "survival_curves_matched.png"))

# Enhanced risk curves with confidence bands
tryCatch({
  plot_survival_curves_ci(result_surv_tmle,
                          file.path(cfg$paths$results, "risk_curves_full_ci.png"))
  plot_survival_curves_ci(result_surv_tmle_matched,
                          file.path(cfg$paths$results, "risk_curves_matched_ci.png"))
  cr_log("Saved risk curve plots with CI bands")
}, error = function(e) cr_log(paste("CI band plots skipped:", e$message)))

# --- Kaplan-Meier for visual reference ---
cr_log("Generating Kaplan-Meier curves...")
km_data <- data.frame(
  time  = time_var[!is.na(time_var) & !is.na(event_var) & !is.na(A)],
  event = event_var[!is.na(time_var) & !is.na(event_var) & !is.na(A)],
  group = factor(A[!is.na(time_var) & !is.na(event_var) & !is.na(A)],
                  labels = c("Control", "Rescue.Co"))
)

if (nrow(km_data) > 0) {
  km_fit <- survfit(Surv(time, event) ~ group, data = km_data)

  png(file.path(cfg$paths$results, "km_curves.png"), width = 800, height = 500)
  plot(km_fit, col = c("steelblue", "coral"), lwd = 2,
       xlab = "Time (days)", ylab = "Survival probability",
       main = "Kaplan-Meier Survival Curves")
  legend("bottomleft", legend = c("Control", "Rescue.Co"),
         col = c("steelblue", "coral"), lwd = 2)
  dev.off()
  cr_log("Saved KM curves")
}

# --- cleanTMLE survival estimators ---
if (!is.null(lock)) {
  cr_log("Running cleanTMLE survival estimators...")

  # Build cr_spec for time-to-event analysis
  surv_data <- data.frame(
    A = A, time = time_var, event = event_var,
    W_matrix[, seq_len(min(5, ncol(W_matrix))), drop = FALSE]
  )
  surv_data <- surv_data[complete.cases(surv_data), ]

  ct_surv_spec <- tryCatch({
    spec <- cleanTMLE::specify_models(data = surv_data)
    spec <- cleanTMLE::identify_outcome(spec, event, type = "time_to_event")
    spec <- cleanTMLE::identify_treatment(spec, A,
      formula = as.formula(paste("~", paste(colnames(W_matrix)[seq_len(min(5, ncol(W_matrix)))], collapse = "+"))))
    spec
  }, error = function(e) {
    cr_log(paste("cleanTMLE cr_spec creation failed:", e$message))
    NULL
  })

  # IPW risk curves
  ct_ipw_risk <- tryCatch({
    cleanTMLE::estimate_ipwrisk(ct_surv_spec,
      risk_time = target_times,
      weight_type = "iptw", nboot = 20, seed = cfg$seed)
  }, error = function(e) {
    cr_log(paste("cleanTMLE estimate_ipwrisk failed:", e$message))
    NULL
  })

  if (!is.null(ct_ipw_risk)) {
    cr_log("cleanTMLE IPW risk estimates:")
    print(ct_ipw_risk)

    tryCatch({
      ipw_plot <- plot(ct_ipw_risk)
      ggsave(file.path(cfg$paths$results, "cleanTMLE_ipw_risk_curves.png"),
             ipw_plot, width = 8, height = 5)
    }, error = function(e) cr_log(paste("IPW plot failed:", e$message)))

    # Weight summary and extreme weights from IPW
    tryCatch({
      wt_tbl <- cleanTMLE::make_wt_summary_table(ct_ipw_risk)
      cr_log("cleanTMLE IPW weight summary:")
      print(wt_tbl)
    }, error = function(e) cr_log(paste("IPW weight summary failed:", e$message)))
  }

  # IPW hazard ratio
  ct_ipwhr <- tryCatch(
    cleanTMLE::estimate_ipwhr(ct_surv_spec),
    error = function(e) {
      cr_log(paste("cleanTMLE estimate_ipwhr failed:", e$message))
      NULL
    }
  )
  if (!is.null(ct_ipwhr)) {
    cr_log("cleanTMLE IPW HR:")
    print(ct_ipwhr)
  }

  # Survival TMLE
  ct_surv_tmle <- tryCatch(
    cleanTMLE::estimate_surv_tmle(
      data        = surv_data,
      treatment   = "A",
      time        = "time",
      event       = "event",
      covariates  = colnames(W_matrix)[seq_len(min(5, ncol(W_matrix)))],
      target_times = target_times,
      sl_library  = cfg$superlearner$candidate_learners
    ),
    error = function(e) {
      cr_log(paste("cleanTMLE survival TMLE failed:", e$message))
      NULL
    }
  )
  if (!is.null(ct_surv_tmle)) {
    cr_log("cleanTMLE survival TMLE:")
    print(ct_surv_tmle)
  }

  audit <- cleanTMLE::record_stage(audit, "Stage 5", "Survival analysis complete")
  save_stage_output(lock, "stage5_lock.rds")
  save_stage_output(audit, "stage5_audit.rds")
}

# --- Save ---
save_stage_output(surv_results, "stage5_survival_results.rds")
save_stage_output(surv_comparison, "stage5_survival_comparison.rds")
save_stage_output(decisions, "stage5_decisions.rds")

if (!is.null(surv_comparison)) {
  write.csv(surv_comparison, file.path(cfg$paths$results, "survival_outcome_comparison.csv"),
            row.names = FALSE)
}

cr_log("Stage 5 complete. Survival outcomes unlocked and analyzed.")
