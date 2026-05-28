#!/usr/bin/env Rscript
# ============================================================================
# Clean-Room TMLE Simulation Study
# ============================================================================
#
# Demonstrates the full cleanTMLE staged workflow:
#
#   * Stage 1a:  create_analysis_lock + attach_estimand
#                + declare_sensitivity_plan + define_negative_control
#   * Stage 1b:  checkpoint_cohort_adequacy + estimate_design_precision
#   * Outcome blinding: mask_outcome (Stage 4 estimators are now refused
#                       on the masked lock)
#   * Stage 2a:  fit_ps_glm + compute_ps_diagnostics + checkpoint_balance
#   * Stage 2b:  run_plasmode_feasibility + run_plasmode_dq_stress
#                + select_tmle_candidate + lock_primary_tmle_spec
#                + gate_check
#   * Stage 3:   run_residual_confounding_stage on the negative-control
#                outcome -> checkpoint_residual_bias
#   * Pre-outcome gate: gate_all + authorize_outcome_analysis
#   * Stage 4:   unmask_outcome -> Crude / IPTW / PS Match / TMLE /
#                TMLE_CF / Match_TMLE
#   * Audit + decision log persisted as RDS artifacts
#
# Three scenarios:
#   A. Good overlap          (overlap_strength = 0.5, no U)
#   B. Marginal overlap      (overlap_strength = 1.5, no U)
#   C. Unmeasured confounding (overlap_strength = 0.5, U with OR = 2)
#
# Six estimators per scenario:
#   Crude, IPTW, PS Match, TMLE, TMLE_CF (cross-fitted), Match_TMLE
#
# Results are saved incrementally so progress is monitorable.
# ============================================================================

# Force line-buffered output when running non-interactively (e.g., nohup).
.flush <- function() if (!interactive()) flush(stdout())

# Install cleanTMLE from local source if needed.
if (!requireNamespace("cleanTMLE", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::install(file.path(getwd(), "cleanTMLE"), quick = TRUE, quiet = TRUE)
  } else {
    stop("Please install cleanTMLE first: devtools::install('cleanTMLE')")
  }
}
library(cleanTMLE)

# ── Configuration ────────────────────────────────────────────────────────────

config <- list(
  # Sample size and rep counts
  n_obs           = 2000L,
  n_reps          = 200L,    # MC replicates per scenario
  n_truth         = 100000L, # large-sample ground truth
  plasmode_reps   = 30L,     # plasmode reps for candidate selection
  dq_reps         = 30L,     # plasmode reps per DQ scenario level

  # SuperLearner library used for: per-replicate TMLE_CF, plasmode-Q, etc.
  sl_library      = c("SL.glm", "SL.glmnet", "SL.mean"),

  # Cross-fitting folds for TMLE_CF
  n_folds         = 2L,

  # Gate thresholds (consumed by checkpoint and gate_check).
  cohort_min_n        = 500L,
  cohort_min_per_arm  = 200L,
  cohort_min_events   = 50L,
  balance_max_smd     = 0.10,
  balance_min_ess_pct = 50,
  gate_targets        = list(
    max_abs_bias = 0.02,
    min_coverage = 0.88,
    se_sd_low    = 0.7,
    se_sd_high   = 1.3
  ),

  results_dir = "results_new",
  seed        = 2026L
)

if (!dir.exists(config$results_dir)) dir.create(config$results_dir, recursive = TRUE)

cat("============================================================\n")
cat("  Clean-Room TMLE Simulation Study\n")
cat("============================================================\n")
cat(sprintf("  N = %d   MC reps = %d   plasmode reps = %d   DQ reps = %d\n",
            config$n_obs, config$n_reps, config$plasmode_reps, config$dq_reps))
cat(sprintf("  SL library: %s\n", paste(config$sl_library, collapse = ", ")))
cat(sprintf("  Cross-fitting: %d folds\n", config$n_folds))
cat("============================================================\n\n")


# ── DGP Functions ────────────────────────────────────────────────────────────

generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL,
                          U_prevalence = 0, U_trt_OR = 1, U_out_OR = 1) {
  if (!is.null(seed)) set.seed(seed)

  age         <- rnorm(n, mean = 55, sd = 10)
  sex         <- rbinom(n, 1, 0.55)
  biomarker   <- rnorm(n, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n, 1, 0.12)

  # Unmeasured confounder (not returned in the data frame the analyst sees).
  U <- if (U_prevalence > 0) rbinom(n, 1, U_prevalence) else rep(0L, n)

  lp_trt <- -0.5 +
    overlap_strength * (0.03 * (age - 55) +
                        0.8 * sex +
                        0.6 * biomarker +
                        0.5 * ckd +
                        0.3 * comorbidity) +
    log(U_trt_OR) * U
  ps_true <- plogis(lp_trt)
  treatment <- rbinom(n, 1, ps_true)

  lp_out <- -2.5 +
    0.015 * (age - 55) +
    0.3 * sex +
    0.2 * biomarker +
    0.6 * ckd +
    0.25 * comorbidity +
    effect_size / 0.15 * treatment +
    log(U_out_OR) * U

  event_24 <- rbinom(n, 1, plogis(lp_out))

  # Negative-control outcome: depends on covariates only, not treatment.
  lp_nc <- -1.0 + 0.01 * (age - 55) + 0.1 * sex + 0.15 * biomarker
  nc_outcome <- rbinom(n, 1, plogis(lp_nc))

  data.frame(
    age         = round(age, 1),
    sex         = sex,
    biomarker   = round(biomarker, 3),
    comorbidity = comorbidity,
    ckd         = ckd,
    treatment   = treatment,
    event_24    = event_24,
    nc_outcome  = nc_outcome,
    stringsAsFactors = FALSE
  )
}


compute_truth <- function(n_truth, overlap_strength, effect_size, seed,
                          U_prevalence = 0, U_out_OR = 1) {
  set.seed(seed)
  age         <- rnorm(n_truth, mean = 55, sd = 10)
  sex         <- rbinom(n_truth, 1, 0.55)
  biomarker   <- rnorm(n_truth, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n_truth, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n_truth, 1, 0.12)
  U           <- if (U_prevalence > 0) rbinom(n_truth, 1, U_prevalence) else rep(0L, n_truth)

  lp_out_base <- -2.5 +
    0.015 * (age - 55) +
    0.3 * sex +
    0.2 * biomarker +
    0.6 * ckd +
    0.25 * comorbidity +
    log(U_out_OR) * U

  coef_trt <- effect_size / 0.15
  risk_1 <- mean(plogis(lp_out_base + coef_trt))
  risk_0 <- mean(plogis(lp_out_base))

  list(risk_1 = risk_1, risk_0 = risk_0, RD = risk_1 - risk_0)
}


# ── Scenarios ────────────────────────────────────────────────────────────────

scenarios <- list(
  good_overlap = list(
    label            = "Scenario A: Good Overlap",
    overlap_strength = 0.5,
    effect_size      = -0.05,
    description      = "Moderate confounding, good PS overlap."
  ),
  marginal_overlap = list(
    label            = "Scenario B: Marginal Overlap",
    overlap_strength = 1.5,
    effect_size      = -0.05,
    description      = "Strong confounding, marginal PS overlap."
  ),
  unmeasured_conf = list(
    label            = "Scenario C: Unmeasured Confounding",
    overlap_strength = 0.5,
    effect_size      = -0.05,
    U_prevalence     = 0.20,
    U_trt_OR         = 2.0,
    U_out_OR         = 2.0,
    description      = "Good measured overlap, but binary U (prev=0.20, OR=2.0) confounds A and Y."
  )
)


# ── Compute Ground Truth ─────────────────────────────────────────────────────

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("Computing ground truth...\n"); .flush()
truths <- lapply(scenarios, function(sc) {
  compute_truth(config$n_truth, sc$overlap_strength, sc$effect_size,
                seed = config$seed,
                U_prevalence = sc$U_prevalence %||% 0,
                U_out_OR     = sc$U_out_OR %||% 1)
})
for (nm in names(truths)) {
  cat(sprintf("  %s: true RD = %.5f\n", scenarios[[nm]]$label, truths[[nm]]$RD))
}
cat("\n"); .flush()


# ── TMLE Candidate Grid ──────────────────────────────────────────────────────
# Vary both library (GLM vs glmnet) and truncation. The grid produces
# meaningfully different candidates under marginal overlap and avoids the
# "all candidates produce identical metrics" tiebreak we observed when the
# grid varied truncation alone.

# NOTE: glmnet candidates dropped. SL.glmnet runs away (CPU/memory spin) on the
# near-positivity-violation marginal-overlap synthetic designs in the DQ stress
# test; glm and glmnet are near-equivalent on this DGP (see manuscript
# Sec. plasmode), and GLM PS gives bounded, stable fits.
tmle_candidates <- list(
  tmle_candidate("glm_t01",    "GLM PS, trunc=0.01",
                 g_library = "SL.glm",     truncation = 0.01),
  tmle_candidate("glm_t05",    "GLM PS, trunc=0.05",
                 g_library = "SL.glm",     truncation = 0.05)
)


# ── DQ Scenario Specification ────────────────────────────────────────────────
# Asymmetric (sensitivity, specificity) for treatment AND outcome
# misclassification, plus a per-arm differential outcome scenario.

dq_spec <- list(
  covariate_missingness = list(
    fractions = c(0.05, 0.10, 0.20),
    variables = NULL    # all covariates
  ),
  treatment_misclass = list(
    sensitivity = c(0.95, 0.90),
    specificity = c(0.99, 0.95)
  ),
  outcome_misclass = list(
    sensitivity = c(0.95, 0.90),
    specificity = c(0.99, 0.95)
  ),
  unmeasured_confounding = list(
    U_prevalence   = 0.20,
    U_treatment_OR = c(1.5, 2.0),
    U_outcome_OR   = c(1.5, 2.0)
  )
)


# ── Single-replicate estimation function ─────────────────────────────────────
# Receives the unmasked lock for Stage 4 estimation. Each estimator is
# wrapped in tryCatch so a single failure does not abort the rep.

run_one_replicate <- function(dat, lock, ps_fit, truth_rd, n_folds = 1L) {
  covariates <- lock$covariates

  results <- list()

  results$Crude <- tryCatch({
    cr <- run_crude_workflow(lock)
    list(estimate = cr$estimate, se = cr$se,
         ci_lower = cr$ci_lower, ci_upper = cr$ci_upper)
  }, error = function(e) list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))

  results$IPTW <- tryCatch({
    ip <- run_iptw_workflow(lock, ps_fit)
    list(estimate = ip$estimate, se = ip$se,
         ci_lower = ip$ci_lower, ci_upper = ip$ci_upper)
  }, error = function(e) list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))

  results$`PS Match` <- tryCatch({
    mt <- run_match_workflow(lock, ps_fit)
    list(estimate = mt$estimate, se = mt$se,
         ci_lower = mt$ci_lower, ci_upper = mt$ci_upper)
  }, error = function(e) list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))

  results$TMLE <- tryCatch({
    g_fit    <- fit_tmle_treatment_mechanism(lock, ps_fit)
    Q_fit    <- fit_tmle_outcome_mechanism(lock, g_fit)
    tmle_upd <- run_tmle_targeting_step(g_fit, Q_fit)
    te       <- extract_tmle_estimate(tmle_upd)
    ate      <- te$estimates$ATE
    list(estimate = ate$estimate, se = ate$se,
         ci_lower = ate$ci_lower, ci_upper = ate$ci_upper)
  }, error = function(e) list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))

  if (n_folds > 1L) {
    results$TMLE_CF <- tryCatch({
      cf <- estimate_tmle_risk_point(
        data       = dat,
        treatment  = lock$treatment,
        outcome    = lock$outcome,
        covariates = covariates,
        sl_library = c("SL.glm", "SL.glmnet", "SL.mean"),
        truncate   = if (!is.null(lock$primary_tmle_spec))
                       lock$primary_tmle_spec$truncation else 0.01,
        n_folds    = n_folds
      )
      ate <- cf$estimates$ATE
      list(estimate = ate$estimate, se = ate$se,
           ci_lower = ate$ci_lower, ci_upper = ate$ci_upper)
    }, error = function(e) list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
  }

  results$Match_TMLE <- tryCatch({
    mt <- run_match_workflow(lock, ps_fit)
    matched_idx <- as.integer(rownames(mt$matched_data))
    mt_tmle <- run_matched_tmle(lock, ps_fit, subset_idx = matched_idx)
    ate <- mt_tmle$estimates$ATE
    list(estimate = ate$estimate, se = ate$se,
         ci_lower = ate$ci_lower, ci_upper = ate$ci_upper)
  }, error = function(e) list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))

  rows <- lapply(names(results), function(m) {
    r <- results[[m]]
    covers <- if (!is.na(r$ci_lower) && !is.na(r$ci_upper))
      as.integer(r$ci_lower <= truth_rd & truth_rd <= r$ci_upper) else NA_integer_
    data.frame(
      method   = m,
      estimate = r$estimate,
      se       = r$se,
      ci_lower = r$ci_lower,
      ci_upper = r$ci_upper,
      covers   = covers,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


build_summary_table <- function(results_df, truth_rd) {
  methods <- unique(results_df$method)
  rows <- lapply(methods, function(m) {
    sub   <- results_df[results_df$method == m, ]
    valid <- !is.na(sub$estimate)
    ests  <- sub$estimate[valid]
    ses   <- sub$se[valid]
    n_ok  <- sum(valid)

    if (n_ok == 0) {
      return(data.frame(
        method = m, n_ok = 0, mean_est = NA, bias = NA, abs_bias = NA,
        emp_sd = NA, mean_se = NA, se_sd_ratio = NA, rmse = NA,
        coverage = NA, mc_se_cov = NA, ci_width = NA,
        stringsAsFactors = FALSE
      ))
    }

    bias_val  <- mean(ests) - truth_rd
    emp_sd    <- sd(ests)
    mean_se   <- mean(ses, na.rm = TRUE)
    se_sd_r   <- if (emp_sd > 0) mean_se / emp_sd else NA
    cov_p     <- mean(sub$covers[valid], na.rm = TRUE)
    mc_se_cov <- sqrt(cov_p * (1 - cov_p) / n_ok)

    data.frame(
      method      = m,
      n_ok        = n_ok,
      mean_est    = round(mean(ests), 5),
      bias        = round(bias_val, 5),
      abs_bias    = round(abs(bias_val), 5),
      emp_sd      = round(emp_sd, 5),
      mean_se     = round(mean_se, 5),
      se_sd_ratio = round(se_sd_r, 3),
      rmse        = round(sqrt(mean((ests - truth_rd)^2)), 5),
      coverage    = round(cov_p, 3),
      mc_se_cov   = round(mc_se_cov, 4),
      ci_width    = round(mean(sub$ci_upper[valid] - sub$ci_lower[valid],
                              na.rm = TRUE), 5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


interim_path <- function(sc_name) {
  file.path(config$results_dir,
            sprintf("interim_%s.rds", sc_name))
}

save_interim <- function(sc_name, rep_results, scenario_meta) {
  obj <- list(
    results  = do.call(rbind, Filter(Negate(is.null), rep_results)),
    meta     = scenario_meta,
    n_done   = sum(!vapply(rep_results, is.null, logical(1))),
    n_total  = length(rep_results),
    saved_at = Sys.time()
  )
  saveRDS(obj, interim_path(sc_name))
}


# ── Main Simulation Loop ────────────────────────────────────────────────────

all_results   <- list()
all_summaries <- list()
all_meta      <- list()
all_audit     <- list()

t_start_global <- Sys.time()

for (sc_name in names(scenarios)) {
  sc       <- scenarios[[sc_name]]
  truth_rd <- truths[[sc_name]]$RD

  # ── Resume guard: skip scenarios already completed in a prior (crashed) run.
  # A scenario is "done" when its lock RDS exists and its interim has all reps.
  # Only fires on restart against a partially populated results_dir; on a fresh
  # run no markers exist so this is a no-op. Reconstructs the in-memory state
  # (all_results/summaries/meta/audit) the final assembly needs.
  .lock_file <- file.path(config$results_dir, sprintf("lock_%s.rds", sc_name))
  .intr_file <- interim_path(sc_name)
  if (file.exists(.lock_file) && file.exists(.intr_file)) {
    .intr <- readRDS(.intr_file)
    if (isTRUE(!is.null(.intr$n_done) && !is.null(.intr$n_total) &&
               .intr$n_done == .intr$n_total)) {
      cat(sprintf("\n=== %s: already complete (%d/%d reps) — resuming, skip ===\n",
                  sc$label, .intr$n_done, .intr$n_total)); .flush()
      all_results[[sc_name]]   <- .intr$results
      all_summaries[[sc_name]] <- build_summary_table(.intr$results, truth_rd)
      all_meta[[sc_name]]      <- .intr$meta
      .audit_file <- file.path(config$results_dir, sprintf("audit_%s.rds", sc_name))
      if (file.exists(.audit_file)) all_audit[[sc_name]] <- readRDS(.audit_file)
      next
    }
  }

  cat(sprintf("\n=== %s ===\n", sc$label))
  cat(sprintf("  %s\n", sc$description))
  cat(sprintf("  True RD = %.5f\n\n", truth_rd))
  .flush()

  # ── Stage 1a: lock + estimand + sensitivity plan + NCO declaration ──
  cat("  Stage 1a: lock + estimand + NCO declaration...\n"); .flush()

  ref_dat <- generate_data(config$n_obs, sc$overlap_strength, sc$effect_size,
                           seed = config$seed,
                           U_prevalence = sc$U_prevalence %||% 0,
                           U_trt_OR     = sc$U_trt_OR %||% 1,
                           U_out_OR     = sc$U_out_OR %||% 1)

  lock <- create_analysis_lock(
    data          = ref_dat,
    treatment     = "treatment",
    outcome       = "event_24",
    covariates    = c("age", "sex", "biomarker", "comorbidity", "ckd"),
    sl_library    = config$sl_library,
    plasmode_reps = config$plasmode_reps,
    seed          = config$seed
  )
  lock <- attach_estimand(lock,
    description          = "Effect of A on 24-month event risk",
    population           = "Adults eligible at index",
    treatment_strategies = c("Treatment", "Control"),
    outcome_label        = "Primary event by 24 months",
    followup             = "24 months",
    contrast             = "risk_difference",
    statistical_estimand = "E_W{E[Y|A=1,W] - E[Y|A=0,W]}"
  )
  lock <- declare_sensitivity_plan(lock, "truncation",
    description = "Truncation sensitivity (0.01, 0.05, 0.10)",
    settings    = list(thresholds = c(0.01, 0.05, 0.10))
  )
  lock <- define_negative_control(lock, "nc_outcome",
    description = "Outcome driven by covariates only, no treatment effect"
  )
  validate_analysis_lock(lock)

  audit  <- create_audit_log(lock)
  audit  <- record_stage(audit, "Stage 1a",
              "Lock created; estimand, sensitivity plan, NCO declared.")
  audit  <- record_decision_log_entry(audit, "Stage 1a",
              decision_type = "lock",
              description   = sprintf("Lock %s for %s", lock$lock_hash, sc$label),
              rationale     = "Pre-outcome lock before any estimation.")

  # ── Stage 1b: cohort adequacy + design precision ──────────────────────
  cat("  Stage 1b: cohort adequacy + design precision...\n"); .flush()
  cp1 <- checkpoint_cohort_adequacy(lock,
    min_n_per_arm = config$cohort_min_per_arm,
    min_events    = config$cohort_min_events
  )
  audit <- record_checkpoint(audit, cp1)
  cat(sprintf("    Check Point 1: %s\n", cp1$decision)); .flush()

  dp <- tryCatch(estimate_design_precision(lock, target_mdd = 0.05),
                 error = function(e) NULL)
  if (!is.null(dp)) {
    audit <- record_decision_log_entry(audit, "Stage 1b",
                decision_type = "design_precision",
                description   = "Design-precision estimate computed",
                rationale     = "Pre-outcome MDD diagnostic.",
                metrics       = list(target_mdd = 0.05))
  }

  # ── Stage 2a: PS diagnostics + balance checkpoint ───────────────────────
  # Stage 2a does not consult the outcome; we run it on the unmasked lock
  # and apply the mask after the plasmode (which needs the real Y to fit
  # Q0 and generate synthetic outcomes).
  cat("  Stage 2a: PS diagnostics + balance checkpoint...\n"); .flush()
  ps_fit_ref <- fit_ps_glm(lock)
  ps_diag    <- compute_ps_diagnostics(ps_fit_ref)
  cp2 <- tryCatch(
    checkpoint_balance(ps_diag,
                       max_smd     = config$balance_max_smd,
                       min_ess_pct = config$balance_min_ess_pct,
                       lock_hash   = lock$lock_hash),
    error = function(e) {
      message("    Balance checkpoint failed: ", e$message)
      NULL
    }
  )
  if (!is.null(cp2)) {
    audit <- record_checkpoint(audit, cp2)
    cat(sprintf("    Check Point 2: %s\n", cp2$decision)); .flush()
  }
  cat(sprintf("    PS range: [%.3f, %.3f]  ESS: %.0f / %d (%.0f%%)\n",
              min(ps_fit_ref$ps), max(ps_fit_ref$ps),
              ps_diag$ess$ess[3], ps_diag$ess$n[3], ps_diag$ess$ess_pct[3])); .flush()

  # ── Stage 2b: plasmode candidate selection + DQ stress ──────────────────
  # Plasmode fits a covariate-only Q0 from the real outcome and replaces
  # Y with synthetic Y_sim. The analyst learns nothing about the realised
  # treatment-outcome association, which is the sense in which the
  # plasmode is "outcome-blind" (Dang et al. 2023).
  cat("  Stage 2b: plasmode candidate selection...\n"); .flush()
  plas <- run_plasmode_feasibility(
    lock,
    tmle_candidates = tmle_candidates,
    effect_sizes    = c(0.03, 0.05),
    reps            = config$plasmode_reps,
    verbose         = FALSE
  )
  best <- select_tmle_candidate(plas, rule = "min_rmse")
  cat(sprintf("    Selected: %s (g=%s, trunc = %.2f)\n",
              best$candidate_id, paste(best$g_library, collapse = "+"),
              best$truncation)); .flush()

  gate_2b <- gate_check(plas$metrics, scenario_name = sc$label,
                         targets = config$gate_targets,
                         method  = best$candidate_id)
  cat(sprintf("    Stage 2b gate: %s\n", gate_2b$decision)); .flush()
  audit <- record_decision_log_entry(audit, "Stage 2b",
              decision_type = "candidate_selection",
              description   = sprintf("Selected %s by min_rmse", best$candidate_id),
              rationale     = "Plasmode-selected TMLE specification.",
              metrics       = list(decision = gate_2b$decision))

  cat("  Stage 2b*: DQ stress test...\n"); .flush()
  dq_results <- tryCatch({
    run_plasmode_dq_stress(
      lock,
      tmle_candidates = tmle_candidates,
      effect_sizes    = c(0.05),
      reps            = config$dq_reps,
      data_quality_scenarios = dq_spec,
      verbose         = TRUE
    )
  }, error = function(e) {
    message("    DQ stress test failed: ", e$message)
    NULL
  })
  if (!is.null(dq_results)) {
    cat("\n    DQ degradation summary (selected candidate):\n")
    dq_deg <- summarize_dq_degradation(dq_results)
    sel_deg <- dq_deg[dq_deg$candidate == best$candidate_id, ]
    if (nrow(sel_deg) > 0) {
      print(sel_deg[, c("scenario", "level", "bias_degraded", "rmse_ratio",
                         "cov_degraded", "cov_drop")], row.names = FALSE)
    }
    cat("\n"); .flush()

    saveRDS(dq_results,
            file.path(config$results_dir,
                      sprintf("dq_stress_%s.rds", sc_name)))
    audit <- record_decision_log_entry(audit, "Stage 2b*",
                decision_type = "dq_stress",
                description   = "Plasmode DQ stress test executed.",
                rationale     = "Quantitative SPIFD2 mapping check.")
  }

  # ── Outcome blinding ────────────────────────────────────────────────────
  # Lock the primary TMLE spec, then mask. Stage 4 estimators on the
  # masked lock will refuse to run; only after the gate authorises do we
  # unmask and proceed to the comparative analysis.
  original_lock <- lock_primary_tmle_spec(lock, best)
  masked        <- mask_outcome(original_lock)
  audit         <- record_stage(audit, "Mask",
                    "Outcome column masked; Stage 3 follows on masked lock.",
                    decision = "OUTCOME-BLIND")

  # ── Stage 3: NCO residual confounding ───────────────────────────────────
  # The negative-control outcome (`nc_outcome`) is NOT the masked column;
  # the residual-confounding check legitimately uses it to probe bias.
  cat("  Stage 3: residual-confounding check via NCO...\n"); .flush()
  stage3 <- tryCatch(
    run_residual_confounding_stage(masked, ps_fit_ref),
    error = function(e) {
      message("    Stage 3 failed: ", e$message)
      NULL
    }
  )
  cp3 <- if (!is.null(stage3) && !is.null(stage3$checkpoint)) stage3$checkpoint else NULL
  if (!is.null(cp3)) {
    audit <- record_checkpoint(audit, cp3)
    cat(sprintf("    Check Point 3: %s\n", cp3$decision)); .flush()
  }

  # ── Pre-outcome gate ────────────────────────────────────────────────────
  gate_args <- Filter(Negate(is.null), list(cp1, cp2, cp3))
  gate_overall <- do.call(gate_all, c(gate_args, list(allow_flag = TRUE)))
  cat(sprintf("    Pre-outcome gate: %s\n", gate_overall$decision)); .flush()
  audit <- record_checkpoint(audit, gate_overall)

  # Authorisation gate. authorize_outcome_analysis() returns a checkpoint
  # whose `authorized` slot is TRUE when no required stage has STOP.
  auth <- tryCatch(authorize_outcome_analysis(audit),
                   error = function(e) list(authorized = FALSE,
                                            decision   = "STOP",
                                            rationale  = conditionMessage(e)))
  auth_passed <- isTRUE(auth$authorized) ||
                 (!is.null(auth$decision) && auth$decision != "STOP")
  if (auth_passed) {
    cat(sprintf("    Outcome access AUTHORISED (gate decision: %s).\n",
                auth$decision)); .flush()
    audit <- record_decision_log_entry(audit, "Gate",
                decision_type = "authorize",
                description   = "Outcome unblinding authorised by gate.",
                rationale     = "All required checkpoints passed.")
  } else {
    cat(sprintf("    Outcome access DENIED: %s\n",
                if (is.null(auth$rationale)) "(no reason given)" else auth$rationale)); .flush()
    audit <- record_decision_log_entry(audit, "Gate",
                decision_type = "authorize",
                description   = "Outcome unblinding NOT authorised.",
                rationale     = "Forced unblinding for the simulation only.")
  }

  # ── Outcome unmasking (force, since this is a simulation study) ────────
  unmasked_lock <- unmask_outcome(masked, original_lock)

  scenario_meta <- list(
    label    = sc$label,
    truth_rd = truth_rd,
    selected = best$candidate_id,
    gate     = gate_overall$decision,
    auth     = auth,
    ps_range = range(ps_fit_ref$ps),
    ess_pct  = ps_diag$ess$ess_pct[3],
    plasmode = plas$metrics,
    dq_stress = if (!is.null(dq_results)) dq_results$metrics else NULL,
    cp1 = cp1, cp2 = cp2, cp3 = cp3
  )
  all_meta[[sc_name]]  <- scenario_meta

  # ── Stage 4: Monte Carlo loop on the unmasked workflow ─────────────────
  cat(sprintf("\n  Running %d Monte Carlo replicates...\n", config$n_reps)); .flush()
  rep_results <- vector("list", config$n_reps)
  t_start <- Sys.time()

  for (rep_i in seq_len(config$n_reps)) {
    rep_seed <- config$seed + rep_i * 1000L + match(sc_name, names(scenarios))
    dat_rep  <- generate_data(config$n_obs, sc$overlap_strength,
                              sc$effect_size, seed = rep_seed,
                              U_prevalence = sc$U_prevalence %||% 0,
                              U_trt_OR     = sc$U_trt_OR %||% 1,
                              U_out_OR     = sc$U_out_OR %||% 1)

    lock_rep <- create_analysis_lock(
      data          = dat_rep,
      treatment     = "treatment",
      outcome       = "event_24",
      covariates    = c("age", "sex", "biomarker", "comorbidity", "ckd"),
      sl_library    = config$sl_library,
      plasmode_reps = config$plasmode_reps,
      seed          = config$seed
    )
    lock_rep <- lock_primary_tmle_spec(lock_rep, best)
    ps_rep   <- fit_ps_glm(lock_rep)

    rep_results[[rep_i]] <- tryCatch(
      run_one_replicate(dat_rep, lock_rep, ps_rep, truth_rd,
                        n_folds = config$n_folds),
      error = function(e) {
        message("  Rep ", rep_i, " failed: ", e$message)
        NULL
      }
    )

    if (rep_i %% 10 == 0 || rep_i == config$n_reps) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      rate    <- elapsed / rep_i
      eta     <- rate * (config$n_reps - rep_i)
      cat(sprintf("    rep %d/%d  (%.1f min elapsed, ~%.1f min remaining)\n",
                  rep_i, config$n_reps, elapsed, eta)); .flush()
      save_interim(sc_name, rep_results, scenario_meta)
    }
  }

  results_df <- do.call(rbind, Filter(Negate(is.null), rep_results))
  all_results[[sc_name]] <- results_df

  summary_df <- build_summary_table(results_df, truth_rd)
  all_summaries[[sc_name]] <- summary_df

  cat(sprintf("\n  --- %s: Summary (%d reps) ---\n", sc$label, config$n_reps))
  cat(sprintf("  True RD: %.5f\n\n", truth_rd))
  print(summary_df, row.names = FALSE); cat("\n")

  # Save the audit + decision log for this scenario.
  saveRDS(audit, file.path(config$results_dir,
                            sprintf("audit_%s.rds", sc_name)))
  write.csv(export_audit_trail(audit),
            file.path(config$results_dir,
                      sprintf("audit_%s.csv", sc_name)),
            row.names = FALSE)
  write.csv(export_decision_log(audit),
            file.path(config$results_dir,
                      sprintf("decision_log_%s.csv", sc_name)),
            row.names = FALSE)
  saveRDS(unmasked_lock, file.path(config$results_dir,
                                    sprintf("lock_%s.rds", sc_name)))
  all_audit[[sc_name]] <- audit
}

total_time <- as.numeric(difftime(Sys.time(), t_start_global, units = "mins"))


# ── Final Output ─────────────────────────────────────────────────────────────

cat("\n================================================================\n")
cat("  PLASMODE CANDIDATE SELECTION RESULTS\n")
cat("================================================================\n\n")

for (sc_name in names(all_meta)) {
  meta <- all_meta[[sc_name]]
  cat(sprintf("--- %s ---\n", meta$label))
  cat(sprintf("  Selected: %s   Gate: %s   Authorised: %s\n",
              meta$selected, meta$gate,
              if (isTRUE(meta$auth$authorized)) "TRUE" else "FALSE"))
  cat(sprintf("  PS range: [%.3f, %.3f]   ESS%%: %.1f\n\n",
              meta$ps_range[1], meta$ps_range[2], meta$ess_pct))
  cat("  Plasmode performance (per candidate, averaged over effect sizes):\n")
  print(meta$plasmode, row.names = FALSE); cat("\n\n")
}

cat("================================================================\n")
cat("  FINAL MONTE CARLO COMPARISON\n")
cat("================================================================\n\n")
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  truth_rd <- truths[[sc_name]]$RD
  cat(sprintf("--- %s (true RD = %.5f) ---\n", sc$label, truth_rd))
  print(all_summaries[[sc_name]], row.names = FALSE); cat("\n")
}

cat(sprintf("Total runtime: %.1f minutes\n\n", total_time))


# ── Persist Final Results ────────────────────────────────────────────────────

final_output <- list(
  config    = config,
  scenarios = scenarios,
  truths    = truths,
  meta      = all_meta,
  results   = all_results,
  summaries = all_summaries,
  audit     = all_audit,
  runtime_minutes = total_time
)
saveRDS(final_output, file.path(config$results_dir, "simulation_results.rds"))
saveRDS(config,       file.path(config$results_dir, "simulation_config.rds"))

for (sc_name in names(all_summaries)) {
  write.csv(all_summaries[[sc_name]],
            file.path(config$results_dir, sprintf("summary_%s.csv", sc_name)),
            row.names = FALSE)
}
for (sc_name in names(all_meta)) {
  write.csv(all_meta[[sc_name]]$plasmode,
            file.path(config$results_dir, sprintf("plasmode_%s.csv", sc_name)),
            row.names = FALSE)
}

cat("Results saved to ", config$results_dir, "/\n", sep = "")
cat("Done.\n")
