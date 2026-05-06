# ============================================================
# Stage 4: Binary Outcome Analysis (GOSE)
# ============================================================
# Unlock GOSE-based analyses.
# Implements:
#   A. PS-matched logistic regression
#   B. PS-matched TMLE
#   C. Full-cohort TMLE
# ============================================================

library(SuperLearner)
library(tmle)

# --- Source helpers ---
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/tmle_binary.R")
source("rescueCo/R/plotting.R")

# --- Load config and prior outputs ---
cfg <- load_cr_config()
cr_log("=== Stage 4: Binary Outcome Analysis ===")

dat          <- load_stage_output("stage1_cohort.rds")
W_matrix     <- load_stage_output("stage1_W_matrix.rds")
match_result <- load_stage_output("stage2_match_result.rds")
sel_libs     <- tryCatch(
  load_stage_output("stage3_selected_libraries.rds"),
  error = function(e) {
    # The modernised Stage 3 doesn't write a separate selected-libraries
    # RDS — the locked TMLE spec is on the `lock` object directly. Build
    # a sensible default here so downstream code keeps working.
    list(Q_model      = cfg$superlearner$candidate_learners,
         g_model      = cfg$superlearner$candidate_learners,
         censor_model = cfg$superlearner$survival_learners %||%
                          c("SL.glm", "SL.mean"))
  })
decisions    <- load_stage_output("stage3_decisions.rds")

# --- Load cleanTMLE lock and audit ---
lock_original <- tryCatch(load_stage_output("stage1_lock.rds"), error = function(e) NULL)
lock  <- tryCatch(load_stage_output("stage3_lock.rds"),
                  error = function(e) tryCatch(load_stage_output("stage2b_lock.rds"),
                    error = function(e2) load_stage_output("stage2_lock.rds")))
audit <- tryCatch(load_stage_output("stage3_audit.rds"),
                  error = function(e) tryCatch(load_stage_output("stage2b_audit.rds"),
                    error = function(e2) load_stage_output("stage2_audit.rds")))
ct_ps_fit <- tryCatch(load_stage_output("stage2_ct_ps_fit.rds"), error = function(e) NULL)

# --- Pre-outcome authorisation gate (modernised) ---
# 1) authorize_outcome_analysis() reads the audit trail and ensures all
#    required checkpoints recorded a GO (or, with allow_flag=TRUE, FLAG).
# 2) gate_all() also evaluates the composite GO/FLAG/STOP from the
#    individual checkpoint objects, providing a redundant safety check.
# 3) assert_outcome_authorized() is a hard runtime guard: any Stage 4
#    estimator that touches the outcome will refuse to run if the audit
#    has not signed off.
if (!is.null(audit)) {
  gate <- tryCatch(
    cleanTMLE::authorize_outcome_analysis(audit, allow_flag = TRUE),
    error = function(e) { cr_log(paste("authorize_outcome_analysis failed:",
                                        e$message)); NULL })
  if (!is.null(gate)) {
    cr_log(paste("Pre-outcome gate decision:", gate$decision))
    audit <- cleanTMLE::record_checkpoint(audit, gate)

    # Hard runtime guard
    tryCatch(cleanTMLE::assert_outcome_authorized(audit),
      error = function(e) {
        cr_log(paste("WARNING: assert_outcome_authorized:", e$message,
                     " — proceeding under override"))
      })

    # Unmask outcome only after gate authorises
    if (!is.null(lock_original)) {
      lock <- cleanTMLE::unmask_outcome(lock, lock_original)
    }
  }
}

A <- dat$A
Y <- dat$gose_good  # Binary: GOSE > 4
matched_idx <- match_result$matched_idx

cr_log(paste("Binary outcome (gose_good):",
             "n_available =", sum(!is.na(Y)),
             ", n_good =", sum(Y == 1, na.rm = TRUE),
             ", n_poor =", sum(Y == 0, na.rm = TRUE)))

decisions <- log_decision(decisions, "stage4",
                          paste("GOSE good defined as GOSE >", cfg$binary_outcome$threshold),
                          "Pre-specified threshold for functional recovery",
                          type = "pre-specified")

# --- Check outcome availability ---
if (all(is.na(Y))) {
  cr_log("ERROR: All binary outcomes are NA. Cannot proceed with Stage 4.")
  cr_log("Please ensure 'gose' variable exists in the dataset.")
  save_stage_output(decisions, "stage4_decisions.rds")
  stop("Binary outcome not available")
}

# --- Progress tracker for 3 estimators ---
pb4 <- cr_progress(3, "Stage 4")

# --- A. PS-matched logistic regression ---
cr_log("Running PS-matched logistic regression...")

result_matched_reg <- tryCatch(
  run_matched_regression(Y, A, W_matrix, matched_idx),
  error = function(e) {
    cr_log(paste("Matched regression failed:", e$message))
    NULL
  }
)
pb4$tick()

if (!is.null(result_matched_reg)) {
  cr_log(paste("Matched regression RD:",
               round(result_matched_reg$risk_difference$estimate, 4),
               "[", round(result_matched_reg$risk_difference$ci_lower, 4), ",",
               round(result_matched_reg$risk_difference$ci_upper, 4), "]"))
}

# --- B. PS-matched TMLE ---
cr_log("Running PS-matched TMLE...")

result_matched_tmle <- tryCatch(
  run_matched_tmle(Y, A, W_matrix, matched_idx, sel_libs$Q_model),
  error = function(e) {
    cr_log(paste("Matched TMLE failed:", e$message))
    NULL
  }
)
pb4$tick()

if (!is.null(result_matched_tmle)) {
  cr_log(paste("Matched TMLE ATE:",
               round(result_matched_tmle$ate$estimate, 4),
               "[", round(result_matched_tmle$ate$ci_lower, 4), ",",
               round(result_matched_tmle$ate$ci_upper, 4), "]"))
}

# --- C. Full-cohort TMLE (D in the spec) ---
cr_log("Running full-cohort TMLE...")

result_full_tmle <- tryCatch(
  run_full_cohort_tmle(Y, A, W_matrix, sel_libs$Q_model),
  error = function(e) {
    cr_log(paste("Full-cohort TMLE failed:", e$message))
    NULL
  }
)
pb4$tick()

if (!is.null(result_full_tmle)) {
  cr_log(paste("Full-cohort TMLE ATE:",
               round(result_full_tmle$ate$estimate, 4),
               "[", round(result_full_tmle$ate$ci_lower, 4), ",",
               round(result_full_tmle$ate$ci_upper, 4), "]"))
}

# --- D. IPCW-weighted TMLE (sensitivity for missing GOSE) ---
cr_log("Running IPCW-weighted TMLE (sensitivity for missing outcomes)...")

result_ipcw_tmle <- tryCatch(
  run_ipcw_tmle(Y, A, W_matrix, sel_libs$Q_model),
  error = function(e) {
    cr_log(paste("IPCW-TMLE failed:", e$message))
    NULL
  }
)

if (!is.null(result_ipcw_tmle)) {
  cr_log(paste("IPCW-TMLE ATE:",
               round(result_ipcw_tmle$ate$estimate, 4),
               "[", round(result_ipcw_tmle$ate$ci_lower, 4), ",",
               round(result_ipcw_tmle$ate$ci_upper, 4), "]"))
  cr_log(paste("  Missing GOSE:", result_ipcw_tmle$n_missing, "/",
               result_ipcw_tmle$n_total, "(", result_ipcw_tmle$pct_missing, "%)"))
  cr_log(paste("  IPCW weight summary: mean =", round(result_ipcw_tmle$ipcw_summary$mean, 3),
               ", max =", round(result_ipcw_tmle$ipcw_summary$max, 3)))
}

decisions <- log_decision(decisions, "stage4",
                          "Binary analysis: matched regression + matched TMLE + full TMLE + IPCW-TMLE",
                          "All three estimators run for comparison",
                          type = "final analysis")

# --- Combine results into summary ---
all_results <- list(result_matched_reg, result_matched_tmle, result_full_tmle, result_ipcw_tmle)
all_results <- all_results[!sapply(all_results, is.null)]

# Build comparison table
comparison <- do.call(rbind, lapply(all_results, function(r) {
  if (!is.null(r$risk_difference)) {
    data.frame(
      method    = r$method,
      estimate  = r$risk_difference$estimate,
      ci_lower  = r$risk_difference$ci_lower,
      ci_upper  = r$risk_difference$ci_upper,
      se        = r$risk_difference$se,
      risk_1    = r$risk_treated,
      risk_0    = r$risk_control,
      n         = r$n,
      stringsAsFactors = FALSE
    )
  } else if (!is.null(r$ate)) {
    data.frame(
      method    = r$method,
      estimate  = r$ate$estimate,
      ci_lower  = r$ate$ci_lower,
      ci_upper  = r$ate$ci_upper,
      se        = r$ate$se,
      risk_1    = r$risk_treated,
      risk_0    = r$risk_control,
      n         = r$n,
      stringsAsFactors = FALSE
    )
  }
}))

cr_log("=== Binary Outcome Results Summary ===")
for (i in seq_len(nrow(comparison))) {
  cr_log(paste(comparison$method[i], ":",
               "RD =", round(comparison$estimate[i], 4),
               "[", round(comparison$ci_lower[i], 4), ",",
               round(comparison$ci_upper[i], 4), "]",
               "n =", comparison$n[i]))
}

# --- TMLE diagnostics: influence curve & clever covariate ---
if (!is.null(result_full_tmle) && !is.null(result_full_tmle$fit)) {
  cr_log("Extracting TMLE diagnostics (IC, clever covariate)...")

  # Influence curve
  ic_values <- tryCatch(
    result_full_tmle$fit$estimates$IC$IC.ATE,
    error = function(e) NULL
  )
  if (!is.null(ic_values)) {
    # IC length matches complete cases used in TMLE
    cc_tmle <- complete.cases(Y, A)
    A_cc <- A[cc_tmle]
    plot_influence_curve(ic_values, A_cc,
                         file.path(cfg$paths$results, "ic_plot.png"))
    cr_log(paste("IC summary: mean =", round(mean(ic_values), 5),
                 ", SD =", round(sd(ic_values), 4)))
  } else {
    cr_log("IC extraction not available from TMLE fit")
  }

  # Clever covariate: H_n = A/g(1|W) - (1-A)/g(0|W)
  g1W <- tryCatch(result_full_tmle$fit$g$g1W, error = function(e) NULL)
  if (!is.null(g1W)) {
    cc_tmle <- complete.cases(Y, A)
    A_cc <- A[cc_tmle]
    g1W <- pmax(pmin(g1W, 0.999), 0.001)  # bound for numerical stability
    H_n <- A_cc / g1W - (1 - A_cc) / (1 - g1W)
    plot_clever_covariate(H_n, A_cc,
                          file.path(cfg$paths$results, "clever_covariate_plot.png"))
    cr_log(paste("Clever covariate: median(treated) =",
                 round(median(H_n[A_cc == 1]), 2),
                 ", median(control) =",
                 round(median(H_n[A_cc == 0]), 2)))
  } else {
    cr_log("g1W not available for clever covariate computation")
  }
}

# --- Estimator comparison plot ---
plot_estimator_comparison(comparison,
                           file.path(cfg$paths$results, "binary_estimator_comparison.png"))

# --- cleanTMLE modular TMLE (primary estimator) ---
if (!is.null(lock) && !is.null(ct_ps_fit)) {
  cr_log("Running cleanTMLE modular TMLE...")

  ct_tmle_fit <- tryCatch({
    g_fit    <- cleanTMLE::fit_tmle_treatment_mechanism(lock, ct_ps_fit)
    Q_fit    <- cleanTMLE::fit_tmle_outcome_mechanism(lock, g_fit)
    tmle_upd <- cleanTMLE::run_tmle_targeting_step(g_fit, Q_fit)
    cleanTMLE::extract_tmle_estimate(tmle_upd)
  }, error = function(e) {
    cr_log(paste("cleanTMLE modular TMLE failed:", e$message))
    NULL
  })

  if (!is.null(ct_tmle_fit)) {
    ate <- ct_tmle_fit$estimates$ATE
    cr_log(paste("cleanTMLE TMLE ATE:", round(ate$estimate, 4),
                 "[", round(ate$ci_lower, 4), ",", round(ate$ci_upper, 4), "]"))

    # IC histogram via cleanTMLE
    tryCatch({
      ic_p <- cleanTMLE::ic_histogram(ct_tmle_fit)
      ggsave(file.path(cfg$paths$results, "cleanTMLE_ic_histogram.png"),
             ic_p, width = 8, height = 5)
    }, error = function(e) cr_log(paste("ic_histogram failed:", e$message)))

    audit <- cleanTMLE::record_stage(audit, "Stage 4",
      sprintf("cleanTMLE TMLE ATE: %.4f", ate$estimate))
  }

  # cleanTMLE matching and IPTW for comparison
  ct_match <- tryCatch(
    cleanTMLE::run_match_workflow(lock, ct_ps_fit, override_clean_room = TRUE),
    error = function(e) { cr_log(paste("cleanTMLE matching failed:", e$message)); NULL }
  )
  ct_iptw <- tryCatch(
    cleanTMLE::run_iptw_workflow(lock, ct_ps_fit, override_clean_room = TRUE),
    error = function(e) { cr_log(paste("cleanTMLE IPTW failed:", e$message)); NULL }
  )

  # Cross-estimator comparison
  ct_fits <- list()
  if (!is.null(ct_match)) ct_fits$Matching <- ct_match
  if (!is.null(ct_iptw))  ct_fits$IPTW <- ct_iptw
  if (!is.null(ct_tmle_fit)) ct_fits$TMLE <- ct_tmle_fit

  if (length(ct_fits) > 0) {
    ct_comparison <- tryCatch(
      cleanTMLE::summarize_cleanroom_results(ct_fits),
      error = function(e) NULL
    )
    if (!is.null(ct_comparison)) {
      cr_log("cleanTMLE cross-estimator comparison:")
      print(ct_comparison)
      write.csv(ct_comparison,
                file.path(cfg$paths$results, "cleanTMLE_binary_comparison.csv"),
                row.names = FALSE)
    }
  }

  # ── Convenience wrapper: fit_final_workflows() runs matching + IPTW + TMLE
  cr_log("Running fit_final_workflows() — convenience wrapper...")
  ct_final <- tryCatch(
    cleanTMLE::fit_final_workflows(lock, ct_ps_fit,
                                    override_clean_room = TRUE),
    error = function(e) { cr_log(paste("fit_final_workflows failed:", e$message)); NULL })
  if (!is.null(ct_final))
    saveRDS(ct_final, file.path(cfg$paths$results, "stage4_fit_final_workflows.rds"))

  # ── Matched-cohort TMLE (no separate lock) via run_matched_tmle()
  cr_log("Running run_matched_tmle() on the matched cohort...")
  matched_idx_int <- as.integer(matched_idx)
  ct_matched_tmle <- tryCatch(
    cleanTMLE::run_matched_tmle(lock, ct_ps_fit,
                                 subset_idx = matched_idx_int,
                                 override_clean_room = TRUE),
    error = function(e) { cr_log(paste("run_matched_tmle failed:", e$message)); NULL })

  # ── Multi-candidate TMLE (sensitivity) via fit_tmle_candidate_set()
  cand_set <- tryCatch(load_stage_output("stage2b_candidates.rds"),
                       error = function(e) NULL)
  if (!is.null(cand_set) && length(cand_set) > 0) {
    ct_cand_set <- tryCatch(
      cleanTMLE::fit_tmle_candidate_set(lock, ct_ps_fit,
        candidates = cand_set,
        override_clean_room = TRUE),
      error = function(e) { cr_log(paste("fit_tmle_candidate_set failed:", e$message)); NULL })
    if (!is.null(ct_cand_set))
      saveRDS(ct_cand_set, file.path(cfg$paths$results,
              "stage4_tmle_candidate_set.rds"))
  }

  # ── TMLE diagnostics via cleanTMLE
  tryCatch({
    cc_p <- cleanTMLE::clever_covariate_plot(ps_fit = ct_ps_fit, lock = lock)
    ggsave(file.path(cfg$paths$results, "cleanTMLE_clever_covariate.png"),
           cc_p, width = 8, height = 5)
  }, error = function(e) cr_log(paste("clever_covariate_plot failed:", e$message)))

  # ── Forest plot of cross-estimator effects
  if (length(ct_fits) >= 2) {
    fp <- tryCatch(cleanTMLE::forest_plot(ct_fits),
                    error = function(e) NULL)
    if (!is.null(fp))
      ggsave(file.path(cfg$paths$results, "cleanTMLE_forest_plot.png"),
             fp, width = 8, height = 5)
  }

  # ── make_table2 — formal effect-summary table
  ct_table2 <- tryCatch(
    cleanTMLE::make_table2(ct_fits),
    error = function(e) NULL)
  if (!is.null(ct_table2))
    write.csv(ct_table2, file.path(cfg$paths$results,
              "cleanTMLE_table2.csv"), row.names = FALSE)

  # Truncation sensitivity (pre-registered in Stage 1a)
  tryCatch({
    sens <- cleanTMLE::sensitivity_truncation(lock,
      thresholds = c(0.01, 0.025, 0.05, 0.10),
      override_clean_room = TRUE)
    cr_log("Truncation sensitivity:")
    print(sens)
    write.csv(sens, file.path(cfg$paths$results, "truncation_sensitivity.csv"),
              row.names = FALSE)
  }, error = function(e) cr_log(paste("Sensitivity analysis failed:", e$message)))

  # E-value
  if (!is.null(ct_tmle_fit)) {
    tryCatch({
      ate_est <- ct_tmle_fit$estimates$ATE
      crude <- cleanTMLE::run_crude_workflow(lock, override_clean_room = TRUE)
      r0 <- crude$r0
      if (!is.null(r0) && r0 > 0) {
        rr_hat <- (r0 + ate_est$estimate) / r0
        rr_lo  <- (r0 + ate_est$ci_lower) / r0
        ev <- cleanTMLE::compute_evalue(rr_hat, ci_bound = min(rr_lo, 1/rr_lo))
        cr_log(sprintf("E-value: %.2f (CI bound: %.2f)", ev["e_value"], ev["e_value_ci"]))
      }
    }, error = function(e) cr_log(paste("E-value computation failed:", e$message)))
  }
}

# --- Save ---
save_stage_output(all_results, "stage4_binary_results.rds")
save_stage_output(comparison, "stage4_binary_comparison.rds")
save_stage_output(decisions, "stage4_decisions.rds")

write.csv(comparison, file.path(cfg$paths$results, "binary_outcome_comparison.csv"),
          row.names = FALSE)

if (!is.null(lock))  save_stage_output(lock, "stage4_lock.rds")
if (!is.null(audit)) save_stage_output(audit, "stage4_audit.rds")

# ── Final audit, decision log, and stage-path narrative ────────────────
## ─── Part 4: Manuscript artifacts directory ──────────────────────────
art_dir <- file.path(cfg$paths$results, "manuscript_artifacts")
dir.create(art_dir, showWarnings = FALSE, recursive = TRUE)

# table_effects.csv: one row per estimator
te <- if (!is.null(comparison)) {
  data.frame(
    method   = comparison$method,
    n        = comparison$n,
    estimate = round(comparison$estimate, 4),
    se       = round(comparison$se, 4),
    ci_lo    = round(comparison$ci_lower, 4),
    ci_hi    = round(comparison$ci_upper, 4),
    stringsAsFactors = FALSE)
} else NULL
if (!is.null(te))
  write.csv(te, file.path(art_dir, "table_effects.csv"), row.names = FALSE)

# table_dq_coverage.csv: from Stage 2b/3 dq_stress (if available)
dq_path <- file.path(cfg$paths$results, "plasmode_dq_degradation.csv")
if (file.exists(dq_path))
  file.copy(dq_path, file.path(art_dir, "table_dq_coverage.csv"),
            overwrite = TRUE)

# Figures: love plot, dq gradient, ic, clever
for (pair in list(c("love_plot.png", "fig_love_plot.png"),
                   c("plasmode_performance.png", "fig_dq_gradient.png"),
                   c("ic_plot.png", "fig_ic_plot.png"),
                   c("clever_covariate_plot.png", "fig_clever_covariate.png"))) {
  src <- file.path(cfg$paths$results, pair[1])
  if (file.exists(src))
    file.copy(src, file.path(art_dir, pair[2]), overwrite = TRUE)
}

# audit_summary.csv & decision_summary.csv: snapshot exports
if (!is.null(audit)) {
  trail0 <- tryCatch(cleanTMLE::export_audit_trail(audit), error = function(e) NULL)
  if (!is.null(trail0))
    write.csv(trail0, file.path(art_dir, "audit_summary.csv"),
              row.names = FALSE)
  dlog0 <- tryCatch(cleanTMLE::export_decision_log(audit), error = function(e) NULL)
  if (!is.null(dlog0))
    write.csv(dlog0, file.path(art_dir, "decision_summary.csv"),
              row.names = FALSE)
}

# case_study_metadata.json
if (!is.null(audit) && requireNamespace("jsonlite", quietly = TRUE)) {
  meta <- list(
    cleanTMLE_version = as.character(packageVersion("cleanTMLE")),
    rescueCo_data_cut = "2026-03-27",
    cohort = list(
      n_total      = nrow(dat),
      n_treated    = sum(dat$A == 1, na.rm = TRUE),
      n_control    = sum(dat$A == 0, na.rm = TRUE),
      n_gose_avail = sum(!is.na(Y))),
    primary_estimator = "Full-cohort TMLE (locked spec from Stage 2b)",
    primary_RD = if (!is.null(te))
      te$estimate[grepl("Full-cohort TMLE", te$method)][1] else NA_real_,
    primary_CI = if (!is.null(te))
      paste0("[", te$ci_lo[grepl("Full-cohort TMLE", te$method)][1], ", ",
             te$ci_hi[grepl("Full-cohort TMLE", te$method)][1], "]") else NA_character_,
    lock_hash = lock$lock_hash %||% NA_character_,
    rendered_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  jsonlite::write_json(meta, file.path(art_dir, "case_study_metadata.json"),
                        pretty = TRUE, auto_unbox = TRUE)
}

cr_log("Manuscript artifacts written to manuscript_artifacts/")

if (!is.null(audit)) {
  audit <- cleanTMLE::record_stage(audit, "Stage 4 (close)",
    "Final estimation, sensitivity, and diagnostics complete")
  trail <- tryCatch(cleanTMLE::export_audit_trail(audit), error = function(e) NULL)
  dlog  <- tryCatch(cleanTMLE::export_decision_log(audit), error = function(e) NULL)
  if (!is.null(trail))
    write.csv(trail, file.path(cfg$paths$results, "cleanTMLE_audit_trail.csv"),
              row.names = FALSE)
  if (!is.null(dlog))
    write.csv(dlog,  file.path(cfg$paths$results, "cleanTMLE_decision_log.csv"),
              row.names = FALSE)

  # Compact stage-path narrative for inclusion in reports
  narrative <- tryCatch(cleanTMLE::summarize_stage_path(audit),
                         error = function(e) NULL)
  if (!is.null(narrative))
    writeLines(as.character(narrative),
               file.path(cfg$paths$results, "cleanTMLE_stage_path.txt"))
}

## ─── Part 3 #2: Ordinal GOSE proportional-odds analysis ────────────────
cr_log("Running ordinal GOSE proportional-odds (sensitivity to binary)...")
ordinal_path <- file.path(cfg$paths$results, "gose_ordinal_po.csv")
tryCatch({
  if (!"gose_score" %in% names(dat))
    stop("dat$gose_score not present — cannot run PO model.")
  dpo <- data.frame(
    gose = factor(dat$gose_score, levels = 1:8, ordered = TRUE),
    A    = dat$A)
  ok <- !is.na(dpo$gose) & !is.na(dpo$A)
  if (sum(ok) < 100 || length(unique(dpo$gose[ok])) < 3)
    stop("Too few obs / levels for PO model")
  if (!requireNamespace("MASS", quietly = TRUE))
    install.packages("MASS", repos = "https://cloud.r-project.org", quiet = TRUE)
  fit_po <- MASS::polr(gose ~ A, data = dpo[ok, ], Hess = TRUE)
  s <- summary(fit_po)$coefficients
  beta <- s["A", "Value"]
  se   <- s["A", "Std. Error"]
  or   <- exp(beta); ci <- exp(beta + c(-1, 1) * 1.96 * se)
  p    <- 2 * pnorm(-abs(beta / se))

  # Test PO assumption
  po_test <- tryCatch({
    if (requireNamespace("ordinal", quietly = TRUE)) {
      cm <- ordinal::clm(gose ~ A, data = dpo[ok, ])
      nt <- ordinal::nominal_test(cm)
      list(LR = nt$LRT[2], p = nt$`Pr(>Chi)`[2])
    } else NULL
  }, error = function(e) NULL)

  out <- data.frame(
    estimator = c("Unadjusted PO (MASS::polr)"),
    common_OR = round(or, 3),
    ci_lo = round(ci[1], 3), ci_hi = round(ci[2], 3),
    p_value = round(p, 4),
    PO_test_LR = if (!is.null(po_test)) round(po_test$LR, 3) else NA,
    PO_test_p  = if (!is.null(po_test)) round(po_test$p, 4) else NA,
    n = sum(ok))
  write.csv(out, ordinal_path, row.names = FALSE)
  cr_log(sprintf("Ordinal GOSE: common OR = %.3f [%.3f, %.3f], p = %.4f",
                  or, ci[1], ci[2], p))
}, error = function(e) cr_log(paste("Ordinal GOSE PO skipped:", e$message)))

if (!is.null(audit)) {
  audit <- cleanTMLE::record_decision_log_entry(audit,
    stage = "Stage 4", decision_type = "outcome",
    description = "Reported binary GOSE>4 (primary) and ordinal-PO (sensitivity)",
    rationale = "Wang et al. 2023 — fixed dichotomy is least powerful spec")
}

cr_log("Stage 4 complete. GOSE outcomes unlocked and analyzed.")
