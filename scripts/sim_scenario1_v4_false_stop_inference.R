#!/usr/bin/env Rscript
# ===========================================================================
# Scenario 1 v4: False-Stop Inference with NP Bootstrap
# ===========================================================================
#
# PURPOSE:
#   Demonstrate that a clean-room overlap flag can be a genuine false stop
#   for TMLE *inferential* performance (not just point estimation). Uses
#   full nonparametric bootstrap (resample + refit nuisance + retarget) as
#   the primary TMLE inferential check, with multiplier bootstrap demoted
#   to a secondary diagnostic.
#
#   Two modes of outcome dependence on the overlap-driving covariate W4:
#
#     +-------------------+-----------+-----------+
#     |                   | Weak W4   | Strong W4 |
#     +-------------------+-----------+-----------+
#     | OK  overlap flag  | baseline  | baseline  |
#     | BAD overlap flag  | false stop| true stop |
#     +-------------------+-----------+-----------+
#
#   Key v4 changes from v3:
#     - Mode-specific overlap severity (s_bad_weak < s_bad_strong)
#     - Moderated nonlinear Q (reduced W1^2, |W3| coefficients)
#     - Reduced treatment heterogeneity in weak mode
#     - Full NP bootstrap for TMLE-ML as primary inference method
#     - Gate criterion based on NP bootstrap coverage
#
# CLEAN-ROOM STAGES:
#   Stage 1:  Pre-specification of estimand, estimators, diagnostics
#   Stage 2:  Outcome-blind overlap diagnostics (PS only; stop/go gate)
#   Stage 2b: Outcome-blind simulation-based OC → justify proceed / stop
#   Stage 3:  Estimation with blinded outcomes
#
# ESTIMATORS:
#   1. Crude RD
#   2. Adjusted regression (misspecified GLM g-computation, analytic SE)
#   3. IPTW RD (Hajek, analytic IF-based SE)
#   4. PS-matched regression (MatchIt NN, pairs bootstrap CI)
#   5. TMLE-GLM (misspecified, IC-based CI)
#   6. TMLE-ML (SL library, IC-based + full NP bootstrap CI)
#   7. Cross-fitted TMLE-ML (V-fold, IC-based CI)
#
# PRESETS:
#   "fast"  — N=400, reps=120, B_np=75   (development/tuning)
#   "eval"  — N=600, reps=300, B_np=150  (final evaluation)
#
# USAGE:
#   Rscript scripts/sim_scenario1_v4_false_stop_inference.R [fast|eval]
#
# OUTPUT:
#   results/scenario1_v4/{weak,strong}/ — truth, per-rep, summary, figures
#   results/scenario1_v4/design_parameters_used.json
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Source shared helpers
source(file.path("R", "sim_helpers_v4.R"))

# ===========================================================================
# Preset selection
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args) >= 1 && args[1] %in% c("fast", "eval")) args[1] else "fast"
message("Using preset: ", preset)

preset_config <- list(
  fast = list(n_main = 400, reps = 120, B_npboot_tmle = 0,
              B_matchboot = 100, B_multboot = 200, V_cf = 3),
  eval = list(n_main = 600, reps = 300, B_npboot_tmle = 150,
              B_matchboot = 200, B_multboot = 300, V_cf = 5)
)[[preset]]

# ===========================================================================
# Parameters
# ===========================================================================

params_base <- list(
  # Sample sizes (from preset)
  n_main    = preset_config$n_main,
  n_truth   = 100000,
  reps      = preset_config$reps,
  master_seed = 20260306,
  preset    = preset,

  # Bootstrap / cross-fitting (from preset)
  B_regboot     = 0,                          # analytic SE for adj regression
  B_iptw        = 0,                          # analytic SE for IPTW
  B_matchboot   = preset_config$B_matchboot,  # pairs bootstrap for PS-match
  B_npboot_tmle = preset_config$B_npboot_tmle,# full NP bootstrap for TMLE-ML
  B_multboot    = preset_config$B_multboot,   # Rademacher multiplier (secondary)
  V_cf          = preset_config$V_cf,         # folds for CV-TMLE

  # PS clipping
  ps_clip_true = c(0.005, 0.995),
  ps_clip_hat  = c(0.01, 0.99),

  # Overlap regime (replicate-level mixture)
  p_bad = 0.50,          # P(bad overlap replicate)
  s_ok  = 0.8,           # treatment-strength multiplier for OK overlap
  # s_bad is mode-specific, set below

  # Treatment model coefficients (before multiplier)
  alpha1_g = 1.5,        # W1
  alpha2_g = 0.8,        # W2
  alpha3_g = -0.5,       # W3
  alpha4_g = 3.0,        # W4 — overlap driver
  g_intercept_shift = -0.5,  # shift P(A=1) below 0.5 → more controls

  # Subject-level W4 mixture (off)
  use_subject_mixture = FALSE,
  sigma_lo  = 0.8,
  sigma_hi  = 2.0,
  p_hiVar   = 0.4,

  # Outcome model (v4 moderated nonlinear)
  # h1(W1) = h1_lin*W1 + h1_sq*W1^2
  # h2(W2) = h2_coef*W2
  # h3(W3) = h3_coef*abs(W3)
  # h4(W4) = beta4_Q*W4
  beta0_Q   = -0.5,      # intercept
  betaA_Q   = 0.40,      # treatment main effect (log-odds)
  h1_lin    = 0.50,       # W1 linear (v3 was 0.6)
  h1_sq     = 0.25,       # W1^2 (v3 was 0.4)
  h2_coef   = 0.70,       # W2
  h3_coef   = 0.35,       # |W3| (v3 was 0.5)
  # beta4_Q, tau_int1, tau_int2, tau_int3 are mode-specific

  # Overlap flag thresholds (Stage 2)
  overlap_flag_thresholds = list(
    extreme_ps   = 0.05,
    extreme_prop = 0.30,    # v3 was 0.35 — slightly tighter
    ess_frac     = 0.30,    # v3 was 0.25 — slightly tighter
    max_w        = 40       # v3 was 50 — slightly tighter
  ),

  # Matching
  match_caliper = NULL,   # NULL → 0.2*sd(logit(ps))

  # SuperLearner libraries for TMLE-ML
  Q_SL_library_ml = c("SL.mean", "SL.glmnet", "SL.xgboost"),
  g_SL_library_ml = c("SL.mean", "SL.glmnet", "SL.xgboost"),

  # Output
  output_dir    = "results/scenario1_v4",
  save_per_rep  = TRUE
)

# Two modes: weak vs strong W4 outcome dependence
# Mode-specific overlap severity AND outcome coefficients
mode_params <- list(
  weak = list(
    s_bad     = 1.50,      # moderate overlap violation (v3 used 2.0 for both)
    beta4_Q   = 0.04,      # W4 barely affects Y
    tau_int1  = 0.15,      # reduced A*I(W1>0) (v3 was 0.25)
    tau_int2  = 0.05,      # reduced A*W2 (v3 was 0.10)
    tau_int3  = 0.0        # no A*W4 interaction
  ),
  strong = list(
    s_bad     = 2.0,       # severe overlap violation
    beta4_Q   = 0.35,      # W4 meaningfully affects Y
    tau_int1  = 0.25,      # A*I(W1>0)
    tau_int2  = 0.10,      # A*W2
    tau_int3  = 0.20       # A*W4 interaction present
  )
)


# ===========================================================================
# Run one replicate
# ===========================================================================

run_one_rep_v4 <- function(r, p, truth_rd) {
  rep_seed <- p$master_seed + r

  # Generate data
  set.seed(rep_seed)
  sim <- sim_data_v4(n = p$n_main, seed = rep_seed, p = p)
  dat <- sim$data
  W   <- dat[, c("W1", "W2", "W3", "W4"), drop = FALSE]
  A   <- dat$A; Y <- dat$Y; n <- nrow(dat)

  # ------------------------------------------------------------------
  # Stage 2: Outcome-blind diagnostics
  # ------------------------------------------------------------------
  ps  <- estimate_ps_v4(dat, p)
  diag <- overlap_diagnostics_v4(ps$g_raw, ps$g_bounded, A, p)

  # ------------------------------------------------------------------
  # Stage 3: Estimation
  # ------------------------------------------------------------------

  # (1) Crude
  set.seed(rep_seed + 100000L)
  crude <- est_crude(Y, A)

  # (2) Adjusted regression (misspecified, analytic SE)
  set.seed(rep_seed + 200000L)
  regadj <- est_reg_adj_v4(dat, B_boot = p$B_regboot)

  # (3) IPTW (analytic IF-based SE)
  set.seed(rep_seed + 300000L)
  iptw <- est_iptw_v4(Y, A, ps$g_bounded, dat, B_boot = p$B_iptw, p = p)

  # (4) PS-matched regression
  set.seed(rep_seed + 400000L)
  psmatch <- est_psmatch_reg_v4(dat, ps$g_bounded, B_boot = p$B_matchboot, p = p)

  # (5) TMLE-GLM (IC-based only)
  set.seed(rep_seed + 500000L)
  tmle_glm <- est_tmle_glm(Y, A, W, p)

  # (6) TMLE-ML with full NP bootstrap (PRIMARY inference)
  set.seed(rep_seed + 600000L)
  tmle_ml_np <- est_tmle_ml_npboot(Y, A, W, p, B = p$B_npboot_tmle)

  # (7) Cross-fitted TMLE-ML (IC-based)
  set.seed(rep_seed + 700000L)
  cvtmle <- est_cvtmle_ml(dat, p, V = p$V_cf)

  # ------------------------------------------------------------------
  # Multiplier bootstrap (secondary diagnostic) for TMLE variants
  # ------------------------------------------------------------------
  mb_glm <- mb_ml <- mb_cv <- list(se = NA, ci_low = NA, ci_high = NA,
                                    ci_low_pct = NA, ci_high_pct = NA,
                                    success = FALSE)
  if (tmle_glm$success && p$B_multboot > 0) {
    set.seed(rep_seed + 800000L)
    mb_glm <- multiplier_bootstrap(tmle_glm$estimate, tmle_glm$eif, p$B_multboot)
  }
  if (tmle_ml_np$success && !is.null(tmle_ml_np$eif) && p$B_multboot > 0) {
    set.seed(rep_seed + 810000L)
    mb_ml <- multiplier_bootstrap(tmle_ml_np$estimate, tmle_ml_np$eif, p$B_multboot)
  }
  if (cvtmle$success && p$B_multboot > 0) {
    set.seed(rep_seed + 820000L)
    mb_cv <- multiplier_bootstrap(cvtmle$estimate, cvtmle$eif, p$B_multboot)
  }

  # ------------------------------------------------------------------
  # Collect results
  # ------------------------------------------------------------------
  data.frame(
    rep = r, seed = rep_seed, n = n,
    overlap_regime_true = sim$overlap_regime,
    s_g_used = sim$s_g,

    # Stage 2 diagnostics
    fraction_extreme = diag$fraction_extreme,
    min_g_raw = diag$min_g_raw, max_g_raw = diag$max_g_raw,
    ess_frac = diag$ess_frac, max_w = diag$max_w,
    overlap_flag = diag$overlap_flag,
    event_rate = mean(Y), n_treated = sum(A), n_control = sum(1 - A),

    # Crude
    crude_est = crude$estimate, crude_se = crude$se,
    crude_ci_lo = crude$ci_low, crude_ci_hi = crude$ci_high,
    crude_success = crude$success,

    # Adjusted regression
    regadj_est = regadj$estimate, regadj_se = regadj$se,
    regadj_ci_lo = regadj$ci_low, regadj_ci_hi = regadj$ci_high,
    regadj_success = regadj$success,

    # IPTW
    iptw_est = iptw$estimate, iptw_se = iptw$se,
    iptw_ci_lo = iptw$ci_low, iptw_ci_hi = iptw$ci_high,
    iptw_success = iptw$success,

    # PS-matched regression
    psmatch_est = psmatch$estimate, psmatch_se = psmatch$se,
    psmatch_ci_lo = psmatch$ci_low, psmatch_ci_hi = psmatch$ci_high,
    psmatch_success = psmatch$success,
    psmatch_n_matched = if (!is.null(psmatch$n_matched)) psmatch$n_matched else NA,

    # TMLE-GLM (IC-based)
    tmle_glm_est = tmle_glm$estimate, tmle_glm_se_ic = tmle_glm$se,
    tmle_glm_ci_ic_lo = tmle_glm$ci_low, tmle_glm_ci_ic_hi = tmle_glm$ci_high,
    tmle_glm_success = tmle_glm$success,
    # TMLE-GLM multiplier bootstrap (secondary)
    tmle_glm_se_mb = mb_glm$se,
    tmle_glm_ci_mb_lo = mb_glm$ci_low, tmle_glm_ci_mb_hi = mb_glm$ci_high,

    # TMLE-ML IC-based
    tmle_ml_est = tmle_ml_np$estimate,
    tmle_ml_se_ic = tmle_ml_np$se_ic,
    tmle_ml_ci_ic_lo = tmle_ml_np$ci_ic_lo,
    tmle_ml_ci_ic_hi = tmle_ml_np$ci_ic_hi,
    tmle_ml_success = tmle_ml_np$success,
    # TMLE-ML NP bootstrap (PRIMARY)
    tmle_ml_se_np = tmle_ml_np$se_np,
    tmle_ml_ci_np_lo = tmle_ml_np$ci_np_lo,
    tmle_ml_ci_np_hi = tmle_ml_np$ci_np_hi,
    tmle_ml_ci_np_norm_lo = tmle_ml_np$ci_np_norm_lo,
    tmle_ml_ci_np_norm_hi = tmle_ml_np$ci_np_norm_hi,
    tmle_ml_n_boot_ok = tmle_ml_np$n_boot_ok,
    # TMLE-ML multiplier bootstrap (secondary)
    tmle_ml_se_mb = mb_ml$se,
    tmle_ml_ci_mb_lo = mb_ml$ci_low, tmle_ml_ci_mb_hi = mb_ml$ci_high,

    # CV-TMLE-ML (IC-based)
    cvtmle_est = cvtmle$estimate, cvtmle_se_ic = cvtmle$se,
    cvtmle_ci_ic_lo = cvtmle$ci_low, cvtmle_ci_ic_hi = cvtmle$ci_high,
    cvtmle_success = cvtmle$success,
    # CV-TMLE multiplier bootstrap (secondary)
    cvtmle_se_mb = mb_cv$se,
    cvtmle_ci_mb_lo = mb_cv$ci_low, cvtmle_ci_mb_hi = mb_cv$ci_high,

    stringsAsFactors = FALSE
  )
}


# ===========================================================================
# Run simulation for one mode
# ===========================================================================

run_mode <- function(mode_name, params_base, mode_params) {
  message("\n", paste(rep("=", 70), collapse = ""))
  message("MODE: ", toupper(mode_name), " W4 outcome dependence")
  message(paste(rep("=", 70), collapse = ""))

  p <- params_base
  mp <- mode_params[[mode_name]]
  p$s_bad     <- mp$s_bad
  p$beta4_Q   <- mp$beta4_Q
  p$tau_int1  <- mp$tau_int1
  p$tau_int2  <- mp$tau_int2
  p$tau_int3  <- mp$tau_int3
  mode_dir <- file.path(p$output_dir, mode_name)
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)
  if (p$save_per_rep)
    dir.create(file.path(mode_dir, "per_rep"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(mode_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  # --- Truth ---
  message("Computing truth (depends on Q0 and W distribution, not g0) ...")
  truth <- compute_truth_v4(p)
  saveRDS(truth, file.path(mode_dir, "truth.rds"))
  message(sprintf("  RD_true = %.6f  (EY1=%.4f, EY0=%.4f)",
                  truth$RD_true, truth$EY1, truth$EY0))

  # --- Replicate loop ---
  message(sprintf("Running %d replicates (N=%d, s_bad=%.2f, beta4=%.2f, tau_int3=%.2f, B_np=%d) ...",
                  p$reps, p$n_main, p$s_bad, p$beta4_Q, p$tau_int3, p$B_npboot_tmle))

  t0 <- Sys.time()
  results_list <- vector("list", p$reps)
  for (r in seq_len(p$reps)) {
    t_rep <- Sys.time()
    results_list[[r]] <- run_one_rep_v4(r, p, truth$RD_true)

    if (p$save_per_rep)
      saveRDS(results_list[[r]],
              file.path(mode_dir, "per_rep", sprintf("rep_%03d.rds", r)))

    elapsed_rep <- as.numeric(difftime(Sys.time(), t_rep, units = "secs"))
    if (r %% 10 == 0 || r == 1) {
      rr <- results_list[[r]]
      elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf("  Rep %d/%d  [%.1fs]  regime=%s flag=%s  TMLE-ML=%.4f  NP.boot.ok=%d  (%.1f min elapsed)",
                      r, p$reps, elapsed_rep,
                      rr$overlap_regime_true, rr$overlap_flag,
                      ifelse(is.na(rr$tmle_ml_est), NA, rr$tmle_ml_est),
                      ifelse(is.na(rr$tmle_ml_n_boot_ok), 0, rr$tmle_ml_n_boot_ok),
                      elapsed_total))
    }
  }

  res <- do.call(rbind, results_list)

  # --- Metrics ---
  build_metrics <- function(df, label) {
    rbind(
      compute_metrics_v4(df$crude_est, df$crude_ci_lo, df$crude_ci_hi,
                         df$crude_success, truth$RD_true, df$crude_se,
                         label, "Crude"),
      compute_metrics_v4(df$regadj_est, df$regadj_ci_lo, df$regadj_ci_hi,
                         df$regadj_success, truth$RD_true, df$regadj_se,
                         label, "Adj. Regression (analytic)"),
      compute_metrics_v4(df$iptw_est, df$iptw_ci_lo, df$iptw_ci_hi,
                         df$iptw_success, truth$RD_true, df$iptw_se,
                         label, "IPTW (analytic)"),
      compute_metrics_v4(df$psmatch_est, df$psmatch_ci_lo, df$psmatch_ci_hi,
                         df$psmatch_success, truth$RD_true, df$psmatch_se,
                         label, "PS-Matched (pairs boot)"),
      compute_metrics_v4(df$tmle_glm_est, df$tmle_glm_ci_ic_lo, df$tmle_glm_ci_ic_hi,
                         df$tmle_glm_success, truth$RD_true, df$tmle_glm_se_ic,
                         label, "TMLE-GLM (IC)"),
      compute_metrics_v4(df$tmle_ml_est, df$tmle_ml_ci_ic_lo, df$tmle_ml_ci_ic_hi,
                         df$tmle_ml_success, truth$RD_true, df$tmle_ml_se_ic,
                         label, "TMLE-ML (IC)"),
      compute_metrics_v4(df$tmle_ml_est, df$tmle_ml_ci_np_lo, df$tmle_ml_ci_np_hi,
                         df$tmle_ml_success, truth$RD_true, df$tmle_ml_se_np,
                         label, "TMLE-ML (NP boot pctl)"),
      compute_metrics_v4(df$tmle_ml_est, df$tmle_ml_ci_np_norm_lo, df$tmle_ml_ci_np_norm_hi,
                         df$tmle_ml_success, truth$RD_true, df$tmle_ml_se_np,
                         label, "TMLE-ML (NP boot norm)"),
      compute_metrics_v4(df$tmle_ml_est, df$tmle_ml_ci_mb_lo, df$tmle_ml_ci_mb_hi,
                         df$tmle_ml_success, truth$RD_true, df$tmle_ml_se_mb,
                         label, "TMLE-ML (mult.boot)"),
      compute_metrics_v4(df$cvtmle_est, df$cvtmle_ci_ic_lo, df$cvtmle_ci_ic_hi,
                         df$cvtmle_success, truth$RD_true, df$cvtmle_se_ic,
                         label, "CV-TMLE-ML (IC)"),
      compute_metrics_v4(df$cvtmle_est, df$cvtmle_ci_mb_lo, df$cvtmle_ci_mb_hi,
                         df$cvtmle_success, truth$RD_true, df$cvtmle_se_mb,
                         label, "CV-TMLE-ML (mult.boot)")
    )
  }

  metrics_all <- build_metrics(res, "all")
  res_ok  <- res[res$overlap_flag == "ok", , drop = FALSE]
  res_bad <- res[res$overlap_flag == "bad", , drop = FALSE]
  metrics_ok  <- if (nrow(res_ok)  >= 5) build_metrics(res_ok,  "overlap_ok")  else NULL
  metrics_bad <- if (nrow(res_bad) >= 5) build_metrics(res_bad, "overlap_bad") else NULL

  summary_overall <- rbind(metrics_all, metrics_ok, metrics_bad)

  # --- Stop/go gate (D3) — based on NP bootstrap coverage ---
  gate <- data.frame(overlap_flag = c("ok", "bad"), stringsAsFactors = FALSE)
  for (fl in c("ok", "bad")) {
    sub <- res[res$overlap_flag == fl, , drop = FALSE]
    if (nrow(sub) < 5) {
      gate$n_reps[gate$overlap_flag == fl] <- nrow(sub)
      gate[gate$overlap_flag == fl, c("tmle_ml_bias", "tmle_ml_cov_np",
                                       "tmle_ml_cov_ic", "tmle_ml_rmse",
                                       "acceptable")] <- NA
      next
    }
    ok_ml <- sub[sub$tmle_ml_success == TRUE, , drop = FALSE]
    if (nrow(ok_ml) < 3) {
      gate$n_reps[gate$overlap_flag == fl] <- nrow(sub)
      gate[gate$overlap_flag == fl, c("tmle_ml_bias", "tmle_ml_cov_np",
                                       "tmle_ml_cov_ic", "tmle_ml_rmse",
                                       "acceptable")] <- NA
      next
    }
    bias_ml <- mean(ok_ml$tmle_ml_est, na.rm = TRUE) - truth$RD_true
    rmse_ml <- sqrt(mean((ok_ml$tmle_ml_est - truth$RD_true)^2, na.rm = TRUE))

    # NP bootstrap coverage (primary gate criterion)
    has_np <- !is.na(ok_ml$tmle_ml_ci_np_lo) & !is.na(ok_ml$tmle_ml_ci_np_hi)
    cov_np <- if (sum(has_np) >= 5)
      mean((ok_ml$tmle_ml_ci_np_lo[has_np] <= truth$RD_true) &
           (truth$RD_true <= ok_ml$tmle_ml_ci_np_hi[has_np])) else NA

    # IC coverage (for comparison)
    has_ic <- !is.na(ok_ml$tmle_ml_ci_ic_lo) & !is.na(ok_ml$tmle_ml_ci_ic_hi)
    cov_ic <- if (sum(has_ic) >= 5)
      mean((ok_ml$tmle_ml_ci_ic_lo[has_ic] <= truth$RD_true) &
           (truth$RD_true <= ok_ml$tmle_ml_ci_ic_hi[has_ic])) else NA

    # Gate: use NP bootstrap if available, otherwise fall back to IC
    cov_gate <- if (!is.na(cov_np)) cov_np else cov_ic
    accept <- !is.na(cov_gate) && (abs(bias_ml) < 0.01) && (cov_gate >= 0.93)

    gate$n_reps[gate$overlap_flag == fl] <- nrow(sub)
    gate$tmle_ml_bias[gate$overlap_flag == fl] <- round(bias_ml, 5)
    gate$tmle_ml_cov_np[gate$overlap_flag == fl] <- round(cov_np, 3)
    gate$tmle_ml_cov_ic[gate$overlap_flag == fl] <- round(cov_ic, 3)
    gate$tmle_ml_rmse[gate$overlap_flag == fl] <- round(rmse_ml, 5)
    gate$acceptable[gate$overlap_flag == fl] <- accept
  }

  # --- Print results ---
  message("\n--- ", toupper(mode_name), " W4 mode results ---")
  message(sprintf("True RD: %.6f", truth$RD_true))
  message(sprintf("Overlap flags: %d ok, %d bad (%.0f%% bad)",
                  nrow(res_ok), nrow(res_bad),
                  100 * nrow(res_bad) / nrow(res)))
  message("\nOverall metrics:")
  print(metrics_all)
  message("\nStop/go gate (NP bootstrap-based):")
  print(gate)

  # --- Save ---
  saveRDS(p, file.path(mode_dir, "params_used.rds"))
  saveRDS(res, file.path(mode_dir, "replicate_results.rds"))
  write.csv(summary_overall, file.path(mode_dir, "summary_overall.csv"),
            row.names = FALSE)
  if (!is.null(metrics_ok))
    write.csv(metrics_ok, file.path(mode_dir, "summary_overlap_ok.csv"),
              row.names = FALSE)
  if (!is.null(metrics_bad))
    write.csv(metrics_bad, file.path(mode_dir, "summary_overlap_bad.csv"),
              row.names = FALSE)
  saveRDS(gate, file.path(mode_dir, "gate_results.rds"))

  # --- Plots ---
  fig_dir <- file.path(mode_dir, "figures")

  # PS distribution by overlap flag
  tryCatch({
    ps_df <- data.frame(
      g_hat = unlist(lapply(seq_len(min(10, p$reps)), function(r) {
        sim <- sim_data_v4(p$n_main, seed = p$master_seed + r, p = p)
        estimate_ps_v4(sim$data, p)$g_raw
      })),
      flag = rep(
        unlist(lapply(seq_len(min(10, p$reps)), function(r) {
          sim <- sim_data_v4(p$n_main, seed = p$master_seed + r, p = p)
          ps <- estimate_ps_v4(sim$data, p)
          diag <- overlap_diagnostics_v4(ps$g_raw, ps$g_bounded, sim$data$A, p)
          rep(diag$overlap_flag, p$n_main)
        })), each = 1)
    )
    p_ps <- ggplot(ps_df, aes(x = g_hat, fill = flag)) +
      geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
      scale_fill_manual(values = c(ok = "steelblue", bad = "firebrick")) +
      labs(title = paste0("PS distribution (", mode_name, " W4 mode)"),
           x = "Estimated PS", y = "Count", fill = "Overlap flag") +
      theme_minimal(base_size = 12)
    ggsave(file.path(fig_dir, "ps_distribution.png"), p_ps, width = 8, height = 5)
  }, error = function(e) message("  Could not create PS plot: ", e$message))

  # Sampling distributions
  tryCatch({
    methods <- c("Crude", "Adj.Reg", "IPTW", "PS-Match",
                 "TMLE-GLM", "TMLE-ML", "CV-TMLE")
    ests <- list(res$crude_est, res$regadj_est, res$iptw_est,
                 res$psmatch_est, res$tmle_glm_est, res$tmle_ml_est,
                 res$cvtmle_est)
    est_long <- do.call(rbind, lapply(seq_along(methods), function(i) {
      data.frame(method = methods[i], estimate = ests[[i]],
                 overlap_flag = res$overlap_flag, stringsAsFactors = FALSE)
    }))
    est_long <- est_long[!is.na(est_long$estimate), ]
    est_long$method <- factor(est_long$method, levels = methods)

    p_dist <- ggplot(est_long, aes(x = estimate, fill = overlap_flag)) +
      geom_histogram(bins = 35, alpha = 0.6, position = "identity") +
      geom_vline(xintercept = truth$RD_true, linewidth = 0.8) +
      facet_wrap(~ method, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = c(ok = "steelblue", bad = "firebrick")) +
      labs(title = paste0("Sampling distributions (", mode_name, " W4)"),
           subtitle = paste0("True RD = ", round(truth$RD_true, 4)),
           x = "Estimated RD", y = "Count", fill = "Overlap") +
      theme_minimal(base_size = 11)
    ggsave(file.path(fig_dir, "sampling_distributions.png"), p_dist,
           width = 9, height = 14)
  }, error = function(e) message("  Could not create sampling dist plot: ", e$message))

  # Bias by estimator and overlap
  tryCatch({
    bias_df <- summary_overall[summary_overall$subset %in% c("overlap_ok", "overlap_bad") &
                               !grepl("mult\\.boot|NP boot", summary_overall$method), ]
    if (nrow(bias_df) > 0) {
      bias_df$method <- factor(bias_df$method)
      p_bias <- ggplot(bias_df, aes(x = method, y = bias, fill = subset)) +
        geom_col(position = "dodge", alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        coord_flip() +
        scale_fill_manual(values = c(overlap_ok = "steelblue",
                                     overlap_bad = "firebrick")) +
        labs(title = paste0("Bias by estimator and overlap (", mode_name, " W4)"),
             x = NULL, y = "Bias", fill = "Overlap") +
        theme_minimal(base_size = 12)
      ggsave(file.path(fig_dir, "bias_by_overlap.png"), p_bias, width = 9, height = 6)
    }
  }, error = function(e) message("  Could not create bias plot: ", e$message))

  # Coverage comparison (all inference methods)
  tryCatch({
    cov_df <- summary_overall[summary_overall$subset == "all" &
                              !is.na(summary_overall$coverage), ]
    if (nrow(cov_df) > 0) {
      cov_df$method <- factor(cov_df$method, levels = rev(cov_df$method))
      cov_df$ci_lo <- pmax(cov_df$coverage - 1.96 * cov_df$mcse_coverage, 0)
      cov_df$ci_hi <- pmin(cov_df$coverage + 1.96 * cov_df$mcse_coverage, 1)
      p_cov <- ggplot(cov_df, aes(x = method, y = coverage, fill = method)) +
        geom_col(alpha = 0.8, width = 0.6) +
        geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
        coord_flip() + ylim(0, 1) +
        labs(title = paste0("Coverage by inference method (", mode_name, " W4)"),
             x = NULL, y = "Coverage") +
        guides(fill = "none") +
        theme_minimal(base_size = 11)
      ggsave(file.path(fig_dir, "coverage_overall.png"), p_cov, width = 10, height = 7)
    }
  }, error = function(e) message("  Could not create coverage plot: ", e$message))

  # NP bootstrap coverage by overlap stratum
  tryCatch({
    np_rows <- summary_overall[grepl("NP boot pctl", summary_overall$method) &
                               summary_overall$subset %in% c("overlap_ok", "overlap_bad") &
                               !is.na(summary_overall$coverage), ]
    if (nrow(np_rows) > 0) {
      np_rows$ci_lo <- pmax(np_rows$coverage - 1.96 * np_rows$mcse_coverage, 0)
      np_rows$ci_hi <- pmin(np_rows$coverage + 1.96 * np_rows$mcse_coverage, 1)
      p_np <- ggplot(np_rows, aes(x = subset, y = coverage, fill = subset)) +
        geom_col(alpha = 0.8, width = 0.5) +
        geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15) +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
        ylim(0, 1) +
        scale_fill_manual(values = c(overlap_ok = "steelblue",
                                     overlap_bad = "firebrick")) +
        labs(title = paste0("TMLE-ML NP bootstrap coverage by overlap (", mode_name, " W4)"),
             x = "Overlap stratum", y = "Coverage") +
        guides(fill = "none") +
        theme_minimal(base_size = 12)
      ggsave(file.path(fig_dir, "npboot_coverage_by_overlap.png"), p_np,
             width = 6, height = 5)
    }
  }, error = function(e) message("  Could not create NP boot coverage plot: ", e$message))

  # Session info
  si <- capture.output(sessionInfo())
  writeLines(si, file.path(mode_dir, "sessionInfo.txt"))

  list(truth = truth, results = res, metrics = summary_overall,
       gate = gate, params = p)
}


# ===========================================================================
# Main execution
# ===========================================================================

message("=== Scenario 1 v4: False-Stop Inference (NP Bootstrap) ===")
message("Preset: ", preset)
message("Start: ", Sys.time())
message("ranger available: ", .has_ranger)
message("MatchIt available: ", .has_MatchIt)

dir.create(params_base$output_dir, recursive = TRUE, showWarnings = FALSE)

out_weak   <- run_mode("weak",   params_base, mode_params)
out_strong <- run_mode("strong", params_base, mode_params)


# ===========================================================================
# Cross-mode comparison summary
# ===========================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("CROSS-MODE COMPARISON")
message(paste(rep("=", 70), collapse = ""))

message("\n--- WEAK W4 (false-stop demonstration) ---")
message("  Gate (NP bootstrap-based):")
print(out_weak$gate)

message("\n--- STRONG W4 (true-stop demonstration) ---")
message("  Gate (NP bootstrap-based):")
print(out_strong$gate)

# Combined gate table
gate_combined <- rbind(
  cbind(mode = "weak",   out_weak$gate),
  cbind(mode = "strong", out_strong$gate)
)
write.csv(gate_combined,
          file.path(params_base$output_dir, "gate_combined.csv"),
          row.names = FALSE)
saveRDS(gate_combined,
        file.path(params_base$output_dir, "gate_combined.rds"))

# Design parameters JSON
tryCatch({
  design_json <- list(
    version = "v4",
    preset = preset,
    date = as.character(Sys.time()),
    base_params = list(
      n_main = params_base$n_main,
      n_truth = params_base$n_truth,
      reps = params_base$reps,
      B_npboot_tmle = params_base$B_npboot_tmle,
      B_multboot = params_base$B_multboot,
      B_matchboot = params_base$B_matchboot,
      V_cf = params_base$V_cf,
      h1_lin = params_base$h1_lin,
      h1_sq = params_base$h1_sq,
      h3_coef = params_base$h3_coef,
      g_intercept_shift = params_base$g_intercept_shift,
      alpha4_g = params_base$alpha4_g,
      s_ok = params_base$s_ok,
      ps_clip_hat = params_base$ps_clip_hat
    ),
    mode_weak = mode_params$weak,
    mode_strong = mode_params$strong,
    overlap_thresholds = params_base$overlap_flag_thresholds
  )
  writeLines(jsonlite::toJSON(design_json, auto_unbox = TRUE, pretty = TRUE),
             file.path(params_base$output_dir, "design_parameters_used.json"))
}, error = function(e) {
  # Fallback if jsonlite not available
  cat(
    "version: v4\n",
    "preset:", preset, "\n",
    "n_main:", params_base$n_main, "\n",
    "reps:", params_base$reps, "\n",
    "B_npboot_tmle:", params_base$B_npboot_tmle, "\n",
    "s_bad_weak:", mode_params$weak$s_bad, "\n",
    "s_bad_strong:", mode_params$strong$s_bad, "\n",
    file = file.path(params_base$output_dir, "design_parameters_used.txt"),
    sep = ""
  )
})

message("\nOutputs saved to: ", params_base$output_dir)
message("Done: ", Sys.time())
