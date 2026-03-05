#!/usr/bin/env Rscript
# ===========================================================================
# sim_scenario1_bad_overlap_tmle_ok.R
#
# Scenario 1: Poor propensity score overlap (practical near-positivity
# violation) but TMLE still achieves low bias and nominal coverage because
# the outcome regression Q(A,W) is well-learned, while a conventional
# unadjusted analysis and IPTW fail.
#
# This script is standalone and does NOT depend on the main repo pipeline.
# It implements a "clean-room style" three-stage flow internally:
#   Stage 1 -- Pre-specification (estimand, estimators, thresholds)
#   Stage 2 -- Outcome-blind diagnostics (PS overlap, ESS, weights)
#   Stage 3 -- Estimation + simulation validation (bias, RMSE, coverage)
#
# Usage:
#   Rscript scripts/sim_scenario1_bad_overlap_tmle_ok.R
#
# Tunable parameters are defined immediately below. Adjust N, reps,
# strength_overlap, truth_N, etc. before running.
#
# ===========================================================================
# TUNING GUIDE
# ===========================================================================
# The DGP separates "what drives treatment" from "what drives outcome".
# W4 is a strong predictor of treatment but a WEAK predictor of outcome,
# while W1 moderately predicts both. This is the key structural feature
# that lets TMLE succeed: the outcome model Q(A,W) is smooth and learnable
# even in regions where the PS is extreme, because PS extremity is mostly
# driven by W4 (which barely affects Y).
#
# To make overlap WORSE (more extreme PS):
#   - Increase strength_overlap (try 1.5 -> 2.0 -> 2.5)
#   - Decrease eps_clip (try 0.005 -> 0.002 -> 0.001)
#
# To make IPTW perform WORSE relative to TMLE:
#   - Increase strength_overlap (larger weights, more instability)
#   - Decrease N (less data to stabilize weighted estimator)
#
# To make TMLE coverage CLOSER to 0.95:
#   - Increase N (IC-based SE more accurate in larger samples)
#   - Increase coef_W4_Q toward 0 (less outcome variation at PS boundary)
#   - Ensure Q.SL.library is rich enough to learn Q well
#
# To make crude bias LARGER:
#   - Increase coef_W1_Q and/or coef_W1_g (stronger shared confounding)
#
# For PRECISION of coverage estimates:
#   - reps=50 gives ~+-6pp MC error; reps=200 gives ~+-3pp; reps=500 ~+-2pp
# ===========================================================================

# ---- Tunable parameters (edit these) ---- #
PARAMS <- list(
  # Simulation size
  N             = 2000,     # sample size per replicate
  reps          = 200,      # MC reps (200+ recommended for stable coverage)
  truth_N       = 100000,   # sample size for ground truth (exact plug-in)
  master_seed   = 20260304, # master seed for reproducibility

  # DGP: treatment mechanism (controls overlap)
  strength_overlap = 1.5,   # global multiplier on g-model coefficients
  coef_W1_g     = 1.5,      # W1 effect on treatment (shared confounder)
  coef_W2_g     = 0.8,      # W2 effect on treatment
  coef_W3_g     = -0.5,     # W3 effect on treatment
  coef_W4_g     = 3.5,      # W4 effect on treatment (overlap driver, weak in Q)
  eps_clip      = 0.005,    # clip g_true to [eps, 1-eps]

  # DGP: outcome model (controls what Q must learn)
  effect_A      = 0.5,      # treatment effect on log-odds scale
  alpha_Q       = -2.0,     # intercept (targets ~10-15% marginal event rate)
  coef_W1_Q     = 0.8,      # W1 effect on outcome (shared confounder)
  coef_W2_Q     = 0.4,      # W2 effect on outcome
  coef_W3_Q     = -0.3,     # W3 effect on outcome
  coef_W4_Q     = 0.05,     # W4 effect on outcome (weak! key for TMLE success)
  nonlinear_Q   = TRUE,     # add 0.3*W1^2 to outcome model

  # TMLE / estimation
  g_trunc       = c(0.01, 0.99),  # truncation bounds for PS in TMLE/IPTW
  Q_SL_library  = c("SL.glm", "SL.gam", "SL.glmnet", "SL.mean"),
  g_SL_library  = c("SL.glm", "SL.glmnet", "SL.mean"),
  iptw_boot_B   = 200,      # bootstrap replicates for IPTW CI

  # Stage 2 diagnostic thresholds
  diag_extreme_threshold = 0.30,  # flag "bad" if >30% of g_hat < 0.05 or > 0.95
  diag_ess_threshold     = 0.20,  # flag "bad" if ESS/N < 0.20

  # Output
  output_dir    = "outputs/sim_scenario1"
)

# ===========================================================================
# 0. Setup
# ===========================================================================

message("=== Scenario 1: Bad Overlap, TMLE OK ===")
message("Start: ", Sys.time())

# Check required packages
for (pkg in c("tmle", "SuperLearner")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package '", pkg, "' is not installed. ",
         "Install with: install.packages('", pkg, "')", call. = FALSE)
  }
}

# Filter SL libraries to only those with installed backing packages
filter_sl_libs <- function(libs) {
  keep <- vapply(libs, function(lib) {
    if (lib %in% c("SL.glm", "SL.mean", "SL.step",
                    "SL.step.interaction", "SL.glm.interaction")) {
      return(TRUE)
    }
    pkg <- sub("^SL\\.", "", lib)
    requireNamespace(pkg, quietly = TRUE)
  }, logical(1))
  out <- libs[keep]
  if (length(out) == 0) out <- c("SL.glm", "SL.mean")
  out
}

PARAMS$Q_SL_library <- filter_sl_libs(PARAMS$Q_SL_library)
PARAMS$g_SL_library <- filter_sl_libs(PARAMS$g_SL_library)

message("Q.SL.library: ", paste(PARAMS$Q_SL_library, collapse = ", "))
message("g.SL.library: ", paste(PARAMS$g_SL_library, collapse = ", "))

# Create output directory
dir.create(PARAMS$output_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# STAGE 1: PRE-SPECIFICATION
# ===========================================================================
# Estimand:
#   RD = E[Y^1] - E[Y^0]  (ATE on risk-difference scale)
#   Y is a binary endpoint. No time-to-event structure in this scenario.
#
# Planned estimators:
#   1. TMLE (primary) -- uses both outcome model Q and treatment model g;
#      doubly robust, so good Q can compensate for poor g estimation.
#   2. Crude RD (comparator) -- unadjusted difference in means.
#      Expected to be biased due to confounding.
#   3. IPTW RD (comparator) -- inverse probability weighted.
#      Expected to be unstable/under-covered due to extreme weights
#      from poor overlap.
#
# Truncation: g(W) clipped to [0.01, 0.99] for IPTW and TMLE.
#
# Diagnostic thresholds (Stage 2, outcome-blind):
#   overlap_flag = "bad" if:
#     - fraction of g_hat in (0, 0.05) or (0.95, 1) exceeds 30%, OR
#     - ESS/N under IPTW weights < 0.20
#   This scenario is DESIGNED to trigger "bad" in most replicates.
#
# Key DGP feature:
#   W4 is a strong driver of treatment assignment but has near-zero
#   effect on the outcome. This means the PS is extreme (poor overlap)
#   but Q(A,W) remains smooth and learnable in the region of extremity.
#   TMLE exploits this through its Q-model targeting step.
# ===========================================================================

message("\n--- Stage 1: Pre-specification ---")
message("Estimand: RD = E[Y^1] - E[Y^0]")
message("Primary estimator: TMLE")
message("Comparators: crude RD, IPTW RD")
message("Truncation bounds: [", PARAMS$g_trunc[1], ", ", PARAMS$g_trunc[2], "]")
message("Overlap flag thresholds: extreme > ", PARAMS$diag_extreme_threshold,
        ", ESS/N < ", PARAMS$diag_ess_threshold)


# ===========================================================================
# DGP: Data-generating function
# ===========================================================================

expit <- function(x) 1 / (1 + exp(-x))

#' Generate one replicate of Scenario 1 data
#'
#' @param N sample size
#' @param seed replicate seed
#' @param p list of parameters (PARAMS)
#' @return list with data (data.frame), g_true, Q_true
generate_data <- function(N, seed = NULL, p = PARAMS) {
  if (!is.null(seed)) set.seed(seed)

  # Covariates
  W1 <- rnorm(N, 0, 1)       # shared confounder (moderate in both g and Q)
  W2 <- rbinom(N, 1, 0.5)    # binary confounder
  W3 <- rnorm(N, 0, 1)       # continuous confounder
  W4 <- rnorm(N, 0, 1)       # overlap driver: strong in g, weak in Q

  # ---- Treatment mechanism ----
  # W4 dominates treatment assignment, creating extreme PS values,
  # but W4 barely affects Y, so Q is smooth at the PS boundary.
  lp_g_raw <- p$strength_overlap * (
    p$coef_W1_g * W1 +
    p$coef_W2_g * W2 +
    p$coef_W3_g * W3 +
    p$coef_W4_g * W4
  )
  # Center so marginal P(A=1) ~ 0.5
  alpha_g <- -mean(lp_g_raw)
  lp_g <- alpha_g + lp_g_raw
  g_true <- expit(lp_g)
  # Clip to [eps, 1-eps] -- practical positivity, not structural violation
  g_true <- pmin(pmax(g_true, p$eps_clip), 1 - p$eps_clip)
  A <- rbinom(N, 1, g_true)

  # ---- Outcome model ----
  # Learnable by GLM (plus optional mild nonlinearity).
  # W4 coefficient is near zero: outcome is nearly unaffected by the
  # variable that drives PS extremity.
  lp_Q <- p$alpha_Q + p$effect_A * A +
    p$coef_W1_Q * W1 + p$coef_W2_Q * W2 +
    p$coef_W3_Q * W3 + p$coef_W4_Q * W4
  if (p$nonlinear_Q) {
    lp_Q <- lp_Q + 0.3 * W1^2
  }
  Q_true <- expit(lp_Q)
  Y <- rbinom(N, 1, Q_true)

  dat <- data.frame(W1 = W1, W2 = W2, W3 = W3, W4 = W4, A = A, Y = Y)

  list(data = dat, g_true = g_true, Q_true = Q_true)
}


# ===========================================================================
# GROUND TRUTH: Compute RD_true via exact plug-in risks
# ===========================================================================

message("\n--- Computing ground truth RD ---")

compute_truth <- function(truth_N, seed, p) {
  set.seed(seed)
  W1 <- rnorm(truth_N, 0, 1)
  W2 <- rbinom(truth_N, 1, 0.5)
  W3 <- rnorm(truth_N, 0, 1)
  W4 <- rnorm(truth_N, 0, 1)

  # Exact risks under A=1 and A=0 (no randomness -- use expit directly)
  lp_base <- p$alpha_Q +
    p$coef_W1_Q * W1 + p$coef_W2_Q * W2 +
    p$coef_W3_Q * W3 + p$coef_W4_Q * W4
  if (p$nonlinear_Q) {
    lp_base <- lp_base + 0.3 * W1^2
  }
  lp_Q1 <- lp_base + p$effect_A * 1
  lp_Q0 <- lp_base + p$effect_A * 0
  p1 <- expit(lp_Q1)
  p0 <- expit(lp_Q0)

  list(
    EY1     = mean(p1),
    EY0     = mean(p0),
    RD_true = mean(p1) - mean(p0),
    truth_N = truth_N,
    note    = "Exact plug-in risks via expit (no Monte Carlo noise in Y)"
  )
}

truth <- compute_truth(PARAMS$truth_N, seed = PARAMS$master_seed, params = PARAMS)
message("  RD_true = ", round(truth$RD_true, 6))
message("  E[Y^1]  = ", round(truth$EY1, 4), ",  E[Y^0] = ", round(truth$EY0, 4))


# ===========================================================================
# STAGE 2 + 3: Replicate loop
# ===========================================================================

message("\n--- Running ", PARAMS$reps, " replicates (N=", PARAMS$N, ") ---")

results <- vector("list", PARAMS$reps)

for (r in seq_len(PARAMS$reps)) {

  rep_seed <- PARAMS$master_seed + r

  # Generate data
  sim <- generate_data(N = PARAMS$N, seed = rep_seed, p = PARAMS)
  dat <- sim$data
  W   <- dat[, c("W1", "W2", "W3", "W4"), drop = FALSE]
  A   <- dat$A
  Y   <- dat$Y
  n   <- nrow(dat)

  # ------------------------------------------------------------------
  # Stage 2: Outcome-blind diagnostics (uses only A and W)
  # ------------------------------------------------------------------
  # Fit PS with GLM (fast, adequate for diagnostics)
  g_fit <- glm(A ~ W1 + W2 + W3 + W4, data = dat, family = "binomial")
  g_hat <- fitted(g_fit)

  # Diagnostics
  fraction_extreme <- mean(g_hat < 0.05 | g_hat > 0.95)
  min_g <- min(g_hat)
  max_g <- max(g_hat)

  # IPTW weights (truncated)
  g_hat_t <- pmin(pmax(g_hat, PARAMS$g_trunc[1]), PARAMS$g_trunc[2])
  w_iptw  <- A / g_hat_t + (1 - A) / (1 - g_hat_t)
  ess     <- sum(w_iptw)^2 / sum(w_iptw^2)
  ess_frac <- ess / n
  max_w   <- max(w_iptw)

  overlap_flag <- ifelse(
    fraction_extreme > PARAMS$diag_extreme_threshold ||
      ess_frac < PARAMS$diag_ess_threshold,
    "bad", "ok"
  )

  # ------------------------------------------------------------------
  # Stage 3: Estimation
  # ------------------------------------------------------------------

  # ---- Estimator 1: TMLE ----
  tmle_res <- tryCatch({
    fit <- tmle::tmle(
      Y = Y, A = A, W = W,
      family = "binomial",
      Q.SL.library = PARAMS$Q_SL_library,
      g.SL.library = PARAMS$g_SL_library,
      gbound = PARAMS$g_trunc
    )
    ate <- fit$estimates$ATE
    list(
      est   = ate$psi,
      se    = sqrt(ate$var.psi),
      ci_lo = ate$CI[1],
      ci_hi = ate$CI[2],
      ok    = TRUE
    )
  }, error = function(e) {
    # Fallback: GLM-only TMLE
    tryCatch({
      fit <- tmle::tmle(
        Y = Y, A = A, W = W,
        family = "binomial",
        Q.SL.library = "SL.glm",
        g.SL.library = "SL.glm",
        gbound = PARAMS$g_trunc
      )
      ate <- fit$estimates$ATE
      list(
        est   = ate$psi,
        se    = sqrt(ate$var.psi),
        ci_lo = ate$CI[1],
        ci_hi = ate$CI[2],
        ok    = TRUE
      )
    }, error = function(e2) {
      list(est = NA, se = NA, ci_lo = NA, ci_hi = NA, ok = FALSE)
    })
  })

  # ---- Estimator 2: Crude (unadjusted) RD ----
  n1 <- sum(A == 1)
  n0 <- sum(A == 0)
  if (n1 >= 1 && n0 >= 1) {
    p1_crude <- mean(Y[A == 1])
    p0_crude <- mean(Y[A == 0])
    rd_crude <- p1_crude - p0_crude
    se_crude <- sqrt(p1_crude * (1 - p1_crude) / n1 +
                       p0_crude * (1 - p0_crude) / n0)
    ci_crude <- rd_crude + c(-1.96, 1.96) * se_crude
  } else {
    rd_crude <- NA_real_; se_crude <- NA_real_
    ci_crude <- c(NA_real_, NA_real_)
  }

  # ---- Estimator 3: IPTW RD (Hajek) with bootstrap CI ----
  iptw_hajek <- function(Y, A, w) {
    denom1 <- sum(w * A)
    denom0 <- sum(w * (1 - A))
    if (denom1 < 1e-10 || denom0 < 1e-10) return(NA_real_)
    r1 <- sum(w * A * Y) / denom1
    r0 <- sum(w * (1 - A) * Y) / denom0
    r1 - r0
  }

  rd_iptw <- iptw_hajek(Y, A, w_iptw)

  # Bootstrap for IPTW SE/CI
  set.seed(rep_seed + 1000000L)
  boot_iptw <- numeric(PARAMS$iptw_boot_B)
  for (b in seq_len(PARAMS$iptw_boot_B)) {
    idx <- sample(n, n, replace = TRUE)
    d_b <- dat[idx, , drop = FALSE]
    A_b <- d_b$A
    Y_b <- d_b$Y

    g_b <- tryCatch({
      fit_b <- glm(A ~ W1 + W2 + W3 + W4, data = d_b, family = "binomial")
      fitted(fit_b)
    }, error = function(e) rep(0.5, n))

    g_b <- pmin(pmax(g_b, PARAMS$g_trunc[1]), PARAMS$g_trunc[2])
    w_b <- A_b / g_b + (1 - A_b) / (1 - g_b)
    boot_iptw[b] <- iptw_hajek(Y_b, A_b, w_b)
  }
  boot_iptw_clean <- boot_iptw[is.finite(boot_iptw)]
  se_iptw <- if (length(boot_iptw_clean) >= 2) sd(boot_iptw_clean) else NA_real_
  ci_iptw <- if (length(boot_iptw_clean) >= 20) {
    unname(quantile(boot_iptw_clean, c(0.025, 0.975)))
  } else if (!is.na(se_iptw)) {
    rd_iptw + c(-1.96, 1.96) * se_iptw
  } else {
    c(NA_real_, NA_real_)
  }

  # ---- Collect replicate results ----
  results[[r]] <- data.frame(
    rep              = r,
    seed             = rep_seed,
    # Stage 2 diagnostics
    fraction_extreme = fraction_extreme,
    min_g            = min_g,
    max_g            = max_g,
    ess_frac         = ess_frac,
    max_w            = max_w,
    overlap_flag     = overlap_flag,
    event_rate       = mean(Y),
    n_treated        = n1,
    n_control        = n0,
    # TMLE
    tmle_est         = tmle_res$est,
    tmle_se          = tmle_res$se,
    tmle_ci_lo       = tmle_res$ci_lo,
    tmle_ci_hi       = tmle_res$ci_hi,
    tmle_ok          = tmle_res$ok,
    # Crude
    crude_est        = rd_crude,
    crude_se         = se_crude,
    crude_ci_lo      = ci_crude[1],
    crude_ci_hi      = ci_crude[2],
    # IPTW
    iptw_est         = rd_iptw,
    iptw_se          = se_iptw,
    iptw_ci_lo       = ci_iptw[1],
    iptw_ci_hi       = ci_iptw[2],
    iptw_boot_ok     = length(boot_iptw_clean),
    stringsAsFactors = FALSE
  )

  if (r %% 25 == 0 || r == 1) {
    message("  Replicate ", r, "/", PARAMS$reps,
            "  overlap=", overlap_flag,
            "  TMLE=", if (tmle_res$ok) round(tmle_res$est, 4) else "FAIL",
            "  crude=", round(rd_crude, 4),
            "  IPTW=", round(rd_iptw, 4))
  }
}


# ===========================================================================
# Aggregate results
# ===========================================================================

res_df <- do.call(rbind, results)

compute_metrics <- function(df, true_rd, label = "all") {
  tmle_ok <- df[df$tmle_ok == TRUE, , drop = FALSE]

  make_row <- function(method, est, ci_lo, ci_hi, se_col = NULL) {
    valid <- !is.na(est)
    if (sum(valid) < 3) {
      return(data.frame(
        subset = label, method = method, n_reps = sum(valid),
        bias = NA, rmse = NA, emp_sd = NA, mean_se = NA, coverage = NA,
        stringsAsFactors = FALSE
      ))
    }
    bias <- mean(est[valid]) - true_rd
    rmse <- sqrt(mean((est[valid] - true_rd)^2))
    emp_sd <- sd(est[valid])
    mean_se <- if (!is.null(se_col)) mean(se_col[valid], na.rm = TRUE) else NA
    covers <- (ci_lo[valid] <= true_rd) & (true_rd <= ci_hi[valid])
    coverage <- mean(covers)
    data.frame(
      subset = label, method = method, n_reps = sum(valid),
      bias = round(bias, 6), rmse = round(rmse, 6),
      emp_sd = round(emp_sd, 6), mean_se = round(mean_se, 6),
      coverage = round(coverage, 4),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_row("TMLE",  tmle_ok$tmle_est, tmle_ok$tmle_ci_lo,
             tmle_ok$tmle_ci_hi, tmle_ok$tmle_se),
    make_row("Crude", df$crude_est, df$crude_ci_lo,
             df$crude_ci_hi, df$crude_se),
    make_row("IPTW",  df$iptw_est, df$iptw_ci_lo,
             df$iptw_ci_hi, df$iptw_se)
  )
}

# All replicates
metrics_all <- compute_metrics(res_df, truth$RD_true, label = "all")

# Conditional on overlap_flag == "bad"
res_bad <- res_df[res_df$overlap_flag == "bad", , drop = FALSE]
if (nrow(res_bad) > 0) {
  metrics_bad <- compute_metrics(res_bad, truth$RD_true, label = "overlap_bad")
} else {
  metrics_bad <- data.frame(
    subset = "overlap_bad", method = "N/A", n_reps = 0,
    bias = NA, rmse = NA, emp_sd = NA, mean_se = NA, coverage = NA,
    stringsAsFactors = FALSE
  )
}

summary_metrics <- rbind(metrics_all, metrics_bad)


# ===========================================================================
# Print results
# ===========================================================================

message("\n", paste(rep("=", 65), collapse = ""))
message("SCENARIO 1 RESULTS: Bad Overlap, TMLE OK")
message(paste(rep("=", 65), collapse = ""))
message("\nTrue RD: ", round(truth$RD_true, 6))
message("Replicates: ", PARAMS$reps, " (N=", PARAMS$N, ")")
message("Strength overlap: ", PARAMS$strength_overlap)
message("PS clip: [", PARAMS$eps_clip, ", ", 1 - PARAMS$eps_clip, "]")
message("W4 coef in g: ", PARAMS$coef_W4_g,
        "  |  W4 coef in Q: ", PARAMS$coef_W4_Q,
        "  (overlap driver, weak in outcome)")
message("\nOverlap diagnostics across replicates:")
message("  Mean fraction g < 0.05 or g > 0.95: ",
        round(mean(res_df$fraction_extreme), 3))
message("  Mean ESS/N: ", round(mean(res_df$ess_frac), 3))
message("  Mean max weight: ", round(mean(res_df$max_w), 1))
message("  Fraction flagged 'bad': ",
        round(mean(res_df$overlap_flag == "bad"), 3))

message("\n--- Estimator Performance (all replicates) ---")
print(metrics_all)

if (nrow(res_bad) > 0) {
  message("\n--- Estimator Performance (overlap_flag == 'bad' only) ---")
  print(metrics_bad)
}

# Quick interpretation
message("\n--- Interpretation ---")
tmle_cov <- metrics_all$coverage[metrics_all$method == "TMLE"]
crude_cov <- metrics_all$coverage[metrics_all$method == "Crude"]
iptw_cov <- metrics_all$coverage[metrics_all$method == "IPTW"]
tmle_bias <- metrics_all$bias[metrics_all$method == "TMLE"]
crude_bias <- metrics_all$bias[metrics_all$method == "Crude"]
if (!is.na(tmle_cov) && tmle_cov >= 0.92) {
  message("  TMLE coverage ", tmle_cov, " is near nominal (target: ~0.95)")
} else {
  message("  TMLE coverage ", tmle_cov, " is below target.",
          " Try: increase N, or reduce coef_W4_Q toward 0")
}
if (!is.na(crude_cov) && crude_cov < 0.10) {
  message("  Crude coverage ", crude_cov, " -- confounding bias dominates (good)")
}
if (!is.na(iptw_cov) && iptw_cov < tmle_cov - 0.05) {
  message("  IPTW coverage ", iptw_cov, " -- worse than TMLE by ",
          round(tmle_cov - iptw_cov, 2), " (overlap instability)")
}


# ===========================================================================
# Save outputs
# ===========================================================================

saveRDS(PARAMS,
        file.path(PARAMS$output_dir, "params_used.rds"))
saveRDS(truth,
        file.path(PARAMS$output_dir, "truth.rds"))
saveRDS(res_df,
        file.path(PARAMS$output_dir, "replicate_results.rds"))
write.csv(summary_metrics,
          file.path(PARAMS$output_dir, "summary_metrics.csv"),
          row.names = FALSE)

# Session info
si <- capture.output(sessionInfo())
writeLines(si, file.path(PARAMS$output_dir, "sessionInfo.txt"))

message("\nOutputs saved to: ", PARAMS$output_dir)
message("  params_used.rds, truth.rds, replicate_results.rds, summary_metrics.csv")
message("\nDone: ", Sys.time())
