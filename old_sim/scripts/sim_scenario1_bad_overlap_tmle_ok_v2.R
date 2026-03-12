#!/usr/bin/env Rscript
# ===========================================================================
# Scenario 1 v2: Bad Overlap, TMLE OK — Enhanced Clean-Room Simulation
# ===========================================================================
#
# PURPOSE:
#   Demonstrate TMLE's double-robustness property under practical near-
#   positivity violations in a clean-room staged workflow. The DGP generates
#   extreme propensity scores driven by a covariate (W4) that is weakly
#   related to the outcome, so Q(A,W) remains learnable despite bad overlap.
#
# CLEAN-ROOM STAGES:
#   Stage 1:  Pre-specification of estimand, estimators, diagnostics
#   Stage 2:  Outcome-blind overlap diagnostics (PS only; stop/go gate)
#   Stage 2b: Outcome-blind operating characteristics (simulation OC)
#   Stage 3:  Estimation with blinded outcomes
#
# ESTIMATORS:
#   1. Crude RD (unadjusted difference in proportions)
#   2. Adjusted regression (outcome regression / g-computation, no targeting)
#   3. IPTW RD (Hajek, with bootstrap CI)
#   4. PS-matched regression (MatchIt nearest-neighbor on logit PS)
#   5. TMLE (IC-based inference)
#   6. TMLE targeted bootstrap (EIF resampling, nuisance fits held fixed)
#   7. TMLE nonparametric bootstrap (full refit; OFF by default)
#
# USAGE:
#   Rscript scripts/sim_scenario1_bad_overlap_tmle_ok_v2.R
#
# OUTPUT:
#   results/scenario1_v2/  — truth, per-rep, summary, figures
# ===========================================================================

# ===========================================================================
# A1. Package imports and reproducibility
# ===========================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Optional packages — degrade gracefully
has_tmle      <- requireNamespace("tmle", quietly = TRUE)
has_SL        <- requireNamespace("SuperLearner", quietly = TRUE)
has_MatchIt   <- requireNamespace("MatchIt", quietly = TRUE)
has_sandwich  <- requireNamespace("sandwich", quietly = TRUE)
has_future    <- requireNamespace("future.apply", quietly = TRUE)
has_digest    <- requireNamespace("digest", quietly = TRUE)

if (!has_tmle) stop("Package 'tmle' is required. Install with: install.packages('tmle')")

expit <- function(x) 1 / (1 + exp(-x))

# ===========================================================================
# A2. Parameters
# ===========================================================================

params <- list(
  # Sample sizes
  n_main    = 500,        # sample size per replicate (small for speed)
  n_truth   = 100000,     # sample size for ground truth plug-in
  reps      = 200,        # Monte Carlo replicates
  master_seed = 20260305, # single master seed


  # Bootstrap settings
  B_tboot   = 500,        # targeted bootstrap (EIF resampling) replicates
  B_npboot  = 0,          # nonparametric bootstrap for TMLE (0 = OFF)
  B_iptw    = 200,        # bootstrap for IPTW CI

  # Propensity score clipping
  ps_clip_true = c(0.005, 0.995),  # clip true g when drawing A

  ps_clip_hat  = c(0.01, 0.99),    # clip estimated g before use

  # DGP: treatment mechanism (controls overlap)
  strength_overlap = 1.5,  # global multiplier on g-model linear predictor
  coef_W1_g = 1.5,         # W1: shared confounder
  coef_W2_g = 0.8,         # W2: binary confounder
  coef_W3_g = -0.5,        # W3: continuous confounder
  coef_W4_g = 3.5,         # W4: overlap driver (strong in g, weak in Q)

  # DGP: outcome model
  effect_A  = 0.5,         # treatment effect on log-odds scale
  alpha_Q   = -2.0,        # intercept (targets ~10-15% marginal event rate)
  coef_W1_Q = 0.8,
  coef_W2_Q = 0.4,
  coef_W3_Q = -0.3,
  coef_W4_Q = 0.05,        # NEAR ZERO — key structural feature
  nonlinear_Q = TRUE,      # adds 0.3*W1^2

  # Overlap flag thresholds (Stage 2)
  overlap_flag_thresholds = list(
    extreme_ps   = 0.05,    # boundary for "extreme" PS
    extreme_prop = 0.25,    # flag if > this fraction is extreme
    ess_frac     = 0.2,     # flag if ESS/N < this
    max_w        = 50       # flag if max weight > this
  ),

  # Matching
  match_caliper   = NULL,   # NULL = 0.2*sd(logit(ps)); set numeric to override
  match_replace   = FALSE,  # matching with replacement toggle

  # Regression specification
  regression_formula = Y ~ A + W1 + W2 + W3 + W4 + I(W1^2),  # correct spec
  misspecify_regression = FALSE,  # if TRUE, drop W1^2 and W4

  # SuperLearner libraries
  Q_SL_library = c("SL.glm", "SL.mean"),
  g_SL_library = c("SL.glm", "SL.mean"),

  # Parallelism
  workers = 1,    # default sequential; set > 1 for parallel

  # Output
  output_dir   = "results/scenario1_v2",
  save_per_rep = TRUE
)

# ===========================================================================
# Filter SL libraries to installed packages
# ===========================================================================

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

params$Q_SL_library <- filter_sl_libs(params$Q_SL_library)
params$g_SL_library <- filter_sl_libs(params$g_SL_library)

# ===========================================================================
# A3. Data-generating process
# ===========================================================================

#' Generate one replicate dataset
#'
#' @param n sample size
#' @param seed replicate-specific seed
#' @param p parameter list
#' @return list with $data (data.frame), $g_true, $Q_true
sim_data <- function(n, seed = NULL, p = params) {
  if (!is.null(seed)) set.seed(seed)

  # Covariates
  W1 <- rnorm(n, 0, 1)        # continuous, shared confounder
  W2 <- rbinom(n, 1, 0.5)     # binary confounder
  W3 <- rnorm(n, 0, 1)        # continuous confounder
  W4 <- rnorm(n, 0, 1)        # overlap driver: strong in g, weak in Q

  # Treatment mechanism
  lp_g_raw <- p$strength_overlap * (
    p$coef_W1_g * W1 + p$coef_W2_g * W2 +
    p$coef_W3_g * W3 + p$coef_W4_g * W4
  )
  alpha_g <- -mean(lp_g_raw)   # center so P(A=1) ~ 0.5
  lp_g <- alpha_g + lp_g_raw
  g_true <- expit(lp_g)
  g_true <- pmin(pmax(g_true, p$ps_clip_true[1]), p$ps_clip_true[2])
  A <- rbinom(n, 1, g_true)

  # Outcome model
  lp_Q <- p$alpha_Q + p$effect_A * A +
    p$coef_W1_Q * W1 + p$coef_W2_Q * W2 +
    p$coef_W3_Q * W3 + p$coef_W4_Q * W4
  if (p$nonlinear_Q) lp_Q <- lp_Q + 0.3 * W1^2
  Q_true <- expit(lp_Q)
  Y <- rbinom(n, 1, Q_true)

  dat <- data.frame(W1 = W1, W2 = W2, W3 = W3, W4 = W4, A = A, Y = Y)
  list(data = dat, g_true = g_true, Q_true = Q_true)
}

# ===========================================================================
# A4. Truth computation
# ===========================================================================

compute_truth <- function(n_truth, seed, p) {
  set.seed(seed)
  W1 <- rnorm(n_truth, 0, 1)
  W2 <- rbinom(n_truth, 1, 0.5)
  W3 <- rnorm(n_truth, 0, 1)
  W4 <- rnorm(n_truth, 0, 1)

  lp_base <- p$alpha_Q +
    p$coef_W1_Q * W1 + p$coef_W2_Q * W2 +
    p$coef_W3_Q * W3 + p$coef_W4_Q * W4
  if (p$nonlinear_Q) lp_base <- lp_base + 0.3 * W1^2

  EY1 <- mean(expit(lp_base + p$effect_A))
  EY0 <- mean(expit(lp_base))

  list(EY1 = EY1, EY0 = EY0, RD_true = EY1 - EY0,
       truth_N = n_truth,
       note = "Exact plug-in via expit (no MC noise in Y)")
}

# ===========================================================================
# A5. Stage 2: Outcome-blind overlap diagnostics
# ===========================================================================

#' Estimate propensity scores via logistic regression
#' @param dat data.frame with A, W1..W4
#' @param p parameter list
#' @return list with $g_raw, $g_bounded
estimate_ps <- function(dat, p = params) {
  g_fit <- glm(A ~ W1 + W2 + W3 + W4, data = dat, family = "binomial")
  g_raw <- fitted(g_fit)
  g_bounded <- pmin(pmax(g_raw, p$ps_clip_hat[1]), p$ps_clip_hat[2])
  list(g_raw = g_raw, g_bounded = g_bounded, fit = g_fit)
}

#' Compute overlap diagnostics (outcome-blind)
#' @param g_raw raw PS estimates
#' @param g_bounded truncated PS estimates
#' @param A treatment vector
#' @param p parameter list
#' @return named list of diagnostics
overlap_diagnostics <- function(g_raw, g_bounded, A, p = params) {
  thr <- p$overlap_flag_thresholds
  n <- length(A)

  # Fraction extreme (raw PS)
  extreme_lo <- thr$extreme_ps
  extreme_hi <- 1 - thr$extreme_ps
  frac_extreme <- mean(g_raw < extreme_lo | g_raw > extreme_hi)

  # IPTW weights on bounded PS
  w_iptw <- A / g_bounded + (1 - A) / (1 - g_bounded)
  ess <- sum(w_iptw)^2 / sum(w_iptw^2)
  ess_frac <- ess / n
  max_w <- max(w_iptw)

  # Flag
  flag <- (frac_extreme > thr$extreme_prop) ||
          (ess_frac < thr$ess_frac) ||
          (max_w > thr$max_w)

  list(
    fraction_extreme = frac_extreme,
    min_g_raw  = min(g_raw),
    max_g_raw  = max(g_raw),
    min_g_bounded = min(g_bounded),
    max_g_bounded = max(g_bounded),
    ess        = ess,
    ess_frac   = ess_frac,
    max_w      = max_w,
    overlap_flag = ifelse(flag, "bad", "ok"),
    w_iptw     = w_iptw
  )
}

# ===========================================================================
# A6. Estimator functions
# ===========================================================================

# ---- (1) Crude RD ----
est_crude <- function(Y, A) {
  tryCatch({
    n1 <- sum(A == 1); n0 <- sum(A == 0)
    if (n1 < 2 || n0 < 2) return(list(estimate = NA, se = NA,
      ci_low = NA, ci_high = NA, warnings = "too few in one arm"))
    p1 <- mean(Y[A == 1]); p0 <- mean(Y[A == 0])
    rd <- p1 - p0
    se <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    list(estimate = rd, se = se,
         ci_low = rd - 1.96 * se, ci_high = rd + 1.96 * se,
         warnings = NA_character_)
  }, error = function(e) {
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         warnings = conditionMessage(e))
  })
}

# ---- (2) Adjusted regression (g-computation, no targeting) ----
est_reg_adj <- function(dat, p = params) {
  tryCatch({
    fml <- if (p$misspecify_regression) {
      Y ~ A + W1 + W2 + W3  # drops W4 and W1^2
    } else {
      Y ~ A + W1 + W2 + W3 + W4 + I(W1^2)
    }
    fit <- glm(fml, data = dat, family = "binomial")

    # G-computation: predict under A=1 and A=0
    d1 <- d0 <- dat
    d1$A <- 1; d0$A <- 0
    EY1 <- mean(predict(fit, newdata = d1, type = "response"))
    EY0 <- mean(predict(fit, newdata = d0, type = "response"))
    rd <- EY1 - EY0

    # Sandwich SE for marginal RD via delta method is complex;
    # use simple paired differences approach for approximate SE
    p1_i <- predict(fit, newdata = d1, type = "response")
    p0_i <- predict(fit, newdata = d0, type = "response")
    rd_i <- p1_i - p0_i
    n <- nrow(dat)
    se <- sd(rd_i) / sqrt(n)

    list(estimate = rd, se = se,
         ci_low = rd - 1.96 * se, ci_high = rd + 1.96 * se,
         warnings = NA_character_)
  }, error = function(e) {
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         warnings = conditionMessage(e))
  })
}

# ---- (3) IPTW RD (Hajek) ----
est_iptw <- function(Y, A, g_bounded, B_boot = 200, dat = NULL, p = params) {
  tryCatch({
    w <- A / g_bounded + (1 - A) / (1 - g_bounded)
    denom1 <- sum(w * A); denom0 <- sum(w * (1 - A))
    if (denom1 < 1e-10 || denom0 < 1e-10) {
      return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                  warnings = "zero denominator"))
    }
    r1 <- sum(w * A * Y) / denom1
    r0 <- sum(w * (1 - A) * Y) / denom0
    rd <- r1 - r0

    se <- NA_real_; ci <- c(NA_real_, NA_real_)
    if (B_boot > 0 && !is.null(dat)) {
      n <- nrow(dat)
      boot_rd <- numeric(B_boot)
      for (b in seq_len(B_boot)) {
        idx <- sample(n, n, replace = TRUE)
        d_b <- dat[idx, , drop = FALSE]
        g_b <- tryCatch({
          fit_b <- glm(A ~ W1 + W2 + W3 + W4, data = d_b, family = "binomial")
          pmin(pmax(fitted(fit_b), p$ps_clip_hat[1]), p$ps_clip_hat[2])
        }, error = function(e) rep(0.5, n))
        w_b <- d_b$A / g_b + (1 - d_b$A) / (1 - g_b)
        d1_b <- sum(w_b * d_b$A); d0_b <- sum(w_b * (1 - d_b$A))
        if (d1_b > 1e-10 && d0_b > 1e-10) {
          boot_rd[b] <- sum(w_b * d_b$A * d_b$Y) / d1_b -
                         sum(w_b * (1 - d_b$A) * d_b$Y) / d0_b
        } else {
          boot_rd[b] <- NA
        }
      }
      ok <- boot_rd[is.finite(boot_rd)]
      se <- if (length(ok) >= 2) sd(ok) else NA_real_
      ci <- if (length(ok) >= 20) {
        unname(quantile(ok, c(0.025, 0.975)))
      } else if (!is.na(se)) {
        rd + c(-1.96, 1.96) * se
      } else {
        c(NA_real_, NA_real_)
      }
    }

    list(estimate = rd, se = se, ci_low = ci[1], ci_high = ci[2],
         warnings = NA_character_)
  }, error = function(e) {
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         warnings = conditionMessage(e))
  })
}

# ---- (4) PS-matched regression ----
est_psmatch_reg <- function(dat, g_bounded, p = params) {
  if (!has_MatchIt) {
    return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                n_matched = NA, smd_before = NA, smd_after = NA,
                warnings = "MatchIt not installed"))
  }
  tryCatch({
    dat$ps <- g_bounded
    dat$logit_ps <- log(g_bounded / (1 - g_bounded))

    # Caliper
    cal <- p$match_caliper
    if (is.null(cal)) cal <- 0.2 * sd(dat$logit_ps)

    m <- MatchIt::matchit(
      A ~ W1 + W2 + W3 + W4,
      data = dat,
      method = "nearest",
      distance = dat$logit_ps,
      caliper = cal,
      replace = p$match_replace
    )

    mdat <- MatchIt::match.data(m)
    n_matched <- nrow(mdat)

    # SMD before/after for W1..W4
    covs <- c("W1", "W2", "W3", "W4")
    smd_before <- sapply(covs, function(v) {
      m1 <- mean(dat[[v]][dat$A == 1]); m0 <- mean(dat[[v]][dat$A == 0])
      s <- sqrt((var(dat[[v]][dat$A == 1]) + var(dat[[v]][dat$A == 0])) / 2)
      if (s < 1e-10) return(0)
      (m1 - m0) / s
    })
    smd_after <- sapply(covs, function(v) {
      m1 <- mean(mdat[[v]][mdat$A == 1]); m0 <- mean(mdat[[v]][mdat$A == 0])
      s <- sqrt((var(mdat[[v]][mdat$A == 1]) + var(mdat[[v]][mdat$A == 0])) / 2)
      if (s < 1e-10) return(0)
      (m1 - m0) / s
    })

    # Fit regression on matched data, g-computation for RD
    fml <- if (p$misspecify_regression) {
      Y ~ A + W1 + W2 + W3
    } else {
      Y ~ A + W1 + W2 + W3 + W4 + I(W1^2)
    }
    fit_m <- glm(fml, data = mdat, family = "binomial")
    d1 <- d0 <- mdat
    d1$A <- 1; d0$A <- 0
    EY1 <- mean(predict(fit_m, newdata = d1, type = "response"))
    EY0 <- mean(predict(fit_m, newdata = d0, type = "response"))
    rd <- EY1 - EY0

    # Robust SE on matched data
    p1_i <- predict(fit_m, newdata = d1, type = "response")
    p0_i <- predict(fit_m, newdata = d0, type = "response")
    rd_i <- p1_i - p0_i
    se <- sd(rd_i) / sqrt(n_matched)

    list(estimate = rd, se = se,
         ci_low = rd - 1.96 * se, ci_high = rd + 1.96 * se,
         n_matched = n_matched,
         n_unmatched = nrow(dat) - n_matched,
         smd_before = smd_before,
         smd_after = smd_after,
         warnings = NA_character_)
  }, error = function(e) {
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         n_matched = NA, n_unmatched = NA,
         smd_before = NA, smd_after = NA,
         warnings = conditionMessage(e))
  })
}

# ---- (5) TMLE ----
est_tmle <- function(Y, A, W, p = params) {
  tryCatch({
    fit <- tmle::tmle(
      Y = Y, A = A, W = W,
      family = "binomial",
      Q.SL.library = p$Q_SL_library,
      g.SL.library = p$g_SL_library,
      gbound = p$ps_clip_hat
    )
    ate <- fit$estimates$ATE
    # Extract EIF (influence curve values)
    ic <- fit$estimates$IC$IC.ATE

    list(estimate = ate$psi,
         se = sqrt(ate$var.psi),
         ci_low = ate$CI[1],
         ci_high = ate$CI[2],
         ic = ic,
         ok = TRUE,
         warnings = NA_character_)
  }, error = function(e) {
    # Fallback to GLM-only
    tryCatch({
      fit <- tmle::tmle(
        Y = Y, A = A, W = W,
        family = "binomial",
        Q.SL.library = "SL.glm",
        g.SL.library = "SL.glm",
        gbound = p$ps_clip_hat
      )
      ate <- fit$estimates$ATE
      ic <- fit$estimates$IC$IC.ATE
      list(estimate = ate$psi, se = sqrt(ate$var.psi),
           ci_low = ate$CI[1], ci_high = ate$CI[2],
           ic = ic, ok = TRUE,
           warnings = paste0("SL fallback: ", conditionMessage(e)))
    }, error = function(e2) {
      list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
           ic = NULL, ok = FALSE,
           warnings = conditionMessage(e2))
    })
  })
}

# ---- (6) Targeted bootstrap for TMLE (EIF resampling) ----
#' Targeted bootstrap: resamples EIF with nuisance fits held fixed.
#' Uses Bayesian bootstrap (Dirichlet weights) for stability.
#'
#' @param psi_hat TMLE point estimate
#' @param ic_values EIF values D_i (one per observation)
#' @param B number of bootstrap draws
#' @return list with se, percentile CI, all bootstrap estimates
targeted_bootstrap_tmle <- function(psi_hat, ic_values, B = 500) {
  if (is.null(ic_values) || length(ic_values) < 5) {
    return(list(se = NA, ci_low = NA, ci_high = NA,
                psi_boot = NULL, warnings = "insufficient IC values"))
  }
  n <- length(ic_values)
  psi_boot <- numeric(B)
  for (b in seq_len(B)) {
    # Bayesian bootstrap: exponential(1) weights, normalized
    w_b <- rexp(n, rate = 1)
    w_b <- w_b / sum(w_b)
    # Bootstrap replicate: psi_hat + weighted mean of IC
    psi_boot[b] <- psi_hat + sum(w_b * ic_values)
  }
  se <- sd(psi_boot)
  ci <- unname(quantile(psi_boot, c(0.025, 0.975)))
  list(se = se, ci_low = ci[1], ci_high = ci[2],
       psi_boot = psi_boot, warnings = NA_character_)
}

# ---- (7) Nonparametric bootstrap for TMLE (full refit; OFF by default) ----
npboot_tmle <- function(dat, B, p = params) {
  if (B <= 0) {
    return(list(se = NA, ci_low = NA, ci_high = NA,
                psi_boot = NULL, warnings = "npboot disabled"))
  }
  tryCatch({
    n <- nrow(dat)
    W_cols <- c("W1", "W2", "W3", "W4")
    psi_boot <- numeric(B)
    n_fail <- 0
    for (b in seq_len(B)) {
      idx <- sample(n, n, replace = TRUE)
      d_b <- dat[idx, , drop = FALSE]
      res_b <- tryCatch({
        fit_b <- tmle::tmle(
          Y = d_b$Y, A = d_b$A, W = d_b[, W_cols, drop = FALSE],
          family = "binomial",
          Q.SL.library = "SL.glm",
          g.SL.library = "SL.glm",
          gbound = p$ps_clip_hat
        )
        fit_b$estimates$ATE$psi
      }, error = function(e) NA_real_)
      psi_boot[b] <- res_b
      if (is.na(res_b)) n_fail <- n_fail + 1
    }
    ok <- psi_boot[is.finite(psi_boot)]
    se <- if (length(ok) >= 10) sd(ok) else NA_real_
    ci <- if (length(ok) >= 20) {
      unname(quantile(ok, c(0.025, 0.975)))
    } else {
      c(NA_real_, NA_real_)
    }
    list(se = se, ci_low = ci[1], ci_high = ci[2],
         psi_boot = ok,
         n_fail = n_fail,
         warnings = if (n_fail > 0) paste0(n_fail, "/", B, " bootstrap failures") else NA_character_)
  }, error = function(e) {
    list(se = NA, ci_low = NA, ci_high = NA,
         psi_boot = NULL, n_fail = NA,
         warnings = conditionMessage(e))
  })
}


# ===========================================================================
# 0. Setup
# ===========================================================================

message("=== Scenario 1 v2: Bad Overlap, TMLE OK ===")
message("Start: ", Sys.time())
message("Q.SL.library: ", paste(params$Q_SL_library, collapse = ", "))
message("g.SL.library: ", paste(params$g_SL_library, collapse = ", "))
message("MatchIt available: ", has_MatchIt)
message("Workers: ", params$workers)

dir.create(params$output_dir, recursive = TRUE, showWarnings = FALSE)
if (params$save_per_rep) {
  dir.create(file.path(params$output_dir, "per_rep"),
             recursive = TRUE, showWarnings = FALSE)
}
dir.create(file.path(params$output_dir, "figures"),
           recursive = TRUE, showWarnings = FALSE)


# ===========================================================================
# STAGE 1: PRE-SPECIFICATION
# ===========================================================================

message("\n--- Stage 1: Pre-specification ---")
message("Estimand: RD = E[Y^1] - E[Y^0]")
message("Primary: TMLE (IC + targeted bootstrap)")
message("Comparators: crude, adjusted regression, IPTW, PS-matched regression")
message("PS clipping: [", params$ps_clip_hat[1], ", ", params$ps_clip_hat[2], "]")


# ===========================================================================
# A4. Truth computation (with caching)
# ===========================================================================

message("\n--- Computing ground truth RD ---")

truth_file <- file.path(params$output_dir, "truth.rds")
# Cache: reuse truth if params hash matches
params_hash <- if (has_digest) {
  digest::digest(params[c("n_truth", "master_seed", "alpha_Q", "effect_A",
                           "coef_W1_Q", "coef_W2_Q", "coef_W3_Q", "coef_W4_Q",
                           "nonlinear_Q")])
} else {
  ""
}

if (file.exists(truth_file)) {
  truth_cached <- readRDS(truth_file)
  if (!is.null(truth_cached$params_hash) && truth_cached$params_hash == params_hash && params_hash != "") {
    truth <- truth_cached
    message("  Loaded cached truth: RD_true = ", round(truth$RD_true, 6))
  } else {
    truth <- compute_truth(params$n_truth, seed = params$master_seed, p = params)
    truth$params_hash <- params_hash
    saveRDS(truth, truth_file)
    message("  Computed fresh truth: RD_true = ", round(truth$RD_true, 6))
  }
} else {
  truth <- compute_truth(params$n_truth, seed = params$master_seed, p = params)
  truth$params_hash <- params_hash
  saveRDS(truth, truth_file)
  message("  Computed truth: RD_true = ", round(truth$RD_true, 6))
}
message("  E[Y^1] = ", round(truth$EY1, 4), ",  E[Y^0] = ", round(truth$EY0, 4))


# ===========================================================================
# STAGE 2 + 3: Replicate loop
# ===========================================================================

message("\n--- Running ", params$reps, " replicates (N=", params$n_main, ") ---")

run_one_rep <- function(r, p = params, truth_rd = truth$RD_true) {
  rep_seed <- p$master_seed + r

  # Generate data
  sim <- sim_data(n = p$n_main, seed = rep_seed, p = p)
  dat <- sim$data
  W   <- dat[, c("W1", "W2", "W3", "W4"), drop = FALSE]
  A   <- dat$A
  Y   <- dat$Y
  n   <- nrow(dat)

  # ------------------------------------------------------------------
  # Stage 2: Outcome-blind diagnostics
  # ------------------------------------------------------------------
  ps <- estimate_ps(dat, p)
  diag <- overlap_diagnostics(ps$g_raw, ps$g_bounded, A, p)

  # ------------------------------------------------------------------
  # Stage 3: Estimation
  # ------------------------------------------------------------------

  # (1) Crude RD
  set.seed(rep_seed + 100000L)
  crude <- est_crude(Y, A)

  # (2) Adjusted regression
  set.seed(rep_seed + 200000L)
  reg_adj <- est_reg_adj(dat, p)

  # (3) IPTW
  set.seed(rep_seed + 300000L)
  iptw <- est_iptw(Y, A, ps$g_bounded, B_boot = p$B_iptw, dat = dat, p = p)

  # (4) PS-matched regression
  set.seed(rep_seed + 400000L)
  ps_match <- est_psmatch_reg(dat, ps$g_bounded, p)

  # (5) TMLE
  set.seed(rep_seed + 500000L)
  tmle_res <- est_tmle(Y, A, W, p)

  # (6) Targeted bootstrap for TMLE
  tboot <- list(se = NA, ci_low = NA, ci_high = NA, psi_boot = NULL, warnings = NA)
  if (tmle_res$ok && !is.null(tmle_res$ic) && p$B_tboot > 0) {
    set.seed(rep_seed + 600000L)
    tboot <- targeted_bootstrap_tmle(tmle_res$estimate, tmle_res$ic, B = p$B_tboot)
  }

  # (7) Nonparametric bootstrap for TMLE (OFF by default)
  npboot <- list(se = NA, ci_low = NA, ci_high = NA, psi_boot = NULL, warnings = "disabled")
  if (p$B_npboot > 0 && tmle_res$ok) {
    set.seed(rep_seed + 700000L)
    npboot <- npboot_tmle(dat, B = p$B_npboot, p = p)
  }

  # ------------------------------------------------------------------
  # Collect results
  # ------------------------------------------------------------------
  row <- data.frame(
    rep              = r,
    seed             = rep_seed,
    n                = n,

    # Stage 2 diagnostics
    fraction_extreme = diag$fraction_extreme,
    min_g_raw        = diag$min_g_raw,
    max_g_raw        = diag$max_g_raw,
    ess_frac         = diag$ess_frac,
    max_w            = diag$max_w,
    overlap_flag     = diag$overlap_flag,
    event_rate       = mean(Y),
    n_treated        = sum(A == 1),
    n_control        = sum(A == 0),

    # Crude
    crude_est        = crude$estimate,
    crude_se         = crude$se,
    crude_ci_lo      = crude$ci_low,
    crude_ci_hi      = crude$ci_high,

    # Adjusted regression
    regadj_est       = reg_adj$estimate,
    regadj_se        = reg_adj$se,
    regadj_ci_lo     = reg_adj$ci_low,
    regadj_ci_hi     = reg_adj$ci_high,

    # IPTW
    iptw_est         = iptw$estimate,
    iptw_se          = iptw$se,
    iptw_ci_lo       = iptw$ci_low,
    iptw_ci_hi       = iptw$ci_high,

    # PS-matched regression
    psmatch_est      = ps_match$estimate,
    psmatch_se       = ps_match$se,
    psmatch_ci_lo    = ps_match$ci_low,
    psmatch_ci_hi    = ps_match$ci_high,
    psmatch_n_matched = if (!is.null(ps_match$n_matched)) ps_match$n_matched else NA_integer_,

    # TMLE (IC-based)
    tmle_est         = tmle_res$estimate,
    tmle_se_ic       = tmle_res$se,
    tmle_ci_ic_lo    = tmle_res$ci_low,
    tmle_ci_ic_hi    = tmle_res$ci_high,
    tmle_ok          = tmle_res$ok,

    # TMLE targeted bootstrap
    tmle_se_tboot    = tboot$se,
    tmle_ci_tboot_lo = tboot$ci_low,
    tmle_ci_tboot_hi = tboot$ci_high,

    # TMLE nonparametric bootstrap
    tmle_se_npboot   = npboot$se,
    tmle_ci_npboot_lo = npboot$ci_low,
    tmle_ci_npboot_hi = npboot$ci_high,

    stringsAsFactors = FALSE
  )

  row
}

# Run replicates (sequential or parallel)
if (params$workers > 1 && has_future) {
  future::plan(future::multisession, workers = params$workers)
  results_list <- future.apply::future_lapply(
    seq_len(params$reps), run_one_rep,
    p = params, truth_rd = truth$RD_true,
    future.seed = TRUE
  )
  future::plan(future::sequential)
} else {
  results_list <- vector("list", params$reps)
  for (r in seq_len(params$reps)) {
    results_list[[r]] <- run_one_rep(r, p = params, truth_rd = truth$RD_true)

    if (params$save_per_rep) {
      saveRDS(results_list[[r]],
              file.path(params$output_dir, "per_rep",
                        sprintf("rep_%03d.rds", r)))
    }

    if (r %% 25 == 0 || r == 1) {
      rr <- results_list[[r]]
      message(sprintf("  Rep %d/%d  overlap=%s  TMLE=%.4f  crude=%.4f  regadj=%.4f  IPTW=%.4f  psmatch=%s",
                      r, params$reps, rr$overlap_flag,
                      ifelse(is.na(rr$tmle_est), NA, rr$tmle_est),
                      ifelse(is.na(rr$crude_est), NA, rr$crude_est),
                      ifelse(is.na(rr$regadj_est), NA, rr$regadj_est),
                      ifelse(is.na(rr$iptw_est), NA, rr$iptw_est),
                      ifelse(is.na(rr$psmatch_est), "NA",
                             sprintf("%.4f", rr$psmatch_est))))
    }
  }
}

res_df <- do.call(rbind, results_list)


# ===========================================================================
# Aggregate results
# ===========================================================================

compute_metrics <- function(df, true_rd, label = "all") {
  make_row <- function(method, est, ci_lo, ci_hi, se_col = NULL) {
    valid <- !is.na(est)
    if (sum(valid) < 3) {
      return(data.frame(
        subset = label, method = method, n_reps = sum(valid),
        bias = NA, rmse = NA, emp_sd = NA, mean_se = NA,
        coverage = NA, mcse_coverage = NA,
        stringsAsFactors = FALSE
      ))
    }
    e <- est[valid]
    bias <- mean(e) - true_rd
    rmse <- sqrt(mean((e - true_rd)^2))
    emp_sd <- sd(e)
    mean_se <- if (!is.null(se_col) && any(!is.na(se_col[valid]))) {
      mean(se_col[valid], na.rm = TRUE)
    } else NA
    covers <- !is.na(ci_lo[valid]) & !is.na(ci_hi[valid]) &
              (ci_lo[valid] <= true_rd) & (true_rd <= ci_hi[valid])
    coverage <- mean(covers)
    mcse_cov <- sqrt(coverage * (1 - coverage) / sum(valid))
    data.frame(
      subset = label, method = method, n_reps = sum(valid),
      bias = round(bias, 6), rmse = round(rmse, 6),
      emp_sd = round(emp_sd, 6), mean_se = round(mean_se, 6),
      coverage = round(coverage, 4), mcse_coverage = round(mcse_cov, 4),
      stringsAsFactors = FALSE
    )
  }

  tmle_ok <- df[!is.na(df$tmle_ok) & df$tmle_ok == TRUE, , drop = FALSE]

  rbind(
    make_row("Crude",            df$crude_est,  df$crude_ci_lo,  df$crude_ci_hi,  df$crude_se),
    make_row("Adj. Regression",  df$regadj_est, df$regadj_ci_lo, df$regadj_ci_hi, df$regadj_se),
    make_row("IPTW",             df$iptw_est,   df$iptw_ci_lo,   df$iptw_ci_hi,   df$iptw_se),
    make_row("PS-Matched Reg.",  df$psmatch_est, df$psmatch_ci_lo, df$psmatch_ci_hi, df$psmatch_se),
    make_row("TMLE (IC)",        tmle_ok$tmle_est, tmle_ok$tmle_ci_ic_lo, tmle_ok$tmle_ci_ic_hi, tmle_ok$tmle_se_ic),
    make_row("TMLE (tboot)",     tmle_ok$tmle_est, tmle_ok$tmle_ci_tboot_lo, tmle_ok$tmle_ci_tboot_hi, tmle_ok$tmle_se_tboot),
    make_row("TMLE (npboot)",    tmle_ok$tmle_est, tmle_ok$tmle_ci_npboot_lo, tmle_ok$tmle_ci_npboot_hi, tmle_ok$tmle_se_npboot)
  )
}

# All replicates
metrics_all <- compute_metrics(res_df, truth$RD_true, label = "all")

# Conditional on overlap_flag == "bad"
res_bad <- res_df[res_df$overlap_flag == "bad", , drop = FALSE]
metrics_bad <- if (nrow(res_bad) > 0) {
  compute_metrics(res_bad, truth$RD_true, label = "overlap_bad")
} else {
  NULL
}

summary_metrics <- rbind(metrics_all, metrics_bad)


# ===========================================================================
# A10. Clean-room gating demonstration (stop/go)
# ===========================================================================

message("\n--- Stop/Go Gating Demonstration ---")

# Overall overlap gate
n_flagged <- sum(res_df$overlap_flag == "bad")
pct_flagged <- round(100 * n_flagged / params$reps, 1)
message("  Replicates flagged 'bad' overlap: ", n_flagged, "/", params$reps,
        " (", pct_flagged, "%)")

# Operating characteristics by flag status
gate_results <- data.frame(
  overlap_flag = c("bad", "ok"),
  stringsAsFactors = FALSE
)
for (flag_val in c("bad", "ok")) {
  sub <- res_df[res_df$overlap_flag == flag_val & !is.na(res_df$tmle_ok) & res_df$tmle_ok, ]
  if (nrow(sub) > 5) {
    bias <- mean(sub$tmle_est, na.rm = TRUE) - truth$RD_true
    # Coverage using targeted bootstrap CI if available
    if (any(!is.na(sub$tmle_ci_tboot_lo))) {
      covers_tboot <- !is.na(sub$tmle_ci_tboot_lo) & !is.na(sub$tmle_ci_tboot_hi) &
                      (sub$tmle_ci_tboot_lo <= truth$RD_true) &
                      (truth$RD_true <= sub$tmle_ci_tboot_hi)
      cov_tboot <- mean(covers_tboot)
    } else {
      cov_tboot <- NA
    }
    covers_ic <- !is.na(sub$tmle_ci_ic_lo) & !is.na(sub$tmle_ci_ic_hi) &
                 (sub$tmle_ci_ic_lo <= truth$RD_true) &
                 (truth$RD_true <= sub$tmle_ci_ic_hi)
    cov_ic <- mean(covers_ic)

    acceptable <- (abs(bias) < 0.01) && (!is.na(cov_tboot) && cov_tboot >= 0.93)
    gate_results$n_reps[gate_results$overlap_flag == flag_val] <- nrow(sub)
    gate_results$tmle_bias[gate_results$overlap_flag == flag_val] <- round(bias, 5)
    gate_results$tmle_cov_ic[gate_results$overlap_flag == flag_val] <- round(cov_ic, 3)
    gate_results$tmle_cov_tboot[gate_results$overlap_flag == flag_val] <- round(cov_tboot, 3)
    gate_results$acceptable[gate_results$overlap_flag == flag_val] <- acceptable

    message(sprintf("  [%s] n=%d, TMLE bias=%.5f, IC cov=%.3f, tboot cov=%s, acceptable=%s",
                    flag_val, nrow(sub), bias, cov_ic,
                    ifelse(is.na(cov_tboot), "NA", sprintf("%.3f", cov_tboot)),
                    acceptable))
  } else {
    gate_results$n_reps[gate_results$overlap_flag == flag_val] <- nrow(sub)
    gate_results$tmle_bias[gate_results$overlap_flag == flag_val] <- NA
    gate_results$tmle_cov_ic[gate_results$overlap_flag == flag_val] <- NA
    gate_results$tmle_cov_tboot[gate_results$overlap_flag == flag_val] <- NA
    gate_results$acceptable[gate_results$overlap_flag == flag_val] <- NA
    message(sprintf("  [%s] n=%d — too few replicates", flag_val, nrow(sub)))
  }
}


# ===========================================================================
# Print results
# ===========================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("SCENARIO 1 v2 RESULTS: Bad Overlap, TMLE OK")
message(paste(rep("=", 70), collapse = ""))
message("\nTrue RD: ", round(truth$RD_true, 6))
message("Replicates: ", params$reps, " (N=", params$n_main, ")")
message("Overlap strength: ", params$strength_overlap)
message("W4 coef in g: ", params$coef_W4_g,
        "  |  W4 coef in Q: ", params$coef_W4_Q)
message("\nOverlap diagnostics across replicates:")
message("  Mean fraction extreme PS: ", round(mean(res_df$fraction_extreme), 3))
message("  Mean ESS/N: ", round(mean(res_df$ess_frac), 3))
message("  Mean max weight: ", round(mean(res_df$max_w), 1))
message("  Fraction flagged 'bad': ", round(mean(res_df$overlap_flag == "bad"), 3))

message("\n--- Estimator Performance (all replicates) ---")
print(metrics_all)

if (!is.null(metrics_bad) && nrow(metrics_bad) > 0) {
  message("\n--- Estimator Performance (overlap_flag == 'bad' only) ---")
  print(metrics_bad)
}


# ===========================================================================
# A8. Plots
# ===========================================================================

message("\n--- Generating figures ---")
fig_dir <- file.path(params$output_dir, "figures")

# -- PS distribution by treatment group --
tryCatch({
  ps_data <- do.call(rbind, lapply(seq_len(min(5, params$reps)), function(r) {
    sim <- sim_data(n = params$n_main, seed = params$master_seed + r, p = params)
    ps <- estimate_ps(sim$data, params)
    data.frame(g_hat = ps$g_raw, g_bounded = ps$g_bounded,
               A = factor(sim$data$A, labels = c("Control", "Treated")))
  }))
  p_ps <- ggplot(ps_data, aes(x = g_hat, fill = A)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    geom_vline(xintercept = c(params$ps_clip_hat[1], params$ps_clip_hat[2]),
               linetype = "dashed", color = "red") +
    labs(title = "Estimated propensity score distribution by treatment",
         subtitle = paste0("First ", min(5, params$reps), " replicates pooled"),
         x = "Estimated PS", y = "Count", fill = "Group") +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "ps_distribution.png"), p_ps, width = 8, height = 5)
}, error = function(e) message("  Warning: could not create PS plot: ", e$message))

# -- Sampling distributions --
tryCatch({
  est_long <- data.frame(
    method = rep(c("Crude", "Adj. Regression", "IPTW", "PS-Matched", "TMLE"),
                 each = nrow(res_df)),
    estimate = c(res_df$crude_est, res_df$regadj_est, res_df$iptw_est,
                 res_df$psmatch_est, res_df$tmle_est),
    stringsAsFactors = FALSE
  )
  est_long <- est_long[!is.na(est_long$estimate), ]
  est_long$method <- factor(est_long$method,
    levels = c("TMLE", "Adj. Regression", "PS-Matched", "IPTW", "Crude"))

  p_dist <- ggplot(est_long, aes(x = estimate, fill = method)) +
    geom_histogram(bins = 40, alpha = 0.7) +
    geom_vline(xintercept = truth$RD_true, linewidth = 0.8) +
    facet_wrap(~ method, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Sampling distributions of RD estimators",
         subtitle = paste0("True RD = ", round(truth$RD_true, 4)),
         x = "Estimated RD", y = "Count") +
    guides(fill = "none") +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "sampling_distributions.png"), p_dist,
         width = 8, height = 10)
}, error = function(e) message("  Warning: could not create sampling dist plot: ", e$message))

# -- Overlap severity vs error --
tryCatch({
  tmle_ok <- res_df[!is.na(res_df$tmle_ok) & res_df$tmle_ok, ]
  if (nrow(tmle_ok) > 10) {
    err_df <- data.frame(
      fraction_extreme = rep(tmle_ok$fraction_extreme, 3),
      abs_error = c(
        abs(tmle_ok$tmle_est - truth$RD_true),
        abs(tmle_ok$iptw_est - truth$RD_true),
        abs(tmle_ok$crude_est - truth$RD_true)
      ),
      method = rep(c("TMLE", "IPTW", "Crude"), each = nrow(tmle_ok)),
      stringsAsFactors = FALSE
    )
    err_df$method <- factor(err_df$method, levels = c("TMLE", "IPTW", "Crude"))

    p_err <- ggplot(err_df, aes(x = fraction_extreme, y = abs_error, color = method)) +
      geom_point(alpha = 0.3, size = 1.5) +
      geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "Absolute error vs overlap severity",
           x = "Fraction of observations with extreme PS",
           y = "|Estimate - Truth|",
           color = "Estimator") +
      theme_minimal(base_size = 12)
    ggsave(file.path(fig_dir, "overlap_vs_error.png"), p_err, width = 8, height = 5)
  }
}, error = function(e) message("  Warning: could not create overlap-error plot: ", e$message))

# -- Coverage vs fraction extreme (TMLE IC vs tboot) --
tryCatch({
  tmle_ok <- res_df[!is.na(res_df$tmle_ok) & res_df$tmle_ok, ]
  if (nrow(tmle_ok) > 20) {
    tmle_ok$covers_ic <- (!is.na(tmle_ok$tmle_ci_ic_lo) &
                          tmle_ok$tmle_ci_ic_lo <= truth$RD_true &
                          truth$RD_true <= tmle_ok$tmle_ci_ic_hi)
    tmle_ok$covers_tboot <- (!is.na(tmle_ok$tmle_ci_tboot_lo) &
                             tmle_ok$tmle_ci_tboot_lo <= truth$RD_true &
                             truth$RD_true <= tmle_ok$tmle_ci_tboot_hi)

    # Bin by deciles of fraction_extreme
    tmle_ok$fe_bin <- cut(tmle_ok$fraction_extreme,
                          breaks = quantile(tmle_ok$fraction_extreme,
                                            probs = seq(0, 1, 0.2), na.rm = TRUE),
                          include.lowest = TRUE)
    if (any(!is.na(tmle_ok$fe_bin))) {
      cov_by_bin <- aggregate(
        cbind(covers_ic, covers_tboot) ~ fe_bin,
        data = tmle_ok, FUN = mean, na.rm = TRUE
      )
      cov_long <- data.frame(
        fe_bin = rep(cov_by_bin$fe_bin, 2),
        coverage = c(cov_by_bin$covers_ic, cov_by_bin$covers_tboot),
        CI_type = rep(c("IC Wald", "Targeted Bootstrap"), each = nrow(cov_by_bin)),
        stringsAsFactors = FALSE
      )
      p_cov <- ggplot(cov_long, aes(x = fe_bin, y = coverage,
                                     fill = CI_type, group = CI_type)) +
        geom_col(position = "dodge", alpha = 0.8) +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
        scale_fill_manual(values = c("IC Wald" = "steelblue",
                                     "Targeted Bootstrap" = "darkorange")) +
        labs(title = "TMLE coverage by overlap severity",
             subtitle = "IC-based vs targeted bootstrap CIs",
             x = "PS extremity quintile", y = "Coverage",
             fill = "CI method") +
        ylim(0, 1) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(fig_dir, "coverage_by_overlap.png"), p_cov,
             width = 8, height = 5)
    }
  }
}, error = function(e) message("  Warning: could not create coverage plot: ", e$message))


# ===========================================================================
# Save outputs
# ===========================================================================

message("\n--- Saving outputs ---")

saveRDS(params, file.path(params$output_dir, "params_used.rds"))
saveRDS(res_df, file.path(params$output_dir, "replicate_results.rds"))
write.csv(summary_metrics, file.path(params$output_dir, "summary_metrics.csv"),
          row.names = FALSE)
saveRDS(gate_results, file.path(params$output_dir, "gate_results.rds"))

# Session info
si <- capture.output(sessionInfo())
writeLines(si, file.path(params$output_dir, "sessionInfo.txt"))

# Results summary markdown
summary_md <- file.path(params$output_dir, "results_summary.md")
sink(summary_md)
cat("# Scenario 1 v2: Bad Overlap, TMLE OK — Results Summary\n\n")
cat(sprintf("**True RD:** %.6f\n\n", truth$RD_true))
cat(sprintf("**Replicates:** %d (N = %d)\n\n", params$reps, params$n_main))
cat("## Estimator Performance (all replicates)\n\n")
cat("| Method | Bias | RMSE | Emp SD | Mean SE | Coverage | MCSE Cov |\n")
cat("|--------|------|------|--------|---------|----------|----------|\n")
for (i in seq_len(nrow(metrics_all))) {
  m <- metrics_all[i, ]
  cat(sprintf("| %s | %.4f | %.4f | %.4f | %s | %s | %s |\n",
              m$method, m$bias, m$rmse, m$emp_sd,
              ifelse(is.na(m$mean_se), "—", sprintf("%.4f", m$mean_se)),
              ifelse(is.na(m$coverage), "—", sprintf("%.2f%%", m$coverage * 100)),
              ifelse(is.na(m$mcse_coverage), "—", sprintf("%.4f", m$mcse_coverage))))
}
cat("\n## Stop/Go Gate\n\n")
cat(sprintf("Replicates flagged 'bad' overlap: %d/%d (%.1f%%)\n\n",
            n_flagged, params$reps, pct_flagged))
cat("| Flag | N | TMLE Bias | IC Cov | tboot Cov | Acceptable |\n")
cat("|------|---|-----------|--------|-----------|------------|\n")
for (i in seq_len(nrow(gate_results))) {
  g <- gate_results[i, ]
  cat(sprintf("| %s | %s | %s | %s | %s | %s |\n",
              g$overlap_flag,
              ifelse(is.na(g$n_reps), "—", as.character(g$n_reps)),
              ifelse(is.na(g$tmle_bias), "—", sprintf("%.5f", g$tmle_bias)),
              ifelse(is.na(g$tmle_cov_ic), "—", sprintf("%.3f", g$tmle_cov_ic)),
              ifelse(is.na(g$tmle_cov_tboot), "—", sprintf("%.3f", g$tmle_cov_tboot)),
              ifelse(is.na(g$acceptable), "—", as.character(g$acceptable))))
}
sink()

message("\nOutputs saved to: ", params$output_dir)
message("  truth.rds, params_used.rds, replicate_results.rds")
message("  summary_metrics.csv, gate_results.rds, results_summary.md")
message("  figures/: ps_distribution, sampling_distributions, overlap_vs_error, coverage_by_overlap")
message("\nDone: ", Sys.time())
