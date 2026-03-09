# ===========================================================================
# R/sim_helpers_v3.R
# Shared helpers for Scenario 1 v3: Mixed Overlap, Nonlinear Q, CV-TMLE
# ===========================================================================
#
# Provides: DGP, truth computation, overlap diagnostics, and all estimator
# wrappers for the v3 simulation. Sourced by the main driver script.
#
# Estimator return contract:
#   list(estimate, se, ci_low, ci_high, success, warnings, ...)
#   - success: TRUE if estimate and CI are usable; FALSE otherwise
#   - CI may be NA even when success=TRUE (e.g., bootstrap disabled)
# ===========================================================================

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

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

# Check for optional packages once
.has_tmle    <- requireNamespace("tmle", quietly = TRUE)
.has_SL      <- requireNamespace("SuperLearner", quietly = TRUE)
.has_MatchIt <- requireNamespace("MatchIt", quietly = TRUE)
.has_ranger  <- requireNamespace("ranger", quietly = TRUE)
.has_gam     <- requireNamespace("gam", quietly = TRUE)
.has_glmnet  <- requireNamespace("glmnet", quietly = TRUE)
.has_digest  <- requireNamespace("digest", quietly = TRUE)

# ---------------------------------------------------------------------------
# A1–A3. Data-generating process
# ---------------------------------------------------------------------------

#' Generate one replicate of Scenario 1 v3 data
#'
#' Supports replicate-level overlap regime (ok/bad via Z_r) and optional
#' subject-level W4 mixture, with nonlinear outcome model + interactions.
#'
#' @param n sample size
#' @param seed replicate-specific seed
#' @param p parameter list (must contain DGP coefficients)
#' @param overlap_regime "ok" or "bad" (overrides Z_r draw if supplied)
#' @return list with $data, $g_true, $Q_true, $overlap_regime
sim_data_v3 <- function(n, seed = NULL, p, overlap_regime = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # --- Overlap regime ---
  if (is.null(overlap_regime)) {
    Z_r <- rbinom(1, 1, p$p_bad)
    overlap_regime <- ifelse(Z_r == 1, "bad", "ok")
  }
  s_g <- ifelse(overlap_regime == "bad", p$s_bad, p$s_ok)

  # --- Covariates ---
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rnorm(n, 0, 1)

  # W4: optional subject-level mixture
 if (p$use_subject_mixture) {
    S_i <- rbinom(n, 1, p$p_hiVar)
    W4 <- ifelse(S_i == 1,
                 rnorm(n, 0, p$sigma_hi),
                 rnorm(n, 0, p$sigma_lo))
  } else {
    W4 <- rnorm(n, 0, 1)
  }

  # --- Treatment mechanism ---
  lp_g_raw <- s_g * (
    p$alpha1_g * W1 + p$alpha2_g * W2 +
    p$alpha3_g * W3 + p$alpha4_g * W4
  )
  alpha0_g <- -mean(lp_g_raw) + p$g_intercept_shift  # shift P(A=1)
  lp_g <- alpha0_g + lp_g_raw
  g_true <- expit(lp_g)
  g_true <- pmin(pmax(g_true, p$ps_clip_true[1]), p$ps_clip_true[2])
  A <- rbinom(n, 1, g_true)

  # --- Outcome model (nonlinear + interactions) ---
  h1 <- 0.6 * W1 + 0.4 * W1^2
  h2 <- 0.7 * W2                       # binary, so W2>0 ≡ W2
  h3 <- 0.5 * abs(W3)
  h4 <- p$beta4_Q * W4

  lp_Q <- p$beta0_Q + p$betaA_Q * A +
    h1 + h2 + h3 + h4 +
    p$tau_int1 * A * (W1 > 0) +
    p$tau_int2 * A * W2 +              # W2^2 = W2 for binary
    p$tau_int3 * A * W4

  Q_true <- expit(lp_Q)
  Y <- rbinom(n, 1, Q_true)

  dat <- data.frame(W1 = W1, W2 = W2, W3 = W3, W4 = W4, A = A, Y = Y)
  list(data = dat, g_true = g_true, Q_true = Q_true,
       overlap_regime = overlap_regime, s_g = s_g)
}

# ---------------------------------------------------------------------------
# B. Truth computation
# ---------------------------------------------------------------------------

#' Compute RD truth via exact plug-in on large sample.
#' Truth depends only on Q0 and the marginal W distribution, not on g0.
compute_truth_v3 <- function(p) {
  set.seed(p$master_seed)
  n <- p$n_truth

  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rnorm(n, 0, 1)

  if (p$use_subject_mixture) {
    S_i <- rbinom(n, 1, p$p_hiVar)
    W4 <- ifelse(S_i == 1,
                 rnorm(n, 0, p$sigma_hi),
                 rnorm(n, 0, p$sigma_lo))
  } else {
    W4 <- rnorm(n, 0, 1)
  }

  h1 <- 0.6 * W1 + 0.4 * W1^2
  h2 <- 0.7 * W2
  h3 <- 0.5 * abs(W3)
  h4 <- p$beta4_Q * W4

  lp_base <- p$beta0_Q + h1 + h2 + h3 + h4

  lp_Q1 <- lp_base + p$betaA_Q +
    p$tau_int1 * (W1 > 0) + p$tau_int2 * W2 + p$tau_int3 * W4
  lp_Q0 <- lp_base  # A=0 → no A terms

  EY1 <- mean(expit(lp_Q1))
  EY0 <- mean(expit(lp_Q0))

  list(EY1 = EY1, EY0 = EY0, RD_true = EY1 - EY0, truth_N = n,
       note = "Exact plug-in (no MC noise in Y); depends on Q0 and W dist, not g0")
}

# ---------------------------------------------------------------------------
# D. PS estimation and overlap diagnostics
# ---------------------------------------------------------------------------

estimate_ps_v3 <- function(dat, p) {
  g_fit <- glm(A ~ W1 + W2 + W3 + W4, data = dat, family = "binomial")
  g_raw <- fitted(g_fit)
  g_bounded <- pmin(pmax(g_raw, p$ps_clip_hat[1]), p$ps_clip_hat[2])
  list(g_raw = g_raw, g_bounded = g_bounded, fit = g_fit)
}

overlap_diagnostics_v3 <- function(g_raw, g_bounded, A, p) {
  thr <- p$overlap_flag_thresholds
  n <- length(A)

  frac_extreme <- mean(g_raw < thr$extreme_ps | g_raw > (1 - thr$extreme_ps))
  w_iptw <- A / g_bounded + (1 - A) / (1 - g_bounded)
  ess <- sum(w_iptw)^2 / sum(w_iptw^2)
  ess_frac <- ess / n
  max_w <- max(w_iptw)

  flag <- (frac_extreme > thr$extreme_prop) ||
          (ess_frac < thr$ess_frac) ||
          (max_w > thr$max_w)

  list(fraction_extreme = frac_extreme,
       min_g_raw = min(g_raw), max_g_raw = max(g_raw),
       min_g_bounded = min(g_bounded), max_g_bounded = max(g_bounded),
       ess = ess, ess_frac = ess_frac, max_w = max_w,
       overlap_flag = ifelse(flag, "bad", "ok"),
       w_iptw = w_iptw)
}

# ---------------------------------------------------------------------------
# C1. Crude RD
# ---------------------------------------------------------------------------

est_crude <- function(Y, A) {
  tryCatch({
    n1 <- sum(A == 1); n0 <- sum(A == 0)
    if (n1 < 2 || n0 < 2)
      return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                  success = FALSE, warnings = "too few in one arm"))
    p1 <- mean(Y[A == 1]); p0 <- mean(Y[A == 0])
    rd <- p1 - p0
    se <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    list(estimate = rd, se = se,
         ci_low = rd - 1.96 * se, ci_high = rd + 1.96 * se,
         success = TRUE, warnings = NA_character_)
  }, error = function(e)
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         success = FALSE, warnings = conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# C1. Adjusted regression (misspecified GLM, bootstrap CI)
# ---------------------------------------------------------------------------

est_reg_adj_v3 <- function(dat, B_boot = 0) {
  tryCatch({
    # Intentionally misspecified: no W1^2, abs(W3), or interactions
    fit <- glm(Y ~ A + W1 + W2 + W3 + W4, data = dat, family = "binomial")
    d1 <- d0 <- dat; d1$A <- 1; d0$A <- 0
    p1_i <- predict(fit, d1, type = "response")
    p0_i <- predict(fit, d0, type = "response")
    rd <- mean(p1_i) - mean(p0_i)

    se <- NA_real_; ci <- c(NA_real_, NA_real_)
    if (B_boot > 0) {
      # Nonparametric bootstrap (optional, slow)
      n <- nrow(dat)
      boot_rd <- numeric(B_boot)
      for (b in seq_len(B_boot)) {
        idx <- sample(n, n, replace = TRUE)
        d_b <- dat[idx, , drop = FALSE]
        fit_b <- tryCatch(
          glm(Y ~ A + W1 + W2 + W3 + W4, data = d_b, family = "binomial"),
          error = function(e) NULL)
        if (is.null(fit_b)) { boot_rd[b] <- NA; next }
        d1b <- d0b <- d_b; d1b$A <- 1; d0b$A <- 0
        boot_rd[b] <- mean(predict(fit_b, d1b, type = "response")) -
                      mean(predict(fit_b, d0b, type = "response"))
      }
      ok <- boot_rd[is.finite(boot_rd)]
      se <- if (length(ok) >= 10) sd(ok) else NA_real_
      ci <- if (length(ok) >= 20) unname(quantile(ok, c(0.025, 0.975)))
            else c(NA_real_, NA_real_)
    } else {
      # Analytic SE: SD of unit-level g-computation RD / sqrt(n)
      n <- nrow(dat)
      rd_i <- p1_i - p0_i
      se <- sd(rd_i) / sqrt(n)
      ci <- rd + c(-1.96, 1.96) * se
    }
    list(estimate = rd, se = se, ci_low = ci[1], ci_high = ci[2],
         success = TRUE, warnings = NA_character_)
  }, error = function(e)
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         success = FALSE, warnings = conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# C1. IPTW (Hajek, bootstrap CI)
# ---------------------------------------------------------------------------

est_iptw_v3 <- function(Y, A, g_bounded, dat, B_boot = 0, p = NULL) {
  tryCatch({
    n <- length(Y)
    w <- A / g_bounded + (1 - A) / (1 - g_bounded)
    d1 <- sum(w * A); d0 <- sum(w * (1 - A))
    if (d1 < 1e-10 || d0 < 1e-10)
      return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                  success = FALSE, warnings = "zero denominator"))
    mu1 <- sum(w * A * Y) / d1
    mu0 <- sum(w * (1 - A) * Y) / d0
    rd <- mu1 - mu0

    se <- NA_real_; ci <- c(NA_real_, NA_real_)
    if (B_boot > 0 && !is.null(p)) {
      # Nonparametric bootstrap (optional, slow)
      boot_rd <- numeric(B_boot)
      for (b in seq_len(B_boot)) {
        idx <- sample(n, n, replace = TRUE)
        d_b <- dat[idx, , drop = FALSE]
        g_b <- tryCatch({
          pmin(pmax(fitted(glm(A ~ W1 + W2 + W3 + W4, data = d_b,
                               family = "binomial")),
                    p$ps_clip_hat[1]), p$ps_clip_hat[2])
        }, error = function(e) rep(0.5, n))
        wb <- d_b$A / g_b + (1 - d_b$A) / (1 - g_b)
        d1b <- sum(wb * d_b$A); d0b <- sum(wb * (1 - d_b$A))
        boot_rd[b] <- if (d1b > 1e-10 && d0b > 1e-10)
          sum(wb * d_b$A * d_b$Y) / d1b - sum(wb * (1 - d_b$A) * d_b$Y) / d0b
        else NA
      }
      ok <- boot_rd[is.finite(boot_rd)]
      se <- if (length(ok) >= 10) sd(ok) else NA_real_
      ci <- if (length(ok) >= 20) unname(quantile(ok, c(0.025, 0.975)))
            else c(NA_real_, NA_real_)
    } else {
      # Analytic Hajek SE: influence-function based
      phi1 <- w * A * (Y - mu1) / d1 * n
      phi0 <- w * (1 - A) * (Y - mu0) / d0 * n
      phi  <- phi1 - phi0
      se <- sqrt(sum(phi^2) / n^2)
      ci <- rd + c(-1.96, 1.96) * se
    }
    list(estimate = rd, se = se, ci_low = ci[1], ci_high = ci[2],
         success = TRUE, warnings = NA_character_)
  }, error = function(e)
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         success = FALSE, warnings = conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# C2. PS-matched regression (pairs bootstrap)
# ---------------------------------------------------------------------------

est_psmatch_reg_v3 <- function(dat, g_bounded, B_boot = 200, p = NULL) {
  if (!.has_MatchIt)
    return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                success = FALSE, n_matched = NA, smd_before = NA,
                smd_after = NA, warnings = "MatchIt not installed"))
  tryCatch({
    dat$logit_ps <- logit(g_bounded)
    cal <- if (!is.null(p$match_caliper)) p$match_caliper
           else 0.2 * sd(dat$logit_ps)

    m <- MatchIt::matchit(A ~ W1 + W2 + W3 + W4, data = dat,
                          method = "nearest", distance = dat$logit_ps,
                          caliper = cal, replace = FALSE)
    mdat <- MatchIt::match.data(m)
    n_matched <- nrow(mdat)

    if (n_matched < 10)
      return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                  success = FALSE, n_matched = n_matched,
                  smd_before = NA, smd_after = NA,
                  warnings = "too few matched"))

    # SMD before/after
    covs <- c("W1", "W2", "W3", "W4")
    smd_before <- sapply(covs, function(v) {
      m1 <- mean(dat[[v]][dat$A == 1]); m0 <- mean(dat[[v]][dat$A == 0])
      s <- sqrt((var(dat[[v]][dat$A == 1]) + var(dat[[v]][dat$A == 0])) / 2)
      if (s < 1e-10) 0 else (m1 - m0) / s
    })
    smd_after <- sapply(covs, function(v) {
      m1 <- mean(mdat[[v]][mdat$A == 1]); m0 <- mean(mdat[[v]][mdat$A == 0])
      s <- sqrt((var(mdat[[v]][mdat$A == 1]) + var(mdat[[v]][mdat$A == 0])) / 2)
      if (s < 1e-10) 0 else (m1 - m0) / s
    })

    # G-computation on matched data (misspecified model)
    fit_m <- glm(Y ~ A + W1 + W2 + W3 + W4, data = mdat, family = "binomial")
    d1 <- d0 <- mdat; d1$A <- 1; d0$A <- 0
    rd <- mean(predict(fit_m, d1, type = "response")) -
          mean(predict(fit_m, d0, type = "response"))

    # Pairs bootstrap
    se <- NA_real_; ci <- c(NA_real_, NA_real_)
    if (B_boot > 0 && "subclass" %in% names(mdat)) {
      pairs <- unique(mdat$subclass)
      n_pairs <- length(pairs)
      if (n_pairs >= 10) {
        boot_rd <- numeric(B_boot)
        for (b in seq_len(B_boot)) {
          bp <- sample(pairs, n_pairs, replace = TRUE)
          bdat <- do.call(rbind, lapply(bp, function(pr)
            mdat[mdat$subclass == pr, , drop = FALSE]))
          fit_b <- tryCatch(
            glm(Y ~ A + W1 + W2 + W3 + W4, data = bdat, family = "binomial"),
            error = function(e) NULL)
          if (is.null(fit_b)) { boot_rd[b] <- NA; next }
          d1b <- d0b <- bdat; d1b$A <- 1; d0b$A <- 0
          boot_rd[b] <- mean(predict(fit_b, d1b, type = "response")) -
                        mean(predict(fit_b, d0b, type = "response"))
        }
        ok <- boot_rd[is.finite(boot_rd)]
        se <- if (length(ok) >= 10) sd(ok) else NA_real_
        ci <- if (length(ok) >= 20) unname(quantile(ok, c(0.025, 0.975)))
              else c(NA_real_, NA_real_)
      }
    }

    list(estimate = rd, se = se, ci_low = ci[1], ci_high = ci[2],
         success = TRUE, n_matched = n_matched,
         n_unmatched = nrow(dat) - n_matched,
         smd_before = smd_before, smd_after = smd_after,
         warnings = NA_character_)
  }, error = function(e)
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         success = FALSE, n_matched = NA, smd_before = NA,
         smd_after = NA, warnings = conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# C3. TMLE-GLM (misspecified Q and g, via tmle package)
# ---------------------------------------------------------------------------

est_tmle_glm <- function(Y, A, W, p) {
  if (!.has_tmle)
    return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                eif = NULL, success = FALSE, warnings = "tmle not installed"))
  tryCatch({
    fit <- tmle::tmle(Y = Y, A = A, W = W, family = "binomial",
                      Q.SL.library = "SL.glm",
                      g.SL.library = "SL.glm",
                      gbound = p$ps_clip_hat)
    ate <- fit$estimates$ATE
    eif <- fit$estimates$IC$IC.ATE

    list(estimate = ate$psi, se = sqrt(ate$var.psi),
         ci_low = ate$CI[1], ci_high = ate$CI[2],
         eif = eif, success = TRUE, warnings = NA_character_)
  }, error = function(e)
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         eif = NULL, success = FALSE, warnings = conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# C3. TMLE-ML (richer SuperLearner library)
# ---------------------------------------------------------------------------

est_tmle_ml <- function(Y, A, W, p) {
  if (!.has_tmle)
    return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                eif = NULL, success = FALSE, warnings = "tmle not installed"))
  tryCatch({
    Q_lib <- filter_sl_libs(p$Q_SL_library_ml)
    g_lib <- filter_sl_libs(p$g_SL_library_ml)

    fit <- tmle::tmle(Y = Y, A = A, W = W, family = "binomial",
                      Q.SL.library = Q_lib, g.SL.library = g_lib,
                      gbound = p$ps_clip_hat)
    ate <- fit$estimates$ATE
    eif <- fit$estimates$IC$IC.ATE

    list(estimate = ate$psi, se = sqrt(ate$var.psi),
         ci_low = ate$CI[1], ci_high = ate$CI[2],
         eif = eif, success = TRUE,
         Q_lib_used = Q_lib, g_lib_used = g_lib,
         warnings = NA_character_)
  }, error = function(e) {
    # Fallback to SL.glm only
    tryCatch({
      fit <- tmle::tmle(Y = Y, A = A, W = W, family = "binomial",
                        Q.SL.library = "SL.glm",
                        g.SL.library = "SL.glm",
                        gbound = p$ps_clip_hat)
      ate <- fit$estimates$ATE
      eif <- fit$estimates$IC$IC.ATE
      list(estimate = ate$psi, se = sqrt(ate$var.psi),
           ci_low = ate$CI[1], ci_high = ate$CI[2],
           eif = eif, success = TRUE,
           warnings = paste0("ML fallback to GLM: ", conditionMessage(e)))
    }, error = function(e2)
      list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
           eif = NULL, success = FALSE, warnings = conditionMessage(e2)))
  })
}

# ---------------------------------------------------------------------------
# C4. Cross-fitted TMLE with ML learners
# ---------------------------------------------------------------------------

#' Fit Q model (classification) on training data
#' Uses ranger if available, otherwise enriched GLM
.fit_Q_ml <- function(train_dat) {
  if (.has_ranger) {
    train_dat$Y_f <- factor(train_dat$Y)
    fit <- ranger::ranger(Y_f ~ A + W1 + W2 + W3 + W4,
                          data = train_dat, probability = TRUE,
                          num.trees = 100, min.node.size = 20)
    return(list(type = "ranger", fit = fit))
  }
  # Fallback: enriched GLM with polynomial/interaction terms
  fit <- glm(Y ~ A * (W1 + W2 + W3 + W4) + I(W1^2) + I(abs(W3)),
             data = train_dat, family = "binomial")
  list(type = "glm_enriched", fit = fit)
}

.predict_Q_ml <- function(model, newdata) {
  if (model$type == "ranger") {
    newdata$Y_f <- factor(0, levels = levels(model$fit$forest$levels))
    pred <- predict(model$fit, data = newdata)$predictions[, "1"]
  } else {
    pred <- predict(model$fit, newdata = newdata, type = "response")
  }
  pmin(pmax(pred, 1e-5), 1 - 1e-5)
}

.fit_g_ml <- function(train_dat) {
  if (.has_ranger) {
    train_dat$A_f <- factor(train_dat$A)
    fit <- ranger::ranger(A_f ~ W1 + W2 + W3 + W4,
                          data = train_dat, probability = TRUE,
                          num.trees = 100, min.node.size = 20)
    return(list(type = "ranger", fit = fit))
  }
  fit <- glm(A ~ W1 + W2 + W3 + W4 + I(W1^2) + I(W3^2),
             data = train_dat, family = "binomial")
  list(type = "glm_enriched", fit = fit)
}

.predict_g_ml <- function(model, newdata, clip) {
  if (model$type == "ranger") {
    newdata$A_f <- factor(0, levels = levels(model$fit$forest$levels))
    pred <- predict(model$fit, data = newdata)$predictions[, "1"]
  } else {
    pred <- predict(model$fit, newdata = newdata, type = "response")
  }
  pmin(pmax(pred, clip[1]), clip[2])
}

est_cvtmle_ml <- function(dat, p, V = 5) {
  tryCatch({
    n <- nrow(dat)
    W_cols <- c("W1", "W2", "W3", "W4")

    # Create folds
    folds <- sample(rep(1:V, length.out = n))

    Qbar_A <- numeric(n)     # Q(A_i, W_i) cross-fitted
    Qbar_1 <- numeric(n)     # Q(1, W_i) cross-fitted
    Qbar_0 <- numeric(n)     # Q(0, W_i) cross-fitted
    ghat   <- numeric(n)     # g(W_i) cross-fitted

    for (v in 1:V) {
      train_idx <- which(folds != v)
      valid_idx <- which(folds == v)
      train_d <- dat[train_idx, , drop = FALSE]
      valid_d <- dat[valid_idx, , drop = FALSE]

      # Fit Q on training
      Q_model <- .fit_Q_ml(train_d)

      # Predict on validation: Q(A,W), Q(1,W), Q(0,W)
      Qbar_A[valid_idx] <- .predict_Q_ml(Q_model, valid_d)
      d1v <- d0v <- valid_d; d1v$A <- 1; d0v$A <- 0
      Qbar_1[valid_idx] <- .predict_Q_ml(Q_model, d1v)
      Qbar_0[valid_idx] <- .predict_Q_ml(Q_model, d0v)

      # Fit g on training, predict on validation
      g_model <- .fit_g_ml(train_d)
      ghat[valid_idx] <- .predict_g_ml(g_model, valid_d, p$ps_clip_hat)
    }

    # Targeting step
    A <- dat$A; Y <- dat$Y
    H <- A / ghat - (1 - A) / (1 - ghat)
    H1 <- 1 / ghat                        # clever covariate under A=1
    H0 <- -1 / (1 - ghat)                 # clever covariate under A=0

    offset_A <- logit(Qbar_A)
    eps_fit <- glm(Y ~ -1 + H, offset = offset_A, family = "binomial")
    eps <- coef(eps_fit)

    # Updated predictions
    Qbar_1_star <- expit(logit(Qbar_1) + eps * H1)
    Qbar_0_star <- expit(logit(Qbar_0) + eps * H0)

    psi <- mean(Qbar_1_star - Qbar_0_star)

    # EIF
    Qbar_A_star <- expit(offset_A + eps * H)
    eif <- H * (Y - Qbar_A_star) + (Qbar_1_star - Qbar_0_star) - psi

    se <- sqrt(var(eif) / n)
    ci <- psi + c(-1.96, 1.96) * se

    list(estimate = psi, se = se, ci_low = ci[1], ci_high = ci[2],
         eif = eif, success = TRUE,
         Q_type = if (.has_ranger) "ranger" else "glm_enriched",
         warnings = NA_character_)
  }, error = function(e)
    list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
         eif = NULL, success = FALSE, warnings = conditionMessage(e)))
}

# ---------------------------------------------------------------------------
# C5. Multiplier bootstrap (Rademacher) for TMLE inference
# ---------------------------------------------------------------------------

#' Multiplier bootstrap for TMLE variance estimation.
#' Uses Rademacher multipliers: xi_i iid ∈ {-1, +1}, E=0, Var=1.
#' psi_b = psi_hat + (1/n) * sum(xi_i * EIF_i)
#'
#' @param psi_hat TMLE point estimate
#' @param eif EIF values D_i (length n)
#' @param B number of bootstrap draws
#' @return list with se, ci_low, ci_high (normal), ci_low_pct, ci_high_pct (percentile)
multiplier_bootstrap <- function(psi_hat, eif, B = 500) {
  if (is.null(eif) || length(eif) < 5)
    return(list(se = NA, ci_low = NA, ci_high = NA,
                ci_low_pct = NA, ci_high_pct = NA,
                success = FALSE, warnings = "insufficient EIF values"))

  n <- length(eif)
  psi_boot <- numeric(B)
  for (b in seq_len(B)) {
    xi <- sample(c(-1L, 1L), n, replace = TRUE)   # Rademacher
    psi_boot[b] <- psi_hat + mean(xi * eif)
  }
  se <- sd(psi_boot)
  ci_norm <- psi_hat + c(-1.96, 1.96) * se
  ci_pct  <- unname(quantile(psi_boot, c(0.025, 0.975)))

  list(se = se,
       ci_low = ci_norm[1], ci_high = ci_norm[2],        # normal
       ci_low_pct = ci_pct[1], ci_high_pct = ci_pct[2],  # percentile
       success = TRUE, warnings = NA_character_)
}

# ---------------------------------------------------------------------------
# E. Metrics computation (correct coverage accounting)
# ---------------------------------------------------------------------------

#' Compute summary metrics for one estimator.
#' Coverage is computed over successes with non-NA CIs only.
#' n_success and n_ci are reported separately.
compute_metrics_v3 <- function(est, ci_lo, ci_hi, success, true_rd,
                               se_col = NULL, label = "all", method = "") {
  n_total <- length(est)
  n_success <- sum(success, na.rm = TRUE)
  has_ci <- success & !is.na(ci_lo) & !is.na(ci_hi)
  n_ci <- sum(has_ci, na.rm = TRUE)

  if (n_success < 3) {
    return(data.frame(
      subset = label, method = method, n_total = n_total,
      n_success = n_success, n_ci = n_ci,
      bias = NA, rmse = NA, emp_sd = NA, mean_se = NA,
      coverage = NA, mcse_coverage = NA,
      stringsAsFactors = FALSE))
  }

  e <- est[success]
  bias <- mean(e) - true_rd
  rmse <- sqrt(mean((e - true_rd)^2))
  emp_sd <- sd(e)
  mean_se <- if (!is.null(se_col) && any(!is.na(se_col[success])))
    mean(se_col[success], na.rm = TRUE) else NA

  coverage <- NA_real_; mcse_cov <- NA_real_
  if (n_ci >= 5) {
    covers <- (ci_lo[has_ci] <= true_rd) & (true_rd <= ci_hi[has_ci])
    coverage <- mean(covers)
    mcse_cov <- sqrt(coverage * (1 - coverage) / n_ci)
  }

  data.frame(subset = label, method = method,
             n_total = n_total, n_success = n_success, n_ci = n_ci,
             bias = round(bias, 6), rmse = round(rmse, 6),
             emp_sd = round(emp_sd, 6),
             mean_se = round(mean_se, 6),
             coverage = round(coverage, 4),
             mcse_coverage = round(mcse_cov, 4),
             stringsAsFactors = FALSE)
}
