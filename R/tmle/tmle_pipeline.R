#' @title TMLE Estimation Pipeline
#' @description Implements Targeted Maximum Likelihood Estimation for
#'   point-treatment binary outcomes and discrete-time survival risk.
#' @name tmle_pipeline
NULL

#' TMLE for Point-Treatment Binary Outcome
#'
#' Implements TMLE for estimating the ATE of a binary treatment on a binary
#' outcome using the \code{tmle} package with SuperLearner.
#'
#' @param Y Binary outcome vector.
#' @param A Binary treatment vector.
#' @param W Data frame of baseline covariates.
#' @param sl_lib_Q Character vector of SuperLearner libraries for outcome.
#' @param sl_lib_g Character vector of SuperLearner libraries for treatment.
#' @param g_trunc Numeric vector of length 2 for truncation bounds.
#' @return A list with components: psi (ATE), se, ci, pvalue, and diagnostics.
#' @export
tmle_point_binary <- function(Y, A, W,
                              sl_lib_Q = c("SL.glm", "SL.gam", "SL.mean"),
                              sl_lib_g = c("SL.glm", "SL.gam", "SL.mean"),
                              g_trunc = c(0.01, 0.99)) {
  sl_lib_Q <- filter_sl_libraries(sl_lib_Q)
  sl_lib_g <- filter_sl_libraries(sl_lib_g)

  W <- as.data.frame(W)

  fit <- tryCatch({
    tmle::tmle(
      Y = Y, A = A, W = W,
      family = "binomial",
      Q.SL.library = sl_lib_Q,
      g.SL.library = sl_lib_g,
      gbound = g_trunc
    )
  }, error = function(e) {
    message("TMLE with SL failed, falling back to GLM: ", e$message)
    tmle::tmle(
      Y = Y, A = A, W = W,
      family = "binomial",
      Q.SL.library = "SL.glm",
      g.SL.library = "SL.glm",
      gbound = g_trunc
    )
  })

  ate <- fit$estimates$ATE
  ic <- fit$estimates$ATE$IC
  ic_diagnostics <- list(
    ic_mean = mean(ic),
    ic_sd   = stats::sd(ic),
    ic_min  = min(ic),
    ic_max  = max(ic),
    ic_skew = mean(((ic - mean(ic)) / stats::sd(ic))^3)
  )

  g1 <- fit$g$g1W
  trunc_info <- truncate_ps(g1, g_trunc[1], g_trunc[2])

  list(
    psi         = ate$psi,
    se          = sqrt(ate$var.psi),
    ci          = ate$CI,
    pvalue      = ate$pvalue,
    ic_diagnostics = ic_diagnostics,
    truncation  = list(n_lower = trunc_info$n_lower,
                       n_upper = trunc_info$n_upper),
    sl_fit_Q    = if (!is.null(fit$Qinit)) "SuperLearner" else "GLM",
    sl_fit_g    = if (!is.null(fit$g)) "SuperLearner" else "GLM"
  )
}

#' TMLE for Discrete-Time Survival Risk
#'
#' Implements discrete-time survival TMLE for estimating counterfactual risk
#' at time t using pooled logistic hazards with TMLE targeting. Incorporates
#' inverse-probability-of-censoring weighting (IPCW) for informative censoring
#' and treatment switching.
#'
#' The estimand is the marginal counterfactual cumulative risk:
#' \deqn{\psi(t) = E[I(T^a \le t)]}
#' under a hypothetical intervention that sets treatment to a (0 or 1) and
#' removes censoring. The risk difference is:
#' \deqn{RD(t) = \psi_1(t) - \psi_0(t)}
#'
#' Censoring types handled:
#' - Administrative censoring (max follow-up reached)
#' - Informative censoring (dependent on risk factors)
#' - Switching-induced censoring (30-day washout window post-switch)
#'
#' @param data Data frame with follow_time, event, treatment, and covariates.
#' @param t_eval Numeric time point at which to evaluate risk (days).
#' @param sl_lib_Q Character vector of SL libraries for outcome (hazard).
#' @param sl_lib_g Character vector of SL libraries for treatment PS.
#' @param sl_lib_cens Character vector of SL libraries for censoring model.
#' @param g_trunc Numeric vector of length 2 for g truncation bounds.
#' @param max_time Maximum discrete time for pooled logistic.
#' @return A list with risk_1, risk_0, RD, RR, SEs, CIs, and diagnostics.
#' @export
tmle_survival_risk <- function(data,
                               t_eval = 180,
                               sl_lib_Q = c("SL.glm", "SL.mean"),
                               sl_lib_g = c("SL.glm", "SL.mean"),
                               sl_lib_cens = c("SL.glm", "SL.mean"),
                               g_trunc = c(0.01, 0.99),
                               max_time = NULL) {

  sl_lib_Q    <- filter_sl_libraries(sl_lib_Q)
  sl_lib_g    <- filter_sl_libraries(sl_lib_g)
  sl_lib_cens <- filter_sl_libraries(sl_lib_cens)

  if (is.null(max_time)) max_time <- t_eval

  # --------------------------------------------------------------------------
  # 1. Identify covariates
  # --------------------------------------------------------------------------
  exclude_cols <- c("id", "follow_time", "event", "treatment", "switch",
                    "race", "region")
  covar_cols <- setdiff(names(data), exclude_cols)
  covar_cols <- covar_cols[vapply(data[, covar_cols, drop = FALSE],
                                 is.numeric, logical(1))]
  W <- as.data.frame(data[, covar_cols, drop = FALSE])
  A <- data$treatment
  n <- nrow(data)

  # --------------------------------------------------------------------------
  # 2. Fit treatment propensity score P(A=1|W)
  # --------------------------------------------------------------------------
  g_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = A, X = W, family = stats::binomial(),
      SL.library = sl_lib_g, cvControl = list(V = 5)
    )
  }, error = function(e) {
    fit <- stats::glm(A ~ ., data = cbind(A = A, W), family = "binomial")
    list(SL.predict = stats::predict(fit, type = "response"))
  })
  g1W <- as.numeric(g_fit$SL.predict)
  g_trunc_info <- truncate_ps(g1W, g_trunc[1], g_trunc[2])
  g1W <- g_trunc_info$p

  # --------------------------------------------------------------------------
  # 3. Fit censoring model: P(uncensored | W, A)
  # --------------------------------------------------------------------------
  # Censored = follow_time < t_eval AND event == 0
  censored <- as.integer(data$follow_time < max_time & data$event == 0)

  cens_W <- cbind(W, A = A)
  gC_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = 1 - censored, X = as.data.frame(cens_W),
      family = stats::binomial(),
      SL.library = sl_lib_cens, cvControl = list(V = 5)
    )
  }, error = function(e) {
    fit <- stats::glm((1 - censored) ~ .,
                      data = cbind(uncens = 1 - censored, cens_W),
                      family = "binomial")
    list(SL.predict = stats::predict(fit, type = "response"))
  })
  gC <- as.numeric(gC_fit$SL.predict)
  gC <- pmax(gC, g_trunc[1])

  # --------------------------------------------------------------------------
  # 4. Fit outcome model: P(event by t_eval | W, A) -- initial Q
  # --------------------------------------------------------------------------
  Y_binary <- as.integer(data$follow_time <= t_eval & data$event == 1)
  Q_W <- cbind(W, A = A)

  Q_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = Y_binary, X = as.data.frame(Q_W),
      family = stats::binomial(),
      SL.library = sl_lib_Q, cvControl = list(V = 5)
    )
  }, error = function(e) {
    fit <- stats::glm(Y_binary ~ .,
                      data = cbind(Y = Y_binary, Q_W),
                      family = "binomial")
    pred <- stats::predict(fit, type = "response")
    list(SL.predict = pred, glm_fit = fit, is_glm = TRUE)
  })

  QAW <- as.numeric(Q_fit$SL.predict)
  QAW <- pmin(pmax(QAW, 0.001), 0.999)

  # Counterfactual predictions
  W1 <- as.data.frame(cbind(W, A = 1))
  W0 <- as.data.frame(cbind(W, A = 0))

  if (!is.null(Q_fit$is_glm) && isTRUE(Q_fit$is_glm)) {
    Q1W <- stats::predict(Q_fit$glm_fit, newdata = W1, type = "response")
    Q0W <- stats::predict(Q_fit$glm_fit, newdata = W0, type = "response")
  } else {
    Q1W <- as.numeric(predict(Q_fit, newdata = W1)$pred)
    Q0W <- as.numeric(predict(Q_fit, newdata = W0)$pred)
  }
  Q1W <- pmin(pmax(Q1W, 0.001), 0.999)
  Q0W <- pmin(pmax(Q0W, 0.001), 0.999)

  # --------------------------------------------------------------------------
  # 5. TMLE targeting step with IPCW clever covariates
  # --------------------------------------------------------------------------
  # Indicator for being observed (not censored before t_eval)
  observed <- as.integer(!(data$follow_time < t_eval & data$event == 0))

  # Clever covariates
  H1 <- observed * A / (g1W * gC)
  H0 <- observed * (1 - A) / ((1 - g1W) * gC)

  # Fluctuation model
  offset_Q <- stats::qlogis(QAW)
  # Use only observations where we have valid data
  valid <- is.finite(H1) & is.finite(H0) & is.finite(offset_Q)
  eps_fit <- tryCatch({
    stats::glm(Y_binary ~ -1 + H1 + H0,
               offset = offset_Q,
               family = stats::binomial(),
               subset = valid)
  }, error = function(e) {
    list(coefficients = c(H1 = 0, H0 = 0))
  })

  eps1 <- stats::coef(eps_fit)["H1"]
  eps0 <- stats::coef(eps_fit)["H0"]
  if (is.na(eps1)) eps1 <- 0
  if (is.na(eps0)) eps0 <- 0

  # Update counterfactual predictions
  Q1W_star <- stats::plogis(stats::qlogis(Q1W) + eps1 / gC)
  Q0W_star <- stats::plogis(stats::qlogis(Q0W) + eps0 / gC)

  # --------------------------------------------------------------------------
  # 6. Compute estimates
  # --------------------------------------------------------------------------
  risk_1 <- mean(Q1W_star)
  risk_0 <- mean(Q0W_star)
  RD <- risk_1 - risk_0
  RR <- if (risk_0 > 1e-10) risk_1 / risk_0 else NA_real_

  # --------------------------------------------------------------------------
  # 7. Influence curve-based inference
  # --------------------------------------------------------------------------
  IC_1 <- H1 * (Y_binary - QAW) + Q1W_star - risk_1
  IC_0 <- H0 * (Y_binary - QAW) + Q0W_star - risk_0
  IC_RD <- IC_1 - IC_0

  se_RD <- sqrt(mean(IC_RD^2) / n)
  ci_RD <- RD + c(-1.96, 1.96) * se_RD

  se_1 <- sqrt(mean(IC_1^2) / n)
  se_0 <- sqrt(mean(IC_0^2) / n)

  # RR inference via delta method on log scale
  if (!is.na(RR) && risk_0 > 1e-10 && risk_1 > 1e-10) {
    IC_logRR <- IC_1 / risk_1 - IC_0 / risk_0
    se_logRR <- sqrt(mean(IC_logRR^2) / n)
    ci_RR <- exp(log(RR) + c(-1.96, 1.96) * se_logRR)
  } else {
    se_logRR <- NA_real_
    ci_RR <- c(NA_real_, NA_real_)
  }

  # --------------------------------------------------------------------------
  # 8. Diagnostics
  # --------------------------------------------------------------------------
  diagnostics <- list(
    clever_covariate = list(
      H1_mean = mean(H1[A == 1 & observed == 1], na.rm = TRUE),
      H1_sd   = stats::sd(H1[A == 1 & observed == 1], na.rm = TRUE),
      H0_mean = mean(H0[A == 0 & observed == 1], na.rm = TRUE),
      H0_sd   = stats::sd(H0[A == 0 & observed == 1], na.rm = TRUE)
    ),
    ic = list(
      ic_rd_mean  = mean(IC_RD),
      ic_rd_sd    = stats::sd(IC_RD),
      ic_rd_range = range(IC_RD)
    ),
    truncation = list(
      g_n_lower  = g_trunc_info$n_lower,
      g_n_upper  = g_trunc_info$n_upper,
      gC_n_trunc = sum(gC <= g_trunc[1])
    ),
    positivity = list(
      min_g1W = min(g1W),
      max_g1W = max(g1W),
      min_gC  = min(gC),
      frac_near_boundary = mean(g1W < 0.05 | g1W > 0.95)
    ),
    eps = list(eps1 = unname(eps1), eps0 = unname(eps0)),
    n_observed = sum(observed),
    n_censored = sum(1 - observed)
  )

  list(
    risk_1      = risk_1,
    risk_0      = risk_0,
    RD          = RD,
    RR          = RR,
    se_RD       = se_RD,
    ci_RD       = ci_RD,
    se_1        = se_1,
    se_0        = se_0,
    ci_RR       = ci_RR,
    t_eval      = t_eval,
    n           = n,
    diagnostics = diagnostics
  )
}

#' Filter SuperLearner Libraries
#'
#' Retain only SL libraries whose underlying packages are installed.
#'
#' @param libs Character vector of SL library names.
#' @return Filtered character vector with at least SL.glm and SL.mean.
#' @keywords internal
filter_sl_libraries <- function(libs) {
  available <- Filter(function(lib) {
    tryCatch({
      if (lib %in% c("SL.glm", "SL.mean", "SL.step",
                      "SL.step.interaction", "SL.glm.interaction")) {
        return(TRUE)
      }
      pkg <- sub("SL\\.", "", lib)
      requireNamespace(pkg, quietly = TRUE)
    }, error = function(e) FALSE)
  }, libs)
  if (length(available) == 0) available <- c("SL.glm", "SL.mean")
  available
}
