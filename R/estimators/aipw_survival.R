#' @title Augmented Inverse Probability of Treatment Weighting (AIPW)
#' @description Doubly-robust AIPW estimator for survival risk combining
#'   an outcome regression with inverse probability weighting, using
#'   SuperLearner for both nuisance models.
#' @name aipw
NULL

#' AIPW Survival Risk Estimator
#'
#' Estimates counterfactual risk at time t using augmented inverse probability
#' weighting (AIPW / doubly-robust). Combines an outcome regression estimate
#' with a bias-correction term based on inverse propensity and censoring
#' weights.
#'
#' Unlike TMLE, AIPW applies the augmentation additively rather than through
#' a targeted fluctuation step. It is doubly robust (consistent if either the
#' outcome or treatment model is correct) but not substitution-based.
#'
#' @param data Data frame with follow_time, event, treatment, and covariates.
#' @param t_eval Numeric time point for risk evaluation.
#' @param sl_lib_Q Character vector of SL libraries for outcome model.
#' @param sl_lib_g Character vector of SL libraries for propensity model.
#' @param sl_lib_cens Character vector of SL libraries for censoring model.
#' @param g_trunc Numeric vector of length 2 for truncation bounds.
#' @return A list with risk_1, risk_0, RD, RR, SE, CI, and diagnostics.
#' @export
aipw_survival <- function(data,
                          t_eval = 180,
                          sl_lib_Q = c("SL.glm", "SL.mean"),
                          sl_lib_g = c("SL.glm", "SL.mean"),
                          sl_lib_cens = c("SL.glm", "SL.mean"),
                          g_trunc = c(0.01, 0.99)) {

  sl_lib_Q    <- filter_sl_libraries(sl_lib_Q)
  sl_lib_g    <- filter_sl_libraries(sl_lib_g)
  sl_lib_cens <- filter_sl_libraries(sl_lib_cens)

  # Identify covariates
  exclude_cols <- c("id", "follow_time", "event", "treatment", "switch",
                    "race", "region")
  covar_cols <- setdiff(names(data), exclude_cols)
  covar_cols <- covar_cols[vapply(data[, covar_cols, drop = FALSE],
                                 is.numeric, logical(1))]

  W <- as.data.frame(data[, covar_cols, drop = FALSE])
  A <- data$treatment
  n <- nrow(data)

  # Binary outcome: event by t_eval
  Y <- as.integer(data$follow_time <= t_eval & data$event == 1)

  # --------------------------------------------------------------------------
  # 1. Fit propensity score P(A=1|W)
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
  trunc_info <- truncate_ps(g1W, g_trunc[1], g_trunc[2])
  g1W <- trunc_info$p

  # --------------------------------------------------------------------------
  # 2. Fit censoring model P(uncensored by t_eval | W, A)
  # --------------------------------------------------------------------------
  censored <- as.integer(data$follow_time < t_eval & data$event == 0)
  cens_W <- cbind(W, A = A)
  gC_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = 1 - censored, X = as.data.frame(cens_W),
      family = stats::binomial(),
      SL.library = sl_lib_cens, cvControl = list(V = 5)
    )
  }, error = function(e) {
    uncens <- 1 - censored
    fit <- stats::glm(uncens ~ ., data = as.data.frame(cens_W),
                      family = "binomial")
    list(SL.predict = stats::predict(fit, type = "response"))
  })
  gC <- pmax(as.numeric(gC_fit$SL.predict), g_trunc[1])

  # --------------------------------------------------------------------------
  # 3. Fit outcome model P(Y=1 | W, A) -- Q-bar
  # --------------------------------------------------------------------------
  Q_W <- cbind(W, A = A)
  Q_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = Y, X = as.data.frame(Q_W),
      family = stats::binomial(),
      SL.library = sl_lib_Q, cvControl = list(V = 5)
    )
  }, error = function(e) {
    fit <- stats::glm(Y ~ ., data = cbind(Y = Y, Q_W), family = "binomial")
    list(SL.predict = stats::predict(fit, type = "response"),
         glm_fit = fit, is_glm = TRUE)
  })

  QAW <- as.numeric(Q_fit$SL.predict)
  QAW <- pmin(pmax(QAW, 0.001), 0.999)

  # Counterfactual predictions Q(1, W) and Q(0, W)
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
  # 4. AIPW estimator: augmented IPW (no targeting step)
  # --------------------------------------------------------------------------
  observed <- as.integer(!(data$follow_time < t_eval & data$event == 0))

  # AIPW components for each potential outcome
  # \hat\psi_1 = (1/n) \sum [ Q1W + (A * observed / (g1W * gC)) * (Y - QAW) ]
  psi_1_i <- Q1W + (A * observed / (g1W * gC)) * (Y - QAW)
  psi_0_i <- Q0W + ((1 - A) * observed / ((1 - g1W) * gC)) * (Y - QAW)

  risk_1 <- mean(psi_1_i)
  risk_0 <- mean(psi_0_i)
  RD <- risk_1 - risk_0
  RR <- if (risk_0 > 1e-10) risk_1 / risk_0 else NA_real_

  # --------------------------------------------------------------------------
  # 5. Influence-curve based inference
  # --------------------------------------------------------------------------
  IC_1 <- psi_1_i - risk_1
  IC_0 <- psi_0_i - risk_0
  IC_RD <- IC_1 - IC_0

  se_RD <- sqrt(mean(IC_RD^2) / n)
  ci_RD <- RD + c(-1.96, 1.96) * se_RD

  se_1 <- sqrt(mean(IC_1^2) / n)
  se_0 <- sqrt(mean(IC_0^2) / n)

  if (!is.na(RR) && risk_0 > 1e-10 && risk_1 > 1e-10) {
    IC_logRR <- IC_1 / risk_1 - IC_0 / risk_0
    se_logRR <- sqrt(mean(IC_logRR^2) / n)
    ci_RR <- exp(log(RR) + c(-1.96, 1.96) * se_logRR)
  } else {
    ci_RR <- c(NA_real_, NA_real_)
  }

  # --------------------------------------------------------------------------
  # 6. Diagnostics
  # --------------------------------------------------------------------------
  diagnostics <- list(
    ic = list(
      ic_rd_mean  = mean(IC_RD),
      ic_rd_sd    = stats::sd(IC_RD)
    ),
    truncation = list(
      g_n_lower  = trunc_info$n_lower,
      g_n_upper  = trunc_info$n_upper,
      gC_n_trunc = sum(gC <= g_trunc[1])
    ),
    positivity = list(
      min_g1W = min(g1W),
      max_g1W = max(g1W),
      frac_near_boundary = mean(g1W < 0.05 | g1W > 0.95)
    ),
    n_observed = sum(observed),
    n_censored = sum(1 - observed)
  )

  list(
    method      = "AIPW",
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
