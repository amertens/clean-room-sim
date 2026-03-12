#' @title IPTW Survival Estimator
#' @description Inverse probability of treatment weighted estimator for
#'   survival risk using weighted Kaplan-Meier or weighted Cox.
#' @name iptw_survival
NULL

#' IPTW Survival Risk Estimator
#'
#' Estimates counterfactual risk at time t using inverse probability of
#' treatment weighting with a weighted Kaplan-Meier estimator.
#'
#' @param data Data frame with follow_time, event, treatment, and covariates.
#' @param t_eval Numeric time point for risk evaluation.
#' @param sl_lib Character vector of SL libraries for PS estimation.
#' @param g_trunc Numeric vector of length 2 for truncation bounds.
#' @param use_stabilized Logical; if TRUE, use stabilized weights.
#' @return A list with risk_1, risk_0, RD, RR, SE, CI, and weight diagnostics.
#' @export
iptw_survival <- function(data,
                          t_eval = 180,
                          sl_lib = c("SL.glm", "SL.mean"),
                          g_trunc = c(0.01, 0.99),
                          use_stabilized = TRUE) {

  sl_lib <- filter_sl_libraries(sl_lib)

  exclude_cols <- c("id", "follow_time", "event", "treatment", "switch",
                    "race", "region")
  covar_cols <- setdiff(names(data), exclude_cols)
  covar_cols <- covar_cols[vapply(data[, covar_cols, drop = FALSE],
                                 is.numeric, logical(1))]

  W <- as.data.frame(data[, covar_cols, drop = FALSE])
  A <- data$treatment
  n <- nrow(data)

  # Fit PS model
  g_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = A, X = W, family = stats::binomial(),
      SL.library = sl_lib, cvControl = list(V = 5)
    )
  }, error = function(e) {
    fit <- stats::glm(A ~ ., data = cbind(A = A, W), family = "binomial")
    list(SL.predict = stats::predict(fit, type = "response"))
  })

  ps <- as.numeric(g_fit$SL.predict)
  trunc_info <- truncate_ps(ps, g_trunc[1], g_trunc[2])
  ps <- trunc_info$p

  # Compute IPW weights
  if (use_stabilized) {
    p_A <- mean(A)
    w <- ifelse(A == 1, p_A / ps, (1 - p_A) / (1 - ps))
  } else {
    w <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
  }

  # Weighted KM risk at t_eval for treated and control groups
  risk_1 <- weighted_km_risk(
    time   = data$follow_time[A == 1],
    event  = data$event[A == 1],
    weights = w[A == 1],
    t_eval = t_eval
  )

  risk_0 <- weighted_km_risk(
    time   = data$follow_time[A == 0],
    event  = data$event[A == 0],
    weights = w[A == 0],
    t_eval = t_eval
  )

  RD <- risk_1 - risk_0
  RR <- if (risk_0 > 1e-10) risk_1 / risk_0 else NA_real_

  # Bootstrap SE (quick: use influence-function approximation)
  # Hajek-style SE approximation for weighted estimator
  Y1 <- as.integer(data$follow_time[A == 1] <= t_eval &
                      data$event[A == 1] == 1)
  Y0 <- as.integer(data$follow_time[A == 0] <= t_eval &
                      data$event[A == 0] == 1)

  IF_1 <- w[A == 1] * (Y1 - risk_1)
  IF_0 <- w[A == 0] * (Y0 - risk_0)

  # Pad to full length
  IC_full <- numeric(n)
  IC_full[A == 1] <- IF_1
  IC_full[A == 0] <- -IF_0

  se_RD <- sqrt(sum(IC_full^2) / n^2)
  ci_RD <- RD + c(-1.96, 1.96) * se_RD

  # Weight diagnostics
  weight_diag <- list(
    mean_w_treated = mean(w[A == 1]),
    mean_w_control = mean(w[A == 0]),
    max_w          = max(w),
    ess_treated    = effective_ss(w[A == 1]),
    ess_control    = effective_ss(w[A == 0]),
    n_truncated_lower = trunc_info$n_lower,
    n_truncated_upper = trunc_info$n_upper
  )

  list(
    method      = "IPTW",
    risk_1      = risk_1,
    risk_0      = risk_0,
    RD          = RD,
    RR          = RR,
    se_RD       = se_RD,
    ci_RD       = ci_RD,
    t_eval      = t_eval,
    n           = n,
    weight_diag = weight_diag
  )
}
