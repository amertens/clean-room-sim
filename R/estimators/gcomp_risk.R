#' @title G-Computation Risk Estimator
#' @description Outcome regression plug-in estimator for counterfactual risk
#'   at time t using SuperLearner.
#' @name gcomp
NULL

#' G-Computation Risk Estimator
#'
#' Estimates counterfactual risk at time t using outcome regression
#' (g-computation / plug-in estimator). Fits a model for P(event by t | W, A)
#' and standardizes over the covariate distribution.
#'
#' @param data Data frame with follow_time, event, treatment, and covariates.
#' @param t_eval Numeric time point for risk evaluation.
#' @param sl_lib Character vector of SL libraries for outcome model.
#' @return A list with risk_1, risk_0, RD, RR, and bootstrap SE/CI.
#' @export
gcomp_risk <- function(data,
                       t_eval = 180,
                       sl_lib = c("SL.glm", "SL.mean")) {

  sl_lib <- filter_sl_libraries(sl_lib)

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

  Q_W <- cbind(W, A = A)

  # Fit outcome model
  Q_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = Y, X = as.data.frame(Q_W),
      family = stats::binomial(),
      SL.library = sl_lib, cvControl = list(V = 5)
    )
  }, error = function(e) {
    fit <- stats::glm(Y ~ ., data = cbind(Y = Y, Q_W), family = "binomial")
    list(SL.predict = stats::predict(fit, type = "response"),
         glm_fit = fit, is_glm = TRUE)
  })

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

  risk_1 <- mean(Q1W)
  risk_0 <- mean(Q0W)
  RD <- risk_1 - risk_0
  RR <- if (risk_0 > 1e-10) risk_1 / risk_0 else NA_real_

  # Simple nonparametric bootstrap for SE
  n_boot <- 200
  boot_rd <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    d_b <- data[idx, , drop = FALSE]
    Y_b <- as.integer(d_b$follow_time <= t_eval & d_b$event == 1)
    W_b <- as.data.frame(d_b[, covar_cols, drop = FALSE])
    A_b <- d_b$treatment

    fit_b <- tryCatch({
      stats::glm(Y_b ~ ., data = cbind(Y = Y_b, W_b, A = A_b),
                 family = "binomial")
    }, error = function(e) NULL)

    if (is.null(fit_b)) { boot_rd[b] <- NA; next }

    W1b <- cbind(W_b, A = 1)
    W0b <- cbind(W_b, A = 0)
    q1 <- stats::predict(fit_b, newdata = W1b, type = "response")
    q0 <- stats::predict(fit_b, newdata = W0b, type = "response")
    boot_rd[b] <- mean(q1) - mean(q0)
  }

  se_RD <- stats::sd(boot_rd, na.rm = TRUE)
  ci_RD <- RD + c(-1.96, 1.96) * se_RD

  list(
    method   = "G-computation",
    risk_1   = risk_1,
    risk_0   = risk_0,
    RD       = RD,
    RR       = RR,
    se_RD    = se_RD,
    ci_RD    = ci_RD,
    t_eval   = t_eval,
    n        = n,
    n_boot   = n_boot
  )
}
