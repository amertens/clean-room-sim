#' @title Cox Proportional Hazards Estimator
#' @description Cox PH model for time-to-event data with PH diagnostics.
#' @name cox_ph
NULL

#' Cox PH Estimator
#'
#' Fits a Cox proportional hazards model and returns the hazard ratio for
#' treatment with robust standard errors. Includes Schoenfeld residual test
#' for PH assumption.
#'
#' @param data Data frame with follow_time, event, treatment, and covariates.
#' @param covar_cols Character vector of covariate column names. If NULL,
#'   all numeric columns excluding id, follow_time, event, treatment, switch
#'   are used.
#' @return A list with HR, CI, p-value, PH test results, and risk estimates.
#' @export
cox_ph_estimator <- function(data, covar_cols = NULL) {
  if (is.null(covar_cols)) {
    exclude <- c("id", "follow_time", "event", "treatment", "switch",
                 "race", "region")
    covar_cols <- setdiff(names(data), exclude)
    covar_cols <- covar_cols[vapply(data[, covar_cols, drop = FALSE],
                                   is.numeric, logical(1))]
  }

  # Build formula
  covars <- paste(covar_cols, collapse = " + ")
  fmla <- stats::as.formula(
    paste0("survival::Surv(follow_time, event) ~ treatment + ", covars)
  )

  # Fit Cox model with robust SE
  cox_fit <- survival::coxph(fmla, data = data, robust = TRUE)
  cox_summary <- summary(cox_fit)

  # Extract treatment HR
  trt_row <- which(rownames(cox_summary$coefficients) == "treatment")
  if (length(trt_row) == 0) {
    stop("Treatment variable not found in Cox model.")
  }

  hr       <- exp(cox_summary$coefficients[trt_row, "coef"])
  se_coef  <- cox_summary$coefficients[trt_row, "robust se"]
  ci_hr    <- exp(cox_summary$coefficients[trt_row, "coef"] +
                    c(-1.96, 1.96) * se_coef)
  pvalue   <- cox_summary$coefficients[trt_row, "Pr(>|z|)"]

  # PH diagnostic: Schoenfeld test
  ph_test <- tryCatch({
    zt <- survival::cox.zph(cox_fit)
    trt_ph <- zt$table["treatment", ]
    list(
      chisq  = trt_ph["chisq"],
      df     = trt_ph["df"],
      pvalue = trt_ph["p"],
      ph_violated = trt_ph["p"] < 0.05
    )
  }, error = function(e) {
    list(chisq = NA, df = NA, pvalue = NA, ph_violated = NA)
  })

  list(
    method    = "Cox PH",
    HR        = hr,
    log_HR    = log(hr),
    se_logHR  = se_coef,
    ci_HR     = ci_hr,
    pvalue    = pvalue,
    ph_test   = ph_test,
    n         = nrow(data),
    n_events  = sum(data$event),
    cox_fit   = cox_fit
  )
}

#' Cox PH Risk Estimate at Time t
#'
#' Uses the fitted Cox model to predict marginal risk at time t under
#' treatment and control via standardization (g-computation with Cox).
#'
#' @param cox_result Output from \code{cox_ph_estimator}.
#' @param data Original data frame.
#' @param t_eval Time at which to evaluate risk.
#' @return A list with risk_1, risk_0, RD, and RR.
#' @export
cox_risk_at_t <- function(cox_result, data, t_eval = 180) {
  cox_fit <- cox_result$cox_fit
  # Predict under all-treated and all-control
  d1 <- d0 <- data
  d1$treatment <- 1L
  d0$treatment <- 0L

  # Get survival predictions
  surv1 <- tryCatch({
    sf1 <- survival::survfit(cox_fit, newdata = d1)
    # Find the survival at t_eval
    idx <- max(which(sf1$time <= t_eval), 1)
    # Average marginal survival across individuals (approximation)
    mean(1 - summary(sf1, times = t_eval)$surv)
  }, error = function(e) {
    # Fallback: use baseline cumulative hazard
    bh <- survival::basehaz(cox_fit, centered = FALSE)
    idx <- max(which(bh$time <= t_eval), 1)
    H0t <- bh$hazard[idx]
    lp1 <- predict(cox_fit, newdata = d1, type = "lp")
    mean(1 - exp(-H0t * exp(lp1)))
  })

  surv0 <- tryCatch({
    sf0 <- survival::survfit(cox_fit, newdata = d0)
    mean(1 - summary(sf0, times = t_eval)$surv)
  }, error = function(e) {
    bh <- survival::basehaz(cox_fit, centered = FALSE)
    idx <- max(which(bh$time <= t_eval), 1)
    H0t <- bh$hazard[idx]
    lp0 <- predict(cox_fit, newdata = d0, type = "lp")
    mean(1 - exp(-H0t * exp(lp0)))
  })

  risk_1 <- surv1
  risk_0 <- surv0
  RD <- risk_1 - risk_0
  RR <- if (risk_0 > 1e-10) risk_1 / risk_0 else NA_real_

  list(risk_1 = risk_1, risk_0 = risk_0, RD = RD, RR = RR, t_eval = t_eval)
}
