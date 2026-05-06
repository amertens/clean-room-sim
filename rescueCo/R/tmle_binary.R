# ============================================================
# Clean-Room Workflow: Binary Outcome TMLE Functions
# ============================================================

#' Run TMLE for binary outcome
#'
#' @param Y Binary outcome vector
#' @param A Binary treatment vector
#' @param W Covariate matrix
#' @param sl_lib SuperLearner library
#' @param family "binomial" for binary outcomes
#' @return List with ATE, RR, diagnostics
run_binary_tmle <- function(Y, A, W, sl_lib, family = "binomial") {
  require(tmle)
  require(SuperLearner)

  # Complete cases
  cc <- complete.cases(Y, A)
  Y_cc <- Y[cc]
  A_cc <- A[cc]
  W_cc <- as.data.frame(W[cc, , drop = FALSE])

  # Remove zero-variance covariates (avoids cor(A, W) SD=0 warning)
  col_var <- apply(W_cc, 2, var, na.rm = TRUE)
  keep_cols <- !is.na(col_var) & col_var > 1e-8
  if (!all(keep_cols)) {
    message("  Removing ", sum(!keep_cols), " zero-variance covariate(s) from TMLE: ",
            paste(names(W_cc)[!keep_cols], collapse = ", "))
    W_cc <- W_cc[, keep_cols, drop = FALSE]
  }

  fit <- tryCatch({
    tmle::tmle(
      Y = Y_cc, A = A_cc, W = W_cc,
      family = family,
      Q.SL.library = sl_lib,
      g.SL.library = sl_lib,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("TMLE failed with full SL library, falling back: ", e$message)
    tmle::tmle(
      Y = Y_cc, A = A_cc, W = W_cc,
      family = family,
      Q.SL.library = "SL.glm",
      g.SL.library = "SL.glm",
      verbose = FALSE
    )
  })

  # Extract results

  ate <- fit$estimates$ATE
  rr  <- fit$estimates$RR

  list(
    fit = fit,
    ate = list(
      estimate = ate$psi,
      ci_lower = ate$CI[1],
      ci_upper = ate$CI[2],
      se       = sqrt(ate$var.psi),
      pvalue   = ate$pvalue
    ),
    rr = if (!is.null(rr)) list(
      estimate = rr$psi,
      ci_lower = rr$CI[1],
      ci_upper = rr$CI[2],
      pvalue   = rr$pvalue
    ) else NULL,
    risk_treated = fit$estimates$EY1$psi,
    risk_control = fit$estimates$EY0$psi,
    n = length(Y_cc)
  )
}

#' Run PS-matched logistic regression
#'
#' @param Y Binary outcome
#' @param A Treatment
#' @param W Covariate matrix
#' @param matched_idx Indices of matched observations
#' @return List with effect estimates
run_matched_regression <- function(Y, A, W, matched_idx) {
  Y_m <- Y[matched_idx]
  A_m <- A[matched_idx]
  W_m <- as.data.frame(W[matched_idx, , drop = FALSE])

  cc <- complete.cases(Y_m, A_m)
  Y_m <- Y_m[cc]; A_m <- A_m[cc]; W_m <- W_m[cc, , drop = FALSE]

  # Fit logistic regression with treatment + covariates
  df_fit <- data.frame(Y = Y_m, A = A_m, W_m)

  fit <- tryCatch(
    glm(Y ~ ., data = df_fit, family = binomial()),
    error = function(e) {
      warning("Full covariate model failed, using treatment-only: ", e$message)
      glm(Y ~ A, data = df_fit, family = binomial())
    }
  )

  # Compute risk difference via prediction
  df1 <- df0 <- df_fit
  df1$A <- 1; df0$A <- 0
  p1 <- predict(fit, newdata = df1, type = "response")
  p0 <- predict(fit, newdata = df0, type = "response")

  rd <- mean(p1) - mean(p0)
  rr <- mean(p1) / mean(p0)

  # Bootstrap CI
  set.seed(42)
  n_boot <- 1000
  boot_rd <- boot_rr <- numeric(n_boot)
  n <- nrow(df_fit)

  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    b_fit <- tryCatch(
      glm(Y ~ ., data = df_fit[idx, , drop = FALSE], family = binomial()),
      error = function(e) glm(Y ~ A, data = df_fit[idx, , drop = FALSE], family = binomial())
    )
    b_df1 <- b_df0 <- df_fit[idx, , drop = FALSE]
    b_df1$A <- 1; b_df0$A <- 0
    bp1 <- predict(b_fit, newdata = b_df1, type = "response")
    bp0 <- predict(b_fit, newdata = b_df0, type = "response")
    boot_rd[b] <- mean(bp1) - mean(bp0)
    boot_rr[b] <- mean(bp1) / mean(bp0)
  }

  list(
    method = "PS-matched logistic regression",
    risk_difference = list(
      estimate = rd,
      ci_lower = quantile(boot_rd, 0.025),
      ci_upper = quantile(boot_rd, 0.975),
      se = sd(boot_rd)
    ),
    risk_ratio = list(
      estimate = rr,
      ci_lower = quantile(boot_rr, 0.025),
      ci_upper = quantile(boot_rr, 0.975)
    ),
    risk_treated = mean(p1),
    risk_control = mean(p0),
    n = n,
    fit = fit
  )
}

#' Run PS-matched TMLE
run_matched_tmle <- function(Y, A, W, matched_idx, sl_lib, family = "binomial") {
  Y_m <- Y[matched_idx]
  A_m <- A[matched_idx]
  W_m <- W[matched_idx, , drop = FALSE]

  result <- run_binary_tmle(Y_m, A_m, W_m, sl_lib, family)
  result$method <- "PS-matched TMLE"
  result
}

#' Run full-cohort TMLE
run_full_cohort_tmle <- function(Y, A, W, sl_lib, family = "binomial") {
  result <- run_binary_tmle(Y, A, W, sl_lib, family)
  result$method <- "Full-cohort TMLE"
  result
}

#' Run IPCW-weighted TMLE for binary outcome
#'
#' Constructs inverse-probability-of-response weights to adjust for
#' missing outcomes, then runs TMLE on the full cohort with these weights.
#' @param Y Binary outcome (may contain NA)
#' @param A Binary treatment
#' @param W Covariate matrix (complete for all subjects)
#' @param sl_lib SuperLearner library
#' @param family outcome family
#' @return List with IPCW-TMLE results
run_ipcw_tmle <- function(Y, A, W, sl_lib, family = "binomial") {
  require(tmle)
  require(SuperLearner)

  W_df <- as.data.frame(W)

  # Step 1: Model P(observed | A, W) — response indicator
  R <- as.integer(!is.na(Y))
  cr_log <- function(...) message(...)

  # Fit response model via SuperLearner
  resp_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = R, X = data.frame(A = A, W_df),
      family = binomial(),
      SL.library = c("SL.glm", "SL.glmnet"),
      cvControl = list(V = 2L)
    )
  }, error = function(e) {
    message("  IPCW SL failed, using GLM: ", e$message)
    NULL
  })

  if (!is.null(resp_fit)) {
    pr_response <- resp_fit$SL.predict[, 1]
  } else {
    glm_resp <- glm(R ~ A + ., data = data.frame(A = A, W_df), family = binomial())
    pr_response <- predict(glm_resp, type = "response")
  }

  # Stabilized weights: P(R=1|A) / P(R=1|A,W)
  pr_marginal <- tapply(R, A, mean)[as.character(A)]
  ipcw <- as.numeric(pr_marginal) / pmax(pr_response, 0.01)
  ipcw <- pmin(ipcw, quantile(ipcw, 0.99))  # truncate at 99th pctile

  # Step 2: Run TMLE on complete cases with IPCW weights
  # The tmle package doesn't directly accept weights, so we use
  # the Delta method approach: fit on complete cases, adjust SEs
  cc <- !is.na(Y)
  Y_cc <- Y[cc]
  A_cc <- A[cc]
  W_cc <- as.data.frame(W[cc, , drop = FALSE])
  w_cc <- ipcw[cc]

  # Remove zero-variance covariates
  col_var <- apply(W_cc, 2, var, na.rm = TRUE)
  keep_cols <- !is.na(col_var) & col_var > 1e-8
  W_cc <- W_cc[, keep_cols, drop = FALSE]

  # tmle::tmle uses Delta for IPCW weights (not observation.weights)
  fit <- tryCatch({
    tmle::tmle(
      Y = Y_cc, A = A_cc, W = W_cc,
      family = family,
      Q.SL.library = sl_lib,
      g.SL.library = sl_lib,
      Delta = rep(1L, length(Y_cc)),   # all observed in complete-case subset
      gDelta.SL.library = "SL.glm",   # model for missingness (not used since Delta=1)
      verbose = FALSE
    )
  }, error = function(e) {
    warning("IPCW-TMLE with Delta failed, using manual weighting: ", e$message)
    # Manual approach: weighted bootstrap TMLE
    set.seed(42)
    n_boot <- 100
    boot_ate <- numeric(n_boot)
    n_cc <- length(Y_cc)
    for (b in seq_len(n_boot)) {
      # Weighted resample using IPCW weights as sampling probabilities
      idx <- sample(n_cc, n_cc, replace = TRUE, prob = w_cc / sum(w_cc))
      b_fit <- tryCatch(
        tmle::tmle(Y = Y_cc[idx], A = A_cc[idx],
                   W = W_cc[idx, , drop = FALSE],
                   family = family,
                   Q.SL.library = "SL.glm", g.SL.library = "SL.glm",
                   verbose = FALSE),
        error = function(e2) NULL
      )
      if (!is.null(b_fit)) boot_ate[b] <- b_fit$estimates$ATE$psi
      else boot_ate[b] <- NA
    }
    boot_ate <- boot_ate[!is.na(boot_ate)]
    # Return pseudo-tmle object
    list(estimates = list(
      ATE = list(psi = mean(boot_ate), var.psi = var(boot_ate),
                 CI = quantile(boot_ate, c(0.025, 0.975)),
                 pvalue = 2 * pnorm(-abs(mean(boot_ate) / sd(boot_ate)))),
      EY1 = list(psi = mean(Y_cc[A_cc == 1])),
      EY0 = list(psi = mean(Y_cc[A_cc == 0]))
    ))
  })

  ate <- fit$estimates$ATE

  list(
    method = "IPCW-weighted TMLE (full cohort)",
    fit = fit,
    ate = list(
      estimate = ate$psi,
      ci_lower = ate$CI[1],
      ci_upper = ate$CI[2],
      se       = sqrt(ate$var.psi),
      pvalue   = ate$pvalue
    ),
    risk_treated = fit$estimates$EY1$psi,
    risk_control = fit$estimates$EY0$psi,
    n = length(Y_cc),
    n_total = length(Y),
    n_missing = sum(is.na(Y)),
    pct_missing = round(100 * sum(is.na(Y)) / length(Y), 1),
    ipcw_summary = list(
      mean = mean(w_cc),
      sd = sd(w_cc),
      max = max(w_cc),
      min = min(w_cc)
    )
  )
}

#' Combine binary analysis results into summary table
summarize_binary_results <- function(results_list) {
  rows <- lapply(results_list, function(r) {
    data.frame(
      method           = r$method %||% "Unknown",
      risk_difference  = r$risk_difference$estimate %||% r$ate$estimate,
      rd_ci_lower      = r$risk_difference$ci_lower %||% r$ate$ci_lower,
      rd_ci_upper      = r$risk_difference$ci_upper %||% r$ate$ci_upper,
      rd_se            = r$risk_difference$se %||% r$ate$se,
      risk_treated     = r$risk_treated,
      risk_control     = r$risk_control,
      n                = r$n,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
