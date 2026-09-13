# ============================================================
# Clean-Room Workflow: Survival TMLE Functions
# Discrete-time survival TMLE implementation
# ============================================================

#' Expand data to person-time (long) format for discrete-time survival
#'
#' @param time_var Vector of event/censor times (integer days or time periods)
#' @param event_var Binary event indicator (1 = event, 0 = censored)
#' @param A Treatment vector
#' @param W Covariate matrix
#' @param max_time Maximum time horizon
#' @param time_grid Optional vector of time points (defaults to 1:max_time)
#' @return Data frame in person-time format
expand_person_time <- function(time_var, event_var, A, W,
                               max_time = NULL, time_grid = NULL) {
  if (is.null(max_time)) max_time <- max(time_var, na.rm = TRUE)
  if (is.null(time_grid)) time_grid <- seq_len(max_time)

  W_df <- as.data.frame(W)
  n <- length(time_var)
  records <- vector("list", n)

  for (i in seq_len(n)) {
    t_i <- min(time_var[i], max_time)
    n_periods <- sum(time_grid <= t_i)
    if (n_periods < 1) next

    periods <- time_grid[seq_len(n_periods)]

    rec <- data.frame(
      id   = i,
      time = periods,
      A    = A[i],
      stringsAsFactors = FALSE
    )

    # Hazard indicator: 1 only at the event time
    rec$Y_hazard <- 0L
    if (event_var[i] == 1 && time_var[i] <= max_time) {
      last_row <- which(rec$time == max(periods[periods <= time_var[i]]))
      if (length(last_row) > 0) rec$Y_hazard[max(last_row)] <- 1L
    }

    # Censoring indicator: 1 if censored at this time
    rec$C <- 0L
    if (event_var[i] == 0 && time_var[i] <= max_time) {
      last_row <- nrow(rec)
      rec$C[last_row] <- 1L
    }

    # Add covariates
    rec <- cbind(rec, W_df[rep(i, nrow(rec)), , drop = FALSE])
    rownames(rec) <- NULL
    records[[i]] <- rec
  }

  do.call(rbind, records)
}

#' Internal: targeted discrete-time survival TMLE core (single failure type)
#'
#' Fits initial main-terms logistic working models for the event hazard,
#' censoring hazard and treatment mechanism, then TARGETS the event hazard by
#' iterating a logistic fluctuation along the clever covariate
#'   H_a(t) = -I(A = a) / (g_a(W) * G_c(t- | a, W)) * S_a(t0 | W) / S_a(t | W)
#' until the fluctuation coefficient converges to 0, so the plug-in solves the
#' efficient influence function for the counterfactual risk F_a(t0) = 1 - S_a(t0).
#' Standard errors are the influence-curve-based (analytic) SEs.
#'
#' Matrices are oriented [time-grid x subject]; cumulative products run down
#' columns (over time).
#'
#' @keywords internal
.discrete_surv_tmle_core <- function(time_b, event_c, A_c, W_mat, grid, target_times,
                                     gbound = 0.025, cbound = 0.02,
                                     max_iter = 50, tol = 1e-6) {
  n <- length(time_b)
  K <- length(grid)
  tix <- match(time_b, grid)                       # observed grid index per subject
  if (anyNA(tix)) return(NULL)

  # Observed person-time (long): subject i contributes rows j = 1..tix[i]
  ri <- rep(seq_len(n), tix)
  rj <- sequence(tix)
  yhaz  <- as.integer(rj == tix[ri] & event_c[ri] == 1L)  # event at this period
  cflag <- as.integer(rj == tix[ri] & event_c[ri] == 0L)  # censored at this period
  day   <- grid[rj]

  # --- Initial working models (main-terms logistic GLM) --------------------
  # Parametric GLMs (not SuperLearner): with only a few dozen events on the
  # discrete grid the SL library (esp. SL.glmnet) is unstable and could shrink
  # the treatment coefficient to exactly 0, which previously collapsed the
  # matched-cohort risk difference to a spurious 0. The primary survival
  # estimator (survtmle) retains the SuperLearner fit.
  keep <- cflag == 0L                              # event-hazard risk set excludes censor rows
  hfit <- suppressWarnings(glm(y ~ ., family = binomial(),
    data = data.frame(y = yhaz[keep], time = day[keep], A = A_c[ri][keep],
                      W_mat[ri, , drop = FALSE][keep, , drop = FALSE])))
  cfit <- suppressWarnings(glm(y ~ ., family = binomial(),
    data = data.frame(y = cflag, time = day, A = A_c[ri], W_mat[ri, , drop = FALSE])))
  gfit <- suppressWarnings(glm(A ~ ., family = binomial(),
    data = data.frame(A = A_c, W_mat)))
  g1 <- pmin(pmax(as.numeric(predict(gfit, type = "response")), gbound), 1 - gbound)

  # Counterfactual full-grid hazard predictions (K x n matrix per arm)
  predmat <- function(fit, a) {
    fi <- rep(seq_len(n), each = K)
    dd <- rep(grid, n)
    px <- data.frame(time = dd, A = a, W_mat[fi, , drop = FALSE])
    matrix(pmin(pmax(as.numeric(predict(fit, newdata = px, type = "response")),
                     1e-6), 1 - 1e-6), nrow = K)
  }
  arms <- list()
  for (a in c(0, 1)) {
    lam  <- predmat(hfit, a)                        # event hazard
    lamc <- predmat(cfit, a)                        # censoring hazard
    GcS    <- apply(1 - lamc, 2, cumprod)
    Gcprev <- rbind(1, GcS[-K, , drop = FALSE])     # censoring survival just before t
    arms[[as.character(a)]] <- list(lam = lam, Gcprev = pmax(Gcprev, cbound))
  }

  est <- vector("list", length(target_times))
  for (ti in seq_along(target_times)) {
    tt <- target_times[ti]
    m0 <- match(tt, grid)
    risks <- c()
    ics   <- list()

    for (a in c(1, 0)) {
      L <- arms[[as.character(a)]]
      lam_t <- L$lam
      Gcp   <- L$Gcprev
      pia   <- if (a == 1) g1 else 1 - g1

      # Observed at-risk uncensored score rows for this arm/horizon
      mask <- (A_c[ri] == a) & (rj <= m0) & (cflag == 0L)
      ri_s <- ri[mask]
      rj_s <- rj[mask]
      y_s  <- yhaz[mask]

      # Clever covariate matrix (m0 x n): H_a(t,i)
      clever <- function(St) {
        Sm0 <- St[m0, ]
        -sweep(sweep(1 / St[seq_len(m0), , drop = FALSE], 2, Sm0, "*") /
                 Gcp[seq_len(m0), , drop = FALSE], 2, pia, "/")
      }

      # --- Targeting: iterate the logistic fluctuation to convergence ---
      for (iter in seq_len(max_iter)) {
        St  <- apply(1 - lam_t, 2, cumprod)
        h   <- clever(St)
        off <- qlogis(pmin(pmax(lam_t[cbind(rj_s, ri_s)], 1e-6), 1 - 1e-6))
        hh  <- h[cbind(rj_s, ri_s)]
        eps <- 0
        if (length(y_s) > 2 && sd(hh) > 0) {
          fe <- suppressWarnings(try(
            glm(y_s ~ -1 + hh + offset(off), family = binomial()), silent = TRUE))
          if (!inherits(fe, "try-error")) {
            eps <- as.numeric(coef(fe)[1])
            if (is.na(eps)) eps <- 0
          }
        }
        lam_t[seq_len(m0), ] <- plogis(
          qlogis(pmin(pmax(lam_t[seq_len(m0), , drop = FALSE], 1e-6), 1 - 1e-6)) + eps * h)
        if (abs(eps) < tol) break
      }

      St     <- apply(1 - lam_t, 2, cumprod)
      Sm0t   <- St[m0, ]
      risk_a <- mean(1 - Sm0t)

      # Influence curve for F_a(t0): score part + W part, centred
      h       <- clever(St)
      resid   <- y_s - lam_t[cbind(rj_s, ri_s)]
      contrib <- h[cbind(rj_s, ri_s)] * resid
      Dout <- numeric(n)
      if (length(ri_s) > 0) {
        agg <- rowsum(contrib, ri_s)
        Dout[as.integer(rownames(agg))] <- agg[, 1]
      }
      ics[[as.character(a)]] <- Dout + (1 - Sm0t) - risk_a
      risks[as.character(a)] <- risk_a
    }

    r1 <- risks["1"]
    r0 <- risks["0"]
    icRD <- ics[["1"]] - ics[["0"]]
    se <- sqrt(var(icRD) / n)
    est[[ti]] <- data.frame(
      time            = tt,
      risk_1          = r1,
      risk_0          = r0,
      risk_difference = r1 - r0,
      risk_ratio      = if (isTRUE(r0 > 0)) r1 / r0 else NA_real_,
      rd_se           = se,
      rd_ci_lower     = (r1 - r0) - 1.96 * se,
      rd_ci_upper     = (r1 - r0) + 1.96 * se,
      row.names       = NULL,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, est)
}

#' Discrete-time survival TMLE (targeted)
#'
#' Targeted maximum-likelihood estimator of the counterfactual cumulative
#' incidence (risk) at each target time, using discrete-time hazards. See
#' \code{.discrete_surv_tmle_core} for the targeting step and influence-curve
#' inference. Unlike a plug-in g-computation, the event hazard is fluctuated
#' toward the efficient influence function using the treatment and censoring
#' mechanisms, so the estimate cannot collapse to a spurious null when the
#' outcome model happens to drop the treatment term.
#'
#' @param time_var Observed event/censor times
#' @param event_var Binary event indicator (1 = event, 0 = censored)
#' @param A Treatment vector
#' @param W Covariate matrix
#' @param target_times Vector of time points at which to estimate risk
#' @param sl_lib,sl_lib_censor Retained for API compatibility; unused. The
#'   discrete-time estimator uses parametric working models (see Details); the
#'   primary survtmle estimator is where the SuperLearner library applies.
#' @param seed Random seed (the estimator is deterministic; kept for API parity)
#' @param cv_folds,n_boot Retained for API compatibility; unused (analytic
#'   influence-curve SEs replace the former bootstrap).
#' @return List with counterfactual risk estimates and IC-based SEs
run_survival_tmle <- function(time_var, event_var, A, W,
                              target_times = c(30, 90, 180),
                              sl_lib = c("SL.glm", "SL.glmnet", "SL.mean"),
                              sl_lib_censor = c("SL.glm", "SL.mean"),
                              seed = 42, cv_folds = 2, n_boot = 20) {
  set.seed(seed)

  # Complete cases
  cc <- complete.cases(time_var, event_var, A)
  time_cc  <- time_var[cc]
  event_cc <- as.integer(event_var[cc])
  A_cc     <- as.integer(A[cc])
  W_cc     <- W[cc, , drop = FALSE]

  # Remove zero-variance covariates
  col_var <- apply(W_cc, 2, var, na.rm = TRUE)
  keep_cols <- !is.na(col_var) & col_var > 1e-8
  if (!all(keep_cols)) {
    message("  Removing ", sum(!keep_cols), " zero-variance covariate(s) from survival TMLE: ",
            paste(colnames(W_cc)[!keep_cols], collapse = ", "))
    W_cc <- W_cc[, keep_cols, drop = FALSE]
  }
  W_mat <- as.matrix(W_cc)
  n <- sum(cc)
  max_time <- max(target_times)

  # Discretise time to a weekly grid (plus the target horizons)
  time_grid <- sort(unique(c(seq(1, max_time, by = 7), target_times)))
  time_grid <- time_grid[time_grid <= max_time]

  # Bin observed times to the grid
  time_binned <- sapply(time_cc, function(t)
    max(time_grid[time_grid <= max(t, min(time_grid))]))

  estimates <- tryCatch(
    .discrete_surv_tmle_core(time_binned, event_cc, A_cc, W_mat,
                             time_grid, target_times),
    error = function(e) {
      warning("Discrete-time survival TMLE failed: ", e$message)
      NULL
    })
  if (is.null(estimates)) return(NULL)

  list(
    method    = "Survival TMLE (discrete-time)",
    estimates = estimates,
    n         = n,
    max_time  = max_time,
    time_grid = time_grid
  )
}

#' Run Cox proportional hazards regression
#'
#' @param time_var Survival time
#' @param event_var Event indicator
#' @param A Treatment
#' @param W Covariate matrix
#' @param matched_idx Optional matched indices
#' @return List with Cox model results
run_cox_regression <- function(time_var, event_var, A, W, matched_idx = NULL) {
  require(survival)

  if (!is.null(matched_idx)) {
    time_var  <- time_var[matched_idx]
    event_var <- event_var[matched_idx]
    A         <- A[matched_idx]
    W         <- W[matched_idx, , drop = FALSE]
  }

  cc <- complete.cases(time_var, event_var, A)
  time_cc  <- time_var[cc]
  event_cc <- event_var[cc]
  A_cc     <- A[cc]
  W_cc     <- as.data.frame(W[cc, , drop = FALSE])

  df <- data.frame(time = time_cc, event = event_cc, A = A_cc, W_cc)

  # Fit Cox model
  cox_fit <- tryCatch(
    coxph(Surv(time, event) ~ ., data = df),
    error = function(e) {
      warning("Full Cox model failed, using treatment-only: ", e$message)
      coxph(Surv(time, event) ~ A, data = df)
    }
  )

  # Extract HR for treatment
  coefs <- summary(cox_fit)$coefficients
  hr_row <- if ("A" %in% rownames(coefs)) coefs["A", ] else coefs[1, ]

  # PH test
  ph_test <- tryCatch(cox.zph(cox_fit), error = function(e) NULL)

  list(
    method    = if (is.null(matched_idx)) "Cox PH regression" else "Matched Cox PH regression",
    fit       = cox_fit,
    hr        = exp(hr_row["coef"]),
    hr_ci     = exp(confint(cox_fit)["A", ]),
    log_hr    = hr_row["coef"],
    log_hr_se = hr_row["se(coef)"],
    pvalue    = hr_row["Pr(>|z|)"],
    ph_test   = ph_test,
    n         = sum(cc),
    n_events  = sum(event_cc)
  )
}

#' Run matched survival TMLE
run_matched_survival_tmle <- function(time_var, event_var, A, W,
                                      matched_idx, target_times,
                                      sl_lib, sl_lib_censor, seed = 42,
                                      cv_folds = 2, n_boot = 20) {
  result <- run_survival_tmle(
    time_var  = time_var[matched_idx],
    event_var = event_var[matched_idx],
    A         = A[matched_idx],
    W         = W[matched_idx, , drop = FALSE],
    target_times  = target_times,
    sl_lib        = sl_lib,
    sl_lib_censor = sl_lib_censor,
    seed          = seed,
    cv_folds      = cv_folds,
    n_boot        = n_boot
  )
  if (!is.null(result)) {
    result$method <- "Matched Survival TMLE (discrete-time)"
  }
  result
}

#' Summarize survival analysis results
summarize_survival_results <- function(results_list) {
  rows <- list()
  for (r in results_list) {
    if (is.null(r)) next
    if (!is.null(r$estimates) && is.data.frame(r$estimates) && nrow(r$estimates) > 0) {
      # Survival TMLE results (both custom and survtmle)
      est <- r$estimates
      n_val <- if (!is.null(r$n)) r$n
               else if ("n" %in% names(est)) est$n[1]
               else NA
      for (i in seq_len(nrow(est))) {
        rows[[length(rows) + 1]] <- data.frame(
          method  = if (!is.null(est$method)) est$method[i] else r$method,
          time    = est$time[i],
          risk_1  = est$risk_1[i],
          risk_0  = est$risk_0[i],
          risk_difference = est$risk_difference[i],
          risk_ratio = est$risk_ratio[i],
          rd_se   = if ("rd_se" %in% names(est)) est$rd_se[i] else NA,
          n       = n_val,
          stringsAsFactors = FALSE
        )
      }
    } else if (!is.null(r$hr)) {
      # Cox model results
      rows[[length(rows) + 1]] <- data.frame(
        method  = r$method,
        time    = NA,
        risk_1  = NA,
        risk_0  = NA,
        risk_difference = NA,
        risk_ratio = r$hr,
        rd_se   = r$log_hr_se,
        n       = r$n,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}
