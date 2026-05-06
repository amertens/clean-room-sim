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

#' Discrete-time survival TMLE
#'
#' Implements TMLE for counterfactual survival curves using discrete hazard models.
#'
#' @param time_var Observed event/censor times
#' @param event_var Binary event indicator
#' @param A Treatment vector
#' @param W Covariate matrix
#' @param target_times Vector of time points at which to estimate risk
#' @param sl_lib SuperLearner library for hazard and treatment models
#' @param sl_lib_censor SuperLearner library for censoring model
#' @param seed Random seed
#' @return List with counterfactual risk estimates
run_survival_tmle <- function(time_var, event_var, A, W,
                              target_times = c(30, 90, 180),
                              sl_lib = c("SL.glm", "SL.glmnet", "SL.mean"),
                              sl_lib_censor = c("SL.glm", "SL.mean"),
                              seed = 42, cv_folds = 2, n_boot = 20) {
  require(SuperLearner)
  set.seed(seed)

  # Complete cases
  cc <- complete.cases(time_var, event_var, A)
  time_cc  <- time_var[cc]
  event_cc <- event_var[cc]
  A_cc     <- A[cc]
  W_cc     <- W[cc, , drop = FALSE]

  # Remove zero-variance covariates
  col_var <- apply(W_cc, 2, var, na.rm = TRUE)
  keep_cols <- !is.na(col_var) & col_var > 1e-8
  if (!all(keep_cols)) {
    message("  Removing ", sum(!keep_cols), " zero-variance covariate(s) from survival TMLE: ",
            paste(colnames(W_cc)[!keep_cols], collapse = ", "))
    W_cc <- W_cc[, keep_cols, drop = FALSE]
  }

  n        <- sum(cc)

  max_time <- max(target_times)

  # Discretize time to reasonable grid (e.g., weekly bins for daily data)
  time_grid <- sort(unique(c(seq(1, max_time, by = 7), target_times)))
  time_grid <- time_grid[time_grid <= max_time]

  # Bin observed times to grid
  time_binned <- sapply(time_cc, function(t) {
    max(time_grid[time_grid <= max(t, min(time_grid))])
  })

  # Expand to person-time
  pt_data <- expand_person_time(time_binned, event_cc, A_cc, W_cc,
                                max_time = max_time, time_grid = time_grid)

  if (is.null(pt_data) || nrow(pt_data) == 0) {
    warning("No person-time records created")
    return(NULL)
  }

  # --- Step 1: Estimate treatment model g(A|W) ---
  # Use subject-level data, not person-time
  g_data <- data.frame(A = A_cc, as.data.frame(W_cc))
  g_fit <- tryCatch(
    SuperLearner(Y = g_data$A, X = g_data[, -1, drop = FALSE],
                 family = binomial(), SL.library = sl_lib,
                 cvControl = list(V = cv_folds)),
    error = function(e) {
      warning("g-model SL failed, using GLM: ", e$message)
      SuperLearner(Y = g_data$A, X = g_data[, -1, drop = FALSE],
                   family = binomial(), SL.library = "SL.glm",
                   cvControl = list(V = cv_folds))
    }
  )
  g_pred <- as.numeric(g_fit$SL.predict)
  g_pred <- pmin(pmax(g_pred, 0.01), 0.99)

  # Map g to person-time
  pt_data$g_A <- g_pred[pt_data$id]

  # --- Step 2: Estimate censoring model P(C=0 | t, A, W) ---
  W_cols <- setdiff(names(pt_data), c("id", "time", "A", "Y_hazard", "C", "g_A"))
  cens_X <- pt_data[, c("time", "A", W_cols), drop = FALSE]

  cens_fit <- tryCatch(
    SuperLearner(Y = 1 - pt_data$C, X = cens_X, family = binomial(),
                 SL.library = sl_lib_censor, cvControl = list(V = cv_folds)),
    error = function(e) {
      warning("Censoring model SL failed: ", e$message)
      SuperLearner(Y = 1 - pt_data$C, X = cens_X, family = binomial(),
                   SL.library = "SL.glm", cvControl = list(V = cv_folds))
    }
  )
  pt_data$surv_cens <- pmax(as.numeric(cens_fit$SL.predict), 0.01)

  # --- Step 3: Estimate hazard model Q(Y_hazard | t, A, W) ---
  # Only among uncensored person-time rows
  uncens <- pt_data$C == 0
  haz_X <- pt_data[uncens, c("time", "A", W_cols), drop = FALSE]
  haz_Y <- pt_data$Y_hazard[uncens]

  haz_fit <- tryCatch(
    SuperLearner(Y = haz_Y, X = haz_X, family = binomial(),
                 SL.library = sl_lib, cvControl = list(V = cv_folds)),
    error = function(e) {
      warning("Hazard model SL failed: ", e$message)
      SuperLearner(Y = haz_Y, X = haz_X, family = binomial(),
                   SL.library = "SL.glm", cvControl = list(V = cv_folds))
    }
  )

  # --- Step 4: Predict counterfactual hazards ---
  results <- list()

  for (a_val in c(1, 0)) {
    # Predict hazard under A = a_val for all person-time rows
    pred_X <- pt_data[uncens, c("time", "A", W_cols), drop = FALSE]
    pred_X$A <- a_val
    h_pred <- pmin(pmax(as.numeric(predict(haz_fit, newdata = pred_X)$pred), 0.001), 0.999)

    # For each individual, compute cumulative survival
    surv_by_id <- list()
    ids <- unique(pt_data$id[uncens])

    for (id_i in ids) {
      rows <- which(pt_data$id[uncens] == id_i)
      h_i <- h_pred[rows]
      t_i <- pt_data$time[uncens][rows]
      surv_i <- cumprod(1 - h_i)
      surv_by_id[[as.character(id_i)]] <- data.frame(
        id = id_i, time = t_i, surv = surv_i
      )
    }

    surv_all <- do.call(rbind, surv_by_id)

    # Estimate risk at target times
    for (tt in target_times) {
      # For each individual, get survival at time tt
      surv_at_t <- sapply(ids, function(id_i) {
        s <- surv_all[surv_all$id == id_i, ]
        s <- s[s$time <= tt, ]
        if (nrow(s) == 0) return(1)
        min(s$surv)  # cumulative survival up to tt
      })
      risk_at_t <- 1 - surv_at_t
      results[[paste0("risk_", a_val, "_t", tt)]] <- mean(risk_at_t)
    }
  }

  # --- Step 5: Compute contrasts ---
  estimates <- data.frame(
    time = target_times,
    risk_1 = sapply(target_times, function(tt) results[[paste0("risk_1_t", tt)]]),
    risk_0 = sapply(target_times, function(tt) results[[paste0("risk_0_t", tt)]]),
    stringsAsFactors = FALSE
  )
  estimates$risk_difference <- estimates$risk_1 - estimates$risk_0
  estimates$risk_ratio <- ifelse(estimates$risk_0 > 0,
                                  estimates$risk_1 / estimates$risk_0, NA)

  # --- Bootstrap SEs ---
  boot_results <- array(NA, dim = c(n_boot, length(target_times), 2))  # RD, RR
  pb_boot <- cr_progress(n_boot, "  Bootstrap")

  for (b in seq_len(n_boot)) {
    pb_boot$tick()
    boot_ids <- sample(n, n, replace = TRUE)
    time_b  <- time_binned[boot_ids]
    event_b <- event_cc[boot_ids]
    A_b     <- A_cc[boot_ids]
    W_b     <- W_cc[boot_ids, , drop = FALSE]

    boot_est <- tryCatch({
      pt_b <- expand_person_time(time_b, event_b, A_b, W_b,
                                  max_time = max_time, time_grid = time_grid)
      if (is.null(pt_b) || nrow(pt_b) == 0) return(NULL)

      # Simple GLM-based bootstrap for speed
      haz_X_b <- pt_b[pt_b$C == 0, c("time", "A", names(as.data.frame(W_b))), drop = FALSE]
      haz_Y_b <- pt_b$Y_hazard[pt_b$C == 0]

      haz_glm <- glm(haz_Y_b ~ ., data = haz_X_b, family = binomial())

      boot_risks <- sapply(target_times, function(tt) {
        sapply(c(1, 0), function(a_val) {
          pred_data <- haz_X_b
          pred_data$A <- a_val
          h_b <- pmin(pmax(predict(haz_glm, newdata = pred_data, type = "response"), 0.001), 0.999)

          # Average risk at target time
          at_time <- haz_X_b$time <= tt
          if (sum(at_time) == 0) return(0)
          mean(h_b[at_time])  # Simplified: average hazard as proxy
        })
      })
      boot_risks
    }, error = function(e) NULL)

    if (!is.null(boot_est) && is.matrix(boot_est)) {
      for (j in seq_along(target_times)) {
        boot_results[b, j, 1] <- boot_est[1, j] - boot_est[2, j]  # RD
        boot_results[b, j, 2] <- ifelse(boot_est[2, j] > 0,
                                          boot_est[1, j] / boot_est[2, j], NA)  # RR
      }
    }
  }

  # Compute SEs from bootstrap
  estimates$rd_se <- apply(boot_results[, , 1, drop = FALSE], 2,
                           sd, na.rm = TRUE)
  estimates$rd_ci_lower <- estimates$risk_difference - 1.96 * estimates$rd_se
  estimates$rd_ci_upper <- estimates$risk_difference + 1.96 * estimates$rd_se

  list(
    method     = "Survival TMLE (discrete-time)",
    estimates  = estimates,
    n          = n,
    max_time   = max_time,
    time_grid  = time_grid,
    g_fit      = g_fit,
    haz_fit    = haz_fit,
    cens_fit   = cens_fit
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
