# ============================================================
# Clean-Room Workflow: Outcome-Blind Simulation Functions
# ============================================================

#' Simulate binary outcomes (logistic model with nonlinear effects)
#'
#' Uses observed covariates and treatment but generates SYNTHETIC outcomes.
#' No real outcome data is used.
#'
#' @param A Treatment vector
#' @param W Covariate matrix
#' @param true_ate True average treatment effect to embed
#' @param seed Random seed
#' @return Simulated binary outcome vector
simulate_binary_outcome <- function(A, W, true_ate = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  W_mat <- as.matrix(W)
  n <- length(A)
  p <- ncol(W_mat)

  # Generate random coefficients for covariates
  set.seed(seed)
  k <- min(p, 10)
  beta <- rnorm(k, 0, 0.3)

  # Use first k covariates — scale to [0,1] so nonlinear terms stay bounded
  W_sub <- W_mat[, seq_len(k), drop = FALSE]
  for (j in seq_len(k)) {
    rng <- range(W_sub[, j], na.rm = TRUE)
    if (diff(rng) > 0) {
      W_sub[, j] <- (W_sub[, j] - rng[1]) / diff(rng)
    } else {
      W_sub[, j] <- 0
    }
  }

  # Linear predictor with nonlinear terms (bounded inputs keep lp finite)
  lp <- -1.5 +
    W_sub %*% beta +
    0.3 * sin(W_sub[, 1]) +
    ifelse(k >= 2, 0.2 * W_sub[, 2]^2, 0) +
    qlogis(0.5 + true_ate / 2) * A

  # Guard against NaN / Inf in lp
  lp[is.na(lp) | is.nan(lp)] <- 0
  lp[lp >  10] <-  10
  lp[lp < -10] <- -10

  prob <- plogis(as.numeric(lp))
  Y <- rbinom(n, 1, prob)

  list(Y = Y, true_ate = true_ate, prob = prob)
}

#' Simulate discrete survival outcomes
#'
#' Generates synthetic survival times from a discrete hazard model.
#'
#' @param A Treatment vector
#' @param W Covariate matrix
#' @param max_time Maximum follow-up time
#' @param true_risk_diff True risk difference at max_time
#' @param censor_rate Probability of censoring per time unit
#' @param seed Random seed
#' @return List with time, event, and true parameters
simulate_survival_outcome <- function(A, W, max_time = 180,
                                       true_risk_diff = 0.03,
                                       censor_rate = 0.005,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  W_mat <- as.matrix(W)
  n <- length(A)
  p <- ncol(W_mat)

  k <- min(p, 10)
  beta <- rnorm(k, 0, 0.1)
  W_sub <- W_mat[, seq_len(k), drop = FALSE]
  # Scale to [0,1] for numerical stability
  for (j in seq_len(k)) {
    rng <- range(W_sub[, j], na.rm = TRUE)
    if (diff(rng) > 0) {
      W_sub[, j] <- (W_sub[, j] - rng[1]) / diff(rng)
    } else {
      W_sub[, j] <- 0
    }
  }

  # Baseline hazard (low, so events accumulate gradually)
  base_haz <- 0.002

  event_time <- rep(max_time, n)
  event_ind  <- rep(0L, n)

  for (i in seq_len(n)) {
    lp_i <- sum(W_sub[i, ] * beta) + log(1 - true_risk_diff * 5) * A[i]
    for (t in seq_len(max_time)) {
      h_it <- plogis(qlogis(base_haz) + lp_i + 0.001 * t)
      if (runif(1) < h_it) {
        event_time[i] <- t
        event_ind[i]  <- 1L
        break
      }
    }
  }

  # Apply random censoring
  censor_time <- rgeom(n, censor_rate) + 1
  censored <- censor_time < event_time
  obs_time  <- pmin(event_time, censor_time)
  obs_event <- ifelse(censored, 0L, event_ind)

  list(
    time           = obs_time,
    event          = obs_event,
    true_risk_diff = true_risk_diff,
    max_time       = max_time,
    n_events       = sum(obs_event),
    n_censored     = sum(censored),
    event_rate     = mean(obs_event)
  )
}

#' Evaluate a single estimator on simulated data
#'
#' @param estimator_fn Function that takes (Y, A, W, ...) and returns estimate
#' @param sim_data Simulated data (list with Y or time/event)
#' @param A Treatment
#' @param W Covariates
#' @param true_param True parameter value
#' @param ... Additional arguments passed to estimator_fn
#' @return List with bias, variance, RMSE, coverage
evaluate_estimator <- function(estimator_fn, sim_data, A, W,
                               true_param, ...) {
  result <- tryCatch(
    estimator_fn(sim_data = sim_data, A = A, W = W, ...),
    error = function(e) {
      list(estimate = NA, ci_lower = NA, ci_upper = NA, se = NA)
    }
  )

  est <- result$estimate
  ci_lo <- result$ci_lower
  ci_hi <- result$ci_upper

  list(
    estimate = est,
    bias     = est - true_param,
    ci_lower = ci_lo,
    ci_upper = ci_hi,
    covered  = !is.na(ci_lo) && !is.na(ci_hi) &&
               ci_lo <= true_param && ci_hi >= true_param,
    se       = result$se
  )
}

#' Run outcome-blind simulation study
#'
#' @param A Observed treatment vector
#' @param W Observed covariate matrix
#' @param matched_idx Matched observation indices
#' @param n_sims Number of simulation repetitions
#' @param sl_lib SuperLearner library
#' @param config Configuration list
#' @param seed Random seed
#' @return List with simulation results for binary and survival
run_simulation_study <- function(A, W, matched_idx, n_sims = 200,
                                  sl_lib, config, seed = 42) {
  require(SuperLearner)
  require(tmle)

  set.seed(seed)
  true_ate_bin  <- config$simulation$true_ate_binary
  true_ate_surv <- config$simulation$true_ate_survival

  # --- Binary outcome simulations ---
  binary_results <- data.frame(
    sim       = integer(0),
    estimator = character(0),
    estimate  = numeric(0),
    bias      = numeric(0),
    covered   = logical(0),
    stringsAsFactors = FALSE
  )

  cc <- !is.na(A)
  A_cc <- A[cc]
  W_cc <- as.matrix(W[cc, , drop = FALSE])
  # Ensure safe column names and no residual NAs
  colnames(W_cc) <- make.names(colnames(W_cc), unique = TRUE)
  W_cc[is.na(W_cc) | is.nan(W_cc) | is.infinite(W_cc)] <- 0

  # Remap matched_idx to complete-case indices
  cc_map <- which(cc)
  matched_cc <- match(matched_idx, cc_map)
  matched_cc <- matched_cc[!is.na(matched_cc)]

  # Pre-build data.frame once (avoids repeated conversion in loop)
  W_df <- as.data.frame(W_cc)

  pb_bin <- cr_progress(n_sims, "Binary sim")
  for (s in seq_len(n_sims)) {
    pb_bin$tick()

    sim <- simulate_binary_outcome(A_cc, W_cc, true_ate = true_ate_bin,
                                    seed = seed + s)
    Y_sim <- sim$Y

    # Estimator 1: PS-matched regression
    est1 <- tryCatch({
      Y_m <- Y_sim[matched_cc]; A_m <- A_cc[matched_cc]
      rd <- mean(Y_m[A_m == 1]) - mean(Y_m[A_m == 0])
      se <- sqrt(var(Y_m[A_m == 1]) / sum(A_m == 1) + var(Y_m[A_m == 0]) / sum(A_m == 0))
      list(estimate = rd, ci_lower = rd - 1.96 * se, ci_upper = rd + 1.96 * se, se = se)
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA, se = NA))

    # Estimator 2: PS-matched TMLE
    est2 <- tryCatch({
      fit <- tmle::tmle(Y = Y_sim[matched_cc], A = A_cc[matched_cc],
                        W = W_df[matched_cc, , drop = FALSE],
                        family = "binomial", Q.SL.library = "SL.glm",
                        g.SL.library = "SL.glm", verbose = FALSE)
      ate <- fit$estimates$ATE
      list(estimate = ate$psi, ci_lower = ate$CI[1], ci_upper = ate$CI[2],
           se = sqrt(ate$var.psi))
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA, se = NA))

    # Estimator 3: Full-cohort TMLE
    est3 <- tryCatch({
      fit <- tmle::tmle(Y = Y_sim, A = A_cc,
                        W = W_df,
                        family = "binomial", Q.SL.library = "SL.glm",
                        g.SL.library = "SL.glm", verbose = FALSE)
      ate <- fit$estimates$ATE
      list(estimate = ate$psi, ci_lower = ate$CI[1], ci_upper = ate$CI[2],
           se = sqrt(ate$var.psi))
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA, se = NA))

    for (est_info in list(
      list(name = "PS-matched regression", est = est1),
      list(name = "PS-matched TMLE", est = est2),
      list(name = "Full-cohort TMLE", est = est3)
    )) {
      e <- est_info$est
      covered <- !is.na(e$ci_lower) && !is.na(e$ci_upper) &&
                 e$ci_lower <= true_ate_bin && e$ci_upper >= true_ate_bin
      binary_results <- rbind(binary_results, data.frame(
        sim = s, estimator = est_info$name,
        estimate = e$estimate, bias = e$estimate - true_ate_bin,
        covered = covered, stringsAsFactors = FALSE
      ))
    }
  }

  # --- Survival outcome simulations (fewer reps for speed) ---
  n_surv_sims <- min(n_sims, 50)
  survival_results <- data.frame(
    sim = integer(0), estimator = character(0),
    estimate = numeric(0), bias = numeric(0),
    covered = logical(0), stringsAsFactors = FALSE
  )

  pb_surv <- cr_progress(n_surv_sims, "Survival sim")
  for (s in seq_len(n_surv_sims)) {
    pb_surv$tick()

    sim_surv <- simulate_survival_outcome(
      A_cc, W_cc, max_time = 180,
      true_risk_diff = true_ate_surv, seed = seed + 10000 + s
    )

    # Estimator E: Cox PH (full) — use first 5 covariates for tractability
    n_cox_vars <- min(5, ncol(W_df))
    estE <- tryCatch({
      df <- data.frame(time = sim_surv$time, event = sim_surv$event,
                        A = A_cc, W_df[, seq_len(n_cox_vars), drop = FALSE])
      cfit <- survival::coxph(survival::Surv(time, event) ~ ., data = df)
      hr <- exp(coef(cfit)["A"])
      ci <- exp(confint(cfit)["A", ])
      list(estimate = hr, ci_lower = ci[1], ci_upper = ci[2], se = NA)
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA, se = NA))

    # Estimator F: Matched Cox
    estF <- tryCatch({
      df <- data.frame(time = sim_surv$time[matched_cc],
                        event = sim_surv$event[matched_cc],
                        A = A_cc[matched_cc],
                        W_df[matched_cc, seq_len(n_cox_vars), drop = FALSE])
      cfit <- survival::coxph(survival::Surv(time, event) ~ ., data = df)
      hr <- exp(coef(cfit)["A"])
      ci <- exp(confint(cfit)["A", ])
      list(estimate = hr, ci_lower = ci[1], ci_upper = ci[2], se = NA)
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA, se = NA))

    for (est_info in list(
      list(name = "Cox PH regression", est = estE),
      list(name = "Matched Cox PH regression", est = estF)
    )) {
      e <- est_info$est
      survival_results <- rbind(survival_results, data.frame(
        sim = s, estimator = est_info$name,
        estimate = e$estimate, bias = NA,
        covered = NA, stringsAsFactors = FALSE
      ))
    }
  }

  # --- Aggregate ---
  binary_summary <- aggregate(
    cbind(estimate, bias, covered) ~ estimator,
    data = binary_results,
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE))
    }
  )

  survival_summary <- aggregate(
    estimate ~ estimator,
    data = survival_results,
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE))
    }
  )

  # Compute per-estimator metrics
  binary_metrics <- do.call(rbind, lapply(unique(binary_results$estimator), function(est) {
    d <- binary_results[binary_results$estimator == est, ]
    data.frame(
      estimator = est,
      mean_estimate = mean(d$estimate, na.rm = TRUE),
      bias = mean(d$bias, na.rm = TRUE),
      variance = var(d$estimate, na.rm = TRUE),
      rmse = sqrt(mean(d$bias^2, na.rm = TRUE)),
      coverage = mean(d$covered, na.rm = TRUE),
      n_valid = sum(!is.na(d$estimate)),
      stringsAsFactors = FALSE
    )
  }))

  list(
    binary_results    = binary_results,
    binary_metrics    = binary_metrics,
    survival_results  = survival_results,
    survival_summary  = survival_summary,
    true_ate_binary   = true_ate_bin,
    true_ate_survival = true_ate_surv,
    n_sims_binary     = n_sims,
    n_sims_survival   = n_surv_sims
  )
}
