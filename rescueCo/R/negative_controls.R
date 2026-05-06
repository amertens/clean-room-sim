# ============================================================
# Clean-Room Workflow: Negative Control & Plasmode Functions
# ============================================================

#' Run TMLE for a single negative-control outcome
#'
#' @param Y_nc Binary negative-control outcome vector
#' @param A Treatment vector
#' @param W Covariate matrix (should NOT include Y_nc)
#' @param sl_lib SuperLearner library
#' @param nc_name Name of the NC outcome (for labeling)
#' @return List with ATE, CI, p-value, and pass indicator
run_negative_control_tmle <- function(Y_nc, A, W, sl_lib, nc_name) {
  require(tmle)
  require(SuperLearner)

  cc <- complete.cases(Y_nc, A)
  Y_cc <- Y_nc[cc]
  A_cc <- A[cc]
  W_cc <- as.data.frame(W[cc, , drop = FALSE])

  fit <- tryCatch({
    tmle::tmle(
      Y = Y_cc, A = A_cc, W = W_cc,
      family = "binomial",
      Q.SL.library = sl_lib,
      g.SL.library = sl_lib,
      verbose = FALSE
    )
  }, error = function(e) {
    tryCatch({
      tmle::tmle(
        Y = Y_cc, A = A_cc, W = W_cc,
        family = "binomial",
        Q.SL.library = "SL.glm",
        g.SL.library = "SL.glm",
        verbose = FALSE
      )
    }, error = function(e2) NULL)
  })

  if (is.null(fit)) {
    return(list(
      nc_name  = nc_name,
      estimate = NA, ci_lower = NA, ci_upper = NA,
      pvalue   = NA, se = NA,
      pass     = NA, error = "TMLE failed"
    ))
  }

  ate <- fit$estimates$ATE
  ci  <- ate$CI
  pass <- ci[1] <= 0 && ci[2] >= 0

  list(
    nc_name  = nc_name,
    estimate = ate$psi,
    ci_lower = ci[1],
    ci_upper = ci[2],
    pvalue   = ate$pvalue,
    se       = sqrt(ate$var.psi),
    pass     = pass,
    error    = NULL
  )
}

#' Run TMLE on all negative-control outcomes
#'
#' @param W Full covariate matrix (NC variables will be extracted from here)
#' @param A Treatment vector
#' @param nc_vars Character vector of NC outcome column names
#' @param sl_lib SuperLearner library
#' @return data.frame with one row per NC outcome
run_all_negative_controls <- function(W, A, nc_vars, sl_lib) {
  results <- vector("list", length(nc_vars))

  for (i in seq_along(nc_vars)) {
    nc_var <- nc_vars[i]

    if (!nc_var %in% colnames(W)) {
      results[[i]] <- data.frame(
        nc_name = nc_var, estimate = NA, ci_lower = NA, ci_upper = NA,
        pvalue = NA, se = NA, pass = NA, error = "Variable not in W",
        stringsAsFactors = FALSE
      )
      next
    }

    # Extract NC outcome from W, remove it from covariate matrix
    Y_nc <- W[, nc_var]
    W_reduced <- W[, setdiff(colnames(W), nc_var), drop = FALSE]

    res <- run_negative_control_tmle(Y_nc, A, W_reduced, sl_lib, nc_var)

    results[[i]] <- data.frame(
      nc_name  = res$nc_name,
      estimate = res$estimate,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      pvalue   = res$pvalue,
      se       = res$se,
      pass     = res$pass,
      error    = res$error %||% "",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}

#' Plasmode simulation: compare estimators under known treatment effect
#'
#' Uses observed covariates, resamples treatment from fitted PS,
#' generates outcomes under known ATE. Compares TMLE, matched regression,
#' and g-computation.
#'
#' @param A Observed treatment vector
#' @param W Covariate matrix
#' @param ps Fitted propensity scores
#' @param matched_idx Matched observation indices
#' @param sl_lib SuperLearner library (use SL.glm for speed)
#' @param n_sims Number of plasmode simulations
#' @param true_ate Known true treatment effect
#' @param seed Random seed
#' @return data.frame with bias, RMSE, coverage by estimator
run_plasmode_simulation <- function(A, W, ps, matched_idx,
                                    sl_lib = "SL.glm",
                                    n_sims = 50, true_ate = 0,
                                    seed = 42) {
  require(tmle)
  set.seed(seed)

  cc <- !is.na(A) & !is.na(ps)
  A_obs <- A[cc]
  W_cc  <- as.data.frame(W[cc, , drop = FALSE])
  ps_cc <- ps[cc]
  n     <- sum(cc)

  # Remap matched_idx to complete-case indices
  cc_map <- which(cc)
  matched_cc <- match(matched_idx, cc_map)
  matched_cc <- matched_cc[!is.na(matched_cc)]

  colnames(W_cc) <- make.names(colnames(W_cc), unique = TRUE)
  W_cc[is.na(W_cc)] <- 0

  results <- data.frame(
    sim = integer(0), estimator = character(0),
    estimate = numeric(0), bias = numeric(0),
    covered = logical(0), stringsAsFactors = FALSE
  )

  for (s in seq_len(n_sims)) {
    # Generate treatment from fitted PS (plasmode)
    A_sim <- rbinom(n, 1, ps_cc)

    # Generate outcome under known effect
    lp <- -1 + 0.3 * scale(W_cc[, 1]) + true_ate * A_sim
    lp <- pmin(pmax(as.numeric(lp), -10), 10)
    Y_sim <- rbinom(n, 1, plogis(lp))

    # Estimator 1: Full-cohort TMLE
    est1 <- tryCatch({
      fit <- tmle::tmle(Y = Y_sim, A = A_sim, W = W_cc,
                        family = "binomial",
                        Q.SL.library = sl_lib, g.SL.library = sl_lib,
                        verbose = FALSE)
      ate <- fit$estimates$ATE
      list(estimate = ate$psi, ci_lower = ate$CI[1], ci_upper = ate$CI[2])
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA))

    # Estimator 2: Matched difference-in-means (using original match structure)
    est2 <- tryCatch({
      if (length(matched_cc) < 10) stop("too few matched")
      Y_m <- Y_sim[matched_cc]; A_m <- A_sim[matched_cc]
      if (sum(A_m == 1) < 2 || sum(A_m == 0) < 2) stop("sparse groups")
      rd <- mean(Y_m[A_m == 1]) - mean(Y_m[A_m == 0])
      se <- sqrt(var(Y_m[A_m == 1]) / sum(A_m == 1) +
                 var(Y_m[A_m == 0]) / sum(A_m == 0))
      list(estimate = rd, ci_lower = rd - 1.96 * se, ci_upper = rd + 1.96 * se)
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA))

    # Estimator 3: G-computation (parametric outcome model)
    est3 <- tryCatch({
      df <- data.frame(Y = Y_sim, A = A_sim, W_cc)
      gfit <- glm(Y ~ ., data = df, family = binomial())
      df1 <- df0 <- df; df1$A <- 1; df0$A <- 0
      rd <- mean(predict(gfit, df1, type = "response")) -
            mean(predict(gfit, df0, type = "response"))
      # Bootstrap SE
      boot_rd <- replicate(200, {
        idx <- sample(n, n, replace = TRUE)
        bfit <- glm(Y ~ ., data = df[idx, ], family = binomial())
        bdf1 <- bdf0 <- df[idx, ]; bdf1$A <- 1; bdf0$A <- 0
        mean(predict(bfit, bdf1, type = "response")) -
          mean(predict(bfit, bdf0, type = "response"))
      })
      list(estimate = rd,
           ci_lower = quantile(boot_rd, 0.025, names = FALSE),
           ci_upper = quantile(boot_rd, 0.975, names = FALSE))
    }, error = function(e) list(estimate = NA, ci_lower = NA, ci_upper = NA))

    for (info in list(
      list(name = "Full-cohort TMLE", est = est1),
      list(name = "PS-matched regression", est = est2),
      list(name = "G-computation", est = est3)
    )) {
      e <- info$est
      covered <- !is.na(e$ci_lower) && !is.na(e$ci_upper) &&
                 e$ci_lower <= true_ate && e$ci_upper >= true_ate
      results <- rbind(results, data.frame(
        sim = s, estimator = info$name,
        estimate = e$estimate, bias = e$estimate - true_ate,
        covered = covered, stringsAsFactors = FALSE
      ))
    }
  }

  # Aggregate metrics
  metrics <- do.call(rbind, lapply(unique(results$estimator), function(est) {
    d <- results[results$estimator == est, ]
    data.frame(
      estimator     = est,
      mean_estimate = mean(d$estimate, na.rm = TRUE),
      bias          = mean(d$bias, na.rm = TRUE),
      rmse          = sqrt(mean(d$bias^2, na.rm = TRUE)),
      coverage      = mean(d$covered, na.rm = TRUE),
      n_valid       = sum(!is.na(d$estimate)),
      stringsAsFactors = FALSE
    )
  }))

  list(raw_results = results, metrics = metrics)
}

#' Evaluate GO / FLAG / STOP gating decision
#'
#' @param nc_results data.frame from run_all_negative_controls
#' @param plasmode_results list from run_plasmode_simulation
#' @param criteria list with nc_max_failures_go, nc_max_failures_flag,
#'   plasmode_min_coverage_go, plasmode_min_coverage_flag, plasmode_max_bias
#' @return list with decision ("GO", "FLAG", "STOP") and rationale
evaluate_go_flag_stop <- function(nc_results, plasmode_results, criteria) {
  # Count NC failures (CI excludes zero)
  nc_failures <- sum(!nc_results$pass, na.rm = TRUE)
  nc_total    <- sum(!is.na(nc_results$pass))

  # Best plasmode TMLE metrics
  pm <- plasmode_results$metrics
  tmle_row <- pm[pm$estimator == "Full-cohort TMLE", ]
  if (nrow(tmle_row) == 0) tmle_row <- pm[1, ]

  pm_coverage <- tmle_row$coverage
  pm_bias     <- abs(tmle_row$bias)

  # Decision logic
  reasons <- character(0)


  if (nc_failures > criteria$nc_max_failures_flag) {
    decision <- "STOP"
    reasons <- c(reasons,
      paste0(nc_failures, "/", nc_total, " NC outcomes rejected null (max for FLAG: ",
             criteria$nc_max_failures_flag, ")"))
  } else if (pm_coverage < criteria$plasmode_min_coverage_flag) {
    decision <- "STOP"
    reasons <- c(reasons,
      paste0("Plasmode coverage = ", round(pm_coverage, 3),
             " < ", criteria$plasmode_min_coverage_flag))
  } else if (nc_failures > criteria$nc_max_failures_go) {
    decision <- "FLAG"
    reasons <- c(reasons,
      paste0(nc_failures, "/", nc_total, " NC outcomes rejected null"))
  } else if (pm_coverage < criteria$plasmode_min_coverage_go) {
    decision <- "FLAG"
    reasons <- c(reasons,
      paste0("Plasmode coverage = ", round(pm_coverage, 3),
             " < ", criteria$plasmode_min_coverage_go))
  } else if (pm_bias > criteria$plasmode_max_bias) {
    decision <- "FLAG"
    reasons <- c(reasons,
      paste0("Plasmode bias = ", round(pm_bias, 4),
             " > ", criteria$plasmode_max_bias))
  } else {
    decision <- "GO"
    reasons <- c(reasons,
      paste0("All ", nc_total, " NC outcomes pass; plasmode coverage = ",
             round(pm_coverage, 3), ", bias = ", round(pm_bias, 4)))
  }

  list(
    decision    = decision,
    rationale   = paste(reasons, collapse = "; "),
    nc_failures = nc_failures,
    nc_total    = nc_total,
    pm_coverage = pm_coverage,
    pm_bias     = pm_bias
  )
}
