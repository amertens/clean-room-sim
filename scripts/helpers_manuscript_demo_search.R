# ===========================================================================
# helpers_manuscript_demo_search.R
# ===========================================================================
# Helper functions for the manuscript demo search script.
# Provides: data generation, truth computation, PS diagnostics,
#           PS matching, TMLE estimation, and candidate scoring.
# ===========================================================================

suppressPackageStartupMessages({
  library(glmnet)
})

# ---- Utility ----

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1 - p))

# ---- Data generation ----

#' Generate one dataset from the HCV-style DGP.
#'
#' @param n         sample size
#' @param seed      random seed
#' @param alpha_A   treatment intercept (set to target ~0.48-0.50 prevalence)
#' @param b_rr      coefficient of renal_risk in treatment model
#' @param alpha_Y   outcome intercept
#' @param beta_A    treatment main effect (fixed = 0.22)
#' @param beta_age1,beta_age2,beta_age3  age polynomial coefficients
#' @param beta_ckd,beta_cirr,beta_dm,beta_nsaid  covariate outcome coefficients
#' @param beta_rr1,beta_rr2  renal_risk outcome coefficients
#' @param beta_ckdA,beta_ageA,beta_age2A,beta_rrA  treatment-by-covariate interactions
#' @return data.frame with columns: age, age_std, ckd, cirrhosis, diabetes,
#'   nsaid, renal_risk, A, Y, ps_true
generate_data <- function(n, seed,
                          alpha_A, b_rr,
                          alpha_Y,
                          beta_A    = 0.22,
                          beta_age1 = 0.25, beta_age2 = -0.20,
                          beta_age3 = 0.10,
                          beta_ckd  = 0.95, beta_cirr = 0.55,
                          beta_dm   = 0.40, beta_nsaid = 0.25,
                          beta_rr1  = 0.03, beta_rr2  = 0.05,
                          beta_ckdA = 0.30, beta_ageA = 0.18,
                          beta_age2A = 0.10, beta_rrA  = 0.00) {
  set.seed(seed)

  age        <- rnorm(n, 55, 10)
  age_std    <- (age - 55) / 10
  ckd        <- rbinom(n, 1, 0.18)
  cirrhosis  <- rbinom(n, 1, 0.22)
  diabetes   <- rbinom(n, 1, 0.20)
  nsaid      <- rbinom(n, 1, 0.25)
  renal_risk <- rnorm(n, 0, 1)

  # Treatment model
  lp_A <- alpha_A +
    0.25 * age_std + 0.45 * ckd + 0.35 * cirrhosis +
    0.20 * diabetes - 0.15 * nsaid + b_rr * renal_risk
  ps_true <- expit(lp_A)
  ps_true <- pmin(pmax(ps_true, 0.005), 0.995)
  A <- rbinom(n, 1, ps_true)

  # Outcome model
  lp_Y <- alpha_Y +
    beta_A * A +
    beta_age1 * age_std + beta_age2 * age_std^2 + beta_age3 * age_std^3 +
    beta_ckd * ckd + beta_cirr * cirrhosis + beta_dm * diabetes +
    beta_nsaid * nsaid +
    beta_rr1 * renal_risk + beta_rr2 * renal_risk^2 +
    beta_ckdA * A * ckd + beta_ageA * A * age_std +
    beta_age2A * A * age_std^2 + beta_rrA * A * renal_risk
  Y <- rbinom(n, 1, expit(lp_Y))

  data.frame(age = age, age_std = age_std, ckd = ckd, cirrhosis = cirrhosis,
             diabetes = diabetes, nsaid = nsaid, renal_risk = renal_risk,
             A = A, Y = Y, ps_true = ps_true)
}

# ---- Truth computation ----

#' Compute the true ATE (risk difference) by Monte Carlo.
compute_truth <- function(n_mc, seed = 999,
                          alpha_Y,
                          beta_A    = 0.22,
                          beta_age1, beta_age2, beta_age3,
                          beta_ckd  = 0.95, beta_cirr = 0.55,
                          beta_dm   = 0.40, beta_nsaid = 0.25,
                          beta_rr1, beta_rr2,
                          beta_ckdA, beta_ageA, beta_age2A, beta_rrA) {
  set.seed(seed)
  age_std    <- rnorm(n_mc, 0, 1)
  ckd        <- rbinom(n_mc, 1, 0.18)
  cirrhosis  <- rbinom(n_mc, 1, 0.22)
  diabetes   <- rbinom(n_mc, 1, 0.20)
  nsaid      <- rbinom(n_mc, 1, 0.25)
  renal_risk <- rnorm(n_mc, 0, 1)

  lp_base <- alpha_Y +
    beta_age1 * age_std + beta_age2 * age_std^2 + beta_age3 * age_std^3 +
    beta_ckd * ckd + beta_cirr * cirrhosis + beta_dm * diabetes +
    beta_nsaid * nsaid +
    beta_rr1 * renal_risk + beta_rr2 * renal_risk^2

  EY1 <- mean(expit(lp_base + beta_A +
                       beta_ckdA * ckd + beta_ageA * age_std +
                       beta_age2A * age_std^2 + beta_rrA * renal_risk))
  EY0 <- mean(expit(lp_base))
  c(EY1 = EY1, EY0 = EY0, RD = EY1 - EY0)
}

# ---- Calibrate alpha_A to target treatment prevalence ----

#' Find alpha_A that yields approximately target_prev treatment prevalence.
calibrate_alpha_A <- function(b_rr, target_prev = 0.49, n_cal = 50000,
                              seed = 42) {
  set.seed(seed)
  age_std    <- rnorm(n_cal, 0, 1)
  ckd        <- rbinom(n_cal, 1, 0.18)
  cirrhosis  <- rbinom(n_cal, 1, 0.22)
  diabetes   <- rbinom(n_cal, 1, 0.20)
  nsaid      <- rbinom(n_cal, 1, 0.25)
  renal_risk <- rnorm(n_cal, 0, 1)

  lp_no_intercept <- 0.25 * age_std + 0.45 * ckd + 0.35 * cirrhosis +
    0.20 * diabetes - 0.15 * nsaid + b_rr * renal_risk

  f <- function(a) {
    mean(expit(a + lp_no_intercept)) - target_prev
  }
  uniroot(f, interval = c(-10, 10))$root
}

# ---- Calibrate alpha_Y to target outcome prevalence ----

calibrate_alpha_Y <- function(target_prev,
                              beta_A    = 0.22,
                              beta_age1, beta_age2, beta_age3,
                              beta_ckd  = 0.95, beta_cirr = 0.55,
                              beta_dm   = 0.40, beta_nsaid = 0.25,
                              beta_rr1, beta_rr2,
                              beta_ckdA, beta_ageA, beta_age2A, beta_rrA,
                              n_cal = 50000, seed = 43) {
  set.seed(seed)
  age_std    <- rnorm(n_cal, 0, 1)
  ckd        <- rbinom(n_cal, 1, 0.18)
  cirrhosis  <- rbinom(n_cal, 1, 0.22)
  diabetes   <- rbinom(n_cal, 1, 0.20)
  nsaid      <- rbinom(n_cal, 1, 0.25)
  renal_risk <- rnorm(n_cal, 0, 1)
  # Use roughly 49% treatment to calibrate marginal outcome rate
  A          <- rbinom(n_cal, 1, 0.49)

  lp_no_int <- beta_A * A +
    beta_age1 * age_std + beta_age2 * age_std^2 + beta_age3 * age_std^3 +
    beta_ckd * ckd + beta_cirr * cirrhosis + beta_dm * diabetes +
    beta_nsaid * nsaid +
    beta_rr1 * renal_risk + beta_rr2 * renal_risk^2 +
    beta_ckdA * A * ckd + beta_ageA * A * age_std +
    beta_age2A * A * age_std^2 + beta_rrA * A * renal_risk

  f <- function(a) mean(expit(a + lp_no_int)) - target_prev
  uniroot(f, interval = c(-10, 10))$root
}

# ---- Traditional PS workflow ----

#' Fit PS via logistic regression (main effects only).
fit_conventional_ps <- function(dat, ps_clip = c(0.025, 0.975)) {
  g_fit <- glm(A ~ age_std + ckd + cirrhosis + diabetes + nsaid + renal_risk,
               data = dat, family = "binomial")
  ps <- fitted(g_fit)
  ps_trunc <- pmin(pmax(ps, ps_clip[1]), ps_clip[2])
  list(ps = ps, ps_trunc = ps_trunc, g_fit = g_fit)
}

#' 1:1 nearest-neighbor PS matching without replacement.
run_ps_match <- function(dat, ps, ps_clip = c(0.025, 0.975)) {
  lps <- logit(pmin(pmax(ps, ps_clip[1]), ps_clip[2]))
  caliper <- 0.2 * sd(lps)

  tidx <- which(dat$A == 1)
  cidx <- which(dat$A == 0)

  matched_t <- matched_c <- integer(0)
  available <- rep(TRUE, length(cidx))

  for (ii in seq_along(tidx)) {
    i <- tidx[ii]
    if (!any(available)) break
    ac <- cidx[available]
    d_lps <- abs(lps[i] - lps[ac])
    j <- which.min(d_lps)
    if (d_lps[j] <= caliper) {
      matched_t <- c(matched_t, i)
      matched_c <- c(matched_c, ac[j])
      available[which(cidx == ac[j])] <- FALSE
    }
  }

  np <- length(matched_t)
  if (np < 20) {
    return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA,
                n_pairs = np, matched_frac = np / length(tidx),
                max_smd = NA))
  }

  rd <- mean(dat$Y[matched_t] - dat$Y[matched_c])
  se <- sd(dat$Y[matched_t] - dat$Y[matched_c]) / sqrt(np)

  covars <- c("age_std", "ckd", "cirrhosis", "diabetes", "nsaid", "renal_risk")
  smds <- sapply(covars, function(v) {
    xt <- dat[[v]][matched_t]; xc <- dat[[v]][matched_c]
    s <- sqrt((var(xt) + var(xc)) / 2)
    if (s < 1e-10) 0 else (mean(xt) - mean(xc)) / s
  })

  list(estimate = rd, se = se,
       ci_low = rd - 1.96 * se, ci_high = rd + 1.96 * se,
       n_pairs = np, matched_frac = np / length(tidx),
       max_smd = max(abs(smds)))
}

#' Compute Stage-2 diagnostics for the traditional workflow.
compute_stage2_diagnostics <- function(dat, ps, ps_trunc) {
  A <- dat$A
  frac_extreme <- mean(ps < 0.05 | ps > 0.95)
  w <- A / ps_trunc + (1 - A) / (1 - ps_trunc)
  ess <- sum(w)^2 / sum(w^2)
  ess_frac <- ess / nrow(dat)
  max_weight <- max(w)

  m <- run_ps_match(dat, ps)

  list(frac_extreme = frac_extreme,
       ess_frac = ess_frac,
       max_weight = max_weight,
       matched_frac = ifelse(is.null(m$matched_frac), NA, m$matched_frac),
       max_smd_match = ifelse(is.null(m$max_smd), NA, m$max_smd),
       match_result = m)
}

#' Determine if the traditional workflow would fail / stop.
traditional_workflow_fail <- function(diag) {
  any(c(
    diag$frac_extreme > 0.15,
    diag$ess_frac < 0.50,
    isTRUE(diag$matched_frac < 0.50),
    isTRUE(diag$max_smd_match > 0.10)
  ), na.rm = FALSE)
}

# ---- TMLE workflow ----

#' Build a rich design matrix for Q estimation.
build_Q_design <- function(dat) {
  A <- dat$A; age <- dat$age_std; ckd <- dat$ckd; cirr <- dat$cirrhosis
  dm <- dat$diabetes; nsaid <- dat$nsaid; rr <- dat$renal_risk
  cbind(
    A         = A,
    age       = age,
    age2      = age^2,
    age3      = age^3,
    ckd       = ckd,
    cirrhosis = cirr,
    diabetes  = dm,
    nsaid     = nsaid,
    rr        = rr,
    rr2       = rr^2,
    A_ckd     = A * ckd,
    A_age     = A * age,
    A_age2    = A * age^2,
    A_rr      = A * rr,
    A_cirr    = A * cirr,
    A_dm      = A * dm,
    age_ckd   = age * ckd,
    age_rr    = age * rr,
    ckd_rr    = ckd * rr
  )
}

#' Fit g via cv.glmnet.
fit_tmle_g <- function(dat, ps_clip = c(0.025, 0.975)) {
  X <- as.matrix(dat[, c("age_std", "ckd", "cirrhosis", "diabetes",
                          "nsaid", "renal_risk")])
  g_fit <- cv.glmnet(X, dat$A, family = "binomial", alpha = 0.5,
                     nfolds = 5, type.measure = "deviance")
  g_hat <- as.numeric(predict(g_fit, X, s = "lambda.min", type = "response"))
  g_hat <- pmin(pmax(g_hat, ps_clip[1]), ps_clip[2])
  list(g_hat = g_hat, g_fit = g_fit)
}

#' Fit Q via cv.glmnet with rich design matrix.
fit_tmle_Q <- function(dat, ps_clip = c(0.025, 0.975)) {
  X_Q <- build_Q_design(dat)
  Q_fit <- cv.glmnet(X_Q, dat$Y, family = "binomial", alpha = 0.5,
                     nfolds = 5, type.measure = "deviance")
  Q_A <- as.numeric(predict(Q_fit, X_Q, s = "lambda.min", type = "response"))

  dat1 <- dat; dat1$A <- 1
  dat0 <- dat; dat0$A <- 0
  Q_1 <- as.numeric(predict(Q_fit, build_Q_design(dat1),
                            s = "lambda.min", type = "response"))
  Q_0 <- as.numeric(predict(Q_fit, build_Q_design(dat0),
                            s = "lambda.min", type = "response"))

  Q_A <- pmin(pmax(Q_A, 1e-5), 1 - 1e-5)
  Q_1 <- pmin(pmax(Q_1, 1e-5), 1 - 1e-5)
  Q_0 <- pmin(pmax(Q_0, 1e-5), 1 - 1e-5)

  list(Q_A = Q_A, Q_1 = Q_1, Q_0 = Q_0, Q_fit = Q_fit)
}

#' Run one TMLE analysis (targeting step + inference).
run_tmle <- function(dat, ps_clip = c(0.025, 0.975)) {
  n <- nrow(dat)

  g <- fit_tmle_g(dat, ps_clip)
  g_hat <- g$g_hat

  q <- fit_tmle_Q(dat, ps_clip)

  # Targeting step
  H  <- dat$A / g_hat - (1 - dat$A) / (1 - g_hat)
  H1 <- 1 / g_hat
  H0 <- -1 / (1 - g_hat)

  eps_fit <- glm(dat$Y ~ -1 + H, offset = logit(q$Q_A), family = "binomial")
  eps <- coef(eps_fit)

  Q_1_star <- expit(logit(q$Q_1) + eps * H1)
  Q_0_star <- expit(logit(q$Q_0) + eps * H0)
  Q_A_star <- expit(logit(q$Q_A) + eps * H)

  psi <- mean(Q_1_star - Q_0_star)

  # EIF-based inference
  eif <- H * (dat$Y - Q_A_star) + (Q_1_star - Q_0_star) - psi
  se  <- sqrt(var(eif) / n)
  ci  <- psi + c(-1.96, 1.96) * se

  list(estimate = psi, se = se, ci_low = ci[1], ci_high = ci[2])
}

# ---- SuperLearner TMLE (optional) ----

#' Run TMLE using SuperLearner for g and Q.
#' Only called if SuperLearner package is available.
run_tmle_sl <- function(dat, ps_clip = c(0.025, 0.975)) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    return(list(estimate = NA, se = NA, ci_low = NA, ci_high = NA))
  }

  n <- nrow(dat)
  X_g <- as.data.frame(dat[, c("age_std", "ckd", "cirrhosis", "diabetes",
                                "nsaid", "renal_risk")])
  X_Q <- as.data.frame(build_Q_design(dat))

  # g via SuperLearner
  sl_g <- SuperLearner::SuperLearner(
    Y = dat$A, X = X_g, family = binomial(),
    SL.library = c("SL.mean", "SL.glm", "SL.glmnet"),
    cvControl = list(V = 5)
  )
  g_hat <- as.numeric(sl_g$SL.predict)
  g_hat <- pmin(pmax(g_hat, ps_clip[1]), ps_clip[2])

  # Q via SuperLearner
  sl_libs <- c("SL.mean", "SL.glm", "SL.glmnet")
  if (requireNamespace("gam", quietly = TRUE)) {
    sl_libs <- c(sl_libs, "SL.gam")
  }
  sl_Q <- SuperLearner::SuperLearner(
    Y = dat$Y, X = X_Q, family = binomial(),
    SL.library = sl_libs,
    cvControl = list(V = 5)
  )

  Q_A <- as.numeric(sl_Q$SL.predict)
  dat1 <- dat; dat1$A <- 1
  dat0 <- dat; dat0$A <- 0
  Q_1 <- as.numeric(predict(sl_Q, newdata = as.data.frame(build_Q_design(dat1)))$pred)
  Q_0 <- as.numeric(predict(sl_Q, newdata = as.data.frame(build_Q_design(dat0)))$pred)

  Q_A <- pmin(pmax(Q_A, 1e-5), 1 - 1e-5)
  Q_1 <- pmin(pmax(Q_1, 1e-5), 1 - 1e-5)
  Q_0 <- pmin(pmax(Q_0, 1e-5), 1 - 1e-5)

  # Targeting step
  H  <- dat$A / g_hat - (1 - dat$A) / (1 - g_hat)
  H1 <- 1 / g_hat
  H0 <- -1 / (1 - g_hat)

  eps_fit <- glm(dat$Y ~ -1 + H, offset = logit(Q_A), family = "binomial")
  eps <- coef(eps_fit)

  Q_1_star <- expit(logit(Q_1) + eps * H1)
  Q_0_star <- expit(logit(Q_0) + eps * H0)
  Q_A_star <- expit(logit(Q_A) + eps * H)

  psi <- mean(Q_1_star - Q_0_star)

  eif <- H * (dat$Y - Q_A_star) + (Q_1_star - Q_0_star) - psi
  se  <- sqrt(var(eif) / n)
  ci  <- psi + c(-1.96, 1.96) * se

  list(estimate = psi, se = se, ci_low = ci[1], ci_high = ci[2])
}

# ---- Candidate evaluation ----

#' Run one replicate: both workflows on the same dataset.
run_one_replicate <- function(dat, true_rd, ps_clip = c(0.025, 0.975)) {
  # -- Traditional PS workflow --
  ps_fit <- fit_conventional_ps(dat, ps_clip)
  diag   <- compute_stage2_diagnostics(dat, ps_fit$ps, ps_fit$ps_trunc)
  m      <- diag$match_result  # already computed inside diagnostics

  psm_est <- m$estimate
  psm_se  <- m$se
  psm_ci_low  <- m$ci_low
  psm_ci_high <- m$ci_high

  # -- TMLE workflow --
  tmle_res <- tryCatch(
    run_tmle(dat, ps_clip),
    error = function(e) list(estimate = NA, se = NA, ci_low = NA, ci_high = NA)
  )

  list(
    # Stage-2 diagnostics
    frac_extreme = diag$frac_extreme,
    ess_frac     = diag$ess_frac,
    max_weight   = diag$max_weight,
    matched_frac = diag$matched_frac,
    max_smd_match = diag$max_smd_match,
    trad_fail    = traditional_workflow_fail(diag),
    # PS-matched results
    psm_est      = psm_est,
    psm_se       = psm_se,
    psm_ci_low   = psm_ci_low,
    psm_ci_high  = psm_ci_high,
    # TMLE results
    tmle_est     = tmle_res$estimate,
    tmle_se      = tmle_res$se,
    tmle_ci_low  = tmle_res$ci_low,
    tmle_ci_high = tmle_res$ci_high
  )
}

#' Summarize replicate-level results for one candidate configuration.
summarize_candidate <- function(rep_results, true_rd) {
  rr <- do.call(rbind, lapply(rep_results, function(x) as.data.frame(x)))

  # Stage-2 diagnostics (median across reps)
  med_frac_extreme <- median(rr$frac_extreme, na.rm = TRUE)
  med_ess_frac     <- median(rr$ess_frac, na.rm = TRUE)
  med_matched_frac <- median(rr$matched_frac, na.rm = TRUE)
  med_max_smd      <- median(rr$max_smd_match, na.rm = TRUE)
  frac_trad_fail   <- mean(rr$trad_fail, na.rm = TRUE)

  # PS-matched performance
  psm_ok <- !is.na(rr$psm_est)
  if (sum(psm_ok) >= 5) {
    psm_bias     <- mean(rr$psm_est[psm_ok]) - true_rd
    psm_rmse     <- sqrt(mean((rr$psm_est[psm_ok] - true_rd)^2))
    psm_emp_sd   <- sd(rr$psm_est[psm_ok])
    psm_mean_se  <- mean(rr$psm_se[psm_ok], na.rm = TRUE)
    psm_has_ci   <- psm_ok & !is.na(rr$psm_ci_low)
    psm_coverage <- if (sum(psm_has_ci) > 0) {
      mean((rr$psm_ci_low[psm_has_ci] <= true_rd) &
             (true_rd <= rr$psm_ci_high[psm_has_ci]))
    } else NA
    psm_n_ok <- sum(psm_ok)
  } else {
    psm_bias <- psm_rmse <- psm_emp_sd <- psm_mean_se <- psm_coverage <- NA
    psm_n_ok <- sum(psm_ok)
  }

  # TMLE performance
  tmle_ok <- !is.na(rr$tmle_est)
  if (sum(tmle_ok) >= 5) {
    tmle_bias     <- mean(rr$tmle_est[tmle_ok]) - true_rd
    tmle_rmse     <- sqrt(mean((rr$tmle_est[tmle_ok] - true_rd)^2))
    tmle_emp_sd   <- sd(rr$tmle_est[tmle_ok])
    tmle_mean_se  <- mean(rr$tmle_se[tmle_ok], na.rm = TRUE)
    tmle_has_ci   <- tmle_ok & !is.na(rr$tmle_ci_low)
    tmle_coverage <- if (sum(tmle_has_ci) > 0) {
      mean((rr$tmle_ci_low[tmle_has_ci] <= true_rd) &
             (true_rd <= rr$tmle_ci_high[tmle_has_ci]))
    } else NA
    tmle_se_sd_ratio <- if (!is.na(tmle_mean_se) && tmle_emp_sd > 0) {
      tmle_mean_se / tmle_emp_sd
    } else NA
    tmle_n_ok <- sum(tmle_ok)
  } else {
    tmle_bias <- tmle_rmse <- tmle_emp_sd <- tmle_mean_se <- NA
    tmle_coverage <- tmle_se_sd_ratio <- NA
    tmle_n_ok <- sum(tmle_ok)
  }

  # Flags
  trad_fail_flag <- (frac_trad_fail >= 0.5)  # majority of reps fail

  tmle_pass <- !is.na(tmle_bias) && !is.na(tmle_coverage) &&
    !is.na(tmle_se_sd_ratio) &&
    abs(tmle_bias) < 0.01 &&
    tmle_coverage >= 0.93 &&
    tmle_se_sd_ratio > 0.80 && tmle_se_sd_ratio < 1.20

  overlap_stressed <- med_frac_extreme > 0.10

  tmle_beats_psm <- FALSE
  if (!is.na(tmle_rmse) && !is.na(psm_rmse)) {
    rmse_better <- tmle_rmse < psm_rmse * 0.95
    cov_better  <- !is.na(tmle_coverage) && !is.na(psm_coverage) &&
      tmle_coverage > psm_coverage + 0.02
    bias_better <- !is.na(tmle_bias) && !is.na(psm_bias) &&
      abs(tmle_bias) < abs(psm_bias) * 0.80
    tmle_beats_psm <- rmse_better || cov_better || bias_better
  }

  # Scoring
  score <- 0
  if (trad_fail_flag)   score <- score + 3
  if (tmle_pass)        score <- score + 3
  if (tmle_beats_psm)   score <- score + 2
  if (overlap_stressed) score <- score + 1
  if (!is.na(tmle_coverage) && tmle_coverage >= 0.93) score <- score + 1
  if (!is.na(tmle_rmse) && !is.na(psm_rmse) && tmle_rmse < psm_rmse) {
    score <- score + 1
  }

  data.frame(
    # Stage-2 diagnostics (medians)
    med_frac_extreme = round(med_frac_extreme, 4),
    med_ess_frac     = round(med_ess_frac, 4),
    med_matched_frac = round(med_matched_frac, 4),
    med_max_smd      = round(med_max_smd, 4),
    frac_trad_fail   = round(frac_trad_fail, 4),
    # PS-matched performance
    psm_n_ok   = psm_n_ok,
    psm_bias   = round2(psm_bias),
    psm_rmse   = round2(psm_rmse),
    psm_emp_sd = round2(psm_emp_sd),
    psm_mean_se = round2(psm_mean_se),
    psm_coverage = round2(psm_coverage),
    # TMLE performance
    tmle_n_ok   = tmle_n_ok,
    tmle_bias   = round2(tmle_bias),
    tmle_rmse   = round2(tmle_rmse),
    tmle_emp_sd = round2(tmle_emp_sd),
    tmle_mean_se = round2(tmle_mean_se),
    tmle_coverage = round2(tmle_coverage),
    tmle_se_sd_ratio = round2(tmle_se_sd_ratio),
    # Flags & score
    trad_fail_flag = trad_fail_flag,
    tmle_pass      = tmle_pass,
    overlap_stressed = overlap_stressed,
    tmle_beats_psm = tmle_beats_psm,
    score          = score,
    true_rd        = round(true_rd, 6),
    stringsAsFactors = FALSE
  )
}

# Small utility to handle rounding NAs
round2 <- function(x, digits = 5) {
  if (is.na(x)) NA else round(x, digits)
}
