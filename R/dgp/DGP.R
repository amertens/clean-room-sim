#' @title Data-Generating Process for HCV-AKI Simulation
#' @description Generates realistic simulated data mimicking a claims-based
#'   cohort study of sofosbuvir-containing (SOF) vs non-SOF DAA regimens
#'   and acute kidney injury (AKI) risk.
#' @name dgp
NULL

#' Generate HCV-AKI Simulation Data
#'
#' Creates a simulated cohort with configurable complexity, including
#' nonlinear confounding, informative censoring, treatment switching,
#' and non-proportional hazards.
#'
#' @param N Integer sample size (pre-exclusion).
#' @param p_sof Numeric target SOF treatment probability.
#' @param h0 Numeric baseline hazard rate.
#' @param HR_early Numeric hazard ratio for treatment in early period.
#' @param HR_late Numeric hazard ratio for treatment in late period (non-PH).
#' @param tau Numeric change-point day for non-proportional hazards.
#' @param max_follow Numeric administrative censoring cutoff in days.
#' @param risk_window Numeric at-risk window after treatment switch in days.
#' @param np_hazard Logical; if TRUE, use piecewise hazard (non-PH).
#' @param dep_censor Logical; if TRUE, censoring depends on risk factors.
#' @param complexity Logical; if TRUE, use nonlinear PS and outcome models.
#' @param switch_on Logical; if TRUE, use hazard-based treatment switching.
#' @param lambda_sw0 Numeric baseline switch hazard.
#' @param gamma_A Numeric log-HR for treatment effect on switching.
#' @param gamma_ckd Numeric log-HR for CKD effect on switching.
#' @param censor_base Numeric base administrative censoring rate.
#' @param treat_override Character; one of "simulate", "all_treated",
#'   "all_control" for counterfactual computation.
#' @param add_missing Logical; if TRUE, introduce missingness.
#' @param impute Logical; if TRUE and add_missing is TRUE, impute via
#'   missForest.
#' @param seed Integer seed for reproducibility; NULL to use current state.
#'
#' @return A data.frame (tibble) with columns: id, age, sex_male, race, region,
#'   ckd, heart_failure, sepsis, dehydration, obstruction, cirrhosis,
#'   portal_htn, esld, hiv, diabetes, hypertension, bmi, overweight_obese,
#'   smoking, alcohol, substance_abuse, cancer, chemo, nsaid, acearb, diuretic,
#'   aminoglycoside, contrast, statin, aspirin, beta_blocker, ccb, art,
#'   treatment, switch, follow_time, event.
#'
#' @export
generate_hcv_data <- function(
    N               = 125000,
    p_sof           = 0.36,
    h0              = 5e-5,
    HR_early        = 1.25,
    HR_late         = 0.70,
    tau             = 45,
    max_follow      = 720,
    risk_window     = 30,
    np_hazard       = TRUE,
    dep_censor      = TRUE,
    complexity      = TRUE,
    switch_on       = TRUE,
    lambda_sw0      = 2.5e-5,
    gamma_A         = 0.60,
    gamma_ckd       = 0.40,
    censor_base     = 1/100,
    treat_override  = c("simulate", "all_treated", "all_control"),
    add_missing     = FALSE,
    impute          = FALSE,
    seed            = NULL)
{
  if (!is.null(seed)) set.seed(seed)
  treat_override <- match.arg(treat_override)
  if (impute) requireNamespace("missForest", quietly = TRUE)

  # --------------------------------------------------------------------------
  # 1. Demography
  # --------------------------------------------------------------------------
  raw <- data.frame(
    id          = seq_len(N),
    age         = pmax(stats::rnorm(N, 48, 13), 18),
    sex_male    = stats::rbinom(N, 1, 0.58),
    race        = sample(c("white", "black", "hispanic", "asian", "other"), N,
                         replace = TRUE,
                         prob = c(.48, .14, .06, .02, .30)),
    region      = sample(c("NE", "MW", "S", "W"), N, replace = TRUE,
                         prob = c(.20, .18, .37, .25)),
    enroll_days = stats::rpois(N, 420),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------------------------
  # 2. Clinical history and concomitant medications
  # --------------------------------------------------------------------------
  add_bin <- function(p) stats::rbinom(N, 1, p)
  raw$ckd             <- add_bin(.08)
  raw$prior_aki       <- add_bin(.05)
  raw$heart_failure   <- add_bin(.07)
  raw$sepsis          <- add_bin(.03)
  raw$dehydration     <- add_bin(.06)
  raw$obstruction     <- add_bin(.04)
  raw$cirrhosis       <- add_bin(.18)
  raw$portal_htn      <- add_bin(.04)
  raw$esld            <- add_bin(.02)
  raw$hiv             <- add_bin(.04)
  raw$diabetes        <- add_bin(.20)
  raw$hypertension    <- add_bin(.45)
  raw$bmi             <- stats::rnorm(N, 28, 5)
  raw$overweight_obese <- add_bin(.20)
  raw$smoking         <- add_bin(.40)
  raw$alcohol         <- add_bin(.18)
  raw$substance_abuse <- add_bin(.25)
  raw$cancer          <- add_bin(.08)
  raw$chemo           <- add_bin(.01)
  raw$nsaid           <- add_bin(.25)
  raw$acearb          <- add_bin(.30)
  raw$diuretic        <- add_bin(.22)
  raw$aminoglycoside  <- add_bin(.05)
  raw$contrast        <- add_bin(.08)
  raw$statin          <- add_bin(.15)
  raw$aspirin         <- add_bin(.10)
  raw$beta_blocker    <- add_bin(.14)
  raw$ccb             <- add_bin(.16)
  raw$art             <- add_bin(.05)
  raw$prior_sof       <- add_bin(.05)
  raw$prior_nonsof    <- add_bin(.05)

  # --------------------------------------------------------------------------
  # 3. Baseline exclusions
  # --------------------------------------------------------------------------
  cohort <- raw[raw$enroll_days >= 365 & raw$age >= 18 &
                  raw$prior_aki == 0 &
                  !(raw$prior_sof == 1 | raw$prior_nonsof == 1), ]
  N_c <- nrow(cohort)

  # --------------------------------------------------------------------------
  # 4. Treatment assignment
  # --------------------------------------------------------------------------
  if (treat_override == "simulate") {
    lp0 <- with(cohort,
                0.015 * age + 0.30 * cirrhosis + 0.25 * ckd +
                  0.15 * hiv + 0.10 * diabetes -
                  0.10 * cancer + stats::rnorm(N_c, 0, 0.6))
    if (complexity) {
      lp0 <- with(cohort, lp0 +
                    0.02 * (bmi^2) / 100 - 0.3 * sin(0.1 * bmi) +
                    0.5 * (age / 50)^3 +
                    1.5 * ckd * cancer + 0.8 * hiv * log1p(age))
    }
    alpha0  <- stats::qlogis(p_sof) - mean(lp0)
    p_trt   <- pmin(pmax(stats::plogis(alpha0 + lp0), 0.05), 0.95)
    cohort$treatment <- stats::rbinom(N_c, 1, p_trt)
  } else {
    cohort$treatment <- ifelse(treat_override == "all_treated", 1L, 0L)
  }

  # --------------------------------------------------------------------------
  # 5. Individual baseline hazard
  # --------------------------------------------------------------------------
  if (!complexity) {
    lp_out <- with(cohort,
                   -2.8 + 0.03 * age + 0.7 * ckd + 0.5 * cirrhosis +
                     0.3 * heart_failure + 0.25 * nsaid + 0.20 * contrast)
  } else {
    lp_out <- with(cohort,
                   -2.8 + 0.03 * age + 0.0005 * age^2 + 0.7 * ckd +
                     0.5 * cirrhosis +
                     0.02 * (bmi^2) / 100 - 0.3 * sin(0.1 * bmi) +
                     0.4 * heart_failure * acearb +
                     0.6 * nsaid * treatment + 0.3 * contrast * log1p(age))
  }
  base_rate <- h0 * exp(lp_out)

  # --------------------------------------------------------------------------
  # 6. Event times (AKI)
  # --------------------------------------------------------------------------
  if (!np_hazard) {
    rate <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    cohort$event_time <- stats::rexp(N_c, rate = rate)
  } else {
    rpexp_piece <- function(n, r1, r2, tau_val) {
      u  <- stats::runif(n)
      p1 <- 1 - exp(-r1 * tau_val)
      t_out <- numeric(n)
      e  <- u <= p1
      t_out[e]  <- -log(1 - u[e]) / r1[e]
      t_out[!e] <- tau_val - log((1 - u[!e]) / (1 - p1[!e])) / r2[!e]
      t_out
    }
    r1 <- base_rate * ifelse(cohort$treatment == 1, HR_early, 1)
    r2 <- base_rate * ifelse(cohort$treatment == 1, HR_late,  1)
    cohort$event_time <- rpexp_piece(N_c, r1, r2, tau)
  }

  # --------------------------------------------------------------------------
  # 7. Administrative censoring
  # --------------------------------------------------------------------------
  if (!dep_censor) {
    censor_admin <- stats::rexp(N_c, rate = censor_base)
  } else {
    c_rate <- censor_base * exp(0.04 * lp_out + 0.03 * cohort$treatment)
    censor_admin <- stats::rexp(N_c, rate = c_rate)
  }
  cohort$censor_admin <- pmin(censor_admin, max_follow)

  # --------------------------------------------------------------------------
  # 8. Treatment-switch censoring
  # --------------------------------------------------------------------------
  if (!switch_on) {
    cohort$tx_days <- ifelse(cohort$treatment == 1,
                             stats::rpois(N_c, 84),
                             stats::rpois(N_c, 70))
    cohort$switch  <- stats::rbinom(N_c, 1, 0.03)
    cohort$censor_switch <- ifelse(cohort$switch == 1,
                                   cohort$tx_days + risk_window,
                                   max_follow)
  } else {
    lambda_sw <- lambda_sw0 *
      exp(gamma_A * cohort$treatment + gamma_ckd * cohort$ckd)
    Sw_lat <- stats::rexp(N_c, rate = lambda_sw)
    cohort$tx_days <- Sw_lat
    cohort$switch  <- as.integer(Sw_lat < max_follow)
    cohort$censor_switch <- pmin(Sw_lat + risk_window, max_follow)
  }

  # --------------------------------------------------------------------------
  # 9. Observed follow-up and event indicator
  # --------------------------------------------------------------------------
  cohort$follow_time <- pmin(cohort$event_time,
                             cohort$censor_admin,
                             cohort$censor_switch)
  cohort$event <- as.integer(cohort$event_time <= cohort$follow_time)

  # --------------------------------------------------------------------------
  # 10. Analysis dataset
  # --------------------------------------------------------------------------
  drop_cols <- c("enroll_days", "prior_aki", "prior_sof", "prior_nonsof",
                 "tx_days", "event_time", "censor_admin", "censor_switch")
  ana <- cohort[, !(names(cohort) %in% drop_cols)]

  # --------------------------------------------------------------------------
  # 11. Optional missingness / imputation
  # --------------------------------------------------------------------------
  if (add_missing) {
    n_ana <- nrow(ana)
    ana$region[sample(n_ana, round(0.05 * n_ana))] <- NA
    ana$ckd[sample(n_ana, round(0.10 * n_ana))]    <- NA
    if (impute) {
      imp_vars <- c("age", "race", "region", "ckd", "cirrhosis", "hiv",
                    "diabetes", "hypertension", "bmi")
      imp_in <- ana[, imp_vars]
      imp_in$race   <- as.factor(imp_in$race)
      imp_in$region <- as.factor(imp_in$region)
      ana[, imp_vars] <- missForest::missForest(
        as.data.frame(imp_in), verbose = FALSE
      )$ximp
    }
  }

  ana
}

#' Compute True Counterfactual Risk
#'
#' Generate large counterfactual datasets under all-treated and all-control
#' to compute the true risk difference at a given time point.
#'
#' @param N_truth Integer sample size for truth computation (large).
#' @param t_eval Numeric time point for risk evaluation.
#' @param seed Integer seed.
#' @param ... Additional arguments passed to \code{generate_hcv_data}.
#'
#' @return A list with components \code{risk_1} (risk under treatment),
#'   \code{risk_0} (risk under control), \code{RD} (risk difference),
#'   \code{RR} (risk ratio).
#' @export
compute_true_risk <- function(N_truth = 500000, t_eval = 180, seed = 99, ...) {
  # Generate all-treated counterfactual
  set.seed(seed)
  d1 <- generate_hcv_data(N = N_truth, treat_override = "all_treated", ...)
  risk_1 <- mean(d1$follow_time <= t_eval & d1$event == 1)

  # Generate all-control counterfactual
  set.seed(seed)
  d0 <- generate_hcv_data(N = N_truth, treat_override = "all_control", ...)
  risk_0 <- mean(d0$follow_time <= t_eval & d0$event == 1)

  list(
    risk_1 = risk_1,
    risk_0 = risk_0,
    RD     = risk_1 - risk_0,
    RR     = if (risk_0 > 0) risk_1 / risk_0 else NA_real_,
    t_eval = t_eval
  )
}
