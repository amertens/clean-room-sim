# ============================================================================
# Scenario data-generating process for the clean-room TMLE simulation.
#
# SINGLE SOURCE OF TRUTH for the simulation DGP. run_simulation.R and the
# sandbox verification/calibration scripts all `source()` this file, so the
# realised-data generator and the large-sample truth cannot drift between the
# production run and the scripts that verify it. Base R only; no dependency on
# `config` or any package (functions read all inputs from their arguments).
# ============================================================================

generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL,
                          U_prevalence = 0, U_trt_OR = 1, U_out_OR = 1,
                          misspec = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  age         <- rnorm(n, mean = 55, sd = 10)
  sex         <- rbinom(n, 1, 0.55)
  biomarker   <- rnorm(n, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n, 1, 0.12)

  # Unmeasured confounder (not returned in the data frame the analyst sees).
  U <- if (U_prevalence > 0) rbinom(n, 1, U_prevalence) else rep(0L, n)

  if (misspec) {
    # Misspecified-surface DGP (Scenario D). Both nuisances are nonlinear:
    # a quadratic in (standardised) age and a sex x biomarker interaction, and
    # the treatment effect is modified by sex. A main-effects GLM is therefore
    # wrong for BOTH the PS and the outcome, so double robustness cannot rescue
    # it, while a flexible learner (gam + interactions) can. The PS nonlinearity
    # is the symmetric interaction, which does not pile the PS at 0/1.
    a <- (age - 55) / 10
    lp_trt <- -0.1 + 0.5 * a + 1.0 * sex * biomarker + 0.5 * ckd +
      0.2 * comorbidity + log(U_trt_OR) * U
    ps_true <- plogis(lp_trt)
    treatment <- rbinom(n, 1, ps_true)
    lp_out <- -0.6 + 0.4 * a + 0.7 * a^2 + 1.0 * sex * biomarker +
      0.6 * ckd + 0.3 * comorbidity +
      effect_size / 0.15 * treatment * (1 + 0.4 * sex) +
      log(U_out_OR) * U
  } else {
    lp_trt <- -0.5 +
      overlap_strength * (0.03 * (age - 55) +
                          0.8 * sex +
                          0.6 * biomarker +
                          0.5 * ckd +
                          0.3 * comorbidity) +
      log(U_trt_OR) * U
    ps_true <- plogis(lp_trt)
    treatment <- rbinom(n, 1, ps_true)
    lp_out <- -2.5 +
      0.015 * (age - 55) +
      0.3 * sex +
      0.2 * biomarker +
      0.6 * ckd +
      0.25 * comorbidity +
      effect_size / 0.15 * treatment +
      log(U_out_OR) * U
  }

  event_24 <- rbinom(n, 1, plogis(lp_out))

  # Negative-control outcome: depends on covariates only, not treatment.
  lp_nc <- -1.0 + 0.01 * (age - 55) + 0.1 * sex + 0.15 * biomarker
  nc_outcome <- rbinom(n, 1, plogis(lp_nc))

  data.frame(
    age         = round(age, 1),
    sex         = sex,
    biomarker   = round(biomarker, 3),
    comorbidity = comorbidity,
    ckd         = ckd,
    treatment   = treatment,
    event_24    = event_24,
    nc_outcome  = nc_outcome,
    stringsAsFactors = FALSE
  )
}


compute_truth <- function(n_truth, overlap_strength, effect_size, seed,
                          U_prevalence = 0, U_out_OR = 1, misspec = FALSE) {
  set.seed(seed)
  age         <- rnorm(n_truth, mean = 55, sd = 10)
  sex         <- rbinom(n_truth, 1, 0.55)
  biomarker   <- rnorm(n_truth, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n_truth, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n_truth, 1, 0.12)
  U           <- if (U_prevalence > 0) rbinom(n_truth, 1, U_prevalence) else rep(0L, n_truth)

  coef_trt <- effect_size / 0.15
  if (misspec) {
    a <- (age - 55) / 10
    lp_out_base <- -0.6 + 0.4 * a + 0.7 * a^2 + 1.0 * sex * biomarker +
      0.6 * ckd + 0.3 * comorbidity + log(U_out_OR) * U
    # Treatment effect is modified by sex, so the marginal contrast integrates
    # the arm-specific coefficient over the covariate distribution.
    risk_1 <- mean(plogis(lp_out_base + coef_trt * (1 + 0.4 * sex)))
    risk_0 <- mean(plogis(lp_out_base))
    return(list(risk_1 = risk_1, risk_0 = risk_0, RD = risk_1 - risk_0))
  }

  lp_out_base <- -2.5 +
    0.015 * (age - 55) +
    0.3 * sex +
    0.2 * biomarker +
    0.6 * ckd +
    0.25 * comorbidity +
    log(U_out_OR) * U

  risk_1 <- mean(plogis(lp_out_base + coef_trt))
  risk_0 <- mean(plogis(lp_out_base))

  list(risk_1 = risk_1, risk_0 = risk_0, RD = risk_1 - risk_0)
}
