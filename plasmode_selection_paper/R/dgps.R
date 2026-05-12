# Data-generating processes for the plasmode-selection paper.
#
# Five DGPs sharing a common baseline-covariate distribution but
# varying how a linear / parametric outcome model approximates the
# truth. Each DGP exposes the same (treatment, baseline-covariate)
# structure so that the propensity-score side is held fixed and only
# the outcome side challenges the candidate libraries.
#
# All DGPs target a marginal RD of -0.05 unless overridden.

# Number of always-present covariates (the analyst's adjustment set).
.K_BASE <- 6L

#' Generate a Single Synthetic Dataset
#'
#' @param n integer; sample size.
#' @param dgp one of "linear", "nonlinear_smooth", "interactions",
#'   "sparse", "high_dim_noise".
#' @param true_RD numeric; the constant marginal risk difference.
#' @param noise_p integer; extra noise covariates for "high_dim_noise".
#' @param seed integer or NULL.
#'
#' @return data.frame with columns treatment, event_24, age, sex,
#'   biomarker, comorbidity, bmi_z, smoke, plus noise covariates
#'   (n01, n02, ...) when applicable.
make_data <- function(n, dgp = c("linear", "nonlinear_smooth",
                                  "interactions", "sparse",
                                  "high_dim_noise"),
                      true_RD = -0.05,
                      noise_p = 20L,
                      seed    = NULL) {
  dgp <- match.arg(dgp)
  if (!is.null(seed)) set.seed(seed)

  age          <- stats::rnorm(n, 60, 10)
  sex          <- stats::rbinom(n, 1, 0.5)
  biomarker    <- stats::rnorm(n, 0, 1)
  comorbidity  <- stats::rbinom(n, 1, 0.3)
  bmi_z        <- stats::rnorm(n, 0, 1)
  smoke        <- stats::rbinom(n, 1, 0.25)

  # Propensity-score model (same across DGPs)
  lp_A <- -0.5 +
           0.02 * (age - 60) +
           0.50 * biomarker +
           0.30 * comorbidity +
           0.20 * bmi_z +
           0.40 * smoke +
           0.10 * sex
  p_A  <- stats::plogis(lp_A)
  A    <- stats::rbinom(n, 1, p_A)

  # Outcome model varies by DGP
  lp_Y0 <- switch(dgp,
    linear = {
      -2.0 +
        0.03 * (age - 60) +
        0.20 * sex +
        0.40 * biomarker +
        0.50 * comorbidity +
        0.10 * bmi_z +
        0.40 * smoke
    },
    nonlinear_smooth = {
      # Smooth nonlinear effects in continuous covariates
      -2.0 +
        0.03 * (age - 60) +
        0.20 * sex +
        0.40 * biomarker + 0.40 * biomarker^2 +
        0.50 * comorbidity +
        0.10 * bmi_z + 0.20 * sin(bmi_z) +
        0.40 * smoke
    },
    interactions = {
      # Two strong two-way interactions
      -2.0 +
        0.03 * (age - 60) +
        0.20 * sex +
        0.40 * biomarker +
        0.50 * comorbidity +
        0.10 * bmi_z +
        0.40 * smoke +
        0.40 * biomarker * comorbidity +
        0.30 * bmi_z * smoke
    },
    sparse = {
      # Only biomarker matters; analyst includes all 6 covariates.
      # Intercept and slope tuned so the p_Y1 = p_Y0 + true_RD clamp
      # at 0.01 fires for < 5 percent of observations, keeping the
      # true marginal causal RD within 0.001 of the target.
      -1.2 + 0.6 * biomarker
    },
    high_dim_noise = {
      # Linear in the 6 real covariates; analyst also adjusts for noise_p extra noise covariates
      -2.0 +
        0.03 * (age - 60) +
        0.20 * sex +
        0.40 * biomarker +
        0.50 * comorbidity +
        0.10 * bmi_z +
        0.40 * smoke
    }
  )

  p_Y0 <- stats::plogis(lp_Y0)
  p_Y1 <- pmin(pmax(p_Y0 + true_RD, 0.01), 0.99)
  Y    <- ifelse(A == 1L,
                 stats::rbinom(n, 1, p_Y1),
                 stats::rbinom(n, 1, p_Y0))

  out <- data.frame(
    treatment   = A,
    event_24    = Y,
    age         = age,
    sex         = sex,
    biomarker   = biomarker,
    comorbidity = comorbidity,
    bmi_z       = bmi_z,
    smoke       = smoke,
    stringsAsFactors = FALSE
  )

  if (dgp == "high_dim_noise" && noise_p > 0L) {
    noise <- matrix(stats::rnorm(n * noise_p), n, noise_p)
    colnames(noise) <- sprintf("n%02d", seq_len(noise_p))
    out <- cbind(out, as.data.frame(noise))
  }
  attr(out, "dgp")          <- dgp
  attr(out, "true_RD")      <- true_RD
  attr(out, "baseline_p")   <- .K_BASE
  attr(out, "n_noise")      <- if (dgp == "high_dim_noise") noise_p else 0L
  out
}

#' Recover the True Marginal RD by Monte Carlo Integration
#'
#' Useful only if the DGP changes and a fresh anchor is needed.
#' Not called by the main simulation; kept for diagnostics.
mc_true_RD <- function(dgp, n = 1e6L, true_RD = -0.05, seed = 99L) {
  d  <- make_data(n, dgp = dgp, true_RD = true_RD, seed = seed)
  mean(d$event_24[d$treatment == 1L]) - mean(d$event_24[d$treatment == 0L])
}
