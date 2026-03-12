#' @title Shared Helper Functions
#' @description Utility functions used across the clean-room pipeline.
#' @name helpers
NULL

#' Truncate Probabilities
#'
#' Clips probabilities to \[lower, upper\] and reports how many were truncated.
#'
#' @param p Numeric vector of probabilities.
#' @param lower Lower bound (default 0.01).
#' @param upper Upper bound (default 0.99).
#' @return A list with components \code{p} (truncated), \code{n_lower} and
#'   \code{n_upper} (counts of values clipped at each bound).
#' @export
truncate_ps <- function(p, lower = 0.01, upper = 0.99) {
  n_lower <- sum(p < lower, na.rm = TRUE)
  n_upper <- sum(p > upper, na.rm = TRUE)
  p_trunc <- pmin(pmax(p, lower), upper)
  list(p = p_trunc, n_lower = n_lower, n_upper = n_upper)
}

#' Standardized Mean Difference
#'
#' Compute the SMD between treated and control for a single covariate.
#'
#' @param x Numeric covariate vector.
#' @param trt Binary treatment indicator.
#' @param weights Optional numeric weights.
#' @return Numeric scalar: the SMD.
#' @export
compute_smd <- function(x, trt, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(x))
  w1 <- weights[trt == 1]
  w0 <- weights[trt == 0]
  x1 <- x[trt == 1]
  x0 <- x[trt == 0]

  wmean <- function(v, w) sum(v * w) / sum(w)
  wvar  <- function(v, w) {
    m <- wmean(v, w)
    sum(w * (v - m)^2) / sum(w)
  }

  m1 <- wmean(x1, w1)
  m0 <- wmean(x0, w0)
  v1 <- wvar(x1, w1)
  v0 <- wvar(x0, w0)

  denom <- sqrt((v1 + v0) / 2)
  if (denom < 1e-12) return(0)
  (m1 - m0) / denom
}

#' Effective Sample Size
#'
#' Compute the effective sample size given weights.
#'
#' @param weights Numeric vector of weights.
#' @return Numeric scalar: effective sample size.
#' @export
effective_ss <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}

#' Capture Session Info to File
#'
#' Writes sessionInfo() to a text file for reproducibility.
#'
#' @param filepath Character path for output file.
#' @return Invisibly returns the filepath.
#' @export
capture_session_info <- function(filepath) {
  si <- utils::capture.output(utils::sessionInfo())
  writeLines(si, filepath)
  invisible(filepath)
}

#' Weighted Kaplan-Meier Risk
#'
#' Compute risk at time t from weighted survival data using a step-function
#' approach (for IPTW-based risk estimation).
#'
#' @param time Numeric vector of follow-up times.
#' @param event Binary event indicator.
#' @param weights Numeric vector of weights.
#' @param t_eval Time at which to evaluate risk.
#' @return Numeric scalar: estimated risk at time t_eval.
#' @export
weighted_km_risk <- function(time, event, weights, t_eval) {
  # Sort by time

  ord <- order(time)
  time <- time[ord]
  event <- event[ord]
  weights <- weights[ord]

  unique_times <- sort(unique(time[event == 1 & time <= t_eval]))
  if (length(unique_times) == 0) return(0)

  surv <- 1.0
  for (tj in unique_times) {
    at_risk <- time >= tj
    events_j <- (time == tj) & (event == 1)
    dj <- sum(weights[events_j])
    nj <- sum(weights[at_risk])
    if (nj > 0) {
      surv <- surv * (1 - dj / nj)
    }
  }
  1 - surv
}

#' Set Deterministic Seed
#'
#' Sets seed with L'Ecuyer-CMRG for parallel-safe reproducibility.
#'
#' @param seed Integer seed value.
#' @return Invisibly returns NULL.
#' @export
set_deterministic_seed <- function(seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  invisible(NULL)
}
