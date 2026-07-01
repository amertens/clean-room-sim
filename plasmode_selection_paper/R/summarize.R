# Post-processing for the plasmode-selection paper.
#
# Turns the per-replicate rows saved by run_simulation() into a per-(DGP,
# workflow) summary that carries an explicit Monte Carlo standard error for
# every reported operating characteristic (bias, RMSE, coverage). This is what
# lets the manuscript show whether two workflows actually differ or are within
# Monte Carlo noise (TODO A.6).

#' Summarise a plasmode-selection simulation with Monte Carlo SEs
#'
#' @param results per-replicate data.frame from [run_simulation()] (or a saved
#'   `sim_*.rds`), with at least columns `dgp`, `workflow`, `estimate`,
#'   `ci_lower`, `ci_upper` (and optionally `se`).
#' @param true_RD numeric; the marginal risk difference the DGPs target.
#'
#' @return data.frame, one row per (`dgp`, `workflow`), with `n_rep`, `bias`,
#'   `mc_se_bias`, `rmse`, `mc_se_rmse`, `coverage`, `mc_se_coverage`,
#'   `emp_se`, `mean_se`, `se_cal`, and `mean_ci_width`. The `mc_se_*` columns
#'   are the Monte Carlo standard errors of the corresponding estimates.
summarize_simulation <- function(results, true_RD = -0.05) {
  stopifnot(is.data.frame(results),
            all(c("dgp", "workflow", "estimate", "ci_lower", "ci_upper")
                %in% names(results)))

  key   <- interaction(results$dgp, results$workflow, drop = TRUE, sep = "::")
  parts <- split(results, key)

  rows <- lapply(parts, function(df) {
    ok  <- is.finite(df$estimate)
    est <- df$estimate[ok]
    n   <- length(est)
    if (n < 1L) return(NULL)

    err    <- est - true_RD
    bias   <- mean(err)
    emp_se <- if (n > 1L) stats::sd(est) else NA_real_
    mc_se_bias <- if (n > 1L) emp_se / sqrt(n) else NA_real_

    rmse   <- sqrt(mean(err^2))
    # Monte Carlo SE of RMSE via the delta method on mean(err^2):
    # SE(mean sq) = sd(sq)/sqrt(n); SE(sqrt(.)) = SE(mean sq) / (2 * rmse).
    mc_se_rmse <- if (n > 1L && rmse > 0)
      sqrt(stats::var(err^2) / n) / (2 * rmse) else NA_real_

    lo <- df$ci_lower[ok]; hi <- df$ci_upper[ok]
    ci_ok    <- is.finite(lo) & is.finite(hi)
    covered  <- (lo[ci_ok] <= true_RD) & (hi[ci_ok] >= true_RD)
    ncov     <- length(covered)
    coverage <- if (ncov > 0L) mean(covered) else NA_real_
    mc_se_cov <- if (ncov > 0L)
      sqrt(coverage * (1 - coverage) / ncov) else NA_real_

    mean_se <- if ("se" %in% names(df)) mean(df$se[ok], na.rm = TRUE) else NA_real_

    data.frame(
      dgp            = df$dgp[1L],
      workflow       = df$workflow[1L],
      n_rep          = n,
      bias           = bias,
      mc_se_bias     = mc_se_bias,
      rmse           = rmse,
      mc_se_rmse     = mc_se_rmse,
      coverage       = coverage,
      mc_se_coverage = mc_se_cov,
      emp_se         = emp_se,
      mean_se        = mean_se,
      se_cal         = if (is.finite(emp_se) && emp_se > 0)
                         mean_se / emp_se else NA_real_,
      mean_ci_width  = mean(hi[ci_ok] - lo[ci_ok]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  rownames(out) <- NULL
  out[order(out$dgp, out$workflow), , drop = FALSE]
}
