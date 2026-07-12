#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Verification of the RECALIBRATED Scenario C in run_simulation.R.
#
# Replicates run_simulation.R's OWN 5-covariate linear DGP exactly (generate_data
# / compute_truth sourced from dgp_scenarios.R, the single source shared with
# run_simulation.R) with the NEW unmeasured-confounder parameters, and runs a
# measured-covariate TMLE (the same
# cleanTMLE estimator the primary TMLE uses) to confirm realised bias and Wald
# coverage.
#
# New Scenario-C parameters (must match run_simulation.R lines 250-252):
#   U_prevalence = 0.40 ; U_trt_OR = 3.0 ; U_out_OR = 3.0
#
# Acceptance: realised TMLE bias ~0.03-0.05 ; coverage < 0.9 (materially broken).
# This does NOT run run_simulation.R end-to-end.
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(cleanTMLE)
  library(SuperLearner)
  library(parallel)
}))

## ---- Scenario-C parameters under test ---------------------------------------
SC <- list(U_prev = 0.40, U_trt_OR = 3.0, U_out_OR = 3.0)

N_OBS  <- 2000L      # config$n_obs
N_TRUTH<- 100000L    # config$n_truth
SEED   <- 2026L      # config$seed
REPS   <- as.integer(Sys.getenv("SC_REPS", "100"))

## ---- DGP: single source shared with run_simulation.R ------------------------
# Source the SAME generate_data()/compute_truth() the production run uses, so the
# verification cannot drift from what it verifies. Run from the repository root.
source("dgp_scenarios.R")

covs <- c("age", "sex", "biomarker", "comorbidity", "ckd")

one_rep <- function(rep_seed, pars, truth_rd) {
  dat <- generate_data(2000L, 0.5, -0.05, seed = rep_seed,
                       U_prevalence = pars$U_prev,
                       U_trt_OR = pars$U_trt_OR, U_out_OR = pars$U_out_OR)
  tryCatch({
    fit <- estimate_tmle_risk_point(
      data = dat, treatment = "treatment", outcome = "event_24",
      covariates = covs, sl_library = "SL.glm", truncate = 0.025, n_folds = 1L)
    ate <- fit$estimates$ATE
    c(est = ate$estimate, se = ate$se,
      covers = as.integer(ate$ci_lower <= truth_rd & truth_rd <= ate$ci_upper))
  }, error = function(e) c(est = NA_real_, se = NA_real_, covers = NA_real_))
}

## ---- Run ---------------------------------------------------------------------
truth <- compute_truth(N_TRUTH, 0.5, -0.05, seed = SEED,
                       U_prevalence = SC$U_prev, U_out_OR = SC$U_out_OR)
cat(sprintf("Scenario C recalibrated: U_prev=%.2f  U_trt_OR=%.1f  U_out_OR=%.1f\n",
            SC$U_prev, SC$U_trt_OR, SC$U_out_OR))
cat(sprintf("True RD (n_truth=%d) = %.5f\n", N_TRUTH, truth$RD))
cat(sprintf("Reps = %d, n = %d\n\n", REPS, N_OBS))

cl <- makeCluster(min(22L, detectCores()))
invisible(clusterEvalQ(cl, suppressWarnings(suppressMessages({
  library(cleanTMLE); library(SuperLearner)
}))))
clusterExport(cl, c("generate_data", "covs", "one_rep"), envir = environment())

seeds <- SEED + seq_len(REPS) * 1000L + 3L
res <- parLapply(cl, seeds, one_rep, pars = SC, truth_rd = truth$RD)
stopCluster(cl)

M   <- do.call(rbind, res)
ok  <- !is.na(M[, "est"])
ests<- M[ok, "est"]; ses <- M[ok, "se"]; cov <- M[ok, "covers"]
bias <- mean(ests) - truth$RD
emp_sd <- sd(ests); mean_se <- mean(ses)
coverage <- mean(cov)
mcse_bias <- emp_sd / sqrt(sum(ok))
mcse_cov  <- sqrt(coverage * (1 - coverage) / sum(ok))

cat("---- Realised primary TMLE (measured covariates only) ----\n")
cat(sprintf("  n_ok        : %d / %d\n", sum(ok), REPS))
cat(sprintf("  mean est    : %.5f\n", mean(ests)))
cat(sprintf("  true RD     : %.5f\n", truth$RD))
cat(sprintf("  BIAS        : %.5f  (MCSE %.5f)\n", bias, mcse_bias))
cat(sprintf("  emp SD      : %.5f\n", emp_sd))
cat(sprintf("  mean Wald SE: %.5f   (se/sd ratio %.3f)\n", mean_se, mean_se / emp_sd))
cat(sprintf("  COVERAGE    : %.3f  (MCSE %.3f)\n", coverage, mcse_cov))
cat("\n")
pass_bias <- abs(bias) >= 0.03 && abs(bias) <= 0.05
pass_cov  <- coverage < 0.9
cat(sprintf("  bias in [0.03,0.05]? %s   coverage < 0.9? %s\n",
            pass_bias, pass_cov))
