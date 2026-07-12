#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Calibration sweep for run_simulation.R's Scenario C (unmeasured confounding).
#
# Replicates run_simulation.R's OWN 5-covariate linear DGP exactly (generate_data
# / compute_truth copied verbatim, non-misspec branch) and reuses cleanTMLE's
# measured-covariate TMLE to find (U_prevalence, U_trt_OR, U_out_OR) such that the
# realised primary TMLE has bias ~0.04 and coverage ~0.6-0.8.
#
# NOTE: this does NOT run run_simulation.R end-to-end; it only replicates the
# Scenario-C generator + a single measured-covariate TMLE per rep.
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(cleanTMLE)
  library(SuperLearner)
  library(parallel)
}))

## ---- DGP: single source shared with run_simulation.R ------------------------
# Source the SAME generate_data()/compute_truth() the production run uses, so the
# calibration cannot drift from the production DGP. Run from the repository root.
source("dgp_scenarios.R")

## ---- Single measured-covariate TMLE (mirrors the primary TMLE) ---------------
covs <- c("age", "sex", "biomarker", "comorbidity", "ckd")

one_rep <- function(rep_seed, pars, truth_rd) {
  dat <- generate_data(2000L, overlap_strength = 0.5, effect_size = -0.05,
                       seed = rep_seed,
                       U_prevalence = pars$U_prev,
                       U_trt_OR = pars$U_trt_OR, U_out_OR = pars$U_out_OR)
  out <- tryCatch({
    fit <- estimate_tmle_risk_point(
      data = dat, treatment = "treatment", outcome = "event_24",
      covariates = covs, sl_library = "SL.glm",
      truncate = 0.025, n_folds = 1L)
    ate <- fit$estimates$ATE
    covers <- as.integer(ate$ci_lower <= truth_rd & truth_rd <= ate$ci_upper)
    c(est = ate$estimate, covers = covers)
  }, error = function(e) c(est = NA_real_, covers = NA_real_))
  out
}

## ---- Sweep -------------------------------------------------------------------
grid <- list(
  list(U_prev = 0.30, U_trt_OR = 2.5, U_out_OR = 2.5),
  list(U_prev = 0.40, U_trt_OR = 2.5, U_out_OR = 2.5),
  list(U_prev = 0.30, U_trt_OR = 3.0, U_out_OR = 3.0),
  list(U_prev = 0.40, U_trt_OR = 3.0, U_out_OR = 3.0),
  list(U_prev = 0.50, U_trt_OR = 3.0, U_out_OR = 3.0),
  list(U_prev = 0.35, U_trt_OR = 3.5, U_out_OR = 3.5)
)

REPS <- as.integer(Sys.getenv("SC_REPS", "60"))
base_seed <- 2026L

cl <- makeCluster(min(22L, detectCores()))
clusterEvalQ(cl, suppressWarnings(suppressMessages({
  library(cleanTMLE); library(SuperLearner)
})))
clusterExport(cl, c("generate_data", "estimate_tmle_risk_point", "covs", "one_rep"),
              envir = environment())

cat(sprintf("Calibration sweep: %d reps/cell, n=2000\n\n", REPS))
cat(sprintf("%-32s %8s %8s %9s %9s\n",
            "cell (Uprev/OR_A/OR_Y)", "true_rd", "bias", "coverage", "n_ok"))

for (pars in grid) {
  truth <- compute_truth(100000L, 0.5, -0.05, seed = base_seed,
                         U_prevalence = pars$U_prev, U_out_OR = pars$U_out_OR)
  seeds <- base_seed + seq_len(REPS) * 1000L + 3L  # mimic rep_seed offset shape
  res <- parLapply(cl, seeds, one_rep, pars = pars, truth_rd = truth$RD)
  M <- do.call(rbind, res)
  ok <- !is.na(M[, "est"])
  bias <- mean(M[ok, "est"]) - truth$RD
  cover <- mean(M[ok, "covers"], na.rm = TRUE)
  cat(sprintf("%-32s %8.4f %8.4f %9.3f %9d\n",
              sprintf("%.2f / %.1f / %.1f", pars$U_prev, pars$U_trt_OR, pars$U_out_OR),
              truth$RD, bias, cover, sum(ok)))
}

stopCluster(cl)
cat("\nDone.\n")
