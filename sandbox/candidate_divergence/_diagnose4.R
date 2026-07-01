#!/usr/bin/env Rscript
# Diagnostic 4: FORWARD direction at the ORIGINAL (correctly specified) DGP.
# Hypothesis: at moderate overlap (~1.5) tiny truncation (aggressive) edges out
# heavy truncation on BASELINE RMSE (truncation clips a few PS -> tiny bias),
# while TREATMENT MISCLASSIFICATION (a weight-tail / positivity threat) makes
# tiny truncation explode -> robust wins worst-case. Need aggressive lowest
# baseline AND aggressive highest worst-case AND robust worst-case under thresh.
suppressWarnings(suppressMessages(library(pkgload)))
.this_dir <- dirname(normalizePath(sub("^--file=", "",
              commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkgload::load_all(file.path(repo_root, "cleanTMLE"), quiet = TRUE)

generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL, U_prevalence = 0, U_trt_OR = 1, U_out_OR = 1) {
  if (!is.null(seed)) set.seed(seed)
  age <- rnorm(n, 55, 10); sex <- rbinom(n, 1, 0.55)
  biomarker <- rnorm(n, 0, 1); comorbidity <- sample(0:2, n, TRUE, c(.5,.3,.2))
  ckd <- rbinom(n, 1, 0.12)
  U <- if (U_prevalence > 0) rbinom(n, 1, U_prevalence) else rep(0L, n)
  lp_trt <- -0.5 + overlap_strength * (0.03*(age-55) + 0.8*sex + 0.6*biomarker +
            0.5*ckd + 0.3*comorbidity) + log(U_trt_OR)*U
  treatment <- rbinom(n, 1, plogis(lp_trt))
  lp_out <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*biomarker + 0.6*ckd +
            0.25*comorbidity + effect_size/0.15*treatment + log(U_out_OR)*U
  event_24 <- rbinom(n, 1, plogis(lp_out))
  nc_outcome <- rbinom(n, 1, plogis(-1 + 0.01*(age-55) + 0.1*sex + 0.15*biomarker))
  data.frame(age = round(age,1), sex = sex, biomarker = round(biomarker,3),
             comorbidity = comorbidity, ckd = ckd, treatment = treatment,
             event_24 = event_24, nc_outcome = nc_outcome)
}

cov <- c("age","sex","biomarker","comorbidity","ckd")
SEED <- 20260530L

run_cell <- function(ov, rob_trunc, reps = 60L, with_bias_threats = TRUE) {
  cands <- list(
    tmle_candidate("aggressive", g_library = "SL.glm", truncation = 0.001),
    tmle_candidate("middle",     g_library = "SL.glm", truncation = 0.025),
    tmle_candidate("robust",     g_library = "SL.glm", truncation = rob_trunc)
  )
  dq_spec <- list(
    treatment_misclass = list(sensitivity = c(0.90, 0.80, 0.70),
                              specificity = c(0.90, 0.80, 0.70))
  )
  if (with_bias_threats) {
    dq_spec$unmeasured_confounding <- list(U_prevalence = 0.20,
        U_treatment_OR = c(1.5, 2.0), U_outcome_OR = c(1.5, 2.0))
    dq_spec$covariate_missingness <- list(fractions = c(0.10, 0.20))
  }
  ref <- generate_data(2000L, ov, -0.05, seed = SEED)
  lock <- create_analysis_lock(ref, "treatment", "event_24", cov,
              sl_library = c("SL.glm","SL.mean"), plasmode_reps = reps,
              seed = SEED + 1L)
  plas <- run_plasmode_feasibility(lock, cands, effect_sizes = c(0.05),
              reps = reps, verbose = FALSE)
  dq <- run_plasmode_dq_stress(lock, cands, effect_sizes = c(0.05), reps = reps,
              data_quality_scenarios = dq_spec, fit_timeout = 30, verbose = FALSE)
  list(plas = plas, dq = dq)
}

scan <- expand.grid(ov = c(1.4, 1.5, 1.6, 1.7), rob = c(0.10, 0.15))
for (i in seq_len(nrow(scan))) {
  ov <- scan$ov[i]; rob <- scan$rob[i]
  r <- run_cell(ov, rob, reps = 60L)
  bm <- tapply(r$plas$metrics$rmse, r$plas$metrics$candidate, mean)
  m <- r$dq$metrics; md <- m[m$scenario != "none", ]
  wc <- tapply(md$rmse, md$candidate, function(x) max(x, na.rm = TRUE))
  wmin <- select_tmle_candidate(r$plas, rule = "min_rmse")$candidate_id
  wmm  <- suppressMessages(select_tmle_candidate(r$plas, rule = "min_max_rmse",
                                dq_results = r$dq)$candidate_id)
  cat(sprintf("\n=== ov=%.1f rob_trunc=%.2f ===\n", ov, rob))
  cat(sprintf("  baseline RMSE:  aggr=%.5f  mid=%.5f  rob=%.5f  (aggr lowest: %s)\n",
              bm["aggressive"], bm["middle"], bm["robust"],
              which.min(bm[c("aggressive","middle","robust")]) == 1))
  cat(sprintf("  worstcase RMSE: aggr=%.5f  mid=%.5f  rob=%.5f  (aggr highest: %s)\n",
              wc["aggressive"], wc["middle"], wc["robust"],
              which.max(wc[c("aggressive","middle","robust")]) == 1))
  for (cid in c("aggressive","robust")) {
    sub <- md[md$candidate == cid, ]; w <- sub[which.max(sub$rmse), ]
    cat(sprintf("    %s worst threat: %s[%s] rmse=%.5f\n",
                cid, w$scenario, w$level, w$rmse))
  }
  cat(sprintf("  >> min_rmse=%s   min_max_rmse=%s   %s\n",
              wmin, wmm, if (wmin != wmm) "**DIVERGE**" else "agree"))
}
cat("\nDONE\n")
