#!/usr/bin/env Rscript
# Diagnostic: scan overlap strength and threat menu; print full per-scenario,
# per-candidate RMSE so we can see which knob makes truncation bite.
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

cands <- list(
  tmle_candidate("aggressive", g_library = "SL.glm", truncation = 0.001),
  tmle_candidate("middle",     g_library = "SL.glm", truncation = 0.025),
  tmle_candidate("robust",     g_library = "SL.glm", truncation = 0.20)
)

dq_spec <- list(
  treatment_misclass = list(sensitivity = c(0.90, 0.80), specificity = c(0.95, 0.85)),
  unmeasured_confounding = list(U_prevalence = 0.20,
                                U_treatment_OR = c(2.0, 3.0),
                                U_outcome_OR   = c(2.0, 3.0)),
  covariate_missingness = list(fractions = c(0.10, 0.20))
)

cov <- c("age","sex","biomarker","comorbidity","ckd")
SEED <- 20260530L

for (ov in c(2.0, 2.5, 3.0)) {
  ref <- generate_data(2000L, ov, -0.05, seed = SEED)
  ps <- predict(glm(reformulate(cov, "treatment"), ref, family = binomial()),
                type = "response")
  cat(sprintf("\n######## overlap = %.1f  PS range [%.4f, %.4f]  frac<.10=%.1f%%  frac<.025=%.1f%%  frac<.001=%.2f%%  trt=%.3f\n",
              ov, min(ps), max(ps), 100*mean(ps<.10), 100*mean(ps<.025),
              100*mean(ps<.001), mean(ref$treatment)))
  lock <- create_analysis_lock(ref, "treatment", "event_24", cov,
                               sl_library = c("SL.glm","SL.mean"),
                               plasmode_reps = 30L, seed = SEED + 1L)
  plas <- run_plasmode_feasibility(lock, cands, effect_sizes = c(0.05),
                                   reps = 30L, verbose = FALSE)
  cat("  baseline RMSE: ")
  bm <- tapply(plas$metrics$rmse, plas$metrics$candidate, mean)
  print(round(bm[c("aggressive","middle","robust")], 5))
  dq <- run_plasmode_dq_stress(lock, cands, effect_sizes = c(0.05), reps = 30L,
                               data_quality_scenarios = dq_spec,
                               fit_timeout = 30, verbose = FALSE)
  m <- dq$metrics
  m$key <- paste(m$scenario, m$level)
  wide <- reshape(m[, c("key","candidate","rmse")], idvar = "key",
                  timevar = "candidate", direction = "wide")
  names(wide) <- sub("rmse.", "", names(wide))
  print(wide[, c("key","aggressive","middle","robust")], row.names = FALSE)
}
cat("\nDONE\n")
