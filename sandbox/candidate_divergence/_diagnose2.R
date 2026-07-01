#!/usr/bin/env Rscript
# Diagnostic 2: introduce outcome-model misspecification so PS truncation
# becomes a real bias/variance lever. True outcome surface is nonlinear
# (interactions); plasmode Q0 generator captures it via SL.glm.interaction,
# while candidates keep a main-effects SL.glm Q (misspecified).
suppressWarnings(suppressMessages(library(pkgload)))
.this_dir <- dirname(normalizePath(sub("^--file=", "",
              commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkgload::load_all(file.path(repo_root, "cleanTMLE"), quiet = TRUE)

# generate_data with an interaction-driven (nonlinear) outcome surface.
generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL, U_prevalence = 0, U_trt_OR = 1,
                          U_out_OR = 1, nl = 0) {
  if (!is.null(seed)) set.seed(seed)
  age <- rnorm(n, 55, 10); sex <- rbinom(n, 1, 0.55)
  biomarker <- rnorm(n, 0, 1); comorbidity <- sample(0:2, n, TRUE, c(.5,.3,.2))
  ckd <- rbinom(n, 1, 0.12)
  U <- if (U_prevalence > 0) rbinom(n, 1, U_prevalence) else rep(0L, n)
  lp_trt <- -0.5 + overlap_strength * (0.03*(age-55) + 0.8*sex + 0.6*biomarker +
            0.5*ckd + 0.3*comorbidity) + log(U_trt_OR)*U
  treatment <- rbinom(n, 1, plogis(lp_trt))
  # Nonlinear / interaction outcome surface (strength controlled by nl):
  inter <- nl * ( 1.5*biomarker*ckd + 1.2*sex*comorbidity +
                  0.04*(age-55)*biomarker )
  lp_out <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*biomarker + 0.6*ckd +
            0.25*comorbidity + inter + effect_size/0.15*treatment + log(U_out_OR)*U
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
  unmeasured_confounding = list(U_prevalence = 0.20,
                                U_treatment_OR = c(2.0, 3.0, 4.0),
                                U_outcome_OR   = c(2.0, 3.0, 4.0)),
  covariate_missingness = list(fractions = c(0.10, 0.20))
)
cov <- c("age","sex","biomarker","comorbidity","ckd")
Q0LIB <- "SL.glm.interaction"
SEED <- 20260530L

scan <- expand.grid(ov = c(1.0, 1.5, 2.0), nl = c(1.0, 2.0))
for (i in seq_len(nrow(scan))) {
  ov <- scan$ov[i]; nl <- scan$nl[i]
  ref <- generate_data(2000L, ov, -0.05, seed = SEED, nl = nl)
  ps <- predict(glm(reformulate(cov, "treatment"), ref, family = binomial()),
                type = "response")
  cat(sprintf("\n######## overlap=%.1f nl=%.1f  PS[%.4f,%.4f] frac<.10=%.1f%% trt=%.3f evt=%.3f\n",
              ov, nl, min(ps), max(ps), 100*mean(ps<.10),
              mean(ref$treatment), mean(ref$event_24)))
  lock <- create_analysis_lock(ref, "treatment", "event_24", cov,
                               sl_library = c("SL.glm","SL.mean"),
                               plasmode_reps = 40L, seed = SEED + 1L)
  plas <- run_plasmode_feasibility(lock, cands, effect_sizes = c(0.05),
              reps = 40L, q0_library = Q0LIB, verbose = FALSE)
  bm <- tapply(plas$metrics$rmse, plas$metrics$candidate, mean)
  bb <- tapply(plas$metrics$bias, plas$metrics$candidate, mean)
  cat("  baseline RMSE:", sprintf("aggr=%.5f mid=%.5f rob=%.5f",
      bm["aggressive"], bm["middle"], bm["robust"]),
      " bias:", sprintf("aggr=%.5f mid=%.5f rob=%.5f",
      bb["aggressive"], bb["middle"], bb["robust"]), "\n")
  dq <- run_plasmode_dq_stress(lock, cands, effect_sizes = c(0.05), reps = 40L,
              data_quality_scenarios = dq_spec, q0_library = Q0LIB,
              fit_timeout = 30, verbose = FALSE)
  m <- dq$metrics; m$key <- paste(m$scenario, m$level)
  wide <- reshape(m[, c("key","candidate","rmse")], idvar = "key",
                  timevar = "candidate", direction = "wide")
  names(wide) <- sub("rmse.", "", names(wide))
  print(wide[, c("key","aggressive","middle","robust")], row.names = FALSE)
}
cat("\nDONE\n")
