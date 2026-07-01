#!/usr/bin/env Rscript
# Diagnostic 3: combine (1) outcome misspecification (heavy truncation pays a
# baseline bias cost -> light wins baseline) with (2) treatment misclassification
# as a positivity-variance threat (light truncation explodes -> heavy wins
# worst-case). Goal direction: min_rmse = aggressive(0.001),
# min_max_rmse = robust(heavy).
suppressWarnings(suppressMessages(library(pkgload)))
.this_dir <- dirname(normalizePath(sub("^--file=", "",
              commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkgload::load_all(file.path(repo_root, "cleanTMLE"), quiet = TRUE)

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
  inter <- nl * ( 1.5*biomarker*ckd + 1.2*sex*comorbidity + 0.04*(age-55)*biomarker )
  lp_out <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*biomarker + 0.6*ckd +
            0.25*comorbidity + inter + effect_size/0.15*treatment + log(U_out_OR)*U
  event_24 <- rbinom(n, 1, plogis(lp_out))
  nc_outcome <- rbinom(n, 1, plogis(-1 + 0.01*(age-55) + 0.1*sex + 0.15*biomarker))
  data.frame(age = round(age,1), sex = sex, biomarker = round(biomarker,3),
             comorbidity = comorbidity, ckd = ckd, treatment = treatment,
             event_24 = event_24, nc_outcome = nc_outcome)
}

cov <- c("age","sex","biomarker","comorbidity","ckd")
Q0LIB <- "SL.glm.interaction"
SEED <- 20260530L

run_cell <- function(ov, nl, rob_trunc, reps = 50L) {
  cands <- list(
    tmle_candidate("aggressive", g_library = "SL.glm", truncation = 0.001),
    tmle_candidate("middle",     g_library = "SL.glm", truncation = 0.025),
    tmle_candidate("robust",     g_library = "SL.glm", truncation = rob_trunc)
  )
  dq_spec <- list(
    treatment_misclass = list(sensitivity = c(0.90, 0.85, 0.80),
                              specificity = c(0.85, 0.75, 0.65)),
    unmeasured_confounding = list(U_prevalence = 0.20,
                                  U_treatment_OR = c(2.0, 3.0),
                                  U_outcome_OR   = c(2.0, 3.0)),
    covariate_missingness = list(fractions = c(0.10, 0.20))
  )
  ref <- generate_data(2000L, ov, -0.05, seed = SEED, nl = nl)
  lock <- create_analysis_lock(ref, "treatment", "event_24", cov,
              sl_library = c("SL.glm","SL.mean"), plasmode_reps = reps,
              seed = SEED + 1L)
  plas <- run_plasmode_feasibility(lock, cands, effect_sizes = c(0.05),
              reps = reps, q0_library = Q0LIB, verbose = FALSE)
  dq <- run_plasmode_dq_stress(lock, cands, effect_sizes = c(0.05), reps = reps,
              data_quality_scenarios = dq_spec, q0_library = Q0LIB,
              fit_timeout = 30, verbose = FALSE)
  list(plas = plas, dq = dq, cands = cands)
}

scan <- expand.grid(ov = c(1.3, 1.6), nl = c(2.0, 3.0), rob = c(0.20, 0.35))
for (i in seq_len(nrow(scan))) {
  ov <- scan$ov[i]; nl <- scan$nl[i]; rob <- scan$rob[i]
  r <- run_cell(ov, nl, rob, reps = 50L)
  bm <- tapply(r$plas$metrics$rmse, r$plas$metrics$candidate, mean)
  m <- r$dq$metrics; md <- m[m$scenario != "none", ]
  wc <- tapply(md$rmse, md$candidate, function(x) max(x, na.rm = TRUE))
  wmin <- select_tmle_candidate(r$plas, rule = "min_rmse")$candidate_id
  wmm  <- select_tmle_candidate(r$plas, rule = "min_max_rmse",
                                dq_results = r$dq)$candidate_id
  cat(sprintf("\n=== ov=%.1f nl=%.1f rob_trunc=%.2f ===\n", ov, nl, rob))
  cat(sprintf("  baseline RMSE:  aggr=%.5f  mid=%.5f  rob=%.5f\n",
              bm["aggressive"], bm["middle"], bm["robust"]))
  cat(sprintf("  worstcase RMSE: aggr=%.5f  mid=%.5f  rob=%.5f\n",
              wc["aggressive"], wc["middle"], wc["robust"]))
  # which threat is each candidate's worst
  for (cid in c("aggressive","robust")) {
    sub <- md[md$candidate == cid, ]
    w <- sub[which.max(sub$rmse), ]
    cat(sprintf("    %s worst threat: %s[%s] rmse=%.5f\n",
                cid, w$scenario, w$level, w$rmse))
  }
  cat(sprintf("  >> min_rmse=%s   min_max_rmse=%s   %s\n",
              wmin, wmm, if (wmin != wmm) "**DIVERGE**" else "agree"))
}
cat("\nDONE\n")
