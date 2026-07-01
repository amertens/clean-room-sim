#!/usr/bin/env Rscript
# Tuning 2: widen the truncation gap so robust is more U-fragile vs middle.
# Fit candidates {0.001, 0.025, 0.10, 0.15, 0.20} once over the full threat
# sweep at two overlaps, then evaluate the robust-vs-middle margin for each
# choice of robust truncation.
suppressWarnings(suppressMessages(library(pkgload)))
.this_dir <- dirname(normalizePath(sub("^--file=", "",
              commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkgload::load_all(file.path(repo_root, "cleanTMLE"), quiet = TRUE)
source(file.path(.this_dir, "_tune_helpers.R"), local = TRUE)

SEED <- 20260530L; COV <- c("age","sex","biomarker","comorbidity","ckd"); reps <- 60L
cands <- list(
  tmle_candidate("aggressive", g_library="SL.glm", truncation=0.001),
  tmle_candidate("middle",     g_library="SL.glm", truncation=0.025),
  tmle_candidate("rob10",      g_library="SL.glm", truncation=0.10),
  tmle_candidate("rob15",      g_library="SL.glm", truncation=0.15),
  tmle_candidate("rob20",      g_library="SL.glm", truncation=0.20))

for (OV in c(1.6, 2.0)) {
  ref <- generate_data(2000L, OV, -0.05, seed = SEED)
  ps <- predict(glm(reformulate(COV,"treatment"), ref, family=binomial()), type="response")
  lock <- create_analysis_lock(ref, "treatment", "event_24", COV,
            sl_library=c("SL.glm","SL.mean"), plasmode_reps=reps, seed=SEED+1L)
  plas <- run_plasmode_feasibility(lock, cands, effect_sizes=c(0.05), reps=reps, verbose=FALSE)
  br <- tapply(plas$metrics$rmse, plas$metrics$candidate, mean)
  dq_pkg <- run_plasmode_dq_stress(lock, cands, effect_sizes=c(0.05), reps=reps,
              data_quality_scenarios=list(
                unmeasured_confounding=list(U_prevalence=0.20,
                    U_treatment_OR=c(3,4,5,6,7), U_outcome_OR=c(3,4,5,6,7)),
                covariate_missingness=list(fractions=c(0.10,0.20))),
              fit_timeout=30, verbose=FALSE)
  pos <- near_positivity_stress(lock, cands, slopes=c(2,3,4), reps=reps, effect_sizes=c(0.05))
  allm <- rbind(dq_pkg$metrics[,c("scenario","level","candidate","rmse")],
                pos$metrics[,c("scenario","level","candidate","rmse")])
  allm <- allm[allm$scenario!="none",]
  wc <- tapply(allm$rmse, allm$candidate, max)
  cat(sprintf("\n#### overlap=%.1f  PS[%.4f,%.4f] frac<.10=%.1f%% frac<.025=%.1f%%\n",
      OV, min(ps), max(ps), 100*mean(ps<.10), 100*mean(ps<.025)))
  cat(sprintf("  baseline RMSE: a=%.5f m=%.5f r10=%.5f r15=%.5f r20=%.5f\n",
      br["aggressive"], br["middle"], br["rob10"], br["rob15"], br["rob20"]))
  cat(sprintf("  worstcase RMSE: a=%.5f m=%.5f r10=%.5f r15=%.5f r20=%.5f\n",
      wc["aggressive"], wc["middle"], wc["rob10"], wc["rob15"], wc["rob20"]))
  for (rb in c("rob10","rob15","rob20")) {
    set3 <- allm[allm$candidate %in% c("aggressive","middle",rb),]
    w3 <- tapply(set3$rmse, set3$candidate, max)
    b3 <- br[c("aggressive","middle",rb)]
    minr <- names(b3)[which.min(b3)]; minmax <- names(w3)[which.min(w3)]
    cat(sprintf("    {agg,mid,%s}: min_rmse=%s min_max=%s  margin(%s_wc - mid_wc)=%.5f  %s\n",
        rb, minr, minmax, rb, unname(w3[rb]-w3["middle"]),
        if (minr!=minmax) "DIVERGE" else "agree"))
  }
}
cat("\nDONE\n")
