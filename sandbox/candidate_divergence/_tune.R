#!/usr/bin/env Rscript
# Tuning: run every threat across a severity RANGE once, then evaluate minimax
# winners for candidate threat-set definitions in memory. Goal: min_rmse=robust,
# min_max_rmse=middle with a clear margin.
suppressWarnings(suppressMessages(library(pkgload)))
.this_dir <- dirname(normalizePath(sub("^--file=", "",
              commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkgload::load_all(file.path(repo_root, "cleanTMLE"), quiet = TRUE)
source(file.path(.this_dir, "_tune_helpers.R"), local = TRUE)  # generate_data, near_positivity_stress

SEED <- 20260530L; OVERLAP <- 1.6; COV <- c("age","sex","biomarker","comorbidity","ckd")
ref <- generate_data(2000L, OVERLAP, -0.05, seed = SEED)
cands <- list(
  tmle_candidate("aggressive", g_library="SL.glm", truncation=0.001),
  tmle_candidate("middle",     g_library="SL.glm", truncation=0.025),
  tmle_candidate("robust",     g_library="SL.glm", truncation=0.10))
lock <- create_analysis_lock(ref, "treatment", "event_24", COV,
          sl_library=c("SL.glm","SL.mean"), plasmode_reps=60L, seed=SEED+1L)

reps <- 60L
plas <- run_plasmode_feasibility(lock, cands, effect_sizes=c(0.05), reps=reps, verbose=FALSE)
base_rmse <- tapply(plas$metrics$rmse, plas$metrics$candidate, mean)

dq_pkg <- run_plasmode_dq_stress(lock, cands, effect_sizes=c(0.05), reps=reps,
            data_quality_scenarios=list(
              unmeasured_confounding=list(U_prevalence=0.20,
                  U_treatment_OR=c(2,3,4,5,6), U_outcome_OR=c(2,3,4,5,6)),
              covariate_missingness=list(fractions=c(0.10,0.20))),
            fit_timeout=30, verbose=FALSE)
pos <- near_positivity_stress(lock, cands, slopes=c(2,3,4), reps=reps, effect_sizes=c(0.05))

allm <- rbind(dq_pkg$metrics[, c("scenario","level","candidate","rmse")],
              pos$metrics[, c("scenario","level","candidate","rmse")])
allm <- allm[allm$scenario != "none", ]
cat("baseline RMSE:", sprintf("a=%.5f m=%.5f r=%.5f", base_rmse["aggressive"],
    base_rmse["middle"], base_rmse["robust"]), "\n\n")
w <- reshape(allm, idvar=c("scenario","level"), timevar="candidate", direction="wide")
names(w) <- sub("rmse.","",names(w))
print(w[, c("scenario","level","aggressive","middle","robust")], row.names=FALSE)

minimax <- function(sub) {
  wc <- tapply(sub$rmse, sub$candidate, max)
  names(wc)[which.min(wc)]
}
cat("\nMinimax winner for candidate threat sets:\n")
sets <- list(
  "pos{2,3,4}+U{2..6}+miss" = allm,
  "pos{2,3}+U{3,4,5}+miss"  = allm[!(allm$scenario=="near_positivity" & allm$level=="slope_x4.0") &
                                    !(allm$scenario=="unmeasured_U" & allm$level %in% c("OR_trt2.0_out2.0","OR_trt6.0_out6.0")), ],
  "pos{2,3}+U{4,5,6}+miss"  = allm[!(allm$scenario=="near_positivity" & allm$level=="slope_x4.0") &
                                    !(allm$scenario=="unmeasured_U" & allm$level %in% c("OR_trt2.0_out2.0","OR_trt3.0_out3.0")), ],
  "pos{2,3,4}+U{4,5,6}+miss"= allm[!(allm$scenario=="unmeasured_U" & allm$level %in% c("OR_trt2.0_out2.0","OR_trt3.0_out3.0")), ]
)
for (nm in names(sets)) {
  s <- sets[[nm]]; wc <- tapply(s$rmse, s$candidate, max)
  cat(sprintf("  %-28s wc: a=%.4f m=%.4f r=%.4f -> minmax=%s\n", nm,
      wc["aggressive"], wc["middle"], wc["robust"], names(wc)[which.min(wc)]))
}
cat("\nDONE\n")
