#!/usr/bin/env Rscript
# Verify the near_positivity threat is now a first-class package threat and that
# the three-truncation grid separates under the package's own
# run_plasmode_dq_stress + select_tmle_candidate (no sandbox machinery).
suppressWarnings(suppressMessages(library(pkgload)))
.d <- dirname(normalizePath(sub("^--file=", "",
        commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
pkgload::load_all(file.path(.d, "..", "..", "cleanTMLE"), quiet = TRUE)

generate_data <- function(n, overlap_strength, effect_size = -0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  age <- rnorm(n,55,10); sex <- rbinom(n,1,0.55); biomarker <- rnorm(n,0,1)
  comorbidity <- sample(0:2,n,TRUE,c(.5,.3,.2)); ckd <- rbinom(n,1,0.12)
  lp_trt <- -0.5 + overlap_strength*(0.03*(age-55)+0.8*sex+0.6*biomarker+0.5*ckd+0.3*comorbidity)
  treatment <- rbinom(n,1,plogis(lp_trt))
  lp_out <- -2.5 + 0.015*(age-55)+0.3*sex+0.2*biomarker+0.6*ckd+0.25*comorbidity +
            effect_size/0.15*treatment
  event_24 <- rbinom(n,1,plogis(lp_out))
  nc_outcome <- rbinom(n,1,plogis(-1+0.01*(age-55)+0.1*sex+0.15*biomarker))
  data.frame(age=round(age,1),sex=sex,biomarker=round(biomarker,3),
             comorbidity=comorbidity,ckd=ckd,treatment=treatment,
             event_24=event_24,nc_outcome=nc_outcome)
}

cov <- c("age","sex","biomarker","comorbidity","ckd")
ref <- generate_data(2000L, 1.6, -0.05, seed = 20260530L)
lock <- create_analysis_lock(ref, "treatment", "event_24", cov,
          sl_library=c("SL.glm","SL.mean"), plasmode_reps=40L, seed=20260531L)
cands <- list(
  tmle_candidate("aggressive", g_library="SL.glm", truncation=0.001),
  tmle_candidate("middle",     g_library="SL.glm", truncation=0.025),
  tmle_candidate("robust",     g_library="SL.glm", truncation=0.20))

plas <- run_plasmode_feasibility(lock, cands, effect_sizes=c(0.05), reps=40L, verbose=FALSE)
dq <- run_plasmode_dq_stress(lock, cands, effect_sizes=c(0.05), reps=40L,
        data_quality_scenarios = list(
          near_positivity = list(slopes = c(2,3,4)),
          unmeasured_confounding = list(U_prevalence=0.20,
              U_treatment_OR=c(3,4,5,6,7,8), U_outcome_OR=c(3,4,5,6,7,8)),
          covariate_missingness = list(fractions=c(0.10,0.20))),
        fit_timeout = 30, verbose = FALSE)

cat("Scenarios run:\n"); print(unique(dq$metrics[,c("scenario","level")]))
bm <- tapply(plas$metrics$rmse, plas$metrics$candidate, mean)
md <- dq$metrics[dq$metrics$scenario!="none",]
wc <- tapply(md$rmse, md$candidate, max)
wmin <- suppressMessages(select_tmle_candidate(plas, "min_rmse")$candidate_id)
wmm  <- suppressMessages(select_tmle_candidate(plas, "min_max_rmse", dq_results=dq)$candidate_id)
cat(sprintf("\nbaseline RMSE: aggr=%.5f mid=%.5f rob=%.5f\n", bm["aggressive"],bm["middle"],bm["robust"]))
cat(sprintf("worstcase RMSE: aggr=%.5f mid=%.5f rob=%.5f\n", wc["aggressive"],wc["middle"],wc["robust"]))
cat(sprintf(">> min_rmse=%s  min_max_rmse=%s  %s\n", wmin, wmm,
            if (wmin!=wmm) "**DIVERGE**" else "agree"))
cat("DONE\n")
