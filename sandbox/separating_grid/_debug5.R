suppressPackageStartupMessages({ library(cleanTMLE) })
ROOT <- "C:/Users/andre/OneDrive/Documents/clean-room-sim"
local({ exprs <- parse(file=file.path(ROOT,"run_simulation.R"))
  for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-"))) {
    lhs <- e[[2]]; if (is.name(lhs) && as.character(lhs) %in% c("generate_data","compute_truth")) eval(e, envir=globalenv()) } })
say <- function(...) { cat(..., "\n"); flush(stdout()) }
COVARS <- c("age","sex","biomarker","comorbidity","ckd")
ENS2 <- c("SL.glm","SL.glmnet","SL.glm.interaction")   # gam-free
d <- generate_data(n=4000, overlap_strength=1.5, effect_size=-0.05, misspec=FALSE, seed=20260608L)
lock <- create_simple_lock(data=d, treatment="treatment", outcome="event_24", covariates=COVARS, sl_library="SL.glm", plasmode_reps=3L, seed=20260608L)
cand <- tmle_candidate("ens2", g_library=ENS2, q_library=ENS2, truncation=0.025)
threats <- list(covariate_missingness=list(fractions=c(0.10)),
                unmeasured_confounding=list(U_prevalence=0.20, U_treatment_OR=c(2.0), U_outcome_OR=c(2.0)))
for (i in 1:6) {
  say(sprintf("=== gam-free ensemble call %d/6 ===", i))
  r <- run_plasmode_dq_stress(lock, tmle_candidates=list(cand), effect_sizes=c(0.05), reps=3L,
                              data_quality_scenarios=threats, q0_library=NULL, verbose=FALSE)
  say("  OK rows:", nrow(r$metrics))
}
say("ALL DONE")
