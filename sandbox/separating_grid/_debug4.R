suppressPackageStartupMessages({ library(cleanTMLE) })
ROOT <- "C:/Users/andre/OneDrive/Documents/clean-room-sim"
local({ exprs <- parse(file=file.path(ROOT,"run_simulation.R"))
  for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-"))) {
    lhs <- e[[2]]; if (is.name(lhs) && as.character(lhs) %in% c("generate_data","compute_truth")) eval(e, envir=globalenv()) } })
say <- function(...) { cat(..., "\n"); flush(stdout()) }
COVARS <- c("age","sex","biomarker","comorbidity","ckd"); ENS <- c("SL.glm","SL.glmnet","SL.gam")
d <- generate_data(n=4000, overlap_strength=1.5, effect_size=-0.05, misspec=FALSE, seed=20260608L)
lock <- create_simple_lock(data=d, treatment="treatment", outcome="event_24", covariates=COVARS, sl_library="SL.glm", plasmode_reps=3L, seed=20260608L)
ens025 <- tmle_candidate("ens_t025", g_library=ENS, q_library=ENS, truncation=0.025)
threats <- list(covariate_missingness=list(fractions=c(0.10)))
say("=== ens_t025 ALONE, run #1 ===")
r1 <- run_plasmode_dq_stress(lock, tmle_candidates=list(ens025), effect_sizes=c(0.05), reps=3L, data_quality_scenarios=threats, q0_library=NULL, verbose=FALSE)
say("run #1 OK rows:", nrow(r1$metrics))
say("=== ens_t025 ALONE, run #2 (same session) ===")
r2 <- run_plasmode_dq_stress(lock, tmle_candidates=list(ens025), effect_sizes=c(0.05), reps=3L, data_quality_scenarios=threats, q0_library=NULL, verbose=FALSE)
say("run #2 OK rows:", nrow(r2$metrics))
say("ALL DONE")
