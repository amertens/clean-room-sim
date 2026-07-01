suppressPackageStartupMessages({ library(cleanTMLE) })
ROOT <- "C:/Users/andre/OneDrive/Documents/clean-room-sim"
local({ exprs <- parse(file=file.path(ROOT,"run_simulation.R"))
  for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-"))) {
    lhs <- e[[2]]; if (is.name(lhs) && as.character(lhs) %in% c("generate_data","compute_truth")) eval(e, envir=globalenv()) } })
say <- function(...) { cat(..., "\n"); flush(stdout()) }
COVARS <- c("age","sex","biomarker","comorbidity","ckd")
ENS <- c("SL.glm","SL.glmnet","SL.glm.interaction")
threatsA <- list(covariate_missingness=list(fractions=c(0.05,0.10,0.20)),
                 treatment_misclass=list(sensitivity=c(0.95,0.90), specificity=c(0.99,0.95)),
                 outcome_misclass=list(sensitivity=c(0.95,0.90), specificity=c(0.99,0.95)),
                 unmeasured_confounding=list(U_prevalence=0.20, U_treatment_OR=c(1.5,2.0), U_outcome_OR=c(1.5,2.0)))
d <- generate_data(n=4000, overlap_strength=1.5, effect_size=-0.05, misspec=FALSE, seed=20260608L)
lock <- create_simple_lock(data=d, treatment="treatment", outcome="event_24", covariates=COVARS, sl_library="SL.glm", plasmode_reps=10L, seed=20260608L)
glm_grid <- list(tmle_candidate("glm_t01",g_library="SL.glm",q_library="SL.glm",truncation=0.01),
                 tmle_candidate("glm_t05",g_library="SL.glm",q_library="SL.glm",truncation=0.05))
ens_grid <- list(tmle_candidate("ens_t01",g_library=ENS,q_library=ENS,truncation=0.01))
timeit <- function(lab, expr) { t<-Sys.time(); force(expr); el<-as.numeric(difftime(Sys.time(),t,units="secs")); say(sprintf(">>> %s: %.1f s", lab, el)); el }
say("timing GLM-only grid (2 cands), full Arm-A threats, reps=10 ...")
gA <- timeit("GLM x2, reps10, full threats", { run_plasmode_feasibility(lock, glm_grid, c(0.05), 10L, q0_library=NULL, verbose=FALSE); run_plasmode_dq_stress(lock, glm_grid, c(0.05), 10L, data_quality_scenarios=threatsA, q0_library=NULL, verbose=FALSE) })
say("timing ENSEMBLE grid (1 cand), full Arm-A threats, reps=10 ...")
eA <- timeit("ENS x1, reps10, full threats", { run_plasmode_feasibility(lock, ens_grid, c(0.05), 10L, q0_library=NULL, verbose=FALSE); run_plasmode_dq_stress(lock, ens_grid, c(0.05), 10L, data_quality_scenarios=threatsA, q0_library=NULL, verbose=FALSE) })
say(sprintf("\nPER-REP: GLM(2 cands)=%.2fs/rep  ENS(1 cand)=%.2fs/rep", gA/10, eA/10))
say("DONE")
