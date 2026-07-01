# Pinpoint which candidate crashes run_plasmode_dq_stress (Arm A, reps=3).
suppressPackageStartupMessages({ library(cleanTMLE) })
ROOT <- "C:/Users/andre/OneDrive/Documents/clean-room-sim"
local({
  exprs <- parse(file = file.path(ROOT, "run_simulation.R"))
  for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-"))) {
    lhs <- e[[2]]
    if (is.name(lhs) && as.character(lhs) %in% c("generate_data","compute_truth"))
      eval(e, envir = globalenv())
  }
})
say <- function(...) { cat(..., "\n"); flush(stdout()) }
COVARS <- c("age","sex","biomarker","comorbidity","ckd")
ENS <- c("SL.glm","SL.glmnet","SL.gam")
cands <- list(
  glm_t01    = tmle_candidate("glm_t01", g_library="SL.glm", q_library="SL.glm", truncation=0.01),
  glm_t05    = tmle_candidate("glm_t05", g_library="SL.glm", q_library="SL.glm", truncation=0.05),
  ens_t01    = tmle_candidate("ens_t01", g_library=ENS, q_library=ENS, truncation=0.01),
  ens_t025   = tmle_candidate("ens_t025", g_library=ENS, q_library=ENS, truncation=0.025),
  ens_cf_t01 = tmle_candidate("ens_cf_t01", g_library=ENS, q_library=ENS, truncation=0.01,
                              cv_scheme="cv_tmle", cv_V=5L)
)
d <- generate_data(n = 4000, overlap_strength = 1.5, effect_size = -0.05,
                   misspec = FALSE, seed = 20260608L)
lock <- create_simple_lock(data = d, treatment = "treatment", outcome = "event_24",
                           covariates = COVARS, sl_library = "SL.glm",
                           plasmode_reps = 3L, seed = 20260608L)
threats <- list(
  covariate_missingness  = list(fractions = c(0.10)),
  unmeasured_confounding = list(U_prevalence=0.20, U_treatment_OR=c(2.0), U_outcome_OR=c(2.0))
)
for (nm in names(cands)) {
  say("\n##### candidate:", nm, "#####")
  r <- tryCatch(
    run_plasmode_dq_stress(lock, tmle_candidates = list(cands[[nm]]),
                           effect_sizes = c(0.05), reps = 3L,
                           data_quality_scenarios = threats,
                           q0_library = NULL, verbose = FALSE),
    error = function(e) { say("R-ERROR:", conditionMessage(e)); NULL })
  if (!is.null(r)) say("OK rows:", nrow(r$metrics))
  say("##### survived:", nm, "#####")
}
say("\nALL DONE")
