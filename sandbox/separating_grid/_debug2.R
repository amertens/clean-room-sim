# Isolate which DQ threat breaks run_plasmode_dq_stress on Arm A.
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
GRID <- list(
  tmle_candidate("glm_t01", g_library="SL.glm", q_library="SL.glm", truncation=0.01),
  tmle_candidate("ens_t01", g_library=ENS, q_library=ENS, truncation=0.01)
)
d <- generate_data(n = 4000, overlap_strength = 1.5, effect_size = -0.05,
                   misspec = FALSE, seed = 20260608L)
lock <- create_simple_lock(data = d, treatment = "treatment", outcome = "event_24",
                           covariates = COVARS, sl_library = "SL.glm",
                           plasmode_reps = 1L, seed = 20260608L)
full <- default_dq_scenarios("regulatory_standard")
full$near_positivity <- NULL   # Arm A keeps no near-positivity slope
say("threat names:", paste(names(full), collapse=", "))

for (thr in names(full)) {
  one <- full[thr]
  say("\n--- threat:", thr, "| spec:", paste(deparse(one[[1]]), collapse=" "))
  r <- tryCatch(
    run_plasmode_dq_stress(lock, tmle_candidates = GRID, effect_sizes = c(0.05),
                           reps = 1L, data_quality_scenarios = one,
                           q0_library = NULL, verbose = FALSE),
    error = function(e) { say("ERROR:", conditionMessage(e));
                          say("CALL:", paste(deparse(conditionCall(e)), collapse=" ")); NULL })
  if (!is.null(r)) { say("OK rows:", nrow(r$metrics), "scenarios:",
                         paste(unique(r$metrics$scenario), collapse=",")) }
}
say("\nDONE")
