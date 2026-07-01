# Minimal repro to surface the run_cell error with an explicit message.
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
d <- generate_data(n = 4000, overlap_strength = 1.5, effect_size = -0.05,
                   misspec = FALSE, seed = 20260608L)
say("data ok; dim:", paste(dim(d), collapse="x"), "; Yrate:", round(mean(d$event_24),3))
lock <- create_simple_lock(data = d, treatment = "treatment", outcome = "event_24",
                           covariates = COVARS, sl_library = "SL.glm",
                           plasmode_reps = 1L, seed = 20260608L)
say("lock ok; hash:", substr(lock$lock_hash,1,12))

cands <- list(
  glm_t01    = tmle_candidate("glm_t01", g_library="SL.glm", q_library="SL.glm", truncation=0.01),
  ens_t01    = tmle_candidate("ens_t01", g_library=ENS, q_library=ENS, truncation=0.01),
  ens_cf_t01 = tmle_candidate("ens_cf_t01", g_library=ENS, q_library=ENS, truncation=0.01,
                              cv_scheme="cv_tmle", cv_V=5L)
)

for (nm in names(cands)) {
  say("\n--- feasibility, candidate:", nm, "---")
  r <- tryCatch(
    run_plasmode_feasibility(lock, tmle_candidates = list(cands[[nm]]),
                             effect_sizes = c(0.05), reps = 1L,
                             q0_library = NULL, verbose = FALSE),
    error = function(e) { say("ERROR:", conditionMessage(e)); say("CALL:", deparse(conditionCall(e))); NULL })
  if (!is.null(r)) { say("OK metrics:"); print(r$metrics) }
}
say("\nDONE")
