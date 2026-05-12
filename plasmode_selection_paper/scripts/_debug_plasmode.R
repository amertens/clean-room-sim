setwd(here::here("plasmode_selection_paper"))
source("R/dgps.R")
source("R/candidates.R")
source("R/workflows.R")

covariates <- c("age", "sex", "biomarker", "comorbidity", "bmi_z", "smoke")
candidates <- build_candidate_grid()

# Reproduce rep 2 of the smoke test
sd <- 2026L + 1000L * 1L + 2L
set.seed(sd)
d <- make_data(400L, dgp = "linear", true_RD = -0.05)

cat("=== Try fit_plasmode WITHOUT screener candidates ===\n")
candidates_no_screener <- candidates[c("param_t01","paramplus_t01","rich_t01","rich_t05")]
out <- tryCatch(
  fit_plasmode(d, covariates, candidates_no_screener,
               inner_reps = 5L, inner_effect = -0.05,
               rule = "min_rmse", seed = sd),
  error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
cat("OUT (no screener):\n"); print(out)

cat("\n=== Run run_plasmode_feasibility directly with full grid ===\n")
library(cleanTMLE)
lock <- create_simple_lock(d, "treatment", "event_24", covariates,
                            sl_library = c("SL.glm","SL.glmnet","SL.gam","SL.ranger","SL.mean"),
                            seed = sd)
plas <- tryCatch(
  run_plasmode_feasibility(lock, tmle_candidates = candidates,
                            effect_sizes = -0.05, reps = 3L),
  error = function(e) { cat("ERROR run_plasmode:", conditionMessage(e), "\n"); NULL }
)
cat("metrics rows:", if (is.null(plas$metrics)) 0 else nrow(plas$metrics), "\n")
if (!is.null(plas$metrics)) print(head(plas$metrics, 12))
