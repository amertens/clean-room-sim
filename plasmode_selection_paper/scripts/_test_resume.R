# Test resume / checkpoint logic. Save a partial rds that contains
# all 4 workflow rows for (linear, rep=1) only, then call
# run_simulation with linear x 2 reps + nonlinear_smooth x 1 rep.
# Expected: linear rep 1 is SKIPPED (already complete), linear rep 2
# runs, then nonlinear_smooth rep 1 runs. Total new fits: 2 outer
# reps -> 8 new rows. Final file has 12 rows total.

options(warn = 1L)
setwd(here::here("plasmode_selection_paper"))
source("R/run_simulation.R")

partial <- data.frame(
  dgp = "linear",
  rep = 1L,
  workflow = c("fixed_parametric","fixed_rich",
               "plasmode_min_rmse","plasmode_fiord_two_stage"),
  library_label = "fake",
  truncation = 0.01,
  selected_candidate = NA_character_,
  estimate = -0.05, se = 0.04,
  ci_lower = -0.13, ci_upper = 0.03,
  stringsAsFactors = FALSE
)
dir.create("results", showWarnings = FALSE)
saveRDS(partial, "results/sim_resume_test.rds")
cat("[seed] wrote partial rds with 4 rows (linear rep 1 only)\n")

t0 <- Sys.time()
sim <- run_simulation(
  dgps       = c("linear", "nonlinear_smooth"),
  n_per_rep  = 400L,
  n_reps     = 1L,
  inner_reps = 5L,
  out_path   = "results/sim_resume_test.rds",
  resume     = TRUE,
  log_progress = TRUE
)
cat("Resume test finished in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "min\n")
cat("Total rows:", nrow(sim), "\n")
cat("By dgp x rep:\n")
print(table(sim$dgp, sim$rep))
cat("By workflow:\n")
print(table(sim$workflow))
cat("Library labels (linear rep 1 should still say 'fake'):\n")
print(unique(sim[sim$dgp == "linear" & sim$rep == 1, "library_label"]))
