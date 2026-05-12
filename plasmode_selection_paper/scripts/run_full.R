# Full study run. Designed for an inferentially-defensible budget.
#
# Run time on a typical workstation: many hours to ~1 day.
# Recommended: run on a cluster or overnight.
#
# The budget below is the conservative paper-quality default. Adjust
# n_reps downward to 200 if the wall-clock budget is tight; coverage
# Monte Carlo standard error at n_reps = 200 is ~0.015.

setwd(here::here("plasmode_selection_paper"))
source("R/run_simulation.R")

t0 <- Sys.time()
sim <- run_simulation(
  dgps       = c("linear", "nonlinear_smooth", "interactions",
                 "sparse", "high_dim_noise"),
  n_per_rep  = 1000L,
  n_reps     = 500L,
  inner_reps = 50L,
  noise_p    = 20L,
  out_path   = "results/sim_full.rds",
  log_progress = TRUE
)
cat("Full study finished in",
    round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 2),
    "hours\n")
cat("Rows: ", nrow(sim), "\n")
