# Pilot run: small budget to verify the pipeline works end-to-end
# and to produce a placeholder results file for the manuscript draft.
#
# Run time on a typical workstation: ~20-30 minutes.

setwd(here::here("plasmode_selection_paper"))
source("R/run_simulation.R")

t0 <- Sys.time()
sim <- run_simulation(
  dgps       = c("linear", "nonlinear_smooth", "sparse"),
  n_per_rep  = 600L,
  n_reps     = 40L,        # tiny budget for the pilot
  inner_reps = 20L,
  out_path   = "results/sim_pilot.rds",
  log_progress = TRUE
)
cat("Pilot finished in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "minutes\n")
cat("Rows: ", nrow(sim), "\n")
cat("Workflows: ", paste(unique(sim$workflow), collapse = ", "), "\n")
