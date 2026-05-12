# Pilot run: small budget to verify the pipeline works end-to-end
# and to produce a placeholder results file for the manuscript draft.
#
# Realistic run-time on a single workstation, per the smoke benchmarks:
# one outer rep at n=400 with 5 inner reps takes ~4-5 minutes, so this
# pilot (5 DGPs x 10 outer x 10 inner at n=500) is roughly 5-7 hours.
# Use scripts/_smoke.R for a quick (~15 min) one-DGP sanity check;
# scripts/run_full.R for the paper-quality 5 DGPs x 500 x 50 budget.

setwd(here::here("plasmode_selection_paper"))
source("R/run_simulation.R")

t0 <- Sys.time()
sim <- run_simulation(
  dgps       = c("linear", "nonlinear_smooth", "interactions",
                 "sparse", "high_dim_noise"),
  n_per_rep  = 500L,
  n_reps     = 10L,        # pilot budget; MC SE of coverage ~0.07
  inner_reps = 10L,
  out_path   = "results/sim_pilot.rds",
  log_progress = TRUE
)
cat("Pilot finished in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "minutes\n")
cat("Rows: ", nrow(sim), "\n")
cat("Workflows: ", paste(unique(sim$workflow), collapse = ", "), "\n")
