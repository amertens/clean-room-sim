# Full study run. Designed for an inferentially-defensible budget.
#
# Realistic wall-clock at this budget on a single workstation, from
# the smoke benchmarks: ~5000 CPU-hours (5 DGPs x 500 reps x ~2 hours
# per rep at n=1000 / 50 inner). This MUST be run on a multi-core
# cluster or with parallelism added to run_simulation; do not start
# it on a laptop. A worthwhile interim target is n_reps = 100 and
# inner_reps = 25 (~125 CPU-hours), which still gives coverage MC SE
# around 0.02.

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
