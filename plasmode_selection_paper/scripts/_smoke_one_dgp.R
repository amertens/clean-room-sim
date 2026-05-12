# Single-DGP smoke. Takes DGP name as the first CLI arg.
args <- commandArgs(trailingOnly = TRUE)
dgp <- if (length(args) >= 1L) args[1L] else "nonlinear_smooth"
cat("[smoke_one_dgp] dgp =", dgp, "\n")

setwd(here::here("plasmode_selection_paper"))
source("R/run_simulation.R")
sim <- run_simulation(
  dgps       = dgp,
  n_per_rep  = 400L,
  n_reps     = 1L,
  inner_reps = 5L,
  out_path   = paste0("results/sim_smoke_", dgp, ".rds"),
  log_progress = TRUE
)
cat("DONE.\n")
print(sim)
