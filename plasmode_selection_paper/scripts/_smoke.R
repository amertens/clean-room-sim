# Minimum-budget smoke test for the simulation pipeline.
# Should run in ~2 minutes.
setwd(here::here("plasmode_selection_paper"))
source("R/run_simulation.R")
t0 <- Sys.time()
sim <- run_simulation(
  dgps       = c("linear"),
  n_per_rep  = 400L,
  n_reps     = 3L,
  inner_reps = 5L,
  out_path   = "results/sim_smoke.rds",
  log_progress = TRUE
)
cat("Smoke test finished in",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 0),
    "seconds\n")
cat("Rows: ", nrow(sim), "\n")
cat("Workflows: ", paste(unique(sim$workflow), collapse = ", "), "\n")
print(sim)
