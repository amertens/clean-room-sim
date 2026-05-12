options(error = function() {
  cat("=== TRACEBACK ===\n")
  traceback(max.lines = 80)
  cat("=== last.error ===\n")
  if (!is.null(e <- get0(".Last.error"))) print(e)
  quit(save = "no", status = 2L)
})
setwd(here::here("plasmode_selection_paper"))
cat("[step] sourcing run_simulation.R\n")
source("R/run_simulation.R")
cat("[step] calling run_simulation\n")
sim <- run_simulation(
  dgps       = c("linear"),
  n_per_rep  = 400L,
  n_reps     = 1L,
  inner_reps = 5L,
  out_path   = "results/sim_smoke.rds",
  log_progress = TRUE
)
cat("[step] done\n")
print(sim)
