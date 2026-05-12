# Five-DGP smoke test. Two reps per DGP, n=400, 5 inner reps.
# Goal: exercise every DGP at least once, surface DGP-specific bugs,
# and check whether plasmode candidate selection actually varies
# across DGPs (the central question of the paper).
options(warn = 1L)  # print warnings as they occur, not at end
setwd(here::here("plasmode_selection_paper"))
source("R/run_simulation.R")

t0 <- Sys.time()
sim <- run_simulation(
  dgps       = c("linear", "nonlinear_smooth", "interactions",
                 "sparse", "high_dim_noise"),
  n_per_rep  = 400L,
  n_reps     = 2L,
  inner_reps = 5L,
  out_path   = "results/sim_smoke_all.rds",
  log_progress = TRUE
)
cat("Smoke-all finished in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "min\n")
cat("Rows:", nrow(sim), "\n\n")

cat("=== Estimates by DGP x workflow x rep ===\n")
print(sim[, c("dgp","rep","workflow","selected_candidate",
              "estimate","se")])

cat("\n=== Plasmode candidate-selection cross-tab ===\n")
psel <- sim[grepl("^plasmode", sim$workflow), ]
print(table(psel$dgp, psel$selected_candidate, psel$workflow))

cat("\n=== Mean estimate by DGP x workflow (true RD = -0.05) ===\n")
agg <- aggregate(estimate ~ dgp + workflow, data = sim, FUN = mean)
print(agg)

cat("\n=== Any NA estimates? ===\n")
nas <- sim[is.na(sim$estimate), ]
if (nrow(nas) == 0L) cat("None.\n") else print(nas)
