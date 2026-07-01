bv <- read.csv("results_new/bootstrap_variance.csv")
r <- bv[bv$scenario == "C: Very Good Overlap" & bv$method == "Match_TMLE", ]
cat("coverage_if:", r$coverage_if, " n:", r$n_reps, "\n")
cat("MCSE:", sqrt(r$coverage_if * (1 - r$coverage_if) / r$n_reps), "\n")
cat("boot coverage:", r$coverage_boot, "\n")
