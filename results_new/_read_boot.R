bv <- read.csv("results_new/bootstrap_variance.csv")
for (sc in c("C: Very Good Overlap", "A: Good Overlap", "B: Marginal Overlap")) {
  sub <- bv[bv$scenario == sc, ]
  cat(sprintf("\n=== %s ===\n", sc))
  print(sub[, c("method","n_reps","se_sd_ratio_if","coverage_if","se_sd_ratio_boot","coverage_boot")])
}
