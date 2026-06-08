deg <- read.csv("rescueCo/results/plasmode_dq_degradation.csv")
# Baseline (clean-data) per-candidate metrics
b <- deg[deg$scenario == "cov_miss" & deg$level == 0.1, ]  # any cell; baseline is constant
b <- unique(b[, c("candidate", "bias_baseline", "rmse_baseline", "cov_baseline")])
b <- b[order(b$rmse_baseline), ]
cat("Baseline (clean-data) per-candidate metrics:\n")
print(b, row.names = FALSE, digits = 4)
