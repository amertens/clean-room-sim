deg <- read.csv("rescueCo/results/plasmode_dq_degradation.csv")
cat("Scenarios in current case-study DQ data:\n")
print(unique(deg$scenario))
cat("\nCandidates:\n"); print(unique(deg$candidate))
# Worst-case (max RMSE_ratio) per candidate across all DQ scenarios
w <- aggregate(rmse_ratio ~ candidate, data = deg, FUN = max)
w <- w[order(w$rmse_ratio), ]
cat("\nPer-candidate worst-case RMSE ratio (min_max_rmse pick is row 1):\n")
print(w, row.names = FALSE, digits = 4)
cat(sprintf("\nminimax pick = %s\n", w$candidate[1]))
# Baseline RMSE (use baseline column, one row per candidate is fine)
b <- aggregate(rmse_baseline ~ candidate, data = deg, FUN = mean)
b <- b[order(b$rmse_baseline), ]
cat("\nPer-candidate baseline RMSE (min_rmse pick is row 1):\n")
print(b, row.names = FALSE, digits = 4)
cat(sprintf("min_rmse pick = %s\n", b$candidate[1]))
