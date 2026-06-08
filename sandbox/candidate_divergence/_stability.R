# Reviewer item #3: selection-rule stability from the divergence study's batches
x <- readRDS("sandbox/candidate_divergence/results/candidate_divergence_full.rds")

cat("=== Per-batch winners (5 batches x 200 reps each) ===\n")
print(x$win_tbl, row.names = FALSE)

cat("\nmin_rmse selection frequency:\n");      print(table(x$win_tbl$win_min_rmse))
cat("min_max_rmse selection frequency:\n");    print(table(x$win_tbl$win_min_max_rmse))

cat("\n=== Worst-case RMSE per candidate x batch ===\n")
wm <- x$worst_mat
print(round(wm, 4))

gap <- apply(wm, 2, function(col) { s <- sort(col); s[2] - s[1] })
cat("\nminimax margin (2nd-best minus best worst-RMSE) per batch:\n")
print(round(gap, 4))
cat(sprintf("mean margin = %.4f   combined MC SE = %.4f   margin/SE = %.1f\n",
            mean(gap), x$combined_mcse, mean(gap) / x$combined_mcse))
cat("\nSTABILITY_DONE\n")
