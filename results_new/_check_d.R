res <- readRDS("results_new/simulation_results.rds")
# Scenario D summary
cat("=== Scenario D summary ===\n")
d <- res$summaries$misspecified
print(d[, c("method","bias","coverage","se_sd_ratio","rmse")], row.names=FALSE)
cat("\nTruth:", res$truths$misspecified$RD, "\n")

# Check meta structure  
cat("\n=== Meta names ===\n")
cat(paste(names(res$meta), collapse=", "), "\n")
if (!is.null(res$meta$misspecified)) {
  cat("misspecified meta names:", paste(names(res$meta$misspecified), collapse=", "), "\n")
}
