res <- readRDS("results_new/simulation_results.rds")
cat("Names:", paste(names(res), collapse=", "), "\n")
cat("Scenario names:", paste(names(res$scenarios), collapse=", "), "\n")
for (sc in names(res$scenarios)) {
  cat(sprintf("\n--- %s ---\n", sc))
  cat("Label:", res$scenarios[[sc]]$label, "\n")
  cat("Truth:", res$truths[[sc]]$RD, "\n")
  cat("Summary cols:", paste(names(res$summaries[[sc]]), collapse=", "), "\n")
  if (!is.null(res$results[[sc]])) {
    cat("Results class:", class(res$results[[sc]]), "\n")
    cat("Results cols:", paste(colnames(res$results[[sc]]), collapse=", "), "\n")
    cat("nrow:", nrow(res$results[[sc]]), "\n")
  }
}
