lines <- readLines("manuscript_methods.qmd")
# 1. regime-distinct sentence
idx <- grep("distinct from the main simulation scenarios", lines)
cat("Regime-distinct sentence:\n")
for (i in idx) cat(sprintf("  L%d: %s\n", i, trimws(lines[i])))

# 2. B=200 and B=500 in bootstrap table caption area
idx2 <- grep("B = 200|B = 500", lines)
cat("\nB=200/B=500 lines:\n")
for (i in idx2) cat(sprintf("  L%d: %s\n", i, trimws(lines[i])))

# 3. replicate count consistency
idx3 <- grep("100-replicate|100 replicates|100 MC rep", lines)
cat("\n100-replicate lines:\n")
for (i in idx3) cat(sprintf("  L%d: %s\n", i, trimws(lines[i])))

# 4. n_sims = 3 used only as smoke test?
idx4 <- grep("n_sims = 3|3 replicates per", lines)
cat("\nn_sims=3 lines:\n")
for (i in idx4) cat(sprintf("  L%d: %s\n", i, trimws(lines[i])))
