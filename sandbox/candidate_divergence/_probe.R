cat(R.version.string, "\n")
for (p in c("devtools", "pkgload", "SuperLearner", "callr", "ggplot2")) {
  cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly = TRUE)))
}
