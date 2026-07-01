lines <- readLines("manuscript_methods.qmd")
# em-dash
hits_em <- grep("—", lines)
# contractions
hits_contr <- grep("n’t|n't", lines)
# "not X but" pattern in prose (not in code)
hits_notbut <- grep("\\bnot\\b.*\\bbut\\b", lines, perl=TRUE)

for (h in c(hits_em, hits_contr)) {
  cat(sprintf("L%d: %s\n", h, trimws(lines[h])))
}
cat("\n--- not...but ---\n")
for (h in hits_notbut[1:min(10, length(hits_notbut))]) {
  cat(sprintf("L%d: %s\n", h, trimws(lines[h])))
}
