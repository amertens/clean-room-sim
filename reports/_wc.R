lines <- readLines("manuscript_methods.qmd")
in_yaml <- FALSE; in_chunk <- FALSE
prose <- c()
for (l in lines) {
  if (grepl("^---$", l)) { in_yaml <- !in_yaml; next }
  if (in_yaml) next
  if (grepl("^```", l)) { in_chunk <- !in_chunk; next }
  if (in_chunk) next
  if (grepl("^#\\|", l)) next
  prose <- c(prose, l)
}
txt <- paste(prose, collapse = " ")
# simple word count: split on whitespace
words <- sum(nchar(strsplit(txt, "\\s+")[[1]]) > 0)
cat("Body prose word count:", words, "\n")
cat("20% target cut:       ", round(words * 0.20), "\n")
