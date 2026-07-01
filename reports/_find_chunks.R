lines <- readLines("manuscript_methods.qmd")
# Find chunks that (a) contain the identifying content and (b) DON'T have a #| label line
targets <- list(
  mc_metrics   = "knitr::kable.*summary_df|knitr::kable.*mc.*metric|knitr::kable.*s,.*digits",
  mc_boxplots  = "geom_boxplot.*alpha.*0\\.7",
  coverage_bar = "coverage.*bar|geom_bar.*coverage|Sampling.*coverage|95.*CI.*coverage.*bar",
  se_calib     = "se_sd_ratio|SE.*calibration|Ratio.*empirical.*SD",
  mc_forest    = "Forest plot.*MC|Monte Carlo.*forest|geom_errorbarh|mc.*forest",
  dq_degrad    = "dq.*gradient|degradation.*gradient|Degradation gradient"
)
for (nm in names(targets)) {
  idxs <- grep(targets[[nm]], lines, perl=TRUE, ignore.case=TRUE)
  if (length(idxs) == 0) { cat(sprintf("  %-14s NOT FOUND\n", nm)); next }
  # Find the chunk start for each hit
  for (i in idxs) {
    # Walk backward to the chunk open line
    for (j in seq(i, max(1, i-30), -1)) {
      if (grepl("^```\\{r", lines[j])) {
        # Check if there's already a #| label
        has_label <- any(grepl("^#\\| label:", lines[seq(j, min(length(lines), j+10))]))
        cat(sprintf("  %-14s  L%d: %s  [has_label=%s]\n",
                    nm, j, trimws(lines[j]), has_label))
        break
      }
    }
  }
}
