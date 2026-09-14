#!/usr/bin/env Rscript
# ── Quick check on simulation progress ──────────────────────────────────
# Run this while run_simulation.R is executing to see interim results.

cat("Checking simulation progress...\n\n")

for (sc in c("good_overlap", "marginal_overlap", "unmeasured_conf")) {
  path <- file.path("results", sprintf("interim_%s.rds", sc))
  if (!file.exists(path)) {
    cat(sprintf("  %s: not started yet\n", sc))
    next
  }

  interim <- readRDS(path)
  cat(sprintf("=== %s ===\n", interim$meta$label))
  cat(sprintf("  Completed: %d / %d reps  (saved at %s)\n",
              interim$n_done, interim$n_total,
              format(interim$saved_at, "%H:%M:%S")))
  cat(sprintf("  True RD: %.5f\n", interim$meta$truth_rd))
  cat(sprintf("  Selected candidate: %s  |  Gate: %s\n\n",
              interim$meta$selected, interim$meta$gate))

  if (nrow(interim$results) > 0) {
    df <- interim$results
    truth <- interim$meta$truth_rd
    methods <- unique(df$method)
    rows <- lapply(methods, function(m) {
      sub   <- df[df$method == m, ]
      valid <- !is.na(sub$estimate)
      ests  <- sub$estimate[valid]
      ses   <- sub$se[valid]
      n_ok  <- sum(valid)
      if (n_ok == 0) return(NULL)
      emp_sd  <- sd(ests)
      mean_se <- mean(ses, na.rm = TRUE)
      data.frame(
        method     = m,
        n_ok       = n_ok,
        bias       = round(mean(ests) - truth, 5),
        rmse       = round(sqrt(mean((ests - truth)^2)), 5),
        coverage   = round(mean(sub$covers[valid], na.rm = TRUE), 3),
        se_sd      = round(if (emp_sd > 0) mean_se / emp_sd else NA, 3),
        stringsAsFactors = FALSE
      )
    })
    summary_df <- do.call(rbind, Filter(Negate(is.null), rows))
    print(summary_df, row.names = FALSE)
  }
  cat("\n")
}
