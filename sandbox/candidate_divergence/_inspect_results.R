library(pkgload)
pkgload::load_all("cleanTMLE", quiet = TRUE)

out <- readRDS("sandbox/candidate_divergence/results/candidate_divergence_full.rds")

cat("=== Top-level names ===\n"); print(names(out))

cat("\n=== agg table ===\n"); print(out$agg)

cat("\n=== win_tbl (per-batch winners) ===\n"); print(out$win_tbl)

cat("\n=== dq_metrics columns ===\n"); print(names(out$dq_metrics))

cat("\n=== dq_metrics scenarios/levels ===\n")
print(unique(out$dq_metrics[, c("scenario","level")]))

cat("\n=== plas_metrics ===\n"); print(out$plas_metrics)

cat("\n=== dq_cell_agg (pooled per-cell) ===\n")
print(head(out$dq_cell_agg, 20))

cat("\n=== separation/mcse ===\n")
cat("separation:", out$separation, " mcse:", out$combined_mcse, "\n")
cat("win_min:", out$win_min_final, " win_mm:", out$win_mm_final, "\n")

cat("\n=== degradation (summarize_dq_degradation output) ===\n")
print(head(out$degradation, 15))
