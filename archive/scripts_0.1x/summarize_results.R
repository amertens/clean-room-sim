#!/usr/bin/env Rscript
# ── Post-processing: extract and display all simulation results ─────────
# Run after run_simulation.R completes (or on interim results).

cat("Loading simulation results...\n\n")

final <- readRDS("results/simulation_results.rds")

cat(sprintf("Config: N=%d, Reps=%d, SL=%s, Folds=%d\n",
            final$config$n_obs, final$config$n_reps,
            paste(final$config$sl_library, collapse="+"),
            final$config$n_folds))
cat(sprintf("Runtime: %.1f minutes\n\n", final$runtime_minutes))


# ── Plasmode Candidate Selection ─────────────────────────────────────────

cat("================================================================\n")
cat("  PLASMODE CANDIDATE SELECTION\n")
cat("================================================================\n\n")

for (sc_name in names(final$meta)) {
  meta <- final$meta[[sc_name]]
  cat(sprintf("--- %s ---\n", meta$label))
  cat(sprintf("  Selected: %s  |  Gate: %s\n", meta$selected, meta$gate))
  cat(sprintf("  PS range: [%.3f, %.3f]  |  ESS%%: %.1f\n\n",
              meta$ps_range[1], meta$ps_range[2], meta$ess_pct))
  print(meta$plasmode, row.names = FALSE)
  cat("\n\n")
}


# ── Monte Carlo Results ──────────────────────────────────────────────────

cat("================================================================\n")
cat("  MONTE CARLO ESTIMATOR COMPARISON\n")
cat("================================================================\n\n")

for (sc_name in names(final$summaries)) {
  truth <- final$truths[[sc_name]]$RD
  cat(sprintf("--- %s (true RD = %.5f) ---\n",
              final$scenarios[[sc_name]]$label, truth))
  print(final$summaries[[sc_name]], row.names = FALSE)
  cat("\n")
}


# ── Plasmode vs Monte Carlo comparison ───────────────────────────────────

cat("================================================================\n")
cat("  PLASMODE vs MONTE CARLO: Does plasmode predict MC performance?\n")
cat("================================================================\n\n")

for (sc_name in names(final$meta)) {
  meta    <- final$meta[[sc_name]]
  mc_summ <- final$summaries[[sc_name]]

  cat(sprintf("--- %s ---\n", meta$label))
  cat(sprintf("  Plasmode selected: %s\n\n", meta$selected))

  # The selected candidate's plasmode metrics
  plas <- meta$plasmode
  selected_plas <- plas[plas$candidate == meta$selected, ]
  if (nrow(selected_plas) > 0) {
    cat("  Plasmode (pre-outcome, averaged over effect sizes):\n")
    cat(sprintf("    Bias:     %.5f\n", mean(selected_plas$bias)))
    cat(sprintf("    RMSE:     %.5f\n", mean(selected_plas$rmse)))
    cat(sprintf("    Coverage: %.3f\n", mean(selected_plas$coverage)))
    cat(sprintf("    SE/SD:    %.3f\n", mean(selected_plas$se_cal, na.rm = TRUE)))
  }

  # The TMLE MC metrics
  tmle_mc <- mc_summ[mc_summ$method == "TMLE", ]
  if (nrow(tmle_mc) > 0) {
    cat("\n  Monte Carlo (real outcome, 500 reps):\n")
    cat(sprintf("    Bias:     %.5f\n", tmle_mc$bias))
    cat(sprintf("    RMSE:     %.5f\n", tmle_mc$rmse))
    cat(sprintf("    Coverage: %.3f\n", tmle_mc$coverage))
    cat(sprintf("    SE/SD:    %.3f\n", tmle_mc$se_sd_ratio))
  }
  cat("\n\n")
}

cat("Done.\n")
