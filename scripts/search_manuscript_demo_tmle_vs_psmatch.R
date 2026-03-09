#!/usr/bin/env Rscript
# ===========================================================================
# search_manuscript_demo_tmle_vs_psmatch.R
# ===========================================================================
#
# PURPOSE:
#   Fast grid search to find a single demonstrative simulation setup for the
#   manuscript. The goal is a realistic DGP where:
#
#     1. A traditional PS-diagnostics-gated workflow would fail or stop
#        (poor overlap, low ESS, poor matching feasibility).
#     2. A clean-room TMLE workflow, validated by outcome-blind simulation,
#        would proceed (low bias, good coverage, calibrated SE).
#     3. TMLE outperforms PS-matched regression in the estimator comparison.
#
#   This is NOT a full simulation study. It identifies candidate setups for
#   a worked manuscript example.
#
# TWO WORKFLOWS COMPARED:
#   Workflow 1 (Traditional PS-gated):
#     - PS via logistic regression (main effects)
#     - 1:1 nearest-neighbor matching, caliper = 0.2 * SD(logit PS)
#     - GO/STOP based on overlap diagnostics and matching feasibility
#
#   Workflow 2 (Clean-room TMLE):
#     - g via cv.glmnet, Q via cv.glmnet with rich feature matrix
#     - GO/STOP based on Stage-3 outcome-blind simulation metrics
#     - Can proceed even when PS diagnostics are poor
#
# STRATEGY:
#   Two-phase search to keep runtime manageable:
#   Phase 1 (Screen):  Few reps at one sample size to identify promising DGPs
#   Phase 2 (Refine):  More reps at multiple sample sizes for top candidates
#
# USAGE:
#   Rscript scripts/search_manuscript_demo_tmle_vs_psmatch.R [fast|refine]
#
# OUTPUT:
#   outputs/manuscript_demo_search/grid_results.csv
#   outputs/manuscript_demo_search/top_demo_candidates.csv
#   outputs/manuscript_demo_search/top_candidate_replicates.rds
#   outputs/manuscript_demo_search/search_summary.txt
# ===========================================================================

suppressPackageStartupMessages({
  library(glmnet)
})

# Source helpers
source(file.path("scripts", "helpers_manuscript_demo_search.R"))

# ===========================================================================
# Configuration
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1 && args[1] == "refine") "refine" else "fast"

FAST_MODE          <- (mode == "fast")
USE_SUPERLEARNER   <- FALSE
RERUN_TOP_WITH_SL  <- TRUE

# Phase 1: Screen — cheap, one sample size, few reps
SCREEN_N       <- 1500
SCREEN_REPS    <- 15
# Phase 2: Refine — more reps, multiple sample sizes
if (FAST_MODE) {
  REFINE_N_VALUES <- c(1000, 1500, 2000)
  REFINE_REPS     <- 40
  TRUTH_N         <- 50000
  TOP_K_REFINE    <- 50   # refine top 50 from Phase 1
  message("=== FAST MODE ===")
} else {
  REFINE_N_VALUES <- c(2000, 3000, 4500)
  REFINE_REPS     <- 100
  TRUTH_N         <- 100000
  TOP_K_REFINE    <- 30
  message("=== REFINE MODE ===")
}

PS_CLIP <- c(0.025, 0.975)

# Output directory
out_dir <- file.path("outputs", "manuscript_demo_search")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# Parameter grid
# ===========================================================================
# Design rationale for the grid:
#   - b_rr drives overlap stress (the PS difficulty)
#   - beta_rr1 controls whether renal_risk matters for the outcome
#     (low rr1 = TMLE can compensate; high rr1 = both workflows struggle)
#   - beta_age3, beta_ckdA, beta_ageA, beta_age2A control outcome complexity
#     (the source of TMLE's advantage over PS matching)
#   - beta_rrA controls treatment-by-renal_risk interaction
#
# The grid values below are from the user's specification. We cover the
# full requested ranges but use strategic pruning + two-phase search to
# keep runtime feasible.
# ===========================================================================

# ===========================================================================
# The grid uses pre-specified heterogeneity profiles to keep the search
# manageable while covering the requested parameter ranges.
#
# Fixed coefficients (same across all candidates):
#   beta_A = 0.22, beta_ckd = 0.95, beta_cirr = 0.55,
#   beta_dm = 0.40, beta_nsaid = 0.25
#   beta_age1 = 0.25, beta_age2 = -0.20
#     (these have minimal impact on the TMLE-vs-PS contrast)
#
# Grid dimensions that matter most:
#   b_rr          — overlap stress (PS difficulty)
#   beta_rr1      — renal_risk in outcome (low = TMLE can compensate)
#   heterogeneity — age nonlinearity + treatment-effect modifiers
#   target_y_prev — event rate
# ===========================================================================

# Treatment model: overlap stress levels
b_rr_values <- c(2.2, 2.6, 3.0)

# Renal_risk in outcome: low means TMLE has room to compensate
beta_rr1_values <- c(0.00, 0.03, 0.10, 0.30)
beta_rr2_values <- c(0.05)   # fix quadratic (secondary)
beta_rrA_values <- c(0.00, 0.20)

# Heterogeneity profiles: pre-defined bundles of (age3, ckdA, ageA, age2A)
# Each profile represents a different level/type of outcome complexity
hetero_profiles <- data.frame(
  beta_age3  = c(0.10, 0.15, 0.10, 0.15, 0.10, 0.00, 0.15),
  beta_ckdA  = c(0.30, 0.20, 0.00, 0.30, 0.20, 0.30, 0.30),
  beta_ageA  = c(0.18, 0.10, 0.18, 0.18, 0.00, 0.18, 0.10),
  beta_age2A = c(0.10, 0.05, 0.10, 0.00, 0.10, 0.10, 0.10),
  stringsAsFactors = FALSE
)

# Event rate
target_y_prev_values <- c(0.10, 0.14, 0.18)

# Build grid from components
dgp_grid <- expand.grid(
  b_rr       = b_rr_values,
  beta_rr1   = beta_rr1_values,
  beta_rr2   = beta_rr2_values,
  beta_rrA   = beta_rrA_values,
  hetero_id  = seq_len(nrow(hetero_profiles)),
  target_y_prev = target_y_prev_values,
  stringsAsFactors = FALSE
)

# Merge heterogeneity profiles
dgp_grid <- cbind(dgp_grid, hetero_profiles[dgp_grid$hetero_id, ])
dgp_grid$hetero_id <- NULL

# Add fixed parameters
dgp_grid$beta_age1 <- 0.25
dgp_grid$beta_age2 <- -0.20

total_dgps <- nrow(dgp_grid)
message("Total DGP configurations: ", total_dgps)

# ---------------------------------------------------------------------------
# Pruning
# ---------------------------------------------------------------------------

prune_idx <- with(dgp_grid, {
  # Strong renal_risk in outcome + overlap stress → both workflows fail
  drop_rr_strong <- (beta_rr1 >= 0.30) & (b_rr >= 2.6)
  # renal_risk × treatment when rr1 already large → exacerbates both
  drop_rrA <- (beta_rrA > 0) & (beta_rr1 >= 0.10)
  drop_rr_strong | drop_rrA
})

dgp_grid <- dgp_grid[!prune_idx, ]
dgp_grid <- unique(dgp_grid)
message("After pruning: ", nrow(dgp_grid), " DGPs ",
        "(removed ", total_dgps - nrow(dgp_grid), ")")

# ===========================================================================
# Pre-compute calibrated intercepts and truth
# ===========================================================================

message("\nCalibrating treatment model intercepts...")
alpha_A_cache <- list()
for (b in unique(dgp_grid$b_rr)) {
  key <- as.character(b)
  alpha_A_cache[[key]] <- calibrate_alpha_A(b_rr = b, target_prev = 0.49)
  message("  b_rr=", b, " -> alpha_A=", round(alpha_A_cache[[key]], 4))
}

message("Calibrating outcome model intercepts...")
alpha_Y_cache <- list()
y_param_cols <- c("beta_age1", "beta_age2", "beta_age3",
                   "beta_rr1", "beta_rr2",
                   "beta_ckdA", "beta_ageA", "beta_age2A",
                   "beta_rrA", "target_y_prev")
y_combos <- unique(dgp_grid[, y_param_cols])

for (i in seq_len(nrow(y_combos))) {
  r <- y_combos[i, ]
  key <- paste(r, collapse = "_")
  alpha_Y_cache[[key]] <- calibrate_alpha_Y(
    target_prev = r$target_y_prev,
    beta_age1 = r$beta_age1, beta_age2 = r$beta_age2, beta_age3 = r$beta_age3,
    beta_rr1 = r$beta_rr1, beta_rr2 = r$beta_rr2,
    beta_ckdA = r$beta_ckdA, beta_ageA = r$beta_ageA,
    beta_age2A = r$beta_age2A, beta_rrA = r$beta_rrA
  )
}
message("  Calibrated ", nrow(y_combos), " unique outcome parameter sets.")

message("Computing ground truth (Monte Carlo, n_mc=", TRUTH_N, ")...")
truth_cache <- list()
for (i in seq_len(nrow(y_combos))) {
  r <- y_combos[i, ]
  key <- paste(r, collapse = "_")
  aY <- alpha_Y_cache[[key]]
  truth_cache[[key]] <- compute_truth(
    n_mc = TRUTH_N, seed = 999, alpha_Y = aY,
    beta_age1 = r$beta_age1, beta_age2 = r$beta_age2, beta_age3 = r$beta_age3,
    beta_rr1 = r$beta_rr1, beta_rr2 = r$beta_rr2,
    beta_ckdA = r$beta_ckdA, beta_ageA = r$beta_ageA,
    beta_age2A = r$beta_age2A, beta_rrA = r$beta_rrA
  )
}
message("  Computed ", nrow(y_combos), " unique truth values.\n")

# ===========================================================================
# Helper: evaluate one DGP configuration at given n and reps
# ===========================================================================

evaluate_dgp <- function(dgp_row, n, reps, seed_base, alpha_A_cache,
                          alpha_Y_cache, truth_cache, ps_clip) {
  alpha_A_val <- alpha_A_cache[[as.character(dgp_row$b_rr)]]
  y_key <- paste(dgp_row[, y_param_cols], collapse = "_")
  alpha_Y_val <- alpha_Y_cache[[y_key]]
  true_rd     <- truth_cache[[y_key]]["RD"]

  rep_results <- vector("list", reps)
  for (r in seq_len(reps)) {
    d <- generate_data(
      n = n, seed = seed_base + r,
      alpha_A = alpha_A_val, b_rr = dgp_row$b_rr,
      alpha_Y = alpha_Y_val,
      beta_age1 = dgp_row$beta_age1, beta_age2 = dgp_row$beta_age2,
      beta_age3 = dgp_row$beta_age3,
      beta_rr1 = dgp_row$beta_rr1, beta_rr2 = dgp_row$beta_rr2,
      beta_ckdA = dgp_row$beta_ckdA, beta_ageA = dgp_row$beta_ageA,
      beta_age2A = dgp_row$beta_age2A, beta_rrA = dgp_row$beta_rrA
    )
    rep_results[[r]] <- run_one_replicate(d, true_rd, ps_clip)
  }

  summary_row <- summarize_candidate(rep_results, true_rd)
  cbind(dgp_row, n = n, summary_row)
}

# ===========================================================================
# PHASE 1: Screen all DGPs cheaply
# ===========================================================================

n_dgps <- nrow(dgp_grid)
message("=== PHASE 1: Screening ", n_dgps, " DGPs ===")
message("  n=", SCREEN_N, ", reps=", SCREEN_REPS, "\n")

screen_results <- vector("list", n_dgps)
t_start <- proc.time()

for (idx in seq_len(n_dgps)) {
  screen_results[[idx]] <- evaluate_dgp(
    dgp_grid[idx, ], n = SCREEN_N, reps = SCREEN_REPS,
    seed_base = idx * 10000,
    alpha_A_cache = alpha_A_cache,
    alpha_Y_cache = alpha_Y_cache,
    truth_cache = truth_cache,
    ps_clip = PS_CLIP
  )

  if (idx %% 100 == 0 || idx == n_dgps) {
    elapsed <- (proc.time() - t_start)["elapsed"]
    rate <- idx / elapsed
    eta <- (n_dgps - idx) / rate / 60
    n_hits <- sum(sapply(screen_results[seq_len(idx)], function(x) {
      !is.null(x) && x$trad_fail_flag && x$score >= 4
    }))
    message(sprintf("  [%d/%d] %.0f%% | %.1f/s | ETA: %.1f min | promising: %d",
                    idx, n_dgps, 100 * idx / n_dgps, rate, eta, n_hits))
  }
}

screen_df <- do.call(rbind, screen_results)
screen_elapsed <- (proc.time() - t_start)["elapsed"]
message(sprintf("\nPhase 1 complete: %.1f min, %d DGPs screened",
                screen_elapsed / 60, n_dgps))

# Rank Phase 1 results
screen_df <- screen_df[order(-screen_df$score), ]

# How many look promising?
n_trad_fail_s <- sum(screen_df$trad_fail_flag, na.rm = TRUE)
n_tmle_pass_s <- sum(screen_df$tmle_pass, na.rm = TRUE)
n_promising   <- sum(screen_df$trad_fail_flag & screen_df$score >= 4, na.rm = TRUE)
message(sprintf("  Trad fails: %d | TMLE passes: %d | Promising (score>=4 + trad_fail): %d",
                n_trad_fail_s, n_tmle_pass_s, n_promising))

# ===========================================================================
# PHASE 2: Refine top candidates at multiple sample sizes
# ===========================================================================

# Select top DGPs to refine
refine_mask <- screen_df$trad_fail_flag & (screen_df$score >= 4)
if (sum(refine_mask) < TOP_K_REFINE) {
  # Fall back to top by score
  refine_mask <- seq_len(min(TOP_K_REFINE, nrow(screen_df)))
} else {
  refine_mask <- which(refine_mask)
}
refine_dgps <- screen_df[head(refine_mask, TOP_K_REFINE),
                          names(dgp_grid)]
refine_dgps <- unique(refine_dgps)
n_refine <- nrow(refine_dgps)

message(sprintf("\n=== PHASE 2: Refining %d DGPs ===", n_refine))
message("  n in {", paste(REFINE_N_VALUES, collapse=", "), "}, reps=",
        REFINE_REPS, "\n")

refine_results <- list()
t_start2 <- proc.time()
total_evals <- n_refine * length(REFINE_N_VALUES)
eval_count <- 0

for (idx in seq_len(n_refine)) {
  for (n_val in REFINE_N_VALUES) {
    eval_count <- eval_count + 1
    res <- evaluate_dgp(
      refine_dgps[idx, ], n = n_val, reps = REFINE_REPS,
      seed_base = 500000 + idx * 10000 + n_val,
      alpha_A_cache = alpha_A_cache,
      alpha_Y_cache = alpha_Y_cache,
      truth_cache = truth_cache,
      ps_clip = PS_CLIP
    )
    refine_results <- c(refine_results, list(res))

    if (eval_count %% 10 == 0 || eval_count == total_evals) {
      elapsed <- (proc.time() - t_start2)["elapsed"]
      rate <- eval_count / elapsed
      eta <- (total_evals - eval_count) / rate / 60
      message(sprintf("  [%d/%d] %.0f%% | %.2f/s | ETA: %.1f min",
                      eval_count, total_evals,
                      100 * eval_count / total_evals, rate, eta))
    }
  }
}

refine_df <- do.call(rbind, refine_results)
refine_elapsed <- (proc.time() - t_start2)["elapsed"]
message(sprintf("\nPhase 2 complete: %.1f min, %d evaluations",
                refine_elapsed / 60, total_evals))

# ===========================================================================
# Combine and rank results
# ===========================================================================

# Combine Phase 1 (screen) + Phase 2 (refine) into full results
results_df <- rbind(screen_df, refine_df)
results_df <- results_df[order(-results_df$score), ]

elapsed_total <- screen_elapsed + refine_elapsed

# ===========================================================================
# Identify top candidates
# ===========================================================================

# Primary target: traditional workflow fails AND TMLE passes AND beats PSM
top_mask <- results_df$trad_fail_flag & results_df$tmle_pass &
  results_df$tmle_beats_psm
top_candidates <- results_df[top_mask, ]
top_candidates <- top_candidates[order(-top_candidates$score), ]

# If fewer than 5 with full pattern, relax to trad_fail + tmle_pass
if (nrow(top_candidates) < 5) {
  relaxed_mask <- results_df$trad_fail_flag & results_df$tmle_pass
  relaxed <- results_df[relaxed_mask & !top_mask, ]
  relaxed <- relaxed[order(-relaxed$score), ]
  top_candidates <- rbind(top_candidates,
                          head(relaxed, 5 - nrow(top_candidates)))
}

# ===========================================================================
# Save outputs
# ===========================================================================

message("\n=== Saving outputs ===")

# 1. Full grid results
write.csv(results_df, file.path(out_dir, "grid_results.csv"),
          row.names = FALSE)
message("  grid_results.csv: ", nrow(results_df), " rows")

# 2. Top candidates
write.csv(head(top_candidates, 20), file.path(out_dir, "top_demo_candidates.csv"),
          row.names = FALSE)
message("  top_demo_candidates.csv: ", min(nrow(top_candidates), 20), " rows")

# 3. Detailed replicates for top 5
if (nrow(top_candidates) >= 1) {
  n_save <- min(5, nrow(top_candidates))
  top_replicates <- vector("list", n_save)

  for (k in seq_len(n_save)) {
    tc <- top_candidates[k, ]
    alpha_A_val <- alpha_A_cache[[as.character(tc$b_rr)]]
    y_key <- paste(tc[, y_param_cols], collapse = "_")
    alpha_Y_val <- alpha_Y_cache[[y_key]]
    true_rd     <- truth_cache[[y_key]]["RD"]

    reps <- vector("list", REFINE_REPS)
    for (r in seq_len(REFINE_REPS)) {
      d <- generate_data(
        n = tc$n, seed = 900000 + k * 10000 + r,
        alpha_A = alpha_A_val, b_rr = tc$b_rr,
        alpha_Y = alpha_Y_val,
        beta_age1 = tc$beta_age1, beta_age2 = tc$beta_age2,
        beta_age3 = tc$beta_age3,
        beta_rr1 = tc$beta_rr1, beta_rr2 = tc$beta_rr2,
        beta_ckdA = tc$beta_ckdA, beta_ageA = tc$beta_ageA,
        beta_age2A = tc$beta_age2A, beta_rrA = tc$beta_rrA
      )
      reps[[r]] <- run_one_replicate(d, true_rd, PS_CLIP)
    }
    top_replicates[[k]] <- list(
      params = as.list(tc[, c("n", "b_rr", y_param_cols)]),
      alpha_A = alpha_A_val,
      alpha_Y = alpha_Y_val,
      true_rd = true_rd,
      replicates = reps
    )
  }
  saveRDS(top_replicates, file.path(out_dir, "top_candidate_replicates.rds"))
  message("  top_candidate_replicates.rds: ", n_save, " candidates")
}

# 4. Search summary
n_total <- nrow(results_df)
n_trad_fail <- sum(results_df$trad_fail_flag, na.rm = TRUE)
n_tmle_pass <- sum(results_df$tmle_pass, na.rm = TRUE)
n_both <- sum(results_df$trad_fail_flag & results_df$tmle_pass, na.rm = TRUE)
n_full_pattern <- sum(top_mask, na.rm = TRUE)

summary_lines <- c(
  "=== Manuscript Demo Search Summary ===",
  sprintf("Date: %s", Sys.time()),
  sprintf("Mode: %s", mode),
  sprintf("Phase 1: n=%d, reps=%d, %d DGPs screened",
          SCREEN_N, SCREEN_REPS, nrow(screen_df)),
  sprintf("Phase 2: n in {%s}, reps=%d, %d DGPs refined",
          paste(REFINE_N_VALUES, collapse=","), REFINE_REPS, n_refine),
  sprintf("Total evaluations: %d", n_total),
  sprintf("Elapsed time: %.1f minutes (screen: %.1f, refine: %.1f)",
          elapsed_total / 60, screen_elapsed / 60, refine_elapsed / 60),
  "",
  "--- Pattern counts (across all evaluations) ---",
  sprintf("Traditional workflow fails: %d (%.1f%%)", n_trad_fail,
          100 * n_trad_fail / n_total),
  sprintf("TMLE passes gate:          %d (%.1f%%)", n_tmle_pass,
          100 * n_tmle_pass / n_total),
  sprintf("Both (trad fail + TMLE pass): %d (%.1f%%)", n_both,
          100 * n_both / n_total),
  sprintf("Full pattern (+ beats PSM):   %d (%.1f%%)", n_full_pattern,
          100 * n_full_pattern / n_total),
  ""
)

# Top 5 settings
summary_lines <- c(summary_lines, "--- Top candidates ---")
n_show <- min(5, nrow(top_candidates))
if (n_show > 0) {
  for (k in seq_len(n_show)) {
    tc <- top_candidates[k, ]
    summary_lines <- c(summary_lines, sprintf(
      paste0("#%d [score=%d] n=%d, b_rr=%.1f, rr1=%.2f, rr2=%.2f, ",
             "age1=%.2f, age2=%.2f, age3=%.2f, ckdA=%.2f, ageA=%.2f, ",
             "age2A=%.2f, rrA=%.2f, y_prev=%.2f"),
      k, tc$score, tc$n, tc$b_rr, tc$beta_rr1, tc$beta_rr2,
      tc$beta_age1, tc$beta_age2,
      tc$beta_age3, tc$beta_ckdA, tc$beta_ageA, tc$beta_age2A,
      tc$beta_rrA, tc$target_y_prev
    ))
    summary_lines <- c(summary_lines, sprintf(
      paste0("   PS diag: frac_extreme=%.3f, ess=%.3f, matched=%.3f, ",
             "max_smd=%.3f, trad_fail=%.0f%%"),
      tc$med_frac_extreme, tc$med_ess_frac, tc$med_matched_frac,
      tc$med_max_smd, tc$frac_trad_fail * 100
    ))
    summary_lines <- c(summary_lines, sprintf(
      "   PSM:  bias=%.5f, RMSE=%.5f, cov=%.3f",
      tc$psm_bias, tc$psm_rmse, tc$psm_coverage
    ))
    summary_lines <- c(summary_lines, sprintf(
      "   TMLE: bias=%.5f, RMSE=%.5f, cov=%.3f, SE/SD=%.3f",
      tc$tmle_bias, tc$tmle_rmse, tc$tmle_coverage, tc$tmle_se_sd_ratio
    ))
    summary_lines <- c(summary_lines, sprintf(
      "   True RD=%.5f | tmle_pass=%s | tmle_beats_psm=%s",
      tc$true_rd, tc$tmle_pass, tc$tmle_beats_psm
    ))
    summary_lines <- c(summary_lines, "")
  }
} else {
  summary_lines <- c(summary_lines,
                     "No candidates found matching the target pattern.",
                     "Consider widening the grid or adjusting thresholds.")
}

writeLines(summary_lines, file.path(out_dir, "search_summary.txt"))
message("  search_summary.txt written")

# ===========================================================================
# Console output
# ===========================================================================

cat("\n")
cat("============================================================\n")
cat("  MANUSCRIPT DEMO SEARCH COMPLETE\n")
cat("============================================================\n")
cat(sprintf("  Phase 1 screened:         %d DGPs\n", nrow(screen_df)))
cat(sprintf("  Phase 2 refined:          %d DGPs x %d n-values\n",
            n_refine, length(REFINE_N_VALUES)))
cat(sprintf("  Total evaluations:        %d\n", n_total))
cat(sprintf("  Traditional fails:        %d\n", n_trad_fail))
cat(sprintf("  TMLE passes:              %d\n", n_tmle_pass))
cat(sprintf("  Both (fail + pass):       %d\n", n_both))
cat(sprintf("  Full pattern (+ beats):   %d\n", n_full_pattern))
cat(sprintf("  Elapsed:                  %.1f min\n", elapsed_total / 60))
cat("============================================================\n\n")

if (nrow(top_candidates) > 0) {
  cat("--- Top 10 Candidates ---\n\n")
  n_print <- min(10, nrow(top_candidates))
  for (k in seq_len(n_print)) {
    tc <- top_candidates[k, ]
    cat(sprintf("[#%d | score=%d] n=%d, b_rr=%.1f\n", k, tc$score,
                tc$n, tc$b_rr))
    cat(sprintf("  DGP: rr1=%.2f rr2=%.2f age3=%.2f ckdA=%.2f ageA=%.2f ",
                tc$beta_rr1, tc$beta_rr2, tc$beta_age3,
                tc$beta_ckdA, tc$beta_ageA))
    cat(sprintf("age2A=%.2f rrA=%.2f y_prev=%.2f\n",
                tc$beta_age2A, tc$beta_rrA, tc$target_y_prev))
    cat(sprintf("  PS diag: extreme=%.1f%% ess=%.1f%% matched=%.1f%%\n",
                tc$med_frac_extreme * 100, tc$med_ess_frac * 100,
                tc$med_matched_frac * 100))
    cat(sprintf("  PSM:  RMSE=%.4f cov=%.3f | TMLE: RMSE=%.4f cov=%.3f\n",
                tc$psm_rmse, tc$psm_coverage,
                tc$tmle_rmse, tc$tmle_coverage))

    # One-line interpretation
    if (tc$trad_fail_flag && tc$tmle_pass && tc$tmle_beats_psm) {
      cat("  >> Poor PS diagnostics / matching feasibility, but TMLE passes",
          "with lower RMSE than PS-matched regression.\n")
    } else if (tc$trad_fail_flag && tc$tmle_pass) {
      cat("  >> Poor PS diagnostics, TMLE passes gate but does not clearly",
          "beat PS-matched regression.\n")
    } else {
      cat("  >> Partial pattern match.\n")
    }
    cat("\n")
  }
} else {
  cat("No candidates matched the target pattern.\n")
  cat("Suggestions:\n")
  cat("  - Increase sample sizes (try refine mode)\n")
  cat("  - Increase reps for more stable estimates\n")
  cat("  - Widen the outcome parameter grid\n")
}

# ===========================================================================
# Optional: Rerun top candidates with SuperLearner
# ===========================================================================

if (RERUN_TOP_WITH_SL && USE_SUPERLEARNER) {
  if (requireNamespace("SuperLearner", quietly = TRUE)) {
    message("\n=== Rerunning top candidates with SuperLearner ===")
    n_rerun <- min(3, nrow(top_candidates))
    for (k in seq_len(n_rerun)) {
      tc <- top_candidates[k, ]
      alpha_A_val <- alpha_A_cache[[as.character(tc$b_rr)]]
      y_key <- paste(tc[, y_param_cols], collapse = "_")
      alpha_Y_val <- alpha_Y_cache[[y_key]]
      true_rd     <- truth_cache[[y_key]]["RD"]

      sl_ests <- numeric(REFINE_REPS)
      for (r in seq_len(REFINE_REPS)) {
        d <- generate_data(
          n = tc$n, seed = 990000 + k * 10000 + r,
          alpha_A = alpha_A_val, b_rr = tc$b_rr,
          alpha_Y = alpha_Y_val,
          beta_age1 = tc$beta_age1, beta_age2 = tc$beta_age2,
          beta_age3 = tc$beta_age3,
          beta_rr1 = tc$beta_rr1, beta_rr2 = tc$beta_rr2,
          beta_ckdA = tc$beta_ckdA, beta_ageA = tc$beta_ageA,
          beta_age2A = tc$beta_age2A, beta_rrA = tc$beta_rrA
        )
        sl_res <- tryCatch(run_tmle_sl(d, PS_CLIP),
                           error = function(e) list(estimate = NA))
        sl_ests[r] <- sl_res$estimate
      }
      ok <- !is.na(sl_ests)
      if (sum(ok) >= 5) {
        cat(sprintf(
          "  SL rerun #%d: bias=%.5f, RMSE=%.5f (vs glmnet RMSE=%.5f)\n",
          k, mean(sl_ests[ok]) - true_rd,
          sqrt(mean((sl_ests[ok] - true_rd)^2)), tc$tmle_rmse))
      }
    }
  } else {
    message("SuperLearner package not available — skipping SL rerun.")
  }
}

message("\nDone. Outputs in: ", out_dir)
