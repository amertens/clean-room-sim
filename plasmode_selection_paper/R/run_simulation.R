# Main simulation driver for the plasmode-selection paper.
#
# Loops over DGPs, sample sizes, and outer replicates. For each
# replicate, runs (a) fixed-library TMLE (one or more library tiers),
# and (b) plasmode-selected TMLE on the cleanTMLE candidate grid.
# Saves per-replicate rows so the manuscript can compute any summary.

suppressPackageStartupMessages({
  library(cleanTMLE)
  library(SuperLearner)
})

source_here <- function(f) source(file.path("R", f))

run_simulation <- function(
  dgps         = c("linear", "nonlinear_smooth", "interactions",
                   "sparse", "high_dim_noise", "positivity_strain"),
  n_per_rep    = 1000L,
  n_reps       = 500L,
  inner_reps   = 50L,
  true_RD      = -0.05,
  fixed_libs   = list(
    parametric = c("SL.glm", "SL.mean"),
    rich       = c("SL.glm", "SL.glmnet", "SL.gam",
                   "SL.ranger", "SL.mean")
  ),
  fixed_trunc  = 0.01,
  noise_p      = 20L,
  out_path     = "results/sim_full.rds",
  checkpoint_every = 5L,
  resume       = TRUE,
  seed_base    = 2026L,
  log_progress = TRUE
) {
  source_here("dgps.R")
  source_here("candidates.R")
  source_here("workflows.R")

  covariates_base <- c("age", "sex", "biomarker", "comorbidity",
                       "bmi_z", "smoke")
  candidates <- build_candidate_grid()

  # ---- Resume support ------------------------------------------
  # If an interrupted-run file exists and resume = TRUE, load it and
  # skip (dgp, rep) cells that are already complete. Each cell is
  # counted complete once all expected workflow rows for that cell
  # are present without NA estimates.
  done_cells <- list()
  results    <- list()
  if (isTRUE(resume) && file.exists(out_path)) {
    prev <- tryCatch(readRDS(out_path), error = function(e) NULL)
    if (!is.null(prev) && nrow(prev) > 0L) {
      results <- list(prev)
      expected_workflows <- c(paste0("fixed_", names(fixed_libs)),
                               "plasmode_min_rmse",
                               "plasmode_fiord_two_stage")
      key   <- paste(prev$dgp, prev$rep, sep = "::")
      tab   <- table(key)
      ekey  <- names(tab)[tab >= length(expected_workflows)]
      done_cells <- as.list(ekey)
      if (log_progress)
        message("Resuming from ", out_path, ": ", length(done_cells),
                 " (dgp x rep) cells already complete.")
    }
  }
  is_done <- function(dgp, rep_i) {
    paste(dgp, rep_i, sep = "::") %in% done_cells
  }

  job_i <- 0L
  total <- length(dgps) * n_reps
  t_start <- Sys.time()

  for (dgp in dgps) {
    if (log_progress) message("==> DGP = ", dgp)
    covariates <- if (dgp == "high_dim_noise") {
      c(covariates_base, sprintf("n%02d", seq_len(noise_p)))
    } else covariates_base

    for (rep_i in seq_len(n_reps)) {
      job_i <- job_i + 1L
      if (is_done(dgp, rep_i)) {
        if (log_progress && rep_i %% 25L == 0L)
          message(sprintf("  [%s] rep %d (already complete, skipping)",
                          dgp, rep_i))
        next
      }
      if (log_progress && rep_i %% 5L == 0L) {
        elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
        rate    <- if (job_i > 0L) elapsed / job_i else NA_real_
        eta     <- rate * (total - job_i)
        message(sprintf("  [%s] rep %d / %d  (%d / %d total, %.1f min elapsed, ~%.1f min remaining)",
                         dgp, rep_i, n_reps, job_i, total, elapsed, eta))
      }
      sd <- seed_base + 1000L * match(dgp, dgps) + rep_i

      set.seed(sd)
      d <- make_data(n_per_rep, dgp = dgp, true_RD = true_RD,
                      noise_p = noise_p)

      for (lib_name in names(fixed_libs)) {
        fix <- tryCatch(
          fit_fixed(d, covariates,
                    library    = fixed_libs[[lib_name]],
                    truncation = fixed_trunc,
                    seed       = sd),
          error = function(e) {
            list(estimate = NA_real_, se = NA_real_,
                 ci_lower = NA_real_, ci_upper = NA_real_,
                 workflow = "fixed",
                 library_label = paste(fixed_libs[[lib_name]], collapse = "+"),
                 truncation = fixed_trunc,
                 selected_candidate = NA_character_)
          })
        results[[length(results) + 1L]] <- data.frame(
          dgp           = dgp,
          rep           = rep_i,
          workflow      = paste0("fixed_", lib_name),
          library_label = fix$library_label,
          truncation    = fix$truncation,
          selected_candidate = NA_character_,
          estimate      = fix$estimate, se = fix$se,
          ci_lower      = fix$ci_lower, ci_upper = fix$ci_upper,
          stringsAsFactors = FALSE)
      }

      pls <- tryCatch(
        fit_plasmode(d, covariates, candidates,
                      inner_reps   = inner_reps,
                      inner_effect = true_RD,
                      rule         = "min_rmse",
                      seed         = sd),
        error = function(e) {
          list(estimate = NA_real_, se = NA_real_,
               ci_lower = NA_real_, ci_upper = NA_real_,
               workflow = "plasmode", library_label = NA_character_,
               truncation = NA_real_, selected_candidate = NA_character_)
        })
      results[[length(results) + 1L]] <- data.frame(
        dgp           = dgp,
        rep           = rep_i,
        workflow      = "plasmode_min_rmse",
        library_label = pls$library_label,
        truncation    = pls$truncation,
        selected_candidate = pls$selected_candidate %||% NA_character_,
        estimate      = pls$estimate, se = pls$se,
        ci_lower      = pls$ci_lower, ci_upper = pls$ci_upper,
        stringsAsFactors = FALSE)

      pls2 <- tryCatch(
        fit_plasmode(d, covariates, candidates,
                      inner_reps   = inner_reps,
                      inner_effect = true_RD,
                      rule         = "fiord_two_stage",
                      seed         = sd),
        error = function(e) {
          list(estimate = NA_real_, se = NA_real_,
               ci_lower = NA_real_, ci_upper = NA_real_,
               workflow = "plasmode", library_label = NA_character_,
               truncation = NA_real_, selected_candidate = NA_character_)
        })
      results[[length(results) + 1L]] <- data.frame(
        dgp           = dgp,
        rep           = rep_i,
        workflow      = "plasmode_fiord_two_stage",
        library_label = pls2$library_label,
        truncation    = pls2$truncation,
        selected_candidate = pls2$selected_candidate %||% NA_character_,
        estimate      = pls2$estimate, se = pls2$se,
        ci_lower      = pls2$ci_lower, ci_upper = pls2$ci_upper,
        stringsAsFactors = FALSE)

      # Persist incrementally so a crash mid-run does not lose work.
      # Atomic write via a temp file + rename guards against
      # half-written rds files if R is interrupted mid-save.
      if (rep_i %% checkpoint_every == 0L) {
        out_df <- do.call(rbind, results)
        .atomic_save_rds(out_df, out_path)
      }
    }
    # End-of-DGP checkpoint (always)
    out_df <- do.call(rbind, results)
    .atomic_save_rds(out_df, out_path)
    if (log_progress)
      message("Saved interim results after DGP ", dgp,
              " (", nrow(out_df), " rows so far).")
  }

  out_df <- do.call(rbind, results)
  .atomic_save_rds(out_df, out_path)
  if (log_progress)
    message("Final write to ", out_path, " (", nrow(out_df), " rows).")
  invisible(out_df)
}

# Atomic save: write to a temp file then rename. On Windows the
# rename is not strictly atomic but is far better than a direct
# saveRDS that can be left half-written if the process is killed.
.atomic_save_rds <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  if (file.exists(path)) file.remove(path)
  file.rename(tmp, path)
}

# Local null-coalescing helper
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
