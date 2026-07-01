#!/usr/bin/env Rscript
# =============================================================================
# separating_grid_study.R
# Realistic separating-candidate grid for manuscript Section 9.4.
#
# Reuses the existing DGP (run_simulation.R) and the cleanTMLE public functions
# (create_simple_lock, tmle_candidate, run_plasmode_feasibility,
# run_plasmode_dq_stress, select_tmle_candidate). No parallel simulation engine.
#
# SCOPE (decided 2026-06-08): Arm A only.
#   Arm A  misspec=FALSE, overlap=1.5 -> marginal overlap, linear truth.
#          The truncation lever: GLM propensity models differing only in PS
#          truncation. This is the actual Section 9.4 minimax mechanism, and it
#          is fast and numerically stable (no SuperLearner).
#   Arm B (SuperLearner library lever) was dropped: SL.gam AND SL.glmnet both
#          leak/crash the R session after enough fits, so that arm needs heavy
#          subprocess isolation for little gain on the Section 9.4 claim.
#
# ARCHITECTURE: controller + per-batch workers, for crash isolation and resume.
#   - Controller (default): calibrates near-positivity, writes + fingerprints the
#     locked config, spawns one fresh R worker per batch, then aggregates.
#   - Worker (GRID_ROLE=worker): runs ONE (arm, batch) cell, saves it, exits, so
#     each process frees its memory and a crash costs one batch, not the run.
#
# ANTI-FABRICATION: severities are default_dq_scenarios("regulatory_standard")
# plus the calibration rule, fixed and locked before the sweep. They are never
# tuned to force a divergence. The study reports whatever it computes.
#
# Env:
#   GRID_ROLE   "controller" (default) | "worker"
#   GRID_ARM    arm id for a worker (default "A")
#   GRID_BATCH_SEED  integer seed for a worker's batch
#   GRID_REPS   override reps (worker + smoke)
#   GRID_MODE   "full" (default) | "smoke" (1 batch x 20 reps, no lock written)
# =============================================================================

suppressPackageStartupMessages({ library(cleanTMLE); library(digest) })

ROOT  <- "C:/Users/andre/OneDrive/Documents/clean-room-sim"
OUT   <- file.path(ROOT, "sandbox", "separating_grid")
CELLS <- file.path(OUT, "cells")
dir.create(CELLS, showWarnings = FALSE, recursive = TRUE)

ROLE <- tolower(Sys.getenv("GRID_ROLE", "controller"))
MODE <- tolower(Sys.getenv("GRID_MODE", "full"))

## ---- locked replication design ----
N_OBS        <- 4000L
EFFECT_SIZES <- c(0.05)
BATCHES      <- 5L
REPS         <- 500L
BASE_SEED    <- 20260608L
SEED_STRIDE  <- 1000000L
COVARS       <- c("age", "sex", "biomarker", "comorbidity", "ckd")
TREATMENT    <- "treatment"
OUTCOME      <- "event_24"
TRUNCS       <- c(0.005, 0.01, 0.025, 0.05, 0.10)
RULES        <- c("min_rmse", "min_max_rmse", "fiord_two_stage")
THRESHOLDS   <- list(max_abs_bias = 0.02, min_coverage = 0.85, max_rmse_ratio = 1.5)

if (MODE == "smoke") { BATCHES <- 1L; REPS <- 20L }
if (nzchar(Sys.getenv("GRID_REPS"))) REPS <- as.integer(Sys.getenv("GRID_REPS"))

## ---- reuse generate_data + compute_truth from run_simulation.R (no side effects) ----
local({
  exprs <- parse(file = file.path(ROOT, "run_simulation.R"))
  for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-"))) {
    lhs <- e[[2]]
    if (is.name(lhs) && as.character(lhs) %in% c("generate_data", "compute_truth"))
      eval(e, envir = globalenv())
  }
})
stopifnot(is.function(generate_data), is.function(compute_truth))

## ---- Arm A candidate grid: GLM propensity, truncation lever only ----
glab <- function(t) sprintf("glm_t%s", sub("^0\\.", "", format(t, trim = TRUE)))
build_grid <- function() lapply(TRUNCS, function(t)
  tmle_candidate(glab(t), g_library = "SL.glm", q_library = "SL.glm", truncation = t))
GRID     <- build_grid()
CAND_IDS <- vapply(GRID, function(x) x$candidate_id, character(1))

## ---- arms (Arm A only) ----
ARMS <- list(
  A = list(misspec = FALSE, q0_library = NULL,
           label = "marginal overlap, linear surface (truncation lever)")
)

## ---- near-positivity calibration: keep a slope only if induced PS stays in envelope ----
# Mirrors ps_pos = plogis(mean(lp) + slope*(lp - mean(lp))) (plasmode_dq.R:915-918).
calibrate_nearpos <- function(misspec, slopes = c(1.5, 2.0), n = 8000L, seed = 101L) {
  d  <- generate_data(n = n, overlap_strength = 1.5, effect_size = -0.05,
                      misspec = misspec, seed = seed)
  ps <- stats::predict(stats::glm(stats::reformulate(COVARS, TREATMENT),
                                   data = d, family = stats::binomial()),
                       type = "response")
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  lp <- stats::qlogis(ps); lpbar <- mean(lp)
  rows <- lapply(slopes, function(s) {
    pp <- stats::plogis(lpbar + s * (lp - lpbar)); q <- stats::quantile(pp, c(0.01, 0.99))
    keep <- q[1] >= 0.01 && q[2] <= 0.99 && min(pp) >= 0.002 && max(pp) <= 0.998
    data.frame(slope = s, q01 = round(q[1], 4), q99 = round(q[2], 4),
               min = round(min(pp), 4), max = round(max(pp), 4), keep = keep, row.names = NULL)
  })
  tab <- do.call(rbind, rows)
  list(keep = tab$slope[tab$keep], table = tab)
}

## ---- anchored, locked threat list (regulatory_standard + calibrated near-pos) ----
build_threats <- function(kept_slopes) {
  dqs <- default_dq_scenarios("regulatory_standard")
  if (length(kept_slopes) > 0L) dqs$near_positivity <- list(slopes = kept_slopes)
  else dqs$near_positivity <- NULL
  dqs
}

## ---- one cell: generate, lock, score baseline + threats ----
run_cell <- function(misspec, q0_library, threats, batch_seed, reps) {
  d <- generate_data(n = N_OBS, overlap_strength = 1.5, effect_size = -0.05,
                     misspec = misspec, seed = batch_seed)
  lock <- create_simple_lock(data = d, treatment = TREATMENT, outcome = OUTCOME,
                             covariates = COVARS, sl_library = "SL.glm",
                             plasmode_reps = reps, seed = batch_seed)
  feas <- run_plasmode_feasibility(lock, tmle_candidates = GRID, effect_sizes = EFFECT_SIZES,
                                   reps = reps, q0_library = q0_library, verbose = FALSE)
  dq <- run_plasmode_dq_stress(lock, tmle_candidates = GRID, effect_sizes = EFFECT_SIZES,
                               reps = reps, data_quality_scenarios = threats,
                               q0_library = q0_library, verbose = FALSE)
  list(feas = feas, dq = dq)
}

cell_path <- function(arm, seed) file.path(CELLS, sprintf("cell_%s_%d_r%d.rds", arm, seed, REPS))

# =============================================================================
# WORKER: run exactly one (arm, batch) cell from the locked config and save it.
# =============================================================================
if (ROLE == "worker") {
  arm  <- Sys.getenv("GRID_ARM", "A")
  seed <- as.integer(Sys.getenv("GRID_BATCH_SEED"))
  sl_env <- Sys.getenv("GRID_NEARPOS_SLOPES")   # calibrated kept slopes, comma-separated
  kept   <- if (nzchar(sl_env)) as.numeric(strsplit(sl_env, ",")[[1]]) else numeric(0)
  threats <- build_threats(kept)
  cat(sprintf("[worker] arm=%s seed=%d reps=%d ...\n", arm, seed, REPS)); flush(stdout())
  t0 <- Sys.time()
  cell <- run_cell(ARMS[[arm]]$misspec, ARMS[[arm]]$q0_library, threats, seed, REPS)
  saveRDS(list(arm = arm, seed = seed, reps = REPS,
               feas = cell$feas, dq = cell$dq),
          cell_path(arm, seed))
  cat(sprintf("[worker] arm=%s seed=%d done (%.1f min)\n", arm, seed,
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush(stdout())
  quit(save = "no", status = 0)
}

# =============================================================================
# CONTROLLER
# =============================================================================
mcse  <- function(v) { v <- v[is.finite(v)]; if (length(v) > 1L) sd(v)/sqrt(length(v)) else NA_real_ }
modef <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return(NA_character_)
                       names(sort(table(v), decreasing = TRUE))[1] }

summarise_arm <- function(cells) {
  pb_base <- pb_worst <- pb_bind <- setNames(vector("list", length(CAND_IDS)), CAND_IDS)
  rule_pick <- setNames(vector("list", length(RULES)), RULES)
  for (cl in cells) {
    fm <- cl$feas$metrics; dm <- cl$dq$metrics
    dm_nb <- dm[dm$scenario != "none", , drop = FALSE]
    for (cid in CAND_IDS) {
      pb_base[[cid]] <- c(pb_base[[cid]], mean(fm$rmse[fm$candidate == cid], na.rm = TRUE))
      sub <- dm_nb[dm_nb$candidate == cid, , drop = FALSE]
      if (nrow(sub) && any(is.finite(sub$rmse))) {
        wi <- which.max(sub$rmse)
        pb_worst[[cid]] <- c(pb_worst[[cid]], sub$rmse[wi])
        pb_bind[[cid]]  <- c(pb_bind[[cid]], paste0(sub$scenario[wi], ":", sub$level[wi]))
      } else { pb_worst[[cid]] <- c(pb_worst[[cid]], NA_real_); pb_bind[[cid]] <- c(pb_bind[[cid]], NA_character_) }
    }
    for (r in RULES)
      rule_pick[[r]] <- c(rule_pick[[r]],
        tryCatch(select_tmle_candidate(cl$feas, rule = r, dq_results = cl$dq)$candidate_id,
                 error = function(e) NA_character_))
  }
  per_cand <- do.call(rbind, lapply(CAND_IDS, function(cid) data.frame(
    candidate = cid,
    baseline_rmse  = round(mean(pb_base[[cid]],  na.rm = TRUE), 5), baseline_mcse  = round(mcse(pb_base[[cid]]),  5),
    worstcase_rmse = round(mean(pb_worst[[cid]], na.rm = TRUE), 5), worstcase_mcse = round(mcse(pb_worst[[cid]]), 5),
    binding_threat = modef(pb_bind[[cid]]), row.names = NULL, stringsAsFactors = FALSE)))
  sels <- lapply(rule_pick, function(p) list(per_batch = p, consensus = modef(p),
                 stable = length(unique(p[!is.na(p)])) == 1L))
  list(per_candidate = per_cand, selections = sels)
}

PHASE       <- tolower(Sys.getenv("GRID_PHASE", "lock"))   # "lock" | "aggregate"
batch_seeds <- BASE_SEED + (seq_len(BATCHES) - 1L) * SEED_STRIDE
cat(sprintf("[controller] phase=%s MODE=%s n=%d batches=%d reps=%d candidates=%s\n",
            PHASE, MODE, N_OBS, BATCHES, REPS, paste(CAND_IDS, collapse = ","))); flush(stdout())

## ---- PHASE "lock": calibrate near-positivity, write + fingerprint the config ----
if (PHASE == "lock") {
  cat("\n== near-positivity calibration ==\n")
  calib <- list()
  for (a in names(ARMS)) {
    cb <- calibrate_nearpos(ARMS[[a]]$misspec); calib[[a]] <- cb
    cat(sprintf("Arm %s (%s):\n", a, ARMS[[a]]$label)); print(cb$table)
    cat(sprintf("  kept slopes: %s\n", if (length(cb$keep)) paste(cb$keep, collapse = ", ") else "NONE"))
  }
  config <- list(
    created_for = "manuscript Section 9.4 realistic separating grid (Arm A only)",
    n_obs = N_OBS, effect_sizes = EFFECT_SIZES, batches = BATCHES, reps = REPS,
    base_seed = BASE_SEED, seed_stride = SEED_STRIDE, batch_seeds = batch_seeds,
    covariates = COVARS, treatment = TREATMENT, outcome = OUTCOME, truncations = TRUNCS,
    grid = lapply(GRID, function(c) c[c("candidate_id","g_library","q_library","truncation","cv_scheme")]),
    arms = lapply(names(ARMS), function(a) list(arm = a, misspec = ARMS[[a]]$misspec,
              q0_library = ARMS[[a]]$q0_library, label = ARMS[[a]]$label,
              nearpos_kept_slopes = calib[[a]]$keep, nearpos_calibration = calib[[a]]$table)),
    threat_preset = "regulatory_standard", thresholds = THRESHOLDS, rules = RULES)
  config$fingerprint <- digest::digest(config, algo = "sha256")
  saveRDS(config, file.path(OUT, "locked_config.rds"))
  if (requireNamespace("yaml", quietly = TRUE)) yaml::write_yaml(config, file.path(OUT, "locked_config.yaml"))
  cat(sprintf("\n[controller] LOCKED fingerprint=%s\n", config$fingerprint))
  for (a in names(ARMS))
    cat(sprintf("WORKER_PLAN arm=%s slopes=%s seeds=%s\n",
                a, paste(calib[[a]]$keep, collapse = ","), paste(batch_seeds, collapse = ",")))
  cat("\nLOCK DONE\n"); quit(save = "no", status = 0)
}

## ---- PHASE "aggregate": read the per-batch cell files, summarise, write results ----
if (PHASE == "aggregate") {
  cfg <- if (file.exists(file.path(OUT, "locked_config.rds"))) readRDS(file.path(OUT, "locked_config.rds")) else NULL
  results <- list()
  for (a in names(ARMS)) {
    cells <- Filter(Negate(is.null), lapply(batch_seeds, function(s) {
      p <- cell_path(a, s); if (file.exists(p)) readRDS(p) else NULL }))
    if (!length(cells)) { cat(sprintf("Arm %s: no cells produced.\n", a)); next }
    results[[a]] <- list(summary = summarise_arm(cells),
                         truth = compute_truth(1e5L, 1.5, -0.05, 7L, misspec = ARMS[[a]]$misspec)$RD,
                         n_batches = length(cells))
  }
  saveRDS(list(config = cfg, results = results, reps = REPS),
          file.path(OUT, "separating_grid_results.rds"))
  tidy <- do.call(rbind, lapply(names(results), function(a) {
    pc <- results[[a]]$summary$per_candidate; pc$arm <- a
    pc$true_RD <- round(results[[a]]$truth, 5); pc$n_batches <- results[[a]]$n_batches; pc }))
  if (!is.null(tidy)) write.csv(tidy, file.path(OUT, "separating_grid_results.csv"), row.names = FALSE)
  for (a in names(results)) {
    cat(sprintf("\n== Arm %s (%s)  [%d batches] ==\n", a, ARMS[[a]]$label, results[[a]]$n_batches))
    print(results[[a]]$summary$per_candidate)
    for (r in RULES) { s <- results[[a]]$summary$selections[[r]]
      cat(sprintf("  rule %-16s -> %s (stable across batches: %s)\n", r, s$consensus, s$stable)) }
  }
  cat("\nAGGREGATE DONE\n"); quit(save = "no", status = 0)
}

stop("Unknown GRID_PHASE: ", PHASE)
