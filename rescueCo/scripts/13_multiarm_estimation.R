# ============================================================
# Multi-arm Stage 4: estimation down the declared ladder
# ============================================================
# Joins the outcomes back onto each lock (the first outcome access in the
# pipeline), then estimates, per lock and outcome, the declared estimand
# ladder: the ATE when the support verdict allows it, and every feasible
# fallback (trimmed ATE with refit, complete-case ATT, augmented ATO).
# Every row carries the estimand, the support verdict, the caveat, and the
# implausibility flags. The protocol contrast C1 is the study's point: its
# ATE is expected to be refused by the ladder and estimated only as a
# labelled override row for the reconciliation table.
#
# IPCW via Delta is the primary handling for the follow-up outcomes
# (death_by_6mo, gose, severe_disability, gose_favorable), complete case for
# the rest, matching the main pipeline's delta flags (analysis/42).
#
# Checkpointed per (contrast x outcome); a crash costs one pair.
# Env: MULTIARM_EST_OUTCOMES  comma list (default:
#        death_by_6mo,time_to_arrival_minutes)
#      MULTIARM_EST_ONLY      comma list of contrasts (default: all)
#      TMLE_WORKERS           parallel workers over pairs (default 4)
# ============================================================

suppressMessages(library(cleanTMLE))
suppressMessages(library(SuperLearner))
suppressMessages(library(tmle))

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                 "bootstrap.R"))
source("rescueCo/R/utils.R")
cr_log("=== Multi-arm Stage 4: estimand-ladder estimation ===")

OUT <- "rescueCo/results/multiarm"
locks    <- readRDS(file.path(OUT, "multiarm_locks.rds"))
outcomes <- readRDS(file.path(OUT, "multiarm_outcomes.rds"))

DELTA_OUTCOMES <- c("death_by_6mo", "gose", "severe_disability",
                    "gose_favorable")
FAMILY <- c(death_by_6mo = "binomial", death_in_hospital = "binomial",
            severe_disability = "binomial", gose_favorable = "binomial",
            gose = "gaussian", transferred_out = "binomial",
            transferred_out_higher_care = "binomial",
            time_to_arrival_minutes = "gaussian",
            time_to_def_intervention_minutes = "gaussian",
            time_arrival_to_intervention_minutes = "gaussian")

oc_env <- trimws(strsplit(Sys.getenv(
  "MULTIARM_EST_OUTCOMES", "death_by_6mo,time_to_arrival_minutes"),
  ",")[[1]])
OUTCOMES <- intersect(oc_env[nzchar(oc_env)], names(FAMILY))
only <- trimws(strsplit(Sys.getenv("MULTIARM_EST_ONLY", ""), ",")[[1]])
CONTRASTS <- if (length(only[nzchar(only)])) only[nzchar(only)] else
  names(locks)

SL_LIB <- build_sl_library(role = "Q", n_eff = 3000,
                           preset = "rwe_wide")$library

pairs <- expand.grid(contrast = CONTRASTS, outcome = OUTCOMES,
                     stringsAsFactors = FALSE)
cr_log(paste("Pairs to estimate:", nrow(pairs), "(",
             paste(OUTCOMES, collapse = ", "), ")"))

# Slim per-contrast design extracts: the full design checkpoints carry the
# SuperLearner fit (hundreds of MB); estimation needs only the scores, the
# support assessment, and the feasibility table.
for (cn in CONTRASTS) {
  lite <- file.path(OUT, paste0("design_lite_", cn, ".rds"))
  full <- file.path(OUT, paste0("design_", cn, ".rds"))
  if (!file.exists(lite) && file.exists(full)) {
    st <- readRDS(full)
    saveRDS(list(ps_raw = st$ps$ps_raw %||% st$ps$ps, ps = st$ps$ps,
                 support = st$support, feasibility = st$feasibility),
            lite)
    rm(st); gc()
  }
}

fit_pair <- function(cn, oc) {
  ck <- file.path(OUT, sprintf("ladder_%s_%s.rds", cn, oc))
  if (file.exists(ck) && !nzchar(Sys.getenv("MULTIARM_FORCE")))
    return(readRDS(ck))
  lk <- locks[[cn]]
  lk$sl_library <- SL_LIB   # refits inside the ladder use the study library
  st <- readRDS(file.path(OUT, paste0("design_lite_", cn, ".rds")))

  # Stage 4 outcome join: the one place outcomes reach a lock.
  oc_col <- outcomes[[oc]][match(lk$data$patient_id, outcomes$patient_id)]
  lk$data[[oc]] <- oc_col
  lk$outcome <- oc
  lk$.outcome_masked <- FALSE
  lk <- cleanTMLE:::.log_design_decision(lk, "outcome_access",
    sprintf("Outcome %s joined for Stage 4 estimation (design stage complete).", oc))
  psf <- cleanTMLE:::wrap_ps_fit(lk, ps_scores = st$ps_raw)

  use_ipcw <- oc %in% DELTA_OUTCOMES && anyNA(oc_col)
  res <- tryCatch(run_estimand_ladder(
    lk, psf,
    support = st$support, feasibility = st$feasibility,
    family = unname(FAMILY[[oc]]), use_ipcw = use_ipcw,
    sl_library = SL_LIB, cv_folds = 10L, prescreen_g = FALSE,
    seed = lk$seed, trim_levels = c(0.05, 0.10),
    override_reason = if (cn == "C1")
      paste("Reconciliation row: the main pipeline reports the C1 ATE as a",
            "SEVERE-flagged diagnostic; estimated here under the same",
            "labelling so the two pipelines can be compared.") else NULL,
    verbose = TRUE), error = function(e) {
      cr_log(paste("[", cn, "x", oc, "] ladder failed:",
                   conditionMessage(e)))
      NULL })
  if (!is.null(res)) {
    out <- cbind(contrast = cn, label = lk$contrast$label, outcome = oc,
                 family = unname(FAMILY[[oc]]), used_ipcw = use_ipcw,
                 res$table, stringsAsFactors = FALSE)
    attr(out, "design_log") <- res$design_log
    saveRDS(out, ck)
    return(out)
  }
  NULL
}

nw <- suppressWarnings(as.integer(Sys.getenv("TMLE_WORKERS", "4")))
if (is.na(nw) || nw < 1L) nw <- 4L
todo <- pairs[!file.exists(file.path(
  OUT, sprintf("ladder_%s_%s.rds", pairs$contrast, pairs$outcome))) |
    nzchar(Sys.getenv("MULTIARM_FORCE")), , drop = FALSE]
cr_log(paste("To fit now:", nrow(todo), "pairs;", nw, "workers"))

if (nrow(todo) > 0) {
  if (nw > 1L && nrow(todo) > 1L) {
    run_parallel <- function() {
      cl <- parallel::makeCluster(min(nw, nrow(todo)))
      on.exit(tryCatch(parallel::stopCluster(cl), error = function(e) NULL),
              add = TRUE)
      proj <- getwd()
      parallel::clusterExport(cl, c("proj", "OUT", "locks", "outcomes",
                                    "DELTA_OUTCOMES", "FAMILY", "SL_LIB",
                                    "fit_pair", "todo"), envir = environment())
      parallel::clusterEvalQ(cl, {
        setwd(proj)
        suppressMessages(library(cleanTMLE))
        suppressMessages(library(SuperLearner))
        suppressMessages(library(tmle))
        source("rescueCo/R/utils.R")
        NULL
      })
      parallel::parLapply(cl, seq_len(nrow(todo)), function(i)
        fit_pair(todo$contrast[i], todo$outcome[i]))
    }
    # One retry with a fresh cluster before serial fallback: socket workers
    # do not survive the machine sleeping.
    got <- NULL
    for (attempt in 1:2) {
      got <- tryCatch(run_parallel(), error = function(e) {
        cr_log(paste("parallel attempt", attempt, "failed:",
                     conditionMessage(e)))
        NULL })
      if (!is.null(got)) break
    }
    if (is.null(got)) {
      cr_log("falling back to serial")
      for (i in seq_len(nrow(todo))) fit_pair(todo$contrast[i],
                                              todo$outcome[i])
    }
  } else {
    for (i in seq_len(nrow(todo))) fit_pair(todo$contrast[i],
                                            todo$outcome[i])
  }
}

# Gather everything estimated so far (all checkpoints, not only this run's).
all_ck <- list.files(OUT, pattern = "^ladder_.*\\.rds$", full.names = TRUE)
rows <- lapply(all_ck, readRDS)
est <- do.call(rbind, rows)
utils::write.csv(est, file.path(OUT, "ladder_estimates.csv"),
                 row.names = FALSE)
cr_log(paste("Wrote ladder_estimates.csv:", nrow(est), "rows from",
             length(all_ck), "pairs"))

logs <- lapply(all_ck, function(f) {
  x <- readRDS(f)
  dl <- attr(x, "design_log")
  if (is.null(dl)) return(NULL)
  cbind(contrast = x$contrast[1], outcome = x$outcome[1], dl)
})
logs <- do.call(rbind, Filter(Negate(is.null), logs))
if (!is.null(logs))
  utils::write.csv(logs, file.path(OUT, "ladder_design_log.csv"),
                   row.names = FALSE)
cr_log("Multi-arm Stage 4 complete.")
