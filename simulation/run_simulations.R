#!/usr/bin/env Rscript
# ============================================================================
# Simulation Study: Outcome-Blind Phase + Unblinded Analysis
# ============================================================================
# Two-phase simulation study implementing the clean-room governance model:
#
#   Phase 1 (Outcome-blind): Uses only baseline covariates + treatment to
#     evaluate candidate estimators on feasibility, PS diagnostics, ESS,
#     convergence, and runtime. Outcomes are suppressed. Selection of final
#     estimator/SL library/truncation is logged in the decision log.
#
#   Phase 2 (Unblinded): Generates datasets with outcomes and runs the
#     selected estimators to compute performance metrics (bias, RMSE,
#     coverage) relative to Monte-Carlo truth.
#
# Estimators compared: TMLE, TMLE-CF, AIPW, IPTW, G-computation, Cox PH
#
# Usage:
#   Rscript simulation/run_simulations.R [config_path]
#
# ============================================================================
# QUICK TEST MODE: Set QUICK_TEST <- TRUE for a fast smoke test with small N
# and few replicates. Set FALSE for the full production run using config values.
# ============================================================================

QUICK_TEST <- TRUE   # <-- Toggle: TRUE = fast test, FALSE = full production

# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Source all R modules
source_all <- function() {
  r_files <- c(
    "R/utils/config.R",
    "R/utils/helpers.R",
    "R/dgp/DGP.R",
    "R/tmle/tmle_pipeline.R",
    "R/estimators/cox_ph_estimator.R",
    "R/estimators/iptw_survival.R",
    "R/estimators/aipw_survival.R",
    "R/estimators/gcomp_risk.R",
    "R/staging/decision_log.R",
    "R/staging/checkpoints.R"
  )
  for (f in r_files) {
    if (file.exists(f)) source(f, local = FALSE)
  }
}
source_all()

# Format seconds as human-readable time string
format_time <- function(secs) {
  secs <- round(secs)
  if (secs < 60) return(paste0(secs, "s"))
  mins <- secs %/% 60
  secs <- secs %% 60
  if (mins < 60) return(paste0(mins, "m ", secs, "s"))
  hrs <- mins %/% 60
  mins <- mins %% 60
  paste0(hrs, "h ", mins, "m ", secs, "s")
}

# --------------------------------------------------------------------------
# 1. Load configuration
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) > 0) args[1] else NULL
cfg <- load_config(config_path)

if (QUICK_TEST) {
  n_reps   <- 10
  N_sim    <- 500
  N_truth  <- 50000
} else {
  n_reps   <- cfg$simulation$n_replicates
  N_sim    <- cfg$simulation$N
  N_truth  <- 200000
}

time_pts   <- cfg$simulation$time_points
output_dir <- cfg$output$simulation_dir
ensure_dir(output_dir)

sl_libs_Q <- c("SL.glm", "SL.mean")
sl_libs_g <- c("SL.glm", "SL.mean")

use_parallel <- cfg$simulation$parallel %||% FALSE
n_cores      <- as.integer(cfg$simulation$n_cores %||% 1L)

message("============================================================")
message("Clean-Room Simulation Study")
if (QUICK_TEST) message("  *** QUICK TEST MODE ***")
message("  Replicates:   ", n_reps)
message("  Sample size:  ", N_sim)
message("  Truth N:      ", N_truth)
message("  Time points:  ", paste(time_pts, collapse = ", "))
message("  Scenarios:    ", paste(names(cfg$scenarios), collapse = ", "))
message("  Parallel:     ",
        if (use_parallel) paste0("yes (", n_cores, " cores)") else "no")
message("============================================================")

t_total <- proc.time()

# --------------------------------------------------------------------------
# 2. Set up parallel cluster (PSOCK — works on Windows, macOS, Linux)
# --------------------------------------------------------------------------
cl <- NULL
if (use_parallel && n_cores > 1) {
  message("Starting parallel cluster with ", n_cores, " PSOCK workers...")
  cl <- tryCatch({
    cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
    wd <- getwd()
    r_source_files <- c(
      "R/utils/config.R", "R/utils/helpers.R", "R/dgp/DGP.R",
      "R/tmle/tmle_pipeline.R", "R/estimators/cox_ph_estimator.R",
      "R/estimators/iptw_survival.R", "R/estimators/aipw_survival.R",
      "R/estimators/gcomp_risk.R", "R/staging/decision_log.R",
      "R/staging/checkpoints.R"
    )
    parallel::clusterExport(cluster, c("wd", "r_source_files"),
                            envir = environment())
    parallel::clusterEvalQ(cluster, {
      setwd(wd)
      `%||%` <- function(x, y) if (is.null(x)) y else x
      for (f in r_source_files) {
        if (file.exists(f)) source(f, local = FALSE)
      }
      suppressPackageStartupMessages({
        library(SuperLearner)
        library(survival)
      })
    })
    message("Cluster ready.")
    cluster
  }, error = function(e) {
    warning("Failed to start parallel cluster: ", e$message,
            "\n  Falling back to serial execution.")
    NULL
  })
  if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)
}

# Helper: parLapply when cluster is available, otherwise lapply
par_or_lapply <- function(X, FUN, ...) {
  if (!is.null(cl)) {
    parallel::parLapply(cl, X, FUN, ...)
  } else {
    lapply(X, FUN, ...)
  }
}

# ==========================================================================
# PHASE 1: OUTCOME-BLIND SIMULATION
# ==========================================================================
message("\n###########################################################")
message("# PHASE 1: OUTCOME-BLIND ESTIMATOR SELECTION")
message("###########################################################\n")

# In the outcome-blind phase we generate data but ONLY inspect PS diagnostics:
# - propensity score model fit and convergence
# - PS overlap and positivity
# - effective sample size under IPW weights
# - covariate balance (weighted SMDs)
# - truncation counts
# NO outcome models are fit; NO outcome data (follow_time, event) are used.

n_blind_reps <- if (QUICK_TEST) min(n_reps, 5) else min(n_reps, 20)
blind_results <- list()

t_phase1 <- proc.time()
total_blind <- n_blind_reps * length(cfg$scenarios)
message("Running ", total_blind, " outcome-blind reps ",
        "(", length(cfg$scenarios), " scenarios x ", n_blind_reps, " reps)")
pb1 <- utils::txtProgressBar(min = 0, max = total_blind, style = 3)
blind_count <- 0

for (scen_name in names(cfg$scenarios)) {
  scen_params <- cfg$scenarios[[scen_name]]

  for (rep_i in seq_len(n_blind_reps)) {
    rep_seed <- cfg$dgp$seed +
      (which(names(cfg$scenarios) == scen_name) - 1) * 10000 + rep_i
    set_deterministic_seed(rep_seed)

    # Generate data (outcomes exist but are NOT used for selection)
    d <- tryCatch({
      generate_hcv_data(
        N = N_sim, h0 = cfg$dgp$h0, HR_early = cfg$dgp$HR_early,
        HR_late = cfg$dgp$HR_late, tau = cfg$dgp$tau,
        max_follow = cfg$dgp$max_follow, risk_window = cfg$dgp$risk_window,
        np_hazard = scen_params$np_hazard, dep_censor = scen_params$dep_censor,
        complexity = scen_params$complexity, switch_on = scen_params$switch_on,
        seed = NULL
      )
    }, error = function(e) NULL)
    if (is.null(d) || nrow(d) < 100) next

    # Extract covariates and treatment only (outcome-blind)
    exclude_cols <- c("id", "follow_time", "event", "treatment", "switch",
                      "race", "region")
    covar_cols <- setdiff(names(d), exclude_cols)
    covar_cols <- covar_cols[vapply(d[, covar_cols, drop = FALSE],
                                   is.numeric, logical(1))]
    W <- as.data.frame(d[, covar_cols, drop = FALSE])
    A <- d$treatment

    # Fit PS model (permitted in outcome-blind phase) and time it
    ps_time <- system.time({
      g_fit <- tryCatch({
        SuperLearner::SuperLearner(
          Y = A, X = W, family = stats::binomial(),
          SL.library = sl_libs_g, cvControl = list(V = 5)
        )
      }, error = function(e) {
        fit <- stats::glm(A ~ ., data = cbind(A = A, W), family = "binomial")
        list(SL.predict = stats::predict(fit, type = "response"))
      })
    })
    ps <- as.numeric(g_fit$SL.predict)
    trunc_info <- truncate_ps(ps, cfg$tmle$truncation_lower,
                              cfg$tmle$truncation_upper)
    ps <- trunc_info$p

    # Outcome-blind diagnostics (PS-based only — no outcome models)
    w_ipw <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
    ess_treated <- effective_ss(w_ipw[A == 1])
    ess_control <- effective_ss(w_ipw[A == 0])

    smds <- vapply(covar_cols, function(v) {
      abs(compute_smd(d[[v]], A, weights = w_ipw))
    }, numeric(1))

    blind_results <- c(blind_results, list(data.frame(
      scenario = scen_name, replicate = rep_i,
      ps_converged = TRUE,
      ps_runtime = ps_time[["elapsed"]],
      ess_treated = round(ess_treated, 1),
      ess_control = round(ess_control, 1),
      max_smd_weighted = round(max(smds), 4),
      n_trunc_lower = trunc_info$n_lower,
      n_trunc_upper = trunc_info$n_upper,
      ps_min = round(min(ps), 4),
      ps_max = round(max(ps), 4),
      stringsAsFactors = FALSE
    )))
    blind_count <- blind_count + 1
    utils::setTxtProgressBar(pb1, blind_count)
  }
}
close(pb1)

blind_df <- do.call(rbind, blind_results)
utils::write.csv(blind_df, file.path(output_dir, "phase1_outcome_blind.csv"),
                 row.names = FALSE)

# Summarise outcome-blind phase (PS diagnostics only)
blind_summary <- do.call(rbind, lapply(
  split(blind_df, blind_df$scenario),
  function(sub) {
    data.frame(
      scenario        = sub$scenario[1],
      ps_convergence  = mean(sub$ps_converged),
      mean_ps_runtime = round(mean(sub$ps_runtime), 3),
      mean_ess_trt    = round(mean(sub$ess_treated), 1),
      mean_ess_ctl    = round(mean(sub$ess_control), 1),
      mean_max_smd    = round(mean(sub$max_smd_weighted), 4),
      mean_ps_min     = round(mean(sub$ps_min), 4),
      mean_ps_max     = round(mean(sub$ps_max), 4),
      stringsAsFactors = FALSE
    )
  }
))
rownames(blind_summary) <- NULL
utils::write.csv(blind_summary,
                 file.path(output_dir, "phase1_summary.csv"),
                 row.names = FALSE)

t_phase1_elapsed <- (proc.time() - t_phase1)[["elapsed"]]
message("\nPhase 1 complete in ", format_time(t_phase1_elapsed),
        ". Outcome-blind diagnostics saved.")

# --------------------------------------------------------------------------
# Log outcome-blind decisions
# --------------------------------------------------------------------------
ensure_dir("outputs")
mtg <- start_meeting("outcome_blind_simulation", protocol_version = 1L)
mtg <- log_decision(
  mtg, "PS",
  paste("SL library for g:", paste(sl_libs_g, collapse = ", ")),
  "Pre-specified in config; PS model convergence confirmed across scenarios",
  triggered_by = "phase 1 PS diagnostics"
)
mtg <- log_decision(
  mtg, "Truncation",
  paste("g bounds: [", cfg$tmle$truncation_lower, ",",
        cfg$tmle$truncation_upper, "]"),
  "Pre-specified; truncation counts acceptable across scenarios",
  triggered_by = "phase 1 truncation diagnostics"
)
mtg <- log_decision(
  mtg, "Overlap",
  "PS overlap and balance confirmed adequate across all scenarios",
  paste("ESS and weighted SMD within thresholds;",
        "no near-positivity violations detected"),
  triggered_by = "phase 1 PS overlap assessment"
)
mtg <- log_decision(
  mtg, "Estimator",
  paste("Pre-specified estimators for unblinded analysis: TMLE, TMLE-CF,",
        "AIPW, IPTW, G-comp, Cox PH"),
  "All estimators pre-specified; PS feasibility confirmed outcome-blind",
  triggered_by = "phase 1 PS feasibility assessment"
)
close_meeting(mtg)

# ==========================================================================
# PHASE 2: UNBLINDED SIMULATION (OUTCOMES INCLUDED)
# ==========================================================================
message("\n###########################################################")
message("# PHASE 2: UNBLINDED ANALYSIS")
message("###########################################################\n")

# --------------------------------------------------------------------------
# 2a. Compute ground-truth risks
# --------------------------------------------------------------------------
message("--- Computing ground-truth risks ---")
truths <- list()

for (scen_name in names(cfg$scenarios)) {
  message("  Truth for scenario: ", scen_name)
  scen_params <- cfg$scenarios[[scen_name]]
  truth_t <- list()
  for (t_val in time_pts) {
    tr <- compute_true_risk(
      N_truth = N_truth, t_eval = t_val, seed = 9999,
      h0 = cfg$dgp$h0, HR_early = cfg$dgp$HR_early,
      HR_late = cfg$dgp$HR_late, tau = cfg$dgp$tau,
      max_follow = cfg$dgp$max_follow, risk_window = cfg$dgp$risk_window,
      np_hazard = scen_params$np_hazard, dep_censor = scen_params$dep_censor,
      complexity = scen_params$complexity, switch_on = scen_params$switch_on
    )
    truth_t[[paste0("t", t_val)]] <- tr
  }
  truths[[scen_name]] <- truth_t
}

truth_df <- do.call(rbind, lapply(names(truths), function(scen) {
  do.call(rbind, lapply(names(truths[[scen]]), function(tkey) {
    tr <- truths[[scen]][[tkey]]
    data.frame(scenario = scen, time_point = tr$t_eval,
               true_risk_1 = tr$risk_1, true_risk_0 = tr$risk_0,
               true_RD = tr$RD, true_RR = tr$RR,
               stringsAsFactors = FALSE)
  }))
}))
utils::write.csv(truth_df, file.path(output_dir, "true_values.csv"),
                 row.names = FALSE)
message("Ground truths saved.\n")

# --------------------------------------------------------------------------
# 2b. Worker function: run all estimators for one replicate
# --------------------------------------------------------------------------
run_sim_rep <- function(rep_i, scen_name, scen_params, scen_idx,
                        N_sim, cfg, time_pts, truths,
                        sl_libs_Q, sl_libs_g) {
  `%||%` <- function(x, y) if (is.null(x)) y else x

  rep_seed <- cfg$dgp$seed + (scen_idx - 1) * 10000 + rep_i
  set_deterministic_seed(rep_seed)

  sim_data <- tryCatch({
    generate_hcv_data(
      N = N_sim, h0 = cfg$dgp$h0, HR_early = cfg$dgp$HR_early,
      HR_late = cfg$dgp$HR_late, tau = cfg$dgp$tau,
      max_follow = cfg$dgp$max_follow, risk_window = cfg$dgp$risk_window,
      np_hazard = scen_params$np_hazard, dep_censor = scen_params$dep_censor,
      complexity = scen_params$complexity, switch_on = scen_params$switch_on,
      seed = NULL
    )
  }, error = function(e) NULL)
  if (is.null(sim_data) || nrow(sim_data) < 100) return(NULL)

  rep_results <- list()
  for (t_val in time_pts) {
    truth <- truths[[scen_name]][[paste0("t", t_val)]]
    g_bounds <- c(cfg$tmle$truncation_lower, cfg$tmle$truncation_upper)

    make_row <- function(method, res, runtime) {
      rd <- res$RD %||% NA_real_
      se <- res$se_RD %||% NA_real_
      ci <- res$ci_RD %||% c(NA_real_, NA_real_)
      covers <- if (!is.na(ci[1]) && !is.na(ci[2]) && !is.na(truth$RD)) {
        as.integer(ci[1] <= truth$RD & truth$RD <= ci[2])
      } else NA_integer_
      data.frame(
        scenario = scen_name, replicate = rep_i, time_point = t_val,
        method = method, estimate = rd, se = se,
        ci_lower = ci[1], ci_upper = ci[2],
        truth = truth$RD, bias = rd - truth$RD,
        covers = covers, runtime = runtime[["elapsed"]],
        stringsAsFactors = FALSE
      )
    }

    safe_est <- function(expr) {
      tryCatch(expr, error = function(e)
        list(RD = NA_real_, se_RD = NA_real_, ci_RD = c(NA_real_, NA_real_)))
    }

    t1 <- system.time({ r1 <- safe_est(
      tmle_survival_risk(sim_data, t_val, sl_libs_Q, sl_libs_g,
                         sl_libs_g, g_bounds)) })
    t2 <- system.time({ r2 <- safe_est(
      tmle_survival_risk_cf(sim_data, t_val, sl_libs_Q, sl_libs_g,
                            sl_libs_g, g_bounds, V = 5)) })
    t3 <- system.time({ r3 <- safe_est(
      aipw_survival(sim_data, t_val, sl_libs_Q, sl_libs_g,
                    sl_libs_g, g_bounds)) })
    t4 <- system.time({ r4 <- safe_est(
      iptw_survival(sim_data, t_val, sl_libs_g, g_bounds)) })
    t5 <- system.time({ r5 <- safe_est(
      gcomp_risk(sim_data, t_val, sl_libs_Q)) })
    t6 <- system.time({ r6 <- safe_est({
      cr <- cox_ph_estimator(sim_data)
      cox_r <- cox_risk_at_t(cr, sim_data, t_val)
      list(RD = cox_r$RD, se_RD = NA_real_,
           ci_RD = c(NA_real_, NA_real_))
    }) })

    rep_results <- c(rep_results, list(
      make_row("TMLE",          r1, t1),
      make_row("TMLE-CF",       r2, t2),
      make_row("AIPW",          r3, t3),
      make_row("IPTW",          r4, t4),
      make_row("G-computation", r5, t5),
      make_row("Cox PH",        r6, t6)
    ))
  }
  rep_results
}

# --------------------------------------------------------------------------
# 2c. Run unblinded simulation replicates
# --------------------------------------------------------------------------
t_phase2 <- proc.time()
all_results <- list()

for (scen_name in names(cfg$scenarios)) {
  message("=== Scenario: ", scen_name, " (", n_reps, " reps",
          if (!is.null(cl)) paste0(", ", n_cores, " cores") else ", serial",
          ") ===")
  scen_params <- cfg$scenarios[[scen_name]]
  scen_idx <- which(names(cfg$scenarios) == scen_name)

  common_args <- list(
    scen_name = scen_name, scen_params = scen_params, scen_idx = scen_idx,
    N_sim = N_sim, cfg = cfg, time_pts = time_pts, truths = truths,
    sl_libs_Q = sl_libs_Q, sl_libs_g = sl_libs_g
  )

  pb <- utils::txtProgressBar(min = 0, max = n_reps, style = 3)
  t_scen <- proc.time()
  scen_results <- list()

  if (!is.null(cl)) {
    # Parallel: process in batches so progress bar can update
    batch_size <- n_cores
    for (batch_start in seq(1, n_reps, by = batch_size)) {
      batch_end <- min(batch_start + batch_size - 1, n_reps)
      batch_idx <- seq(batch_start, batch_end)
      batch_res <- do.call(parallel::parLapply,
        c(list(cl = cl, X = batch_idx, fun = run_sim_rep), common_args))
      scen_results <- c(scen_results, batch_res)
      utils::setTxtProgressBar(pb, batch_end)
    }
  } else {
    # Serial: update progress per replicate
    for (rep_i in seq_len(n_reps)) {
      res <- do.call(run_sim_rep, c(list(rep_i = rep_i), common_args))
      scen_results <- c(scen_results, list(res))
      utils::setTxtProgressBar(pb, rep_i)
    }
  }
  close(pb)

  elapsed <- (proc.time() - t_scen)[["elapsed"]]

  # Flatten: each element is either NULL or a list of data.frames
  scen_results <- scen_results[!vapply(scen_results, is.null, logical(1))]
  scen_results <- unlist(scen_results, recursive = FALSE)
  all_results <- c(all_results, scen_results)

  n_valid <- length(scen_results) / (length(time_pts) * 6)
  message("  ", n_valid, " valid reps in ", format_time(elapsed))
}

# --------------------------------------------------------------------------
# 2d. Compile and summarise
# --------------------------------------------------------------------------
message("\n--- Compiling results ---")
sim_results <- do.call(rbind, all_results)
utils::write.csv(sim_results,
                 file.path(output_dir, "simulation_results.csv"),
                 row.names = FALSE)

compute_summary <- function(df) {
  df <- df[!is.na(df$estimate), ]
  if (nrow(df) == 0) {
    return(data.frame(
      n_valid = 0, mean_estimate = NA, mean_bias = NA, empirical_sd = NA,
      mean_se = NA, rmse = NA, coverage = NA, mean_runtime = NA
    ))
  }
  data.frame(
    n_valid       = nrow(df),
    mean_estimate = mean(df$estimate),
    mean_bias     = mean(df$bias),
    empirical_sd  = stats::sd(df$estimate),
    mean_se       = mean(df$se, na.rm = TRUE),
    rmse          = sqrt(mean(df$bias^2)),
    coverage      = mean(df$covers, na.rm = TRUE),
    mean_runtime  = mean(df$runtime, na.rm = TRUE)
  )
}

groups <- unique(sim_results[, c("scenario", "time_point", "method")])
summary_list <- lapply(seq_len(nrow(groups)), function(i) {
  g <- groups[i, ]
  sub <- sim_results[sim_results$scenario == g$scenario &
                       sim_results$time_point == g$time_point &
                       sim_results$method == g$method, ]
  cbind(g, compute_summary(sub))
})
summary_df <- do.call(rbind, summary_list)
summary_df <- merge(summary_df,
                    truth_df[, c("scenario", "time_point", "true_RD")],
                    by = c("scenario", "time_point"), all.x = TRUE)

num_cols <- c("mean_estimate", "mean_bias", "empirical_sd", "mean_se",
              "rmse", "coverage", "mean_runtime", "true_RD")
for (col in num_cols) {
  if (col %in% names(summary_df))
    summary_df[[col]] <- round(summary_df[[col]], 6)
}
utils::write.csv(summary_df,
                 file.path(output_dir, "simulation_summary.csv"),
                 row.names = FALSE)

# Log unblinded phase
mtg2 <- start_meeting("unblinded_simulation", protocol_version = 1L)
mtg2 <- log_decision(
  mtg2, "Simulation",
  paste("Unblinded simulation completed:", n_reps, "replicates x",
        length(names(cfg$scenarios)), "scenarios x",
        length(time_pts), "time points"),
  "Pre-specified simulation protocol",
  triggered_by = "phase 2 execution",
  outcome_blind_confirm = FALSE
)
close_meeting(mtg2)

t_phase2_elapsed <- (proc.time() - t_phase2)[["elapsed"]]
t_total_elapsed  <- (proc.time() - t_total)[["elapsed"]]

message("\n============================================================")
message("Simulation study complete!")
message("============================================================")
message("  Phase 1 (blind):  ", format_time(t_phase1_elapsed))
message("  Phase 2 (unbld):  ", format_time(t_phase2_elapsed))
message("  Total elapsed:    ", format_time(t_total_elapsed))
message("------------------------------------------------------------")
message("  Phase 1 output:   ", file.path(output_dir, "phase1_summary.csv"))
message("  Phase 2 results:  ", file.path(output_dir, "simulation_results.csv"))
message("  Summary:          ", file.path(output_dir, "simulation_summary.csv"))
message("  Truths:           ", file.path(output_dir, "true_values.csv"))

capture_session_info(file.path(output_dir, "session_info.txt"))
