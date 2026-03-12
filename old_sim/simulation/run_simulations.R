#!/usr/bin/env Rscript
# ============================================================================
# Simulation Study: Comparing TMLE, IPTW, G-computation, and Cox PH
# ============================================================================
# Runs across five scenarios, computing bias, RMSE, coverage, and runtime.
#
# Usage:
#   Rscript simulation/run_simulations.R [config_path]
#   Rscript simulation/run_simulations.R config/default.yml
# ============================================================================

# Source all R modules
source_all <- function() {
  r_files <- c(
    "R/utils/config.R",
    "R/utils/helpers.R",
    "R/dgp/DGP.R",
    "R/tmle/tmle_pipeline.R",
    "R/estimators/cox_ph_estimator.R",
    "R/estimators/iptw_survival.R",
    "R/estimators/gcomp_risk.R",
    "R/staging/decision_log.R",
    "R/staging/checkpoints.R"
  )
  for (f in r_files) {
    if (file.exists(f)) source(f, local = FALSE)
  }
}
source_all()

# --------------------------------------------------------------------------
# 1. Load configuration
# --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) > 0) args[1] else NULL
cfg <- load_config(config_path)

n_reps     <- cfg$simulation$n_replicates
N_sim      <- cfg$simulation$N
time_pts   <- cfg$simulation$time_points
output_dir <- cfg$output$simulation_dir
ensure_dir(output_dir)

message("============================================================")
message("Simulation Study Configuration")
message("  Replicates:   ", n_reps)
message("  Sample size:  ", N_sim)
message("  Time points:  ", paste(time_pts, collapse = ", "))
message("  Scenarios:    ", paste(names(cfg$scenarios), collapse = ", "))
message("============================================================")

# --------------------------------------------------------------------------
# 2. Compute truth for each scenario
# --------------------------------------------------------------------------
message("\n--- Computing ground-truth risks ---")

# Use a large sample to compute truth (Monte Carlo)
N_truth <- 200000
truths <- list()

for (scen_name in names(cfg$scenarios)) {
  message("  Truth for scenario: ", scen_name)
  scen_params <- cfg$scenarios[[scen_name]]

  # Build parameter list for DGP (exclude toggle params from base)
  dgp_args <- list(
    N = N_truth,
    h0 = cfg$dgp$h0,
    HR_early = cfg$dgp$HR_early,
    HR_late = cfg$dgp$HR_late,
    tau = cfg$dgp$tau,
    max_follow = cfg$dgp$max_follow,
    risk_window = cfg$dgp$risk_window,
    np_hazard = scen_params$np_hazard,
    dep_censor = scen_params$dep_censor,
    complexity = scen_params$complexity,
    switch_on = scen_params$switch_on,
    seed = 9999
  )

  truth_t <- list()
  for (t_val in time_pts) {
    tr <- compute_true_risk(
      N_truth = N_truth, t_eval = t_val, seed = 9999,
      h0 = dgp_args$h0, HR_early = dgp_args$HR_early,
      HR_late = dgp_args$HR_late, tau = dgp_args$tau,
      max_follow = dgp_args$max_follow,
      risk_window = dgp_args$risk_window,
      np_hazard = dgp_args$np_hazard,
      dep_censor = dgp_args$dep_censor,
      complexity = dgp_args$complexity,
      switch_on = dgp_args$switch_on
    )
    truth_t[[paste0("t", t_val)]] <- tr
  }
  truths[[scen_name]] <- truth_t
}

# Save truths
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
# 3. Run simulations
# --------------------------------------------------------------------------
all_results <- list()

for (scen_name in names(cfg$scenarios)) {
  message("=== Scenario: ", scen_name, " ===")
  scen_params <- cfg$scenarios[[scen_name]]

  for (rep_i in seq_len(n_reps)) {
    rep_seed <- cfg$dgp$seed + (which(names(cfg$scenarios) == scen_name) - 1) *
      10000 + rep_i
    set_deterministic_seed(rep_seed)

    if (rep_i %% 50 == 1 || rep_i == n_reps) {
      message("  Replicate ", rep_i, "/", n_reps)
    }

    # Generate data
    sim_data <- tryCatch({
      generate_hcv_data(
        N = N_sim,
        h0 = cfg$dgp$h0, HR_early = cfg$dgp$HR_early,
        HR_late = cfg$dgp$HR_late, tau = cfg$dgp$tau,
        max_follow = cfg$dgp$max_follow,
        risk_window = cfg$dgp$risk_window,
        np_hazard = scen_params$np_hazard,
        dep_censor = scen_params$dep_censor,
        complexity = scen_params$complexity,
        switch_on = scen_params$switch_on,
        seed = NULL  # seed already set
      )
    }, error = function(e) {
      message("    DGP failed: ", e$message)
      NULL
    })

    if (is.null(sim_data) || nrow(sim_data) < 100) next

    for (t_val in time_pts) {
      truth <- truths[[scen_name]][[paste0("t", t_val)]]

      # --- TMLE ---
      tmle_time <- system.time({
        tmle_res <- tryCatch({
          tmle_survival_risk(
            data = sim_data, t_eval = t_val,
            sl_lib_Q = c("SL.glm", "SL.mean"),
            sl_lib_g = c("SL.glm", "SL.mean"),
            sl_lib_cens = c("SL.glm", "SL.mean"),
            g_trunc = c(cfg$tmle$truncation_lower,
                        cfg$tmle$truncation_upper)
          )
        }, error = function(e) {
          list(RD = NA_real_, se_RD = NA_real_, ci_RD = c(NA_real_, NA_real_))
        })
      })

      # --- IPTW ---
      iptw_time <- system.time({
        iptw_res <- tryCatch({
          iptw_survival(
            data = sim_data, t_eval = t_val,
            sl_lib = c("SL.glm", "SL.mean"),
            g_trunc = c(cfg$tmle$truncation_lower,
                        cfg$tmle$truncation_upper)
          )
        }, error = function(e) {
          list(RD = NA_real_, se_RD = NA_real_, ci_RD = c(NA_real_, NA_real_))
        })
      })

      # --- G-computation ---
      gcomp_time <- system.time({
        gcomp_res <- tryCatch({
          gcomp_risk(data = sim_data, t_eval = t_val,
                     sl_lib = c("SL.glm", "SL.mean"))
        }, error = function(e) {
          list(RD = NA_real_, se_RD = NA_real_, ci_RD = c(NA_real_, NA_real_))
        })
      })

      # --- Cox PH ---
      cox_time <- system.time({
        cox_res <- tryCatch({
          cr <- cox_ph_estimator(sim_data)
          cox_risk <- cox_risk_at_t(cr, sim_data, t_eval = t_val)
          # Use delta-method SE from log(HR) -> approximate RD SE
          list(RD = cox_risk$RD, se_RD = NA_real_,
               ci_RD = c(NA_real_, NA_real_),
               HR = cr$HR, ci_HR = cr$ci_HR,
               ph_pvalue = cr$ph_test$pvalue,
               ph_violated = cr$ph_test$ph_violated)
        }, error = function(e) {
          list(RD = NA_real_, se_RD = NA_real_, ci_RD = c(NA_real_, NA_real_),
               HR = NA_real_, ph_pvalue = NA_real_, ph_violated = NA)
        })
      })

      # Collect results
      add_result <- function(method, res, runtime) {
        rd <- res$RD %||% NA_real_
        se <- res$se_RD %||% NA_real_
        ci <- res$ci_RD %||% c(NA_real_, NA_real_)
        covers <- if (!is.na(ci[1]) && !is.na(ci[2]) && !is.na(truth$RD)) {
          as.integer(ci[1] <= truth$RD & truth$RD <= ci[2])
        } else NA_integer_

        data.frame(
          scenario   = scen_name,
          replicate  = rep_i,
          time_point = t_val,
          method     = method,
          estimate   = rd,
          se         = se,
          ci_lower   = ci[1],
          ci_upper   = ci[2],
          truth      = truth$RD,
          bias       = rd - truth$RD,
          covers     = covers,
          runtime    = runtime[["elapsed"]],
          stringsAsFactors = FALSE
        )
      }

      all_results <- c(all_results, list(
        add_result("TMLE",          tmle_res,  tmle_time),
        add_result("IPTW",          iptw_res,  iptw_time),
        add_result("G-computation", gcomp_res, gcomp_time),
        add_result("Cox PH",        cox_res,   cox_time)
      ))
    }
  }
}

# --------------------------------------------------------------------------
# 4. Compile results
# --------------------------------------------------------------------------
message("\n--- Compiling results ---")
sim_results <- do.call(rbind, all_results)
utils::write.csv(sim_results,
                 file.path(output_dir, "simulation_results.csv"),
                 row.names = FALSE)

# --------------------------------------------------------------------------
# 5. Compute summary statistics
# --------------------------------------------------------------------------
compute_summary <- function(df) {
  # Remove NAs
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

# Group by scenario, time_point, method
groups <- unique(sim_results[, c("scenario", "time_point", "method")])
summary_list <- lapply(seq_len(nrow(groups)), function(i) {
  g <- groups[i, ]
  subset_df <- sim_results[sim_results$scenario == g$scenario &
                             sim_results$time_point == g$time_point &
                             sim_results$method == g$method, ]
  summ <- compute_summary(subset_df)
  cbind(g, summ)
})
summary_df <- do.call(rbind, summary_list)

# Add truth column
summary_df <- merge(summary_df, truth_df[, c("scenario", "time_point",
                                              "true_RD")],
                    by = c("scenario", "time_point"), all.x = TRUE)

# Round for display
num_cols <- c("mean_estimate", "mean_bias", "empirical_sd", "mean_se",
              "rmse", "coverage", "mean_runtime", "true_RD")
for (col in num_cols) {
  if (col %in% names(summary_df)) {
    summary_df[[col]] <- round(summary_df[[col]], 6)
  }
}

utils::write.csv(summary_df,
                 file.path(output_dir, "simulation_summary.csv"),
                 row.names = FALSE)

message("\nSimulation complete!")
message("  Results:  ", file.path(output_dir, "simulation_results.csv"))
message("  Summary:  ", file.path(output_dir, "simulation_summary.csv"))
message("  Truths:   ", file.path(output_dir, "true_values.csv"))

# --------------------------------------------------------------------------
# 6. Session info
# --------------------------------------------------------------------------
capture_session_info(file.path(output_dir, "session_info.txt"))

# Null-coalescing (ensure available)
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x
