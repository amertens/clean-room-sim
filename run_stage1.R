#!/usr/bin/env Rscript
# ============================================================================
# Run Stage 1: Cohort Feasibility and Build
# ============================================================================
# Usage: Rscript run_stage1.R [config_path] [scenario]
# Example: Rscript run_stage1.R config/default.yml simple
# ============================================================================

r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE,
                      full.names = TRUE)
for (f in r_files) source(f, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else NULL
scenario    <- if (length(args) >= 2) args[2] else "simple"

cfg <- load_config(config_path)
params <- get_scenario_params(cfg, scenario)

message("=== Stage 1: Cohort Build (scenario: ", scenario, ") ===")
set_deterministic_seed(params$seed)

data <- do.call(generate_hcv_data, c(list(N = params$N), params))
result <- stage1_build_cohort(data = data, cfg = cfg)

message("Stage 1 complete. Checkpoint: ", result$checkpoint)
