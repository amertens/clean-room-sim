#!/usr/bin/env Rscript
# ============================================================================
# Run Stage 2: Design Adequacy Checks
# ============================================================================
# Usage: Rscript run_stage2.R [config_path] [scenario]
# Requires: Stage 1 checkpoint PASS
# ============================================================================

r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE,
                      full.names = TRUE)
for (f in r_files) source(f, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else NULL
scenario    <- if (length(args) >= 2) args[2] else "simple"

cfg <- load_config(config_path)
params <- get_scenario_params(cfg, scenario)

message("=== Stage 2: Design Checks (scenario: ", scenario, ") ===")
set_deterministic_seed(params$seed)

data <- do.call(generate_hcv_data, params)
result <- stage2_design_checks(cohort = data, cfg = cfg)

message("Stage 2 complete. Checkpoint: ", result$checkpoint)
