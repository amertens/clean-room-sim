#!/usr/bin/env Rscript
# ============================================================================
# Run Stage 3 + Stage 4: Estimation and Reporting
# ============================================================================
# Usage: Rscript run_stage3_report.R [config_path] [scenario]
# Requires: Stage 2 checkpoint PASS
# ============================================================================

r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE,
                      full.names = TRUE)
for (f in r_files) source(f, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else NULL
scenario    <- if (length(args) >= 2) args[2] else "simple"

cfg <- load_config(config_path)
params <- get_scenario_params(cfg, scenario)

# --- Stage 3 ---
message("=== Stage 3: Estimation (scenario: ", scenario, ") ===")
set_deterministic_seed(params$seed)

data <- do.call(generate_hcv_data, c(list(N = params$N), params))
stage3_out <- stage3_estimation(cohort = data, cfg = cfg)

cat("\n--- Results Summary ---\n")
print(stage3_out$summary)

# --- Stage 4 ---
message("\n=== Stage 4: Reporting ===")
stage4_report(cfg = cfg)

message("\nPipeline complete. Report at: ", cfg$output$report_dir)
