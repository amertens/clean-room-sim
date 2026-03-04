#!/usr/bin/env Rscript
# ===========================================================================
# 00_run_all_stages.R
# Single entrypoint that runs Stage 1 -> Stage 2 -> Stage 3 in order,
# writing one RDS file per stage to outputs/.
#
# Usage:
#   Rscript scripts/00_run_all_stages.R               # defaults
#   Rscript scripts/00_run_all_stages.R --quick        # B=200 bootstrap
#   Rscript scripts/00_run_all_stages.R --config path  # custom config
# ===========================================================================

message("=== Clean-Room Staged Analysis Pipeline ===")
message("Start: ", Sys.time())

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args
config_path <- NULL
if ("--config" %in% args) {
  idx <- which(args == "--config")
  if (idx < length(args)) config_path <- args[idx + 1]
}

# Set working directory to project root
if (file.exists("scripts/00_run_all_stages.R")) {
  # Already in project root
} else if (file.exists("../scripts/00_run_all_stages.R")) {
  setwd("..")
}

# Source all R modules
r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE,
                      full.names = TRUE)
for (f in r_files) source(f)

# Load configuration
cfg <- load_config(config_path)

# Record seed
master_seed <- cfg$dgp$seed
message("Master seed: ", master_seed)

# Bootstrap replicates
B <- if (quick_mode) (cfg$bootstrap$B_quick %||% 200) else (cfg$bootstrap$B %||% 500)
message("Bootstrap replicates: ", B)

# Ensure output directories
ensure_dir("outputs")
ensure_dir("outputs/stage1")
ensure_dir("outputs/stage2")
ensure_dir("outputs/stage3")
ensure_dir("logs")

# Remove stale outputs so pipeline runs fresh
unlink("outputs/checkpoint_1.json")
unlink("outputs/checkpoint_2.json")
unlink("outputs/decision_log.csv")

# ===========================================================================
# Stage 1: Cohort Feasibility
# ===========================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STAGE 1: Cohort Feasibility and Build")
message(paste(rep("=", 60), collapse = ""))

set_deterministic_seed(master_seed)
cohort_data <- generate_hcv_data(
  N          = cfg$dgp$N,
  h0         = cfg$dgp$h0,
  HR_early   = cfg$dgp$HR_early,
  HR_late    = cfg$dgp$HR_late,
  tau        = cfg$dgp$tau,
  max_follow = cfg$dgp$max_follow,
  np_hazard  = FALSE,
  dep_censor = FALSE,
  complexity = FALSE,
  switch_on  = FALSE,
  seed       = master_seed
)

stage1_out <- stage1_build_cohort(
  data       = cohort_data,
  cfg        = cfg,
  output_dir = "outputs/stage1"
)

message("Stage 1 checkpoint: ", stage1_out$checkpoint)
message("  N = ", stage1_out$report$N_total,
        ", Events = ", stage1_out$report$overall_N_events,
        ", Event rate = ", round(stage1_out$report$overall_event_rate, 4))

if (stage1_out$checkpoint != "PASS") {
  message("Stage 1 FAILED. Stopping pipeline.")
  saveRDS(list(stage1 = stage1_out), "outputs/stage1_feasibility.rds")
  quit(status = 1)
}

saveRDS(stage1_out, "outputs/stage1_feasibility.rds")

# ===========================================================================
# Stage 2: Design Adequacy
# ===========================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STAGE 2: Design Adequacy Checks")
message(paste(rep("=", 60), collapse = ""))

stage2_out <- stage2_design_checks(
  cohort     = stage1_out$cohort,
  cfg        = cfg,
  output_dir = "outputs/stage2"
)

message("Stage 2 checkpoint: ", stage2_out$checkpoint)
message("  Max weighted SMD = ", stage2_out$diagnostics$max_smd_weighted,
        ", ESS fraction = ", stage2_out$diagnostics$ess_frac)

if (stage2_out$checkpoint != "PASS") {
  message("Stage 2 FAILED. Stopping pipeline.")
  saveRDS(stage2_out, "outputs/stage2_design_adequacy.rds")
  quit(status = 1)
}

saveRDS(stage2_out, "outputs/stage2_design_adequacy.rds")

# ===========================================================================
# Stage 3: Estimation
# ===========================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("STAGE 3: Outcome Modeling and Estimation")
message(paste(rep("=", 60), collapse = ""))

stage3_out <- stage3_estimation(
  cohort        = stage1_out$cohort,
  stage2_result = stage2_out,
  cfg           = cfg,
  output_dir    = "outputs/stage3",
  B             = B
)

message("\nStage 3 summary:")
print(stage3_out$summary)

saveRDS(stage3_out, "outputs/stage3_estimation.rds")

# ===========================================================================
# Save decision log as RDS
# ===========================================================================
if (file.exists("outputs/decision_log.csv")) {
  decision_log <- utils::read.csv("outputs/decision_log.csv",
                                   stringsAsFactors = FALSE)
  saveRDS(decision_log, "outputs/decision_log.rds")
}

# ===========================================================================
# Session info
# ===========================================================================
capture_session_info("logs/session_info.txt")

# Git hash if available
git_hash <- tryCatch(
  system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE),
  error = function(e) "unknown"
)
writeLines(c(
  paste("Pipeline completed:", Sys.time()),
  paste("Git commit:", git_hash),
  paste("Master seed:", master_seed),
  paste("Bootstrap B:", B),
  paste("Config:", config_path %||% "config/default.yml")
), "logs/pipeline_run.log")

message("\n=== Pipeline completed successfully ===")
message("End: ", Sys.time())
message("Outputs saved to outputs/")
