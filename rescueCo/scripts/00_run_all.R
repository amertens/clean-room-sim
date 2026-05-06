# ============================================================
# Clean-Room TMLE Workflow: Master Runner
# RescueCo Kenya Trauma Registry Case Study
# ============================================================
#
# This script runs the full staged clean-room pipeline.
# Each stage saves intermediate RDS outputs so individual
# stages can be re-run without recomputing everything.
#
# Usage:
#   source("rescueCo/scripts/00_run_all.R")
#
# Or run individual stages:
#   source("rescueCo/scripts/01_stage1_build_cohort.R")
#   source("rescueCo/scripts/02_stage2_design_diagnostics.R")
#   ... etc.
# ============================================================

cat("============================================================\n")
cat("  Clean-Room TMLE Workflow: RescueCo Case Study\n")
cat("  Started:", format(Sys.time()), "\n")
cat("============================================================\n\n")

# --- Auto-detect layout ---
# Sources bootstrap.R which sets the working directory to the parent of
# rescueCo/ regardless of how the script is invoked. After bootstrap,
# all subsequent paths use the "rescueCo/..." prefix unambiguously.
boot_path <- if (dir.exists("rescueCo")) "rescueCo/R/bootstrap.R" else
             if (dir.exists("R")) "R/bootstrap.R" else
             stop("Cannot locate bootstrap.R; run from project root or inside rescueCo/.")
source(boot_path)

if (!file.exists("rescueCo/config/clean_room_config.yml")) {
  stop("Bootstrap failed to locate rescueCo/config/clean_room_config.yml")
}

# --- Stage 1: Build Cohort ---
cat(">>> STAGE 1: Building analytic cohort...\n")
source("rescueCo/scripts/01_stage1_build_cohort.R", local = new.env())
cat(">>> Stage 1 complete.\n\n")

# --- Stage 1b: Exploratory Data Analysis (optional render) ---
# Skipped if quarto isn't on PATH — the EDA HTML is decorative and the
# manuscript artifacts come from later stages.
cat(">>> STAGE 1b: EDA render (skipped if quarto not on PATH)...\n")
eda_qmd <- "rescueCo/scripts/01b_stage1b_eda.qmd"
quarto_bin <- Sys.which("quarto")
if (nzchar(quarto_bin) && file.exists(eda_qmd)) {
  eda_rc <- tryCatch(
    system2(quarto_bin, args = c("render", shQuote(normalizePath(eda_qmd))),
             stdout = TRUE, stderr = TRUE),
    error = function(e) { cat("EDA render skipped:", e$message, "\n"); NULL })
  eda_src <- "rescueCo/scripts/01b_stage1b_eda.html"
  if (file.exists(eda_src)) {
    dir.create("rescueCo/reports", showWarnings = FALSE, recursive = TRUE)
    file.rename(eda_src, "rescueCo/reports/stage1b_eda.html")
  }
} else {
  cat(">>> Stage 1b EDA render skipped (quarto not on PATH).\n")
}
cat("\n")

# --- Stage 2: Design Diagnostics ---
cat(">>> STAGE 2: Design adequacy checks...\n")
source("rescueCo/scripts/02_stage2_design_diagnostics.R", local = new.env())
cat(">>> Stage 2 complete.\n\n")

# --- Stage 2b: Negative Control Analysis ---
cat(">>> STAGE 2b: Negative control analysis...\n")
source("rescueCo/scripts/02b_stage2b_negative_controls.R", local = new.env())
cat(">>> Stage 2b complete.\n\n")

# --- Stage 3: Outcome-Blind Simulation ---
cat(">>> STAGE 3: Outcome-blind simulation...\n")
source("rescueCo/scripts/03_stage3_outcome_blind_simulation.R", local = new.env())
cat(">>> Stage 3 complete.\n\n")

# --- Stage 4: Binary Outcome Analysis ---
cat(">>> STAGE 4: Binary outcome analysis (GOSE)...\n")
source("rescueCo/scripts/04_stage4_binary_outcome_analysis.R", local = new.env())
cat(">>> Stage 4 complete.\n\n")

# --- Stage 5: Survival Outcome Analysis ---
cat(">>> STAGE 5: Survival outcome analysis...\n")
source("rescueCo/scripts/05_stage5_survival_outcome_analysis.R", local = new.env())
cat(">>> Stage 5 complete.\n\n")

# --- Stage 6: Render Report ---
cat(">>> STAGE 6: Rendering report...\n")
source("rescueCo/scripts/06_stage6_render_report.R", local = new.env())
cat(">>> Stage 6 complete.\n\n")

cat("============================================================\n")
cat("  Pipeline finished:", format(Sys.time()), "\n")
cat("  Results: rescueCo/results/\n")
cat("  Report:  rescueCo/reports/\n")
cat("============================================================\n")
