# ============================================================
# Stage 6: Render Report & Decision Log
# ============================================================
# Compiles final decision log and renders the Quarto report.
# ============================================================

# --- Auto-detect layout (project root vs inside rescueCo/) ---
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))

# --- Source helpers ---
source("rescueCo/R/utils.R")

# --- Load config ---
cfg <- load_cr_config()
cr_log("=== Stage 6: Render Report ===")

# --- Compile decision log (loads from latest available stage) ---
decisions <- tryCatch(
  load_stage_output("stage5_decisions.rds"),
  error = function(e) tryCatch(
    load_stage_output("stage4_decisions.rds"),
    error = function(e2) tryCatch(
      load_stage_output("stage2b_decisions.rds"),
      error = function(e3) load_stage_output("stage2_decisions.rds")
    )
  )
)

# Save final decision log
write.csv(decisions, file.path(cfg$paths$results, "decision_log.csv"), row.names = FALSE)
cr_log(paste("Decision log saved:", nrow(decisions), "decisions recorded"))

# Print decision summary
cr_log("=== Decision Log Summary ===")
for (i in seq_len(nrow(decisions))) {
  cr_log(paste("[", decisions$type[i], "]",
               decisions$stage[i], ":",
               decisions$decision[i]))
}

# --- Export cleanTMLE audit trail ---
audit <- tryCatch(
  load_stage_output("stage5_audit.rds"),
  error = function(e) tryCatch(
    load_stage_output("stage4_audit.rds"),
    error = function(e2) tryCatch(
      load_stage_output("stage2b_audit.rds"),
      error = function(e3) NULL
    )
  )
)

if (!is.null(audit)) {
  cr_log("Exporting cleanTMLE audit trail...")
  audit_trail <- tryCatch(cleanTMLE::export_audit_trail(audit), error = function(e) NULL)
  if (!is.null(audit_trail)) {
    write.csv(audit_trail, file.path(cfg$paths$results, "cleanTMLE_audit_trail.csv"),
              row.names = FALSE)
    cr_log(paste("cleanTMLE audit trail:", nrow(audit_trail), "entries"))
  }
  audit_dlog <- tryCatch(cleanTMLE::export_decision_log(audit), error = function(e) NULL)
  if (!is.null(audit_dlog)) {
    write.csv(audit_dlog, file.path(cfg$paths$results, "cleanTMLE_decision_log.csv"),
              row.names = FALSE)
    cr_log(paste("cleanTMLE decision log:", nrow(audit_dlog), "entries"))
  }
  tryCatch(cleanroomGov::build_stage_manifest(audit),
           error = function(e) cr_log(paste("Stage manifest failed:", e$message)))
}

# --- Render Quarto report ---
report_path <- "rescueCo/reports/rescueco_clean_room_case_study.qmd"

if (file.exists(report_path)) {
  cr_log("Rendering Quarto report...")
  # Use system2() to render in a separate process — avoids temp-dir issues
  report_abs <- normalizePath(report_path, mustWork = TRUE)
  rc <- system2("quarto", args = c("render", shQuote(report_abs)),
                stdout = TRUE, stderr = TRUE)
  cat(rc, sep = "\n")
  # Check if HTML was produced
  report_html <- sub("\\.qmd$", ".html", report_path)
  if (file.exists(report_html)) {
    cr_log("Report rendered successfully.")
  } else {
    cr_log("Quarto render may have failed. Render manually:")
    cr_log(paste("  quarto render", report_path))
  }
} else {
  cr_log(paste("Report template not found at:", report_path))
}

# --- Final summary ---
results_dir <- cfg$paths$results
result_files <- list.files(results_dir, pattern = "\\.(rds|csv|png)$")

cr_log("=== Pipeline Complete ===")
cr_log(paste("Result files (", length(result_files), "):"))
for (f in result_files) {
  cr_log(paste("  ", f))
}

cr_log("Stage 6 complete.")
