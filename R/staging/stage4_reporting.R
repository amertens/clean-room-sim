#' @title Stage 4: Reporting
#' @description Compiles stage artifacts into a single locked report.
#' @name stage4
NULL

#' Generate Stage 4 Report
#'
#' Compiles results from all stages into a structured report directory
#' with summary tables, decision log, and session info.
#'
#' @param cfg Config list from \code{load_config}.
#' @param output_dir Character path for report outputs.
#' @return Invisibly returns the report directory path.
#' @export
stage4_report <- function(cfg = NULL, output_dir = "outputs/report") {
  if (is.null(cfg)) cfg <- load_config()
  ensure_dir(output_dir)

  # Enforce prerequisite
  require_checkpoint_pass("checkpoint_2",
                          output_dir = file.path(cfg$output$base_dir))

  # Check that Stage 3 outputs exist
  stage3_dir <- cfg$output$stage3_dir
  estimates_file <- file.path(stage3_dir, "stage3_estimates.csv")
  if (!file.exists(estimates_file)) {
    stop("Stage 3 estimates not found at ", estimates_file,
         ". Run stage3_estimation first.", call. = FALSE)
  }

  # --------------------------------------------------------------------------
  # Compile report artifacts
  # --------------------------------------------------------------------------

  # 1. Copy key files to report directory
  files_to_copy <- c(
    file.path(cfg$output$stage1_dir, "stage1_report.json"),
    file.path(cfg$output$stage2_dir, "stage2_diagnostics.json"),
    file.path(cfg$output$stage2_dir, "ess_table.csv"),
    file.path(stage3_dir, "stage3_estimates.csv"),
    file.path(stage3_dir, "cox_hr_secondary.csv"),
    file.path(cfg$output$base_dir, "checkpoint_1.json"),
    file.path(cfg$output$base_dir, "checkpoint_2.json")
  )

  for (f in files_to_copy) {
    if (file.exists(f)) {
      file.copy(f, file.path(output_dir, basename(f)), overwrite = TRUE)
    }
  }

  # Copy diagnostic plots
  plot_files <- list.files(cfg$output$stage2_dir, pattern = "\\.png$",
                           full.names = TRUE)
  for (f in plot_files) {
    file.copy(f, file.path(output_dir, basename(f)), overwrite = TRUE)
  }

  # 2. Copy decision log
  log_file <- file.path(cfg$output$base_dir, "decision_log.csv")
  if (file.exists(log_file)) {
    file.copy(log_file, file.path(output_dir, "decision_log.csv"),
              overwrite = TRUE)
  }

  # 3. Write config used
  yaml::write_yaml(cfg, file.path(output_dir, "config_used.yml"))

  # 4. Capture session info
  capture_session_info(file.path(output_dir, "session_info.txt"))

  # 5. Write report metadata
  metadata <- list(
    report_generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    r_version        = paste(R.version$major, R.version$minor, sep = "."),
    stages_completed = c("stage1", "stage2", "stage3", "stage4"),
    locked           = TRUE
  )
  jsonlite::write_json(metadata, file.path(output_dir, "report_metadata.json"),
                       pretty = TRUE, auto_unbox = TRUE)

  # Log decisions
  mtg <- start_meeting("stage4", protocol_version = 1L)
  mtg <- log_decision(
    mtg, "Report",
    "Final report compiled and locked",
    "All stage checkpoints passed; artifacts compiled",
    triggered_by = "stage completion"
  )
  close_meeting(mtg)

  message("Stage 4 report generated at: ", output_dir)
  invisible(output_dir)
}
