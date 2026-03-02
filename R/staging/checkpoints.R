#' @title Checkpoint Management
#' @description Read and write checkpoint files for stage gating.
#' @name checkpoints
NULL

#' Write Checkpoint
#'
#' Writes a checkpoint JSON file with PASS/FAIL status and criteria details.
#'
#' @param stage Character stage name (e.g., "checkpoint_1").
#' @param status Character "PASS" or "FAIL".
#' @param criteria Named list of criteria and their values.
#' @param output_dir Character path to output directory.
#' @return The file path to the written checkpoint, invisibly.
#' @export
write_checkpoint <- function(stage, status, criteria,
                             output_dir = "outputs") {
  stopifnot(status %in% c("PASS", "FAIL"))
  checkpoint <- list(
    stage      = stage,
    status     = status,
    timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    criteria   = criteria
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fpath <- file.path(output_dir, paste0(stage, ".json"))
  jsonlite::write_json(checkpoint, fpath, pretty = TRUE, auto_unbox = TRUE)
  message(stage, ": ", status)
  invisible(fpath)
}

#' Read Checkpoint
#'
#' Reads a checkpoint JSON file and returns its contents.
#'
#' @param stage Character stage name (e.g., "checkpoint_1").
#' @param output_dir Character path to output directory.
#' @return Named list with checkpoint data, or NULL if not found.
#' @export
read_checkpoint <- function(stage, output_dir = "outputs") {
  fpath <- file.path(output_dir, paste0(stage, ".json"))
  if (!file.exists(fpath)) return(NULL)
  jsonlite::read_json(fpath, simplifyVector = TRUE)
}

#' Require Checkpoint Pass
#'
#' Stops execution if the required checkpoint has not passed.
#'
#' @param stage Character stage name to check.
#' @param output_dir Character path to output directory.
#' @return Invisibly returns TRUE if checkpoint passed.
#' @export
require_checkpoint_pass <- function(stage, output_dir = "outputs") {
  cp <- read_checkpoint(stage, output_dir)
  if (is.null(cp)) {
    stop("Checkpoint '", stage, "' not found. Run the preceding stage first.",
         call. = FALSE)
  }
  if (cp$status != "PASS") {
    stop("Checkpoint '", stage, "' status is FAIL. Cannot proceed. ",
         "Review diagnostics and address issues before continuing.",
         call. = FALSE)
  }
  invisible(TRUE)
}
