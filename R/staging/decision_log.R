#' @title Decision Log for Clean-Room Governance
#' @description Functions for logging analysis decisions with timestamps,
#'   justifications, and outcome-blindness confirmation.
#' @name decision_log
NULL

#' Initialize Decision Log
#'
#' Create or read the decision log CSV file.
#'
#' @param log_path Character path to decision log CSV.
#' @return Data frame of existing log entries (empty if new).
#' @export
init_decision_log <- function(log_path = "outputs/decision_log.csv") {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  cols <- c("meeting_id", "date_time", "stage", "model_component",
            "decision", "justification", "triggered_by",
            "outcome_blind_confirm", "approver", "protocol_version")
  if (file.exists(log_path)) {
    log <- utils::read.csv(log_path, stringsAsFactors = FALSE)
    # Ensure all columns present
    for (col in cols) {
      if (!col %in% names(log)) log[[col]] <- NA
    }
    return(log)
  }
  log <- data.frame(matrix(ncol = length(cols), nrow = 0),
                    stringsAsFactors = FALSE)
  names(log) <- cols
  utils::write.csv(log, log_path, row.names = FALSE)
  log
}

#' Start a Decision Meeting
#'
#' Creates a meeting context with a unique ID and timestamp.
#'
#' @param stage Character stage identifier (e.g., "stage1", "stage2").
#' @param approver Character name of the approver.
#' @param protocol_version Integer protocol version.
#' @return A list representing the meeting context.
#' @export
start_meeting <- function(stage, approver = "analyst",
                          protocol_version = 1L) {
  list(
    meeting_id       = paste0("MTG-", format(Sys.time(), "%Y%m%d%H%M%S"),
                              "-", sample(1000:9999, 1)),
    date_time        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stage            = stage,
    approver         = approver,
    protocol_version = protocol_version,
    decisions        = list()
  )
}

#' Log a Decision
#'
#' Record an analysis decision within a meeting context.
#'
#' @param meeting Meeting context from \code{start_meeting}.
#' @param model_component Character; e.g., "PS", "Q", "Censoring", "Cohort".
#' @param decision Character description of the decision.
#' @param justification Character justification for the decision.
#' @param triggered_by Character description of what diagnostic triggered
#'   the decision.
#' @param outcome_blind_confirm Logical; TRUE if this decision was made
#'   without seeing outcome-by-treatment comparisons.
#' @return Updated meeting context with new decision appended.
#' @export
log_decision <- function(meeting, model_component, decision, justification,
                         triggered_by = "pre-specified",
                         outcome_blind_confirm = TRUE) {
  entry <- data.frame(
    meeting_id           = meeting$meeting_id,
    date_time            = meeting$date_time,
    stage                = meeting$stage,
    model_component      = model_component,
    decision             = decision,
    justification        = justification,
    triggered_by         = triggered_by,
    outcome_blind_confirm = outcome_blind_confirm,
    approver             = meeting$approver,
    protocol_version     = meeting$protocol_version,
    stringsAsFactors     = FALSE
  )

  meeting$decisions <- c(meeting$decisions, list(entry))
  meeting
}

#' Close a Meeting and Write Decisions
#'
#' Appends all logged decisions from a meeting to the decision log file.
#'
#' @param meeting Meeting context with decisions.
#' @param log_path Character path to decision log CSV.
#' @return Invisibly returns the full updated log.
#' @export
close_meeting <- function(meeting, log_path = "outputs/decision_log.csv") {
  if (length(meeting$decisions) == 0) {
    message("No decisions logged in meeting ", meeting$meeting_id)
    return(invisible(init_decision_log(log_path)))
  }
  new_rows <- do.call(rbind, meeting$decisions)
  existing <- init_decision_log(log_path)
  updated  <- rbind(existing, new_rows)
  utils::write.csv(updated, log_path, row.names = FALSE)
  message("Logged ", nrow(new_rows), " decision(s) for meeting ",
          meeting$meeting_id)
  invisible(updated)
}

#' Guard: Enforce Outcome-Blind Stage 2 Decisions
#'
#' Checks that Stage 3 results do not exist when regenerating Stage 2 outputs.
#' If Stage 3 results are found, requires incrementing the protocol version.
#'
#' @param stage3_dir Character path to Stage 3 outputs directory.
#' @param current_protocol_version Integer current protocol version.
#' @param log_path Character path to decision log.
#' @return Invisibly returns TRUE if check passes.
#' @export
guard_outcome_blind <- function(stage3_dir = "outputs/stage3",
                                current_protocol_version = 1L,
                                log_path = "outputs/decision_log.csv") {
  stage3_files <- list.files(stage3_dir, full.names = TRUE)
  if (length(stage3_files) > 0) {
    # Check if protocol version has been incremented
    log <- init_decision_log(log_path)
    if (nrow(log) > 0) {
      max_s3_version <- max(
        log$protocol_version[log$stage == "stage3"], 0,
        na.rm = TRUE
      )
      if (current_protocol_version <= max_s3_version) {
        stop(
          "OUTCOME BLINDNESS VIOLATION: Stage 3 results exist in '",
          stage3_dir, "'. ",
          "To regenerate Stage 2 outputs, you must increment the protocol ",
          "version above ", max_s3_version,
          " and document the justification in the decision log.",
          call. = FALSE
        )
      }
    }
    warning("Stage 3 results exist. Proceeding with incremented protocol ",
            "version ", current_protocol_version, ".")
  }
  invisible(TRUE)
}
