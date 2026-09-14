# Human-readable summaries of a cleanTMLE audit log. These read the audit's
# recorded stage entries and render a manifest / narrative; they compute
# nothing and are read by neither the estimation core nor the gate. Moved out
# of cleanTMLE in phase 2 of the estimation-vs-governance split.

#' Build a Stage Manifest
#'
#' Produces a compact textual summary of the staged analysis path taken, based
#' on the audit log.
#'
#' @param audit A \code{cleantmle_audit} (from
#'   \code{cleanTMLE::create_audit_log()}).
#'
#' @return A character string (printed, and returned invisibly) showing the
#'   stage path.
#'
#' @export
build_stage_manifest <- function(audit) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  if (length(audit$entries) == 0L) {
    msg <- "No stage entries recorded."
    cat(msg, "\n")
    return(invisible(msg))
  }

  lines <- vapply(seq_along(audit$entries), function(i) {
    e <- audit$entries[[i]]
    dec <- if (!is.na(e$decision)) sprintf(" [%s]", e$decision) else ""
    sprintf("  %d. %s: %s%s", i, e$stage, e$action, dec)
  }, character(1))

  header <- sprintf("Clean-Room Stage Path (lock: %s)", audit$lock_hash)
  sep    <- paste(rep("-", nchar(header)), collapse = "")
  msg    <- paste(c(header, sep, lines), collapse = "\n")
  cat(msg, "\n")
  invisible(msg)
}


#' Summarise the Stage Path as a Compact Narrative
#'
#' Generates a compact, human-readable narrative string showing the sequence of
#' stages and checkpoint decisions recorded in the audit log.
#'
#' @param audit A \code{cleantmle_audit}.
#'
#' @return The narrative string, returned invisibly after printing.
#'
#' @export
summarize_stage_path <- function(audit) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  if (length(audit$entries) == 0L) {
    msg <- "(No stages recorded)"
    cat(msg, "\n")
    return(invisible(msg))
  }

  parts <- vapply(audit$entries, function(e) {
    if (!is.na(e$decision) && nzchar(e$decision)) {
      sprintf("%s: %s", e$stage, e$decision)
    } else {
      e$stage
    }
  }, character(1L))

  narrative <- paste(parts, collapse = " -> ")
  cat(narrative, "\n")
  invisible(narrative)
}
