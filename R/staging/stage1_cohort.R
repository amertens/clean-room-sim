#' @title Stage 1: Cohort Feasibility and Build
#' @description Constructs the analysis cohort and produces feasibility
#'   summaries. Outputs are restricted to cohort-level summaries; NO
#'   outcome-by-treatment comparisons are permitted.
#' @name stage1
NULL

#' Build Cohort (Stage 1)
#'
#' Generates the analysis cohort from the DGP and produces a stage 1
#' report with cohort size, attrition, variable missingness, and overall
#' (non-treatment-stratified) event rates.
#'
#' @param data Data frame from \code{generate_hcv_data}, or NULL to generate.
#' @param spec Named list of DGP parameters (used if data is NULL).
#' @param cfg Config list from \code{load_config}.
#' @param output_dir Character path for Stage 1 outputs.
#' @return A list with components:
#'   \describe{
#'     \item{cohort}{The analysis data frame.}
#'     \item{report}{A list of feasibility summaries.}
#'     \item{checkpoint}{PASS/FAIL status.}
#'   }
#' @export
stage1_build_cohort <- function(data = NULL, spec = NULL, cfg = NULL,
                                output_dir = "outputs/stage1") {
  if (is.null(cfg)) cfg <- load_config()
  ensure_dir(output_dir)

  # Generate data if not provided
  if (is.null(data)) {
    if (is.null(spec)) spec <- cfg$dgp
    data <- do.call(generate_hcv_data, spec)
  }

  # --- Feasibility summaries (outcome-blind) ---
  N_total    <- nrow(data)
  N_treated  <- sum(data$treatment == 1)
  N_control  <- sum(data$treatment == 0)
  frac_treated <- N_treated / N_total
  frac_control <- N_control / N_total

  # Overall event rate (NOT stratified by treatment)
  overall_event_rate <- mean(data$event)
  overall_N_events   <- sum(data$event)

  # Follow-up summary (NOT stratified by treatment)
  fu_summary <- summary(data$follow_time)

  # Variable missingness
  miss_summary <- vapply(data, function(x) mean(is.na(x)), numeric(1))
  miss_summary <- miss_summary[miss_summary > 0]
  if (length(miss_summary) == 0) {
    miss_summary <- c(none = 0)
  }

  # Variable type summary
  numeric_vars <- names(data)[vapply(data, is.numeric, logical(1))]
  numeric_vars <- setdiff(numeric_vars, c("id", "treatment", "event",
                                           "follow_time", "switch"))

  # Build report
  report <- list(
    N_total          = N_total,
    N_treated        = N_treated,
    N_control        = N_control,
    frac_treated     = frac_treated,
    frac_control     = frac_control,
    overall_event_rate = overall_event_rate,
    overall_N_events   = overall_N_events,
    followup_summary   = as.list(fu_summary),
    missingness        = as.list(miss_summary),
    n_covariates       = length(numeric_vars),
    covariate_names    = numeric_vars
  )

  # Save report
  jsonlite::write_json(report, file.path(output_dir, "stage1_report.json"),
                       pretty = TRUE, auto_unbox = TRUE)

  # --- Checkpoint evaluation ---
  thresholds <- cfg$stage1
  criteria <- list(
    N_total       = list(value = N_total,
                         threshold = thresholds$min_N,
                         pass = N_total >= thresholds$min_N),
    frac_treated  = list(value = round(frac_treated, 4),
                         threshold = thresholds$min_treated_frac,
                         pass = frac_treated >= thresholds$min_treated_frac),
    frac_control  = list(value = round(frac_control, 4),
                         threshold = thresholds$min_control_frac,
                         pass = frac_control >= thresholds$min_control_frac),
    max_missingness = list(
      value = if (length(miss_summary) > 0)
        round(max(miss_summary), 4) else 0,
      threshold = thresholds$max_missingness,
      pass = all(miss_summary <= thresholds$max_missingness)
    )
  )

  all_pass <- all(vapply(criteria, function(x) x$pass, logical(1)))
  status <- if (all_pass) "PASS" else "FAIL"
  write_checkpoint("checkpoint_1", status, criteria,
                   output_dir = dirname(output_dir))

  # Log initial decisions
  mtg <- start_meeting("stage1", protocol_version = 1L)
  mtg <- log_decision(mtg, "Cohort", "Cohort built and feasibility assessed",
                      paste("N =", N_total, "; events =", overall_N_events),
                      triggered_by = "pre-specified protocol")
  close_meeting(mtg)

  list(cohort = data, report = report, checkpoint = status)
}
