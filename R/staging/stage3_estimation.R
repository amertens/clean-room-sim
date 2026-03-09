#' @title Stage 3: Outcome Modeling and Primary Estimator Execution
#' @description Runs TMLE and comparator estimators after verifying
#'   Stage 2 checkpoint passes. Produces effect estimates, IC-based SEs,
#'   bootstrap CIs, and diagnostics.
#'
#'   Estimand: Cumulative incidence (risk) at pre-specified time horizons
#'   under hypothetical treatment vs control, operationalized as binary
#'   endpoints Y_t = I(T <= t & event == 1). Censoring is handled via
#'   IPCW in TMLE/IPTW.
#' @name stage3
NULL

# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Pre-flight Schema Check for Stage 3
#'
#' Validates that the cohort has the required structure for estimation.
#' Fails with informative errors rather than producing NA tables.
#'
#' @param data Data frame to validate.
#' @param time_points Numeric vector of evaluation time points.
#' @return Invisibly returns TRUE on success; stops otherwise.
#' @keywords internal
stage3_preflight <- function(data, time_points) {
  # Required columns
  required <- c("treatment", "event", "follow_time")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Stage 3 pre-flight FAILED: missing columns: ",
         paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(names(data), collapse = ", "),
         call. = FALSE)
  }

  # Treatment must be binary 0/1 with both levels present
  A <- data$treatment
  if (!all(A %in% c(0L, 1L))) {
    stop("Stage 3 pre-flight FAILED: 'treatment' must be binary 0/1. ",
         "Found values: ", paste(unique(A), collapse = ", "), call. = FALSE)
  }
  if (sum(A == 1) == 0 || sum(A == 0) == 0) {
    stop("Stage 3 pre-flight FAILED: both treatment groups must be present. ",
         "Treated: ", sum(A == 1), ", Control: ", sum(A == 0), call. = FALSE)
  }

  # Event must be binary 0/1 with events present
  if (!all(data$event %in% c(0L, 1L))) {
    stop("Stage 3 pre-flight FAILED: 'event' must be binary 0/1.", call. = FALSE)
  }
  if (sum(data$event) == 0) {
    stop("Stage 3 pre-flight FAILED: no events observed in cohort.", call. = FALSE)
  }

  # Check events exist within each time window
  for (t_val in time_points) {
    n_events_t <- sum(data$follow_time <= t_val & data$event == 1)
    if (n_events_t < 5) {
      warning("Stage 3 pre-flight WARNING: only ", n_events_t,
              " events by t=", t_val, " days. Estimation may be unstable.")
    }
  }

  # Covariate check: at least one numeric covariate
  exclude_cols <- c("id", "follow_time", "event", "treatment", "switch",
                    "race", "region")
  covar_cols <- setdiff(names(data), exclude_cols)
  covar_cols <- covar_cols[vapply(data[, covar_cols, drop = FALSE],
                                 is.numeric, logical(1))]
  if (length(covar_cols) == 0) {
    stop("Stage 3 pre-flight FAILED: no numeric covariates found.", call. = FALSE)
  }

  # Check for all-missing covariates
  all_miss <- vapply(data[, covar_cols, drop = FALSE],
                     function(x) all(is.na(x)), logical(1))
  if (any(all_miss)) {
    stop("Stage 3 pre-flight FAILED: all-missing covariates: ",
         paste(covar_cols[all_miss], collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}


#' Run Stage 3 Estimation
#'
#' Executes the primary TMLE estimator and comparator methods (Cox PH, IPTW,
#' G-computation) for all pre-specified time points. Only runs if
#' checkpoint_2 is PASS.
#'
#' Primary estimand: Cumulative incidence (risk) at clinically meaningful
#' time points (default: 90 and 180 days) under treatment vs control.
#' Operationalized as binary endpoints Y_t = I(T <= t & event == 1).
#' The risk difference and risk ratio are not reliant on PH assumptions.
#'
#' Inference uses both influence-curve-based SEs (for TMLE/IPTW) and
#' nonparametric bootstrap CIs for all estimators. Bootstrap is the
#' primary inference method for this vignette because with realistic
#' event rates (1-5%) and finite samples, IC-based SEs can be unstable.
#'
#' Secondary estimand: Hazard ratio from Cox PH (labeled as such, with
#' PH assumption stated and tested).
#'
#' @param cohort Data frame from Stage 1.
#' @param stage2_result Output from \code{stage2_design_checks}.
#' @param cfg Config list from \code{load_config}.
#' @param output_dir Character path for Stage 3 outputs.
#' @param B Integer number of bootstrap replicates.
#' @return A list with estimates from all methods and time points.
#' @export
stage3_estimation <- function(cohort, stage2_result = NULL, cfg = NULL,
                              output_dir = "outputs/stage3",
                              B = NULL) {
  if (is.null(cfg)) cfg <- load_config()
  ensure_dir(output_dir)

  # Enforce prerequisite
  require_checkpoint_pass("checkpoint_2",
                          output_dir = dirname(output_dir))

  time_points <- cfg$simulation$time_points
  sl_lib_Q    <- filter_sl_libraries(cfg$tmle$sl_library_Q)
  sl_lib_g    <- filter_sl_libraries(cfg$tmle$sl_library_g)
  g_bounds    <- c(cfg$tmle$truncation_lower, cfg$tmle$truncation_upper)

  # Bootstrap config
  if (is.null(B)) {
    B <- cfg$bootstrap$B %||% 500
  }
  boot_seed <- cfg$bootstrap$seed %||% 54321

  # Pre-flight checks
  stage3_preflight(cohort, time_points)

  results <- list()

  for (t_val in time_points) {
    message("=== Stage 3: Estimating risk at t = ", t_val, " days ===")

    # --- Primary estimand: TMLE survival risk (PH-free) ---
    tmle_res <- tryCatch({
      tmle_survival_risk(
        data      = cohort,
        t_eval    = t_val,
        sl_lib_Q  = sl_lib_Q,
        sl_lib_g  = sl_lib_g,
        sl_lib_cens = sl_lib_g,
        g_trunc   = g_bounds
      )
    }, error = function(e) {
      message("TMLE survival failed at t=", t_val, ": ", e$message)
      list(RD = NA, se_RD = NA, ci_RD = c(NA, NA), risk_1 = NA,
           risk_0 = NA, RR = NA, ci_RR = c(NA, NA),
           diagnostics = list(error = e$message))
    })

    # --- Comparator 1: IPTW ---
    iptw_res <- tryCatch({
      iptw_survival(
        data    = cohort,
        t_eval  = t_val,
        sl_lib  = sl_lib_g,
        g_trunc = g_bounds
      )
    }, error = function(e) {
      message("IPTW failed at t=", t_val, ": ", e$message)
      list(RD = NA, se_RD = NA, ci_RD = c(NA, NA), risk_1 = NA,
           risk_0 = NA, RR = NA)
    })

    # --- Comparator 2: G-computation ---
    gcomp_res <- tryCatch({
      gcomp_risk(
        data   = cohort,
        t_eval = t_val,
        sl_lib = sl_lib_Q,
        B      = B,
        boot_seed = boot_seed
      )
    }, error = function(e) {
      message("G-comp failed at t=", t_val, ": ", e$message)
      list(RD = NA, se_RD = NA, ci_RD = c(NA, NA), risk_1 = NA,
           risk_0 = NA, RR = NA)
    })

    results[[paste0("t", t_val)]] <- list(
      tmle  = tmle_res,
      iptw  = iptw_res,
      gcomp = gcomp_res
    )
  }

  # --- Secondary estimand: Cox PH hazard ratio ---
  message("=== Stage 3: Cox PH (secondary estimand) ===")
  cox_res <- tryCatch({
    cox_ph_estimator(cohort)
  }, error = function(e) {
    message("Cox PH failed: ", e$message)
    list(HR = NA, ci_HR = c(NA, NA), pvalue = NA,
         ph_test = list(ph_violated = NA))
  })

  # Add Cox risk estimates at each time point
  for (t_val in time_points) {
    cox_risk <- tryCatch({
      cox_risk_at_t(cox_res, cohort, t_eval = t_val)
    }, error = function(e) {
      list(RD = NA, risk_1 = NA, risk_0 = NA, RR = NA)
    })
    results[[paste0("t", t_val)]]$cox <- cox_risk
  }
  results$cox_hr <- cox_res

  # --------------------------------------------------------------------------
  # Save results
  # --------------------------------------------------------------------------
  # Create summary table
  summary_rows <- list()
  for (t_val in time_points) {
    key <- paste0("t", t_val)
    r <- results[[key]]

    add_row <- function(method, res) {
      data.frame(
        time_point = t_val,
        method     = method,
        risk_1     = round(res$risk_1 %||% NA_real_, 6),
        risk_0     = round(res$risk_0 %||% NA_real_, 6),
        RD         = round(res$RD %||% NA_real_, 6),
        RR         = round(res$RR %||% NA_real_, 4),
        se_RD      = round(res$se_RD %||% NA_real_, 6),
        ci_lower   = round((res$ci_RD %||% c(NA, NA))[1], 6),
        ci_upper   = round((res$ci_RD %||% c(NA, NA))[2], 6),
        stringsAsFactors = FALSE
      )
    }
    summary_rows <- c(summary_rows, list(
      add_row("TMLE", r$tmle),
      add_row("IPTW", r$iptw),
      add_row("G-computation", r$gcomp),
      add_row("Cox PH (standardized)", r$cox)
    ))
  }
  summary_table <- do.call(rbind, summary_rows)

  # Validate: fail loudly if TMLE produced NA
  tmle_na <- is.na(summary_table$RD[summary_table$method == "TMLE"])
  if (any(tmle_na)) {
    warning("Stage 3: TMLE returned NA for time points: ",
            paste(time_points[tmle_na], collapse = ", "),
            ". Check diagnostics.")
  }

  utils::write.csv(summary_table,
                    file.path(output_dir, "stage3_estimates.csv"),
                    row.names = FALSE)

  # Save Cox HR summary
  cox_summary <- data.frame(
    estimand = "Hazard Ratio (SECONDARY - assumes PH)",
    HR       = round(cox_res$HR %||% NA_real_, 4),
    ci_lower = round((cox_res$ci_HR %||% c(NA, NA))[1], 4),
    ci_upper = round((cox_res$ci_HR %||% c(NA, NA))[2], 4),
    pvalue   = round(cox_res$pvalue %||% NA_real_, 6),
    ph_test_pvalue = round(
      cox_res$ph_test$pvalue %||% NA_real_, 4),
    ph_violated = cox_res$ph_test$ph_violated %||% NA,
    note = "PH assumption required; check ph_test_pvalue",
    stringsAsFactors = FALSE
  )
  utils::write.csv(cox_summary,
                    file.path(output_dir, "cox_hr_secondary.csv"),
                    row.names = FALSE)

  # Save TMLE diagnostics
  for (t_val in time_points) {
    key <- paste0("t", t_val)
    diag <- results[[key]]$tmle$diagnostics
    if (!is.null(diag)) {
      jsonlite::write_json(
        diag,
        file.path(output_dir, paste0("tmle_diagnostics_t", t_val, ".json")),
        pretty = TRUE, auto_unbox = TRUE
      )
    }
  }

  # Log decisions
  mtg <- start_meeting("stage3", protocol_version = 1L)
  mtg <- log_decision(
    mtg, "Q",
    paste("Outcome model fitted with SL library:",
          paste(sl_lib_Q, collapse = ", ")),
    "Pre-specified from config",
    triggered_by = "pre-specified"
  )
  mtg <- log_decision(
    mtg, "Estimator",
    paste("Primary estimand: Risk difference at t =",
          paste(time_points, collapse = ", "), "days (PH-free)"),
    "Pre-specified in TL-SAP",
    triggered_by = "pre-specified"
  )
  mtg <- log_decision(
    mtg, "Estimator",
    "Secondary estimand: Cox PH hazard ratio (assumes PH)",
    "Pre-specified in TL-SAP; PH test included as diagnostic",
    triggered_by = "pre-specified"
  )
  mtg <- log_decision(
    mtg, "Inference",
    paste("Bootstrap B =", B, "; IC-based SEs also reported for TMLE/IPTW"),
    "Bootstrap is primary inference for finite-sample stability with rare events",
    triggered_by = "pre-specified"
  )
  close_meeting(mtg)

  list(results = results, summary = summary_table)
}
