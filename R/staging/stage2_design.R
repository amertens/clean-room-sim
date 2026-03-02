#' @title Stage 2: Design Adequacy Checks
#' @description Evaluates propensity score overlap, covariate balance, weight
#'   distributions, and effective sample size. Explicitly does NOT compute
#'   treatment effect estimates or fit outcome models.
#' @name stage2
NULL

#' Run Design Checks (Stage 2)
#'
#' Fits propensity score models, evaluates overlap and balance, and produces
#' diagnostic outputs. This stage is restricted to treatment/covariate
#' information only -- no outcome models or effect estimates.
#'
#' @param cohort Data frame from Stage 1.
#' @param spec Named list of analysis parameters.
#' @param cfg Config list from \code{load_config}.
#' @param output_dir Character path for Stage 2 outputs.
#' @param protocol_version Integer protocol version for outcome-blind guard.
#' @return A list with components:
#'   \describe{
#'     \item{ps}{Fitted propensity scores.}
#'     \item{weights}{IPW weights.}
#'     \item{diagnostics}{SMDs, ESS, truncation info.}
#'     \item{checkpoint}{PASS/FAIL status.}
#'   }
#' @export
stage2_design_checks <- function(cohort, spec = NULL, cfg = NULL,
                                 output_dir = "outputs/stage2",
                                 protocol_version = 1L) {
  if (is.null(cfg)) cfg <- load_config()
  ensure_dir(output_dir)

  # Enforce prerequisites
  require_checkpoint_pass("checkpoint_1",
                          output_dir = dirname(output_dir))
  guard_outcome_blind(file.path(dirname(output_dir), "stage3"),
                      current_protocol_version = protocol_version)

  # --------------------------------------------------------------------------
  # Identify covariates for PS model
  # --------------------------------------------------------------------------
  exclude_cols <- c("id", "treatment", "event", "follow_time", "switch",
                    "race", "region")
  covar_cols <- setdiff(names(cohort), exclude_cols)
  covar_cols <- covar_cols[vapply(cohort[, covar_cols, drop = FALSE],
                                 is.numeric, logical(1))]

  W <- as.data.frame(cohort[, covar_cols, drop = FALSE])
  A <- cohort$treatment

  # --------------------------------------------------------------------------
  # Fit propensity score model using SuperLearner if available
  # --------------------------------------------------------------------------
  sl_libs <- cfg$tmle$sl_library_g
  # Filter to available learners
  available_libs <- Filter(function(lib) {
    tryCatch({
      pkg <- sub("SL\\.", "", lib)
      if (pkg %in% c("glm", "mean", "step")) return(TRUE)
      requireNamespace(pkg, quietly = TRUE)
    }, error = function(e) FALSE)
  }, sl_libs)
  if (length(available_libs) == 0) available_libs <- c("SL.glm", "SL.mean")

  ps_fit <- tryCatch({
    SuperLearner::SuperLearner(
      Y = A, X = W, family = stats::binomial(),
      SL.library = available_libs,
      cvControl = list(V = 5)
    )
  }, error = function(e) {
    message("SuperLearner failed for PS, falling back to GLM: ", e$message)
    fit <- stats::glm(A ~ ., data = cbind(A = A, W), family = "binomial")
    list(SL.predict = stats::predict(fit, type = "response"),
         library.predict = matrix(stats::predict(fit, type = "response"),
                                  ncol = 1),
         libraryNames = "SL.glm",
         fallback = TRUE)
  })

  ps_raw <- as.numeric(ps_fit$SL.predict)

  # --------------------------------------------------------------------------
  # Truncation
  # --------------------------------------------------------------------------
  trunc <- truncate_ps(ps_raw,
                       lower = cfg$tmle$truncation_lower,
                       upper = cfg$tmle$truncation_upper)
  ps <- trunc$p

  # --------------------------------------------------------------------------
  # IPW weights
  # --------------------------------------------------------------------------
  weights <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
  weights_trimmed <- pmin(weights, cfg$stage2$max_weight)

  # --------------------------------------------------------------------------
  # SMDs (unweighted and weighted)
  # --------------------------------------------------------------------------
  smd_unweighted <- vapply(covar_cols, function(v) {
    compute_smd(cohort[[v]], A)
  }, numeric(1))

  smd_weighted <- vapply(covar_cols, function(v) {
    compute_smd(cohort[[v]], A, weights = weights_trimmed)
  }, numeric(1))

  # --------------------------------------------------------------------------
  # Effective sample size
  # --------------------------------------------------------------------------
  ess_treated <- effective_ss(weights_trimmed[A == 1])
  ess_control <- effective_ss(weights_trimmed[A == 0])
  ess_total   <- ess_treated + ess_control
  ess_frac    <- ess_total / length(A)

  # --------------------------------------------------------------------------
  # Overlap diagnostics
  # --------------------------------------------------------------------------
  ps_overlap_region <- mean(ps > cfg$stage2$min_ps_overlap &
                              ps < (1 - cfg$stage2$min_ps_overlap))

  # --------------------------------------------------------------------------
  # Diagnostic plots (saved to files)
  # --------------------------------------------------------------------------
  # PS overlap plot
  grDevices::png(file.path(output_dir, "ps_overlap.png"),
                 width = 800, height = 500)
  plot_data <- data.frame(ps = ps, treatment = factor(A, labels = c("Control",
                                                                     "SOF")))
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ps, fill = treatment)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(title = "Propensity Score Overlap",
                  x = "P(SOF | W)", y = "Density", fill = "Treatment") +
    ggplot2::theme_minimal() +
    ggplot2::geom_vline(xintercept = c(cfg$tmle$truncation_lower,
                                        cfg$tmle$truncation_upper),
                         linetype = "dashed", color = "red")
  print(p1)
  grDevices::dev.off()

  # SMD Love plot
  smd_df <- data.frame(
    variable   = rep(covar_cols, 2),
    smd        = c(smd_unweighted, smd_weighted),
    type       = rep(c("Unweighted", "Weighted"),
                     each = length(covar_cols)),
    stringsAsFactors = FALSE
  )

  grDevices::png(file.path(output_dir, "smd_love_plot.png"),
                 width = 800, height = 600)
  p2 <- ggplot2::ggplot(smd_df, ggplot2::aes(x = abs(smd),
                                               y = stats::reorder(variable,
                                                                   abs(smd)),
                                               color = type)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = cfg$stage2$max_smd,
                         linetype = "dotted", color = "red") +
    ggplot2::labs(title = "Covariate Balance (SMD Love Plot)",
                  x = "|SMD|", y = "Covariate", color = "Weighting") +
    ggplot2::theme_minimal()
  print(p2)
  grDevices::dev.off()

  # Weight distribution histogram
  grDevices::png(file.path(output_dir, "weight_distribution.png"),
                 width = 800, height = 500)
  wt_df <- data.frame(weights = weights_trimmed,
                       treatment = factor(A, labels = c("Control", "SOF")))
  p3 <- ggplot2::ggplot(wt_df, ggplot2::aes(x = weights, fill = treatment)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    ggplot2::labs(title = "IPW Weight Distribution",
                  x = "Weight", y = "Count", fill = "Treatment") +
    ggplot2::theme_minimal()
  print(p3)
  grDevices::dev.off()

  # ESS table
  ess_table <- data.frame(
    Group = c("Treated (SOF)", "Control", "Total"),
    N_actual = c(sum(A == 1), sum(A == 0), length(A)),
    ESS = round(c(ess_treated, ess_control, ess_total), 1),
    ESS_fraction = round(c(ess_treated / sum(A == 1),
                            ess_control / sum(A == 0),
                            ess_frac), 3),
    stringsAsFactors = FALSE
  )
  utils::write.csv(ess_table, file.path(output_dir, "ess_table.csv"),
                    row.names = FALSE)

  # --------------------------------------------------------------------------
  # Save diagnostics
  # --------------------------------------------------------------------------
  diagnostics <- list(
    smd_unweighted = as.list(round(smd_unweighted, 4)),
    smd_weighted   = as.list(round(smd_weighted, 4)),
    max_smd_unweighted = round(max(abs(smd_unweighted)), 4),
    max_smd_weighted   = round(max(abs(smd_weighted)), 4),
    ess_treated    = round(ess_treated, 1),
    ess_control    = round(ess_control, 1),
    ess_frac       = round(ess_frac, 4),
    ps_overlap_frac = round(ps_overlap_region, 4),
    n_truncated_lower = trunc$n_lower,
    n_truncated_upper = trunc$n_upper,
    ps_summary     = as.list(round(summary(ps), 4))
  )

  jsonlite::write_json(diagnostics,
                       file.path(output_dir, "stage2_diagnostics.json"),
                       pretty = TRUE, auto_unbox = TRUE)

  # --------------------------------------------------------------------------
  # Checkpoint evaluation
  # --------------------------------------------------------------------------
  thresholds <- cfg$stage2
  criteria <- list(
    max_smd_weighted = list(
      value = round(max(abs(smd_weighted)), 4),
      threshold = thresholds$max_smd,
      pass = max(abs(smd_weighted)) <= thresholds$max_smd
    ),
    ess_frac = list(
      value = round(ess_frac, 4),
      threshold = thresholds$min_ess_frac,
      pass = ess_frac >= thresholds$min_ess_frac
    ),
    ps_overlap = list(
      value = round(ps_overlap_region, 4),
      threshold = thresholds$min_ps_overlap,
      pass = ps_overlap_region >= thresholds$min_ps_overlap
    )
  )

  all_pass <- all(vapply(criteria, function(x) x$pass, logical(1)))
  status <- if (all_pass) "PASS" else "FAIL"
  write_checkpoint("checkpoint_2", status, criteria,
                   output_dir = dirname(output_dir))

  # Log decisions
  mtg <- start_meeting("stage2", protocol_version = protocol_version)
  mtg <- log_decision(
    mtg, "PS",
    paste("PS model fitted with", paste(available_libs, collapse = ", ")),
    "Pre-specified SL library from config",
    triggered_by = "pre-specified"
  )
  mtg <- log_decision(
    mtg, "PS",
    paste("Truncation at [", cfg$tmle$truncation_lower, ",",
          cfg$tmle$truncation_upper, "]",
          "; truncated:", trunc$n_lower, "lower,", trunc$n_upper, "upper"),
    "Pre-specified truncation bounds from config",
    triggered_by = "pre-specified"
  )
  mtg <- log_decision(
    mtg, "Balance",
    paste("Max weighted SMD =", round(max(abs(smd_weighted)), 4)),
    paste("Threshold =", thresholds$max_smd),
    triggered_by = "SMD diagnostic"
  )
  close_meeting(mtg)

  list(
    ps          = ps,
    weights     = weights_trimmed,
    diagnostics = diagnostics,
    checkpoint  = status,
    ps_fit      = ps_fit,
    covar_cols  = covar_cols
  )
}
