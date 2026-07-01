# Estimand and sensitivity-plan declarations attached to a cleanroom analysis
# lock. These only stash structured metadata on the lock (a list); they do not
# estimate anything and are not read by the cleanTMLE estimation core or the
# pre-outcome gate. Moved out of cleanTMLE in phase 2 of the estimation-vs-
# governance split.

#' Attach an Estimand to an Analysis Lock
#'
#' Enriches a \code{cleanroom_lock} with a structured description of the
#' causal and statistical estimand. The lock hash is unchanged (estimand
#' metadata does not invalidate the analytic specification).
#'
#' @param lock A \code{cleanroom_lock} from
#'   \code{cleanTMLE::create_analysis_lock()}.
#' @param description Character; free-text description of the causal question.
#' @param population Character; description of the target population.
#' @param treatment_strategies Character vector of length 2 describing the
#'   treatment contrast (e.g. \code{c("Drug A", "Placebo")}).
#' @param outcome_label Character; description of the primary outcome.
#' @param followup Character; follow-up horizon description.
#' @param contrast Character; estimand type, e.g. \code{"risk_difference"}.
#' @param statistical_estimand Character; the statistical (observed-data)
#'   estimand under identification assumptions.
#'
#' @return A modified \code{cleanroom_lock} with element \code{estimand}.
#'
#' @export
attach_estimand <- function(lock,
                            description          = NULL,
                            population           = NULL,
                            treatment_strategies = NULL,
                            outcome_label        = NULL,
                            followup             = NULL,
                            contrast             = "risk_difference",
                            statistical_estimand = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  lock$estimand <- list(
    description          = description,
    population           = population,
    treatment_strategies = treatment_strategies,
    outcome_label        = outcome_label,
    followup             = followup,
    contrast             = contrast,
    statistical_estimand = statistical_estimand
  )
  lock
}


#' Declare a Sensitivity Analysis Plan
#'
#' Attaches a named sensitivity analysis plan to the analysis lock. The plan is
#' stored as metadata and does not alter the analytic specification hash.
#' Multiple calls append additional plans.
#'
#' @param lock A \code{cleanroom_lock}.
#' @param label Character; short name for this sensitivity analysis.
#' @param description Character; what the sensitivity analysis evaluates.
#' @param settings A named list of parameters (e.g.
#'   \code{list(truncation = c(0.01, 0.05, 0.10))}).
#'
#' @return Modified \code{cleanroom_lock} with appended sensitivity plan.
#'
#' @export
declare_sensitivity_plan <- function(lock, label, description = NULL,
                                     settings = list()) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!is.character(label) || length(label) != 1L)
    stop("`label` must be a single character string.", call. = FALSE)

  plan <- list(
    label       = label,
    description = description,
    settings    = settings,
    declared_at = Sys.time()
  )

  if (is.null(lock$sensitivity_plans)) lock$sensitivity_plans <- list()
  lock$sensitivity_plans[[label]] <- plan
  lock
}
