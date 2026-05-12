# Workflow definitions: fixed-library TMLE vs plasmode-selected TMLE.
#
# Each workflow takes a synthetic dataset and returns a list with
# the point estimate, the IF-based SE, the 95% CI, and (for the
# plasmode workflow) the selected candidate ID and the inner
# plasmode RMSE table.

library(cleanTMLE)

# Accessor that handles the nested extract_tmle_estimate() output.
.get_est <- function(fit) {
  if (is.null(fit) || is.null(fit$estimates) ||
      is.null(fit$estimates$ATE)) return(list(estimate = NA_real_,
                                              se = NA_real_,
                                              ci_lower = NA_real_,
                                              ci_upper = NA_real_))
  list(
    estimate = unname(fit$estimates$ATE$estimate),
    se       = unname(fit$estimates$ATE$se),
    ci_lower = unname(fit$estimates$ATE$ci_lower),
    ci_upper = unname(fit$estimates$ATE$ci_upper)
  )
}

# ---------------------------------------------------------------
# 1. Fixed-library TMLE
# ---------------------------------------------------------------
#
# Traditional analyst who picks the SuperLearner library and the PS
# truncation level a priori. Default choice is the parametric library
# (Tier 1), which is what a reflex GLM-based analyst would do.

fit_fixed <- function(data, covariates,
                       library = c("SL.glm", "SL.mean"),
                       truncation = 0.01,
                       seed = 1L) {
  lock <- create_simple_lock(
    data       = data,
    treatment  = "treatment",
    outcome    = "event_24",
    covariates = covariates,
    sl_library = library,
    seed       = seed
  )
  spec <- tmle_candidate(
    candidate_id = "fixed_library",
    g_library    = library,
    q_library    = library,
    truncation   = truncation
  )
  lock <- lock_primary_tmle_spec(lock, spec)

  ps   <- fit_ps_superlearner(lock)
  g    <- fit_tmle_treatment_mechanism(lock, ps)
  Qf   <- fit_tmle_outcome_mechanism(lock, g)
  upd  <- run_tmle_targeting_step(g, Qf)
  out  <- .get_est(extract_tmle_estimate(upd))
  out$workflow            <- "fixed"
  out$library_label       <- paste(library, collapse = "+")
  out$truncation          <- truncation
  out$selected_candidate  <- NA_character_
  out
}

# ---------------------------------------------------------------
# 2. Plasmode-selected TMLE
# ---------------------------------------------------------------
#
# cleanTMLE candidate-selection workflow. Runs an outcome-blind
# inner plasmode loop on the outer dataset, picks the candidate
# with the smallest RMSE on that loop, then fits the locked
# candidate on the outer real outcome.

fit_plasmode <- function(data, covariates, candidates,
                          inner_reps    = 30L,
                          inner_effect  = -0.05,
                          rule          = c("min_rmse",
                                            "fiord_two_stage"),
                          seed          = 1L) {
  rule <- match.arg(rule)

  lock <- create_simple_lock(
    data       = data,
    treatment  = "treatment",
    outcome    = "event_24",
    covariates = covariates,
    sl_library = c("SL.glm", "SL.glmnet", "SL.gam",
                   "SL.ranger", "SL.mean"),
    seed       = seed
  )

  # Inner plasmode: candidate selection on outcome-blind synthetic Y
  plas <- tryCatch(
    run_plasmode_feasibility(
      lock,
      tmle_candidates = candidates,
      effect_sizes    = inner_effect,
      reps            = inner_reps
    ),
    error = function(e) NULL
  )

  # Pick a candidate. The inner plasmode can return all-NaN metrics
  # for some seeds (e.g. SL predictions hit the [0,1] boundary and
  # rbinom returns NA, so n_converged = 0 for every candidate). In
  # that case fall back to the default rich candidate so the outer
  # workflow still produces an estimate; flag the fallback via the
  # selected_candidate tag so the manuscript can report frequency.
  default_id <- if ("rich_t01" %in% names(candidates)) "rich_t01"
                else names(candidates)[1L]
  fallback   <- FALSE

  pick_from_agg <- function(m) {
    m <- m[is.finite(m$rmse), , drop = FALSE]
    if (nrow(m) == 0L) return(NA_character_)
    if (rule == "min_rmse") {
      agg <- aggregate(rmse ~ candidate, data = m, FUN = mean)
      as.character(agg$candidate[which.min(agg$rmse)])
    } else {
      agg <- aggregate(cbind(rmse, coverage, mean_se) ~ candidate,
                       data = m, FUN = mean)
      keep <- abs(agg$coverage - 0.95) <= 0.02
      if (any(keep))
        as.character(agg$candidate[keep][which.min(agg$mean_se[keep])])
      else
        as.character(agg$candidate[which.min(abs(agg$coverage - 0.95))])
    }
  }

  selected_id <- tryCatch({
    if (!is.null(plas) && !is.null(plas$metrics) &&
        nrow(plas$metrics) > 0L) pick_from_agg(plas$metrics)
    else NA_character_
  }, error = function(e) NA_character_)

  if (is.na(selected_id) || length(selected_id) == 0L) {
    selected_id <- default_id
    fallback    <- TRUE
  }

  selected <- candidates[[selected_id]]
  lock <- lock_primary_tmle_spec(lock, selected)

  # fit_ps_superlearner uses lock$sl_library directly (it does not
  # honor the locked candidate's g_library). Override the lock's
  # library so the outer PS fit actually matches the selected
  # candidate. Truncation is honored by fit_tmle_treatment_mechanism
  # via the locked spec.
  lock$sl_library <- selected$g_library

  ps   <- fit_ps_superlearner(lock, truncate = selected$truncation)
  g    <- fit_tmle_treatment_mechanism(lock, ps)
  Qf   <- fit_tmle_outcome_mechanism(lock, g)
  upd  <- run_tmle_targeting_step(g, Qf)
  out  <- .get_est(extract_tmle_estimate(upd))
  out$workflow            <- "plasmode"
  out$library_label       <- paste(unlist(selected$g_library), collapse = "+")
  out$truncation          <- selected$truncation
  out$selected_candidate  <- if (fallback)
    paste0(selected_id, "_FALLBACK") else selected_id
  out
}
