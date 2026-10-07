# The answer key on one dataset: the candidates fitted on the real (simulated)
# outcome with the same fitter the map uses, for every estimand and truncation.
# Folds and seed are those of the map run on the same covariates and treatment.

run_key <- function(scenario, n, idx) {
  t0 <- Sys.time()
  d <- draw_data(scenario, n, idx)
  seed <- plan_seed(scenario, n, idx)
  plan <- make_plan(seed, reps = 2L)  # reps is not used here
  X <- cleanTMLE:::.design_matrix(d, COVARIATES)
  folds <- cleanTMLE:::.make_folds(n, plan$K, seed)
  res <- tryCatch(
    cleanTMLE:::.with_seed(seed, cleanTMLE:::.fit_and_target(
      X, d$a, d$y, plan$candidates$library$glm, ESTIMANDS, plan$candidates$truncation,
      plan, folds, seed)),
    error = function(e) cleanTMLE:::.failed_rows(plan, conditionMessage(e)))
  out <- data.frame(scenario = scenario$id, overlap = scenario$overlap, n = n, idx = idx,
                    estimand = res$estimand,
                    truncation = res$truncation,
                    candidate = sprintf("glm_t%s", res$truncation),
                    estimate = res$estimate, se = res$se, ci_lower = res$ci_lower,
                    ci_upper = res$ci_upper, failed = res$failed,
                    stringsAsFactors = FALSE)
  # Claim 3 bookkeeping: the two-sided p-value of each estimate, and, within each
  # candidate, which estimand has the smallest p-value (the outcome-informed choice).
  out$p_value <- 2 * stats::pnorm(-abs(out$estimate / out$se))
  out$min_p <- FALSE
  for (cand in unique(out$candidate)) {
    k <- which(out$candidate == cand & is.finite(out$p_value))
    if (length(k)) out$min_p[k[which.min(out$p_value[k])]] <- TRUE
  }
  out$sec_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out
}
