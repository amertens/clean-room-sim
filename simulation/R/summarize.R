# Summaries behind Table 1 and Figure 1.

# Population truths for every scenario, long format.
truth_table <- function(scenarios, Xp = pop_covariates()) {
  do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    tv <- population_truth(scenarios[i, ], Xp)
    data.frame(scenario = scenarios$id[i], estimand = names(tv), truth = unname(tv),
               stringsAsFactors = FALSE)
  }))
}

# Answer key per (scenario, n, estimand, candidate): bias and coverage against
# the population truth over the key datasets, with Monte Carlo SEs. A candidate
# passes when |bias| <= tol_bias and coverage >= tol_cov on the point estimates,
# with at most 5% failed fits (the package's rule for a usable cell). An
# estimand is estimable when some candidate passes.
answer_key <- function(key_rows, truths, tol_bias = 0.02, tol_cov = 0.90) {
  k <- merge(key_rows, truths, by = c("scenario", "estimand"))
  grp <- unique(k[c("scenario", "n", "estimand", "candidate")])
  cand <- do.call(rbind, lapply(seq_len(nrow(grp)), function(i) {
    g <- grp[i, ]
    s <- k[k$scenario == g$scenario & k$n == g$n & k$estimand == g$estimand &
             k$candidate == g$candidate, ]
    ok <- s$failed %in% FALSE & is.finite(s$estimate) & is.finite(s$ci_lower) &
      is.finite(s$ci_upper)
    R <- sum(ok)
    e <- s$estimate[ok] - s$truth[ok]
    cv <- s$ci_lower[ok] <= s$truth[ok] & s$truth[ok] <= s$ci_upper[ok]
    bias <- if (R) mean(e) else NA_real_
    cov <- if (R) mean(cv) else NA_real_
    data.frame(g, M = nrow(s), failed = nrow(s) - R, truth = s$truth[1],
               bias = bias, bias_mcse = if (R > 1) stats::sd(e) / sqrt(R) else NA_real_,
               coverage = cov,
               coverage_mcse = if (R > 1) sqrt(cov * (1 - cov) / R) else NA_real_,
               halfwidth = if (R) mean((s$ci_upper[ok] - s$ci_lower[ok]) / 2) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  cand$pass <- !is.na(cand$bias) & cand$failed / cand$M <= 0.05 &
    abs(cand$bias) <= tol_bias & cand$coverage >= tol_cov
  est <- stats::aggregate(pass ~ scenario + n + estimand, cand, any)
  names(est)[names(est) == "pass"] <- "estimable"
  list(candidates = cand, estimable = est)
}

# Sensitivity and specificity of a binary verdict (TRUE = feasible) for
# estimable, with binomial Monte Carlo SEs, by group.
sens_spec <- function(verdict, estimable, group) {
  out <- do.call(rbind, lapply(split(seq_along(verdict), group), function(i) {
    v <- verdict[i]; e <- estimable[i]
    se <- function(p, m) if (m > 0) sqrt(p * (1 - p) / m) else NA_real_
    sens <- if (any(e)) mean(v[e]) else NA_real_
    spec <- if (any(!e)) mean(!v[!e]) else NA_real_
    data.frame(n_estimable = sum(e), n_not = sum(!e),
               sensitivity = sens, sensitivity_mcse = se(sens, sum(e)),
               specificity = spec, specificity_mcse = se(spec, sum(!e)))
  }))
  cbind(group = names(split(seq_along(verdict), group)), out, row.names = NULL)
}

# Direction of each diagnostic: TRUE when larger values mean better overlap.
DIAGNOSTICS <- c(share_outside = FALSE, max_weight = FALSE, ess_frac = TRUE, pref_share = TRUE)

# Threshold maximizing Youden's index over every scenario-dataset-estimand.
youden_threshold <- function(x, estimable, higher_better) {
  s <- if (higher_better) x else -x
  cuts <- sort(unique(s))
  J <- vapply(cuts, function(t) mean(s[estimable] >= t) + mean(s[!estimable] < t) - 1,
              numeric(1))
  t <- cuts[which.max(J)]
  list(threshold = if (higher_better) t else -t, youden = max(J))
}

# Agreement of every verdict source with the answer key, by model relation.
# `map_rows` binds the rows of both map stages; the diagnostics are read from
# the generate rows (both variants share the same assess_design() call).
agreement <- function(map_rows, estimable, scenarios) {
  m <- merge(map_rows, estimable, by = c("scenario", "n", "estimand"))
  m$relation <- scenarios$model[match(m$scenario, scenarios$id)]
  out <- list()
  for (v in c("generate", "sample")) {
    d <- m[m$variant == v, ]
    if (nrow(d)) out[[length(out) + 1L]] <- cbind(
      verdict = paste("map", v), threshold = NA_real_,
      sens_spec(d$status == "feasible", d$estimable, d$relation))
  }
  d <- m[m$variant == "generate", ]
  for (dg in names(DIAGNOSTICS)) {
    if (!nrow(d) || length(unique(d$estimable)) < 2L) next
    ok <- is.finite(d[[dg]])
    yt <- youden_threshold(d[[dg]][ok], d$estimable[ok], DIAGNOSTICS[[dg]])
    call_ok <- if (DIAGNOSTICS[[dg]]) d[[dg]] >= yt$threshold else d[[dg]] <= yt$threshold
    out[[length(out) + 1L]] <- cbind(verdict = dg, threshold = yt$threshold,
                                     sens_spec(call_ok %in% TRUE, d$estimable, d$relation))
  }
  do.call(rbind, out)
}
