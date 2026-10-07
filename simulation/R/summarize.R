# Summaries behind Table 1, Figure 1, Table S1 and Figure 2, and the pilot's
# contrast checks. Monte Carlo SEs treat the dataset index as the repetition:
# datasets with different indices are independent, while one map run (one
# overlap level and index) is shared by the six outcome-model scenarios of its
# level, so per-index means are formed first and their spread gives the SE.

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

# Direction of each diagnostic: TRUE when larger values mean better overlap.
DIAGNOSTICS <- c(share_outside = FALSE, max_weight = FALSE, ess_frac = TRUE, pref_share = TRUE)

# Threshold maximizing Youden's index (sensitivity + specificity - 1).
youden_threshold <- function(x, estimable, higher_better) {
  ok <- is.finite(x)
  x <- x[ok]; estimable <- estimable[ok]
  if (!length(x) || all(estimable) || !any(estimable))
    return(list(threshold = NA_real_, youden = NA_real_))
  s <- if (higher_better) x else -x
  cuts <- sort(unique(s))
  J <- vapply(cuts, function(t) mean(s[estimable] >= t) + mean(s[!estimable] < t) - 1,
              numeric(1))
  t <- cuts[which.max(J)]
  list(threshold = if (higher_better) t else -t, youden = max(J))
}

# One row per (scenario, n, idx, estimand) with the answer-key label, both map
# verdicts and the diagnostics of that dataset's assess_design() call.
run_table <- function(map_rows, estimable, scenarios) {
  g <- map_rows[map_rows$variant == "generate", ]
  s <- map_rows[map_rows$variant == "sample", c("overlap", "n", "idx", "estimand", "status")]
  names(s)[5] <- "status_sample"
  m <- merge(g, s, by = c("overlap", "n", "idx", "estimand"), all.x = TRUE)
  m <- merge(scenarios[c("id", "overlap", "model", "effect")], m, by = "overlap")
  names(m)[names(m) == "id"] <- "scenario"
  m <- merge(m, estimable, by = c("scenario", "n", "estimand"))
  m$call_generate <- m$status == "feasible"
  m$call_sample <- m$status_sample == "feasible"
  m
}

# Diagnostic calls at each diagnostic's best threshold, chosen per (n,
# estimand) across every scenario and dataset (all model relations together).
add_diagnostic_calls <- function(rt) {
  thr <- list()
  for (dg in names(DIAGNOSTICS)) rt[[paste0("call_", dg)]] <- NA
  for (nn in unique(rt$n)) for (e in unique(rt$estimand)) {
    k <- rt$n == nn & rt$estimand == e
    for (dg in names(DIAGNOSTICS)) {
      yt <- youden_threshold(rt[[dg]][k], rt$estimable[k], DIAGNOSTICS[[dg]])
      x <- rt[[dg]][k]
      # When every scenario-dataset has the same label there is no threshold to
      # choose; the diagnostic then calls that label everywhere.
      lab <- unique(rt$estimable[k])
      rt[[paste0("call_", dg)]][k] <- if (!is.na(yt$threshold)) {
        if (DIAGNOSTICS[[dg]]) x >= yt$threshold else x <= yt$threshold
      } else if (length(lab) == 1L) lab else NA
      thr[[length(thr) + 1L]] <- data.frame(n = nn, estimand = e, diagnostic = dg,
                                            threshold = yt$threshold, youden = yt$youden)
    }
  }
  attr(rt, "thresholds") <- do.call(rbind, thr)
  rt
}

SOURCES <- c("generate", "sample", names(DIAGNOSTICS))

# Per-index mean of "call correct" over the rows `k` of rt, for one source.
per_idx_correct <- function(rt, k, src) {
  call <- rt[[paste0("call_", src)]][k] %in% TRUE
  correct <- as.numeric(call == rt$estimable[k])
  tapply(correct, rt$idx[k], mean)
}

# Sensitivity and specificity of each source per (n, estimand, relation).
agreement_table <- function(rt) {
  grp <- unique(rt[c("n", "estimand", "model")])
  do.call(rbind, lapply(seq_len(nrow(grp)), function(i) {
    g <- grp[i, ]
    base <- rt$n == g$n & rt$estimand == g$estimand & rt$model == g$model
    do.call(rbind, lapply(SOURCES, function(src) {
      one <- function(lab) {
        k <- base & rt$estimable == lab
        if (!any(k)) return(c(NA_real_, NA_real_, 0))
        v <- per_idx_correct(rt, k, src)
        c(mean(v), stats::sd(v) / sqrt(length(v)), length(unique(rt$scenario[k])))
      }
      se <- one(TRUE); sp <- one(FALSE)
      data.frame(g, source = src, sensitivity = se[1], sensitivity_mcse = se[2],
                 n_scen_estimable = se[3], specificity = sp[1], specificity_mcse = sp[2],
                 n_scen_not = sp[3], stringsAsFactors = FALSE)
    }))
  }))
}

# Pilot check on a per-index difference vector (pilot_check.R's .pilot_verdict).
check_d <- function(d, what, a, b, k = 5) {
  .pilot_verdict(d, what, a, b, k, NULL)
}

# The contrasts Table 1 needs, per (n, estimand) and measure (sensitivity,
# specificity): generate map vs the best diagnostic of that relation, and
# generate vs sample, per relation; and inside vs outside for the generate map.
table1_contrasts <- function(rt, agree, k = 5) {
  out <- list()
  for (nn in unique(rt$n)) for (e in unique(rt$estimand)) for (lab in c(TRUE, FALSE)) {
    measure <- if (lab) "sensitivity" else "specificity"
    for (rel in c("inside", "near", "outside")) {
      kk <- rt$n == nn & rt$estimand == e & rt$model == rel & rt$estimable == lab
      if (!any(kk)) next
      ag <- agree[agree$n == nn & agree$estimand == e & agree$model == rel &
                    agree$source %in% names(DIAGNOSTICS), ]
      score <- ifelse(is.na(ag$sensitivity), 0, ag$sensitivity) +
        ifelse(is.na(ag$specificity), 0, ag$specificity)
      best <- ag$source[which.max(score)]
      vg <- per_idx_correct(rt, kk, "generate")
      for (other in c(best, "sample")) {
        vo <- per_idx_correct(rt, kk, other)
        out[[length(out) + 1L]] <- cbind(
          data.frame(n = nn, estimand = e, relation = rel, measure = measure,
                     comparison = if (other == "sample") "generate vs sample" else
                       "generate vs best diagnostic", other = other,
                     a_value = mean(vg), b_value = mean(vo), stringsAsFactors = FALSE),
          check_d(vg - vo[names(vg)], measure, "generate", other, k))
      }
    }
    ki <- rt$n == nn & rt$estimand == e & rt$model == "inside" & rt$estimable == lab
    ko <- rt$n == nn & rt$estimand == e & rt$model == "outside" & rt$estimable == lab
    if (any(ki) && any(ko)) {
      vi <- per_idx_correct(rt, ki, "generate")
      vo <- per_idx_correct(rt, ko, "generate")
      ids <- intersect(names(vi), names(vo))
      out[[length(out) + 1L]] <- cbind(
        data.frame(n = nn, estimand = e, relation = "inside - outside", measure = measure,
                   comparison = "inside vs outside (generate)", other = "generate",
                   a_value = mean(vi), b_value = mean(vo), stringsAsFactors = FALSE),
        check_d(vi[ids] - vo[ids], measure, "inside", "outside", k))
    }
  }
  do.call(rbind, out)
}

# Claim 3: per scenario and dataset, the prespecified choice (the first
# feasible estimand of the ladder in that dataset's generate map) and the
# outcome-informed choice (smallest p-value), with coverage of the chosen
# estimand's truth and rejection of zero at 0.05. One candidate (glm_t0.01).
claim3_rows <- function(key_rows, map_rows, truths, scenarios) {
  g <- map_rows[map_rows$variant == "generate", ]
  k <- merge(key_rows[key_rows$candidate == sprintf("glm_t%s", TRUNCATIONS[1]), ],
             truths, by = c("scenario", "estimand"))
  k$covers <- k$ci_lower <= k$truth & k$truth <= k$ci_upper
  k$rejects <- k$p_value <= 0.05
  keys <- unique(k[c("scenario", "overlap", "n", "idx")])
  keys <- merge(keys, unique(g[c("overlap", "n", "idx")]), by = c("overlap", "n", "idx"))
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    r <- keys[i, ]
    kr <- k[k$scenario == r$scenario & k$n == r$n & k$idx == r$idx, ]
    gm <- g[g$overlap == r$overlap & g$n == r$n & g$idx == r$idx, ]
    feas <- ESTIMANDS[ESTIMANDS %in% gm$estimand[gm$status == "feasible"]]
    pre <- if (length(feas)) kr[kr$estimand == feas[1], ] else kr[0, ]
    mp <- kr[kr$min_p %in% TRUE, ]
    data.frame(r, effect = scenarios$effect[match(r$scenario, scenarios$id)],
               model = scenarios$model[match(r$scenario, scenarios$id)],
               pre_estimand = if (nrow(pre)) pre$estimand else NA_character_,
               pre_covers = if (nrow(pre)) pre$covers else NA,
               pre_rejects = if (nrow(pre)) pre$rejects else NA,
               minp_estimand = if (nrow(mp)) mp$estimand[1] else NA_character_,
               minp_covers = if (nrow(mp)) mp$covers[1] else NA,
               minp_rejects = if (nrow(mp)) mp$rejects[1] else NA,
               stringsAsFactors = FALSE)
  }))
}

# Claim 3 per (n, overlap, effect): coverage (effect) or rejection (null) of
# each choice, and the pilot check of the per-index difference prespecified
# minus smallest-p (averaged over the three outcome models of the level).
claim3_summary <- function(c3, k = 5) {
  c3 <- c3[!is.na(c3$pre_estimand) & !is.na(c3$minp_estimand), ]
  grp <- unique(c3[c("n", "overlap", "effect")])
  do.call(rbind, lapply(seq_len(nrow(grp)), function(i) {
    g <- grp[i, ]
    s <- c3[c3$n == g$n & c3$overlap == g$overlap & c3$effect == g$effect, ]
    what <- if (g$effect == "null") "rejects" else "covers"
    pre <- s[[paste0("pre_", what)]]; mp <- s[[paste0("minp_", what)]]
    d <- tapply(as.numeric(pre) - as.numeric(mp), s$idx, mean)
    cbind(g, measure = if (what == "covers") "coverage" else "rejection",
          prespecified = mean(pre), smallest_p = mean(mp), datasets = length(pre),
          check_d(d, if (what == "covers") "coverage" else "rejection",
                  "prespecified", "smallest p", k))
  }))
}

# Claim 4: the map's expected half-width for its selected candidate against the
# realized half-width of that candidate on the outcome, where the map said
# feasible.
claim4_rows <- function(key_rows, map_rows, scenarios) {
  g <- map_rows[map_rows$variant == "generate" & map_rows$status == "feasible", ]
  g <- merge(g, scenarios[c("id", "overlap", "model", "effect")], by = "overlap")
  names(g)[names(g) == "id"] <- "scenario"
  k <- key_rows
  k$realized <- (k$ci_upper - k$ci_lower) / 2
  m <- merge(g[c("scenario", "overlap", "model", "effect", "n", "idx", "estimand",
                 "selected", "ci_halfwidth")],
             k[c("scenario", "n", "idx", "estimand", "candidate", "realized")],
             by.x = c("scenario", "n", "idx", "estimand", "selected"),
             by.y = c("scenario", "n", "idx", "estimand", "candidate"))
  m$ratio <- m$realized / m$ci_halfwidth
  m
}

claim4_summary <- function(c4) {
  grp <- unique(c4[c("n", "estimand")])
  do.call(rbind, lapply(seq_len(nrow(grp)), function(i) {
    s <- c4[c4$n == grp$n[i] & c4$estimand == grp$estimand[i] & is.finite(c4$ratio), ]
    v <- tapply(s$ratio, s$idx, mean)
    data.frame(grp[i, ], rows = nrow(s), mean_ratio = mean(s$ratio),
               ratio_mcse = stats::sd(v) / sqrt(length(v)),
               median_ratio = stats::median(s$ratio),
               cor_expected_realized = if (nrow(s) > 2) stats::cor(s$ci_halfwidth, s$realized)
                 else NA_real_)
  }))
}
