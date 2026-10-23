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

# Calls an analyst could make without the answer key of their own study
# (PLAN.md decision rules, claim 1a's second and third comparators): the
# package's overlap grade (good below 5% outside the band, moderate below 20%;
# cleanTMLE assess_design()), read as estimable for the ATE and ATT when good
# and for the trimmed ATE and ATO when good or moderate; and the diagnostic
# that was best at the other sample size, with the threshold tuned there. With
# one sample size in rt, as in the contrast pilot, the carried calls stay NA.
add_analyst_calls <- function(rt) {
  thr <- attr(rt, "thresholds")
  rt$call_grade <- ifelse(rt$estimand %in% c("ATE", "ATT"), rt$share_outside < 0.05,
                          rt$share_outside < 0.20)
  rt$call_carried <- NA
  ns <- unique(rt$n)
  for (nn in ns) for (e in unique(rt$estimand)) {
    other <- setdiff(ns, nn)
    if (!length(other)) next
    t <- thr[thr$n == other[1] & thr$estimand == e & is.finite(thr$threshold), ]
    if (!nrow(t)) next
    t <- t[which.max(t$youden), ]
    k <- rt$n == nn & rt$estimand == e
    x <- rt[[t$diagnostic]][k]
    rt$call_carried[k] <- if (DIAGNOSTICS[[t$diagnostic]]) x >= t$threshold else x <= t$threshold
  }
  attr(rt, "thresholds") <- thr
  rt
}

SOURCES <- c("generate", "sample", names(DIAGNOSTICS), "grade", "carried")
sources_in <- function(rt) SOURCES[paste0("call_", SOURCES) %in% names(rt)]

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
    do.call(rbind, lapply(sources_in(rt), function(src) {
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

# ---------------------------------------------------------------------------
# PLAN.md, "Decision rules, fixed before the contrast pilot is summarized"
# (2026-10-23). Every rule reads a per-index difference d, paired on the same
# datasets, and its 95% interval mean(d) +/- 1.96 sd(d) / sqrt(R). Sizing for
# the full run: k = 3 for each "better" or "falls" comparison (n_sim_safe), and
# the equivalence bound (2 * 2 * sd(d) / margin)^2 of pilot_equivalence() for
# each margin (n_sim_margin).
# ---------------------------------------------------------------------------

MARGIN_CLAIM1 <- 0.05
MARGIN_CLAIM3 <- 0.03
INSIDE_NEAR <- c("inside", "near")
OUTSIDE <- c("outside", "nonlinear")

paired_interval <- function(d) {
  d <- d[is.finite(d)]
  R <- length(d)
  m <- if (R) mean(d) else NA_real_
  se <- if (R > 1) stats::sd(d) / sqrt(R) else NA_real_
  data.frame(datasets = R, diff = m, diff_mcse = se, lower = m - 1.96 * se,
             upper = m + 1.96 * se)
}

sizing <- function(d, measure, a, b, margin, k = 3) {
  pv <- check_d(d, measure, a, b, k)
  data.frame(pilot_z = pv$pilot_z, n_sim_safe = pv$n_sim_safe,
             n_sim_margin = if (is.finite(pv$sd_d) && pv$sd_d > 0)
               ceiling((2 * 2 * pv$sd_d / margin)^2) else NA_real_)
}

# Per-index balanced accuracy of a source over the rows k of rt: the mean of
# its sensitivity and specificity within each dataset index, across the
# overlap levels and truths in k. NULL when k holds one label class only, so
# the estimand tests nothing there.
per_idx_balanced <- function(rt, k, src) {
  if (length(unique(rt$estimable[k])) < 2L) return(NULL)
  sens <- per_idx_correct(rt, k & rt$estimable, src)
  spec <- per_idx_correct(rt, k & !rt$estimable, src)
  ids <- intersect(names(sens), names(spec))
  (sens[ids] + spec[ids]) / 2
}

# Claim 1a, per (n, estimand), over the truths inside or near the family: the
# verdict's balanced accuracy against the best tuned diagnostic (chosen in
# this cell, generous to the comparator), the diagnostic carried over from the
# other sample size, and the package's grade. "better" when the interval's
# lower bound is above 0; against the best tuned diagnostic also "no worse"
# when it is above -margin; otherwise "unresolved".
claim1a_rules <- function(rt, margin = MARGIN_CLAIM1, k = 3) {
  out <- list()
  for (nn in unique(rt$n)) for (e in unique(rt$estimand)) {
    kin <- rt$n == nn & rt$estimand == e & rt$model %in% INSIDE_NEAR
    vg <- per_idx_balanced(rt, kin, "generate")
    if (is.null(vg)) {
      out[[length(out) + 1L]] <- data.frame(n = nn, estimand = e, comparator = NA_character_,
                                            source = NA_character_, verdict_value = NA_real_,
                                            comparator_value = NA_real_,
                                            rule = "uninformative: the label does not vary")
      next
    }
    ba <- vapply(names(DIAGNOSTICS), function(dg) {
      v <- per_idx_balanced(rt, kin, dg)
      if (is.null(v)) NA_real_ else mean(v)
    }, numeric(1))
    comps <- c(best_tuned = names(which.max(ba)), carried = "carried", grade = "grade")
    for (cn in names(comps)) {
      src <- comps[[cn]]
      calls <- rt[[paste0("call_", src)]]
      if (is.null(calls) || all(is.na(calls[kin]))) {
        out[[length(out) + 1L]] <- data.frame(n = nn, estimand = e, comparator = cn, source = src,
                                              verdict_value = mean(vg), comparator_value = NA_real_,
                                              rule = "not available")
        next
      }
      vo <- per_idx_balanced(rt, kin, src)
      ids <- intersect(names(vg), names(vo))
      d <- vg[ids] - vo[ids]
      iv <- paired_interval(d)
      rule <- if (!is.finite(iv$lower)) "unresolved" else if (iv$lower > 0) "better" else
        if (cn == "best_tuned" && iv$lower > -margin) "no worse" else "unresolved"
      out[[length(out) + 1L]] <- cbind(
        data.frame(n = nn, estimand = e, comparator = cn, source = src,
                   verdict_value = mean(vg[ids]), comparator_value = mean(vo[ids]), rule = rule),
        iv, sizing(d, "balanced accuracy", "generate", src, margin, k))
    }
  }
  rx_bind(out)
}

# Claim 1a per sample size: the stronger form needs "better" against all three
# comparators for every informative estimand; the weaker form needs "better"
# or "no worse" against the best tuned diagnostic and "better" against the
# other two. A comparator that is "not available" leaves the form undecided.
claim1a_summary <- function(r1) {
  do.call(rbind, lapply(unique(r1$n), function(nn) {
    s <- r1[r1$n == nn & !startsWith(r1$rule, "uninformative"), ]
    if (!nrow(s)) return(data.frame(n = nn, informative = "", stronger = NA, weaker = NA))
    avail <- !any(s$rule == "not available")
    strong <- all(s$rule == "better")
    weak <- all(s$rule[s$comparator == "best_tuned"] %in% c("better", "no worse")) &&
      all(s$rule[s$comparator != "best_tuned"] == "better")
    data.frame(n = nn, informative = paste(unique(s$estimand), collapse = ", "),
               stronger = if (avail) strong else NA, weaker = if (avail) weak else NA)
  }))
}

# Claim 1b, per (n, estimand): the verdict's balanced accuracy for the truths
# outside the family minus that for the truths inside or near it. "falls" when
# the upper bound is below 0; "not shown to fall" when the interval lies within
# +/- margin; otherwise "unresolved".
claim1b_rule <- function(rt, margin = MARGIN_CLAIM1, k = 3) {
  out <- list()
  for (nn in unique(rt$n)) for (e in unique(rt$estimand)) {
    base <- rt$n == nn & rt$estimand == e
    vin <- per_idx_balanced(rt, base & rt$model %in% INSIDE_NEAR, "generate")
    vout <- per_idx_balanced(rt, base & rt$model %in% OUTSIDE, "generate")
    if (is.null(vin) || is.null(vout)) {
      out[[length(out) + 1L]] <- data.frame(n = nn, estimand = e,
                                            rule = "uninformative: a label does not vary")
      next
    }
    ids <- intersect(names(vin), names(vout))
    d <- vout[ids] - vin[ids]
    iv <- paired_interval(d)
    rule <- if (!is.finite(iv$upper)) "unresolved" else if (iv$upper < 0) "falls" else
      if (iv$lower > -margin && iv$upper < margin) "not shown to fall" else "unresolved"
    out[[length(out) + 1L]] <- cbind(
      data.frame(n = nn, estimand = e, inside_near = mean(vin[ids]), outside = mean(vout[ids]),
                 rule = rule),
      iv, sizing(d, "balanced accuracy", "outside", "inside or near", margin, k))
  }
  rx_bind(out)
}

# Binds data frames whose columns differ, filling the missing ones with NA.
rx_bind <- function(l) {
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
}

# Claim 3: per scenario and dataset, the prespecified choice (the first
# feasible estimand of the ladder in that dataset's generate map) and two
# outcome-informed choices among the same feasible estimands, the smallest
# p-value and the narrowest CI, with coverage of the chosen estimand's truth
# and rejection of zero at 0.05. One candidate (glm_t0.01). key.R's min_p flag
# is not read, because it ranges over all four estimands, feasible or not; the
# choice is made here so that key.R, whose text is in the label fits' cache
# key, stays unchanged.
claim3_rows <- function(key_rows, map_rows, truths, scenarios) {
  g <- map_rows[map_rows$variant == "generate", ]
  k <- merge(key_rows[key_rows$candidate == sprintf("glm_t%s", TRUNCATIONS[1]), ],
             truths, by = c("scenario", "estimand"))
  k <- k[is.finite(k$p_value) & is.finite(k$ci_lower) & is.finite(k$ci_upper), ]
  k$covers <- k$ci_lower <= k$truth & k$truth <= k$ci_upper
  k$rejects <- k$p_value <= 0.05
  keys <- unique(k[c("scenario", "overlap", "n", "idx")])
  keys <- merge(keys, unique(g[c("overlap", "n", "idx")]), by = c("overlap", "n", "idx"))
  pick <- function(x) list(estimand = if (nrow(x)) x$estimand[1] else NA_character_,
                                covers = if (nrow(x)) x$covers[1] else NA,
                                rejects = if (nrow(x)) x$rejects[1] else NA)
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    r <- keys[i, ]
    kr <- k[k$scenario == r$scenario & k$n == r$n & k$idx == r$idx, ]
    gm <- g[g$overlap == r$overlap & g$n == r$n & g$idx == r$idx, ]
    feas <- ESTIMANDS[ESTIMANDS %in% gm$estimand[gm$status == "feasible"]]
    kf <- kr[kr$estimand %in% feas, ]
    kf <- kf[order(match(kf$estimand, ESTIMANDS)), ]   # ties go to the earlier estimand
    pre <- pick(kf[kf$estimand %in% feas[1], ])
    mp <- pick(kf[which.min(kf$p_value), ])
    nw <- pick(kf[which.min(kf$ci_upper - kf$ci_lower), ])
    data.frame(r, effect = scenarios$effect[match(r$scenario, scenarios$id)],
               model = scenarios$model[match(r$scenario, scenarios$id)],
               n_feasible = length(feas),
               pre_estimand = pre$estimand, pre_covers = pre$covers, pre_rejects = pre$rejects,
               minp_estimand = mp$estimand, minp_covers = mp$covers, minp_rejects = mp$rejects,
               narrow_estimand = nw$estimand, narrow_covers = nw$covers, narrow_rejects = nw$rejects,
               stringsAsFactors = FALSE)
  }))
}

# The share of scenario-datasets with 0 to 4 feasible estimands, per (n,
# overlap); the claim 3 rules use those with at least two.
claim3_feasible <- function(c3) {
  tab <- stats::aggregate(list(rows = c3$n_feasible), c3[c("n", "overlap", "n_feasible")], length)
  tab$share <- tab$rows / stats::ave(tab$rows, tab$n, tab$overlap, FUN = sum)
  tab
}

# Claim 3 rules, on the scenario-datasets with at least two feasible
# estimands, per (n, overlap, effect) and pooled over overlap levels per (n,
# effect). Each rule's difference is its rate minus the prespecified choice's,
# averaged over the scenarios of a dataset index. Smallest p-value, null
# scenarios: "inflates the type I error" when the lower bound is above 0.
# Narrowest CI: "within margin" when the interval lies within +/- margin; the
# claim that selecting on width changes neither measure needs "within margin"
# for both coverage (effect) and rejection (null).
claim3_rules <- function(c3, margin = MARGIN_CLAIM3, k = 3) {
  c3 <- c3[c3$n_feasible >= 2 & !is.na(c3$pre_estimand) & !is.na(c3$minp_estimand) &
             !is.na(c3$narrow_estimand), ]
  per_level <- unique(c3[c("n", "overlap", "effect")])
  pooled <- unique(c3[c("n", "effect")]); pooled$overlap <- "pooled"
  grp <- rbind(per_level, pooled[names(per_level)])
  out <- list()
  for (i in seq_len(nrow(grp))) {
    g <- grp[i, ]
    s <- c3[c3$n == g$n & c3$effect == g$effect & (g$overlap == "pooled" | c3$overlap == g$overlap), ]
    if (!nrow(s)) next
    what <- if (g$effect == "null") "rejects" else "covers"
    measure <- if (what == "covers") "coverage" else "rejection"
    pre <- as.numeric(s[[paste0("pre_", what)]])
    for (rl in c("minp", "narrow")) {
      oth <- as.numeric(s[[paste0(rl, "_", what)]])
      d <- tapply(oth - pre, s$idx, mean)
      iv <- paired_interval(d)
      rule <- if (rl == "minp") {
        if (what == "rejects" && is.finite(iv$lower) && iv$lower > 0) "inflates the type I error" else
          if (what == "rejects") "not shown" else "reported"
      } else if (is.finite(iv$lower) && iv$lower > -margin && iv$upper < margin) "within margin" else
        "no claim"
      out[[length(out) + 1L]] <- cbind(
        data.frame(g, measure = measure, choice = if (rl == "minp") "smallest p" else "narrowest CI",
                   prespecified = mean(pre), outcome_informed = mean(oth), rows = nrow(s),
                   rule = rule, stringsAsFactors = FALSE),
        iv, sizing(d, measure, rl, "prespecified", margin, k))
    }
  }
  rx_bind(out)
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
