# The outcome-free support map on one dataset, in two variants:
#   generate: cleanTMLE::simulate_design() as implemented (treatment drawn
#             from the single-fit propensity score in every repetition);
#   sample:   the same loop, but each resampled patient keeps the treatment
#             observed in the dataset.
# Both read the design diagnostics from the same assess_design() call.

COVARIATES <- c("X1", "X2", "X3", "X4")

plan_seed <- function(scenario, n, idx) seed_from("plan", scenario$id, n, idx)

make_plan <- function(seed, reps) {
  cleanTMLE::analysis_plan(
    outcome = "y", estimands = ESTIMANDS, ps_library = "glm",
    candidates = list(library = list(glm = "glm"), truncation = c(0.01, 0.05)),
    trim_band = TRIM_BAND, K = 5L, reps = reps, max_reps = reps, seed = seed)
}

# Local copy of cleanTMLE:::.one_rep (0.4.1) for the sample-treatment variant.
# Without `pairs` (no extension batches run, since max_reps = reps). The single
# difference is marked; the generate-treatment draw is still made and
# discarded, so the resampled rows, folds and outcome draws use the same random
# numbers as repetition r of the generate variant.
one_rep_sample <- function(r, X, g, surfaces, plan, A_obs) {
  pairs <- expand.grid(library = names(plan$candidates$library), surface = names(surfaces),
                       stringsAsFactors = FALSE)[c("surface", "library")]
  n <- nrow(X)
  seed <- plan$seed + r
  draw <- cleanTMLE:::.with_seed(seed, {
    idx <- sample.int(n, n, replace = TRUE)
    a_gen <- stats::rbinom(n, 1L, g[idx])  # drawn and discarded (keeps the RNG stream)
    a <- as.integer(A_obs[idx])            # THE DIFFERENCE: observed treatment kept
    ys <- lapply(surfaces, function(s)
      stats::rbinom(n, 1L, ifelse(a == 1L, s$q1[idx], s$q0[idx])))
    list(idx = idx, a = a, ys = ys)
  })
  Xs <- X[draw$idx, , drop = FALSE]
  folds <- cleanTMLE:::.make_folds(n, plan$K, seed)[draw$idx]
  g_fit <- tryCatch(cleanTMLE:::.fit_g(Xs, draw$a, plan$ps_library, folds, plan$V, seed,
                                       groups = draw$idx),
                    error = function(e) conditionMessage(e))
  rows <- list()
  for (i in seq_len(nrow(pairs))) {
    sn <- pairs$surface[i]
    lib <- pairs$library[i]
    res <- if (is.character(g_fit)) cleanTMLE:::.failed_rows(plan, g_fit) else tryCatch(
      cleanTMLE:::.with_seed(seed, cleanTMLE:::.fit_and_target(
        Xs, draw$a, draw$ys[[sn]], plan$candidates$library[[lib]],
        plan$estimands, plan$candidates$truncation, plan, folds, seed, g_fit,
        groups = draw$idx)),
      error = function(e) cleanTMLE:::.failed_rows(plan, conditionMessage(e)))
    rows[[length(rows) + 1L]] <- cbind(rep = r, surface = sn, library = lib, res)
  }
  do.call(rbind, rows)
}

# The sample-treatment map: the package's truths, metrics and verdict rule on
# the results of one_rep_sample().
simulate_sample <- function(lock, design) {
  p <- lock$plan
  X <- cleanTMLE:::.design_matrix(lock$data, lock$covariates)
  g <- design$g_world
  surfaces <- cleanTMLE:::.surfaces(p, X, g, p$seed)
  truths <- do.call(rbind, lapply(names(surfaces), function(sn) {
    tv <- cleanTMLE:::.true_values(surfaces[[sn]]$q1, surfaces[[sn]]$q0, g, p$trim_band)
    data.frame(surface = sn, estimand = names(tv), truth = unname(tv),
               stringsAsFactors = FALSE)
  }))
  results <- do.call(rbind, lapply(seq_len(p$reps), function(r)
    one_rep_sample(r, X, g, surfaces, p, design$A)))
  metrics <- cleanTMLE:::.sim_metrics(results, truths, p)
  list(verdict = cleanTMLE:::.sim_verdict(metrics, p), metrics = metrics)
}

# Overlap diagnostics per estimand from one assess_design() result.
design_diagnostics <- function(design, n) {
  g <- design$g
  A <- design$A
  P <- mean(A)
  # Walker's preference score: the propensity score with the treated share
  # removed on the logit scale.
  pref <- stats::plogis(stats::qlogis(g) - stats::qlogis(P))
  e <- design$ess
  data.frame(estimand = e$estimand, share_outside = design$share_outside_band,
             max_weight = e$max_weight, ess_frac = (e$ess_treated + e$ess_control) / n,
             pref_share = mean(pref >= 0.3 & pref <= 0.7),
             c_statistic = design$c_statistic, stringsAsFactors = FALSE)
}

run_map <- function(scenario, n, idx, variant = c("generate", "sample"), reps = 200L) {
  variant <- match.arg(variant)
  t0 <- Sys.time()
  d <- draw_data(scenario, n, idx)
  plan <- make_plan(plan_seed(scenario, n, idx), reps)
  lock <- cleanTMLE::create_analysis_lock(d[c("id", "a", COVARIATES)], treatment = "a",
                                          covariates = COVARIATES, plan = plan, id = "id")
  design <- cleanTMLE::assess_design(lock)
  t1 <- Sys.time()
  sim <- if (variant == "generate") {
    s <- suppressMessages(cleanTMLE::simulate_design(lock, design))
    list(verdict = s$verdict, metrics = s$metrics)
  } else simulate_sample(lock, design)
  t2 <- Sys.time()
  v <- sim$verdict
  out <- merge(v[c("estimand", "status", "selected", "ci_halfwidth", "unresolved")],
               design_diagnostics(design, n), by = "estimand", sort = FALSE)
  out <- out[match(ESTIMANDS, out$estimand), ]
  out <- cbind(scenario = scenario$id, n = n, idx = idx, variant = variant, reps = reps,
               out, sec_assess = as.numeric(difftime(t1, t0, units = "secs")),
               sec_total = as.numeric(difftime(t2, t0, units = "secs")),
               stringsAsFactors = FALSE)
  rownames(out) <- NULL
  attr(out, "metrics") <- cbind(scenario = scenario$id, n = n, idx = idx, variant = variant,
                                sim$metrics, stringsAsFactors = FALSE)
  out
}
