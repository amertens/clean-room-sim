#!/usr/bin/env Rscript
# ============================================================================
# Candidate-Divergence Monte Carlo Demonstration  (standalone sandbox study)
# ----------------------------------------------------------------------------
# Shows that select_tmle_candidate(rule = "min_rmse") and
# select_tmle_candidate(rule = "min_max_rmse", dq_results = ...) choose
# DIFFERENT TMLE candidates once a prespecified data-quality threat sweep is
# taken into account.
#
# Candidates differ only in propensity-score truncation:
#   aggressive = 0.001, middle = 0.025, robust = 0.20  (g-library SL.glm).
#
# Prespecified threat sweep (three threats, increasing severity):
#   near_positivity        - covariate->treatment slopes amplified so a
#                            subgroup approaches deterministic treatment
#                            (estimated PS -> 0/1). First-class package threat
#                            in run_plasmode_dq_stress().
#   unmeasured_confounding - latent U shifting A and Y (package threat).
#   covariate_missingness  - MCAR + median imputation (package threat).
#
# The whole sweep runs through run_plasmode_dq_stress() and
# select_tmle_candidate(); this script is a reproducible driver over package
# functions, with no estimator machinery of its own.
#
# Result: heavy truncation (robust) minimises BASELINE RMSE (it sheds IPW
# variance), so min_rmse selects robust. But robust is the most fragile
# candidate under unmeasured confounding (clipping discards the extreme-weight
# observations that carry the residual-confounding signal), while aggressive is
# destroyed by the positivity stress. The MIDDLE truncation hedges both threats
# and has the lowest worst-case RMSE, so min_max_rmse selects middle. The
# minimax rule therefore changes the decision.
#
# Self-contained: touches nothing outside sandbox/candidate_divergence/.
# Usage:  Rscript divergence_study.R [smoke|direction|full]
# ============================================================================

suppressWarnings(suppressMessages(library(pkgload)))

.this_dir <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())

repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkg_dir   <- file.path(repo_root, "cleanTMLE")
res_dir   <- file.path(.this_dir, "results")
fig_dir   <- file.path(.this_dir, "figures")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

pkgload::load_all(pkg_dir, quiet = TRUE)

.flush <- function() if (!interactive()) flush(stdout())
`%||%` <- function(a, b) if (is.null(a)) b else a

MODE <- { a <- commandArgs(TRUE); if (length(a) >= 1L) a[[1]] else "smoke" }

cfg <- switch(MODE,
  smoke     = list(n_obs = 1000L, reps = 10L,  n_batches = 1L, n_truth = 50000L),
  direction = list(n_obs = 2000L, reps = 40L,  n_batches = 1L, n_truth = 100000L),
  full      = list(n_obs = 2000L, reps = 200L, n_batches = 5L, n_truth = 200000L),
  stop("Unknown MODE: ", MODE))

# ── Fixed design ────────────────────────────────────────────────────────────
SEED         <- 20260530L
OVERLAP      <- 1.6      # marginal overlap so truncation engages
EFFECT       <- -0.05
EFFECT_SIZES <- c(0.05)
COVARIATES   <- c("age", "sex", "biomarker", "comorbidity", "ckd")

# Prespecified threat severities.
POS_SLOPES   <- c(2.0, 3.0, 4.0)                     # near-positivity amplification
U_ORS        <- c(3.0, 4.0, 5.0, 6.0, 7.0, 8.0)      # unmeasured-confounding OR sweep
MISS_FRACS   <- c(0.10, 0.20)                        # covariate missingness

# Locked (prespecified) decision thresholds.
LOCK_THRESHOLDS <- list(
  max_abs_bias   = 0.02,
  min_coverage   = 0.88,
  max_rmse_ratio = 2.0    # RMSE under stress must stay <= 2.0x its own baseline
)

# ── DGP (verbatim from run_simulation.R; do NOT source the harness) ─────────
generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL,
                          U_prevalence = 0, U_trt_OR = 1, U_out_OR = 1) {
  if (!is.null(seed)) set.seed(seed)
  age         <- rnorm(n, mean = 55, sd = 10)
  sex         <- rbinom(n, 1, 0.55)
  biomarker   <- rnorm(n, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n, 1, 0.12)
  U <- if (U_prevalence > 0) rbinom(n, 1, U_prevalence) else rep(0L, n)
  lp_trt <- -0.5 +
    overlap_strength * (0.03 * (age - 55) + 0.8 * sex + 0.6 * biomarker +
                        0.5 * ckd + 0.3 * comorbidity) +
    log(U_trt_OR) * U
  treatment <- rbinom(n, 1, plogis(lp_trt))
  lp_out <- -2.5 + 0.015 * (age - 55) + 0.3 * sex + 0.2 * biomarker +
    0.6 * ckd + 0.25 * comorbidity + effect_size / 0.15 * treatment +
    log(U_out_OR) * U
  event_24 <- rbinom(n, 1, plogis(lp_out))
  lp_nc <- -1.0 + 0.01 * (age - 55) + 0.1 * sex + 0.15 * biomarker
  nc_outcome <- rbinom(n, 1, plogis(lp_nc))
  data.frame(age = round(age, 1), sex = sex, biomarker = round(biomarker, 3),
             comorbidity = comorbidity, ckd = ckd, treatment = treatment,
             event_24 = event_24, nc_outcome = nc_outcome,
             stringsAsFactors = FALSE)
}

compute_truth <- function(n_truth, overlap_strength, effect_size, seed,
                          U_prevalence = 0, U_out_OR = 1) {
  set.seed(seed)
  age         <- rnorm(n_truth, mean = 55, sd = 10)
  sex         <- rbinom(n_truth, 1, 0.55)
  biomarker   <- rnorm(n_truth, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n_truth, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n_truth, 1, 0.12)
  U <- if (U_prevalence > 0) rbinom(n_truth, 1, U_prevalence) else rep(0L, n_truth)
  lp_out_base <- -2.5 + 0.015 * (age - 55) + 0.3 * sex + 0.2 * biomarker +
    0.6 * ckd + 0.25 * comorbidity + log(U_out_OR) * U
  coef_trt <- effect_size / 0.15
  list(risk_1 = mean(plogis(lp_out_base + coef_trt)),
       risk_0 = mean(plogis(lp_out_base)),
       RD = mean(plogis(lp_out_base + coef_trt)) - mean(plogis(lp_out_base)))
}

make_candidates <- function() list(
  tmle_candidate("aggressive", "GLM PS, trunc = 0.001 (efficiency-first)",
                 g_library = "SL.glm", truncation = 0.001),
  tmle_candidate("middle",     "GLM PS, trunc = 0.025 (intermediate)",
                 g_library = "SL.glm", truncation = 0.025),
  tmle_candidate("robust",     "GLM PS, trunc = 0.20 (robustness-first)",
                 g_library = "SL.glm", truncation = 0.20))


# ── One batch: feasibility + combined DQ stress on a fresh lock seed ────────
run_one_batch <- function(ref_dat, batch_seed, reps, verbose = FALSE) {
  lock <- create_analysis_lock(
    data = ref_dat, treatment = "treatment", outcome = "event_24",
    covariates = COVARIATES, sl_library = c("SL.glm", "SL.mean"),
    plasmode_reps = reps, seed = batch_seed)
  lock <- attach_estimand(lock,
    description = "Effect of A on 24-month event risk",
    population = "Adults eligible at index",
    treatment_strategies = c("Treatment", "Control"),
    outcome_label = "Primary event by 24 months", followup = "24 months",
    contrast = "risk_difference",
    statistical_estimand = "E_W{E[Y|A=1,W] - E[Y|A=0,W]}")

  cands <- make_candidates()

  plas <- run_plasmode_feasibility(lock, tmle_candidates = cands,
            effect_sizes = EFFECT_SIZES, reps = reps, verbose = FALSE)

  # The near_positivity threat is now a first-class threat in the package, so
  # the whole sweep runs through run_plasmode_dq_stress and select_tmle_candidate
  # with no sandbox-specific machinery.
  dq <- run_plasmode_dq_stress(lock, tmle_candidates = cands,
            effect_sizes = EFFECT_SIZES, reps = reps,
            data_quality_scenarios = list(
              near_positivity = list(slopes = POS_SLOPES),
              unmeasured_confounding = list(U_prevalence = 0.20,
                  U_treatment_OR = U_ORS, U_outcome_OR = U_ORS),
              covariate_missingness = list(fractions = MISS_FRACS)),
            fit_timeout = 30, verbose = verbose)

  list(lock = lock, plas = plas, dq = dq)
}

baseline_rmse <- function(plas) {
  m <- plas$metrics; tapply(m$rmse, m$candidate, mean)
}
worst_rmse <- function(dq) {
  m <- dq$metrics; m <- m[m$scenario != "none", ]
  tapply(m$rmse, m$candidate, function(x) max(x, na.rm = TRUE))
}
worst_threat <- function(dq, cid) {
  m <- dq$metrics; m <- m[m$scenario != "none" & m$candidate == cid, ]
  m[which.max(m$rmse), c("scenario", "level", "rmse")]
}

# ============================================================================
cat("============================================================\n")
cat(sprintf("  Candidate-Divergence Study  [MODE = %s]\n", MODE))
cat(sprintf("  n=%d  reps=%d  batches=%d  overlap=%.2f  seed=%d\n",
            cfg$n_obs, cfg$reps, cfg$n_batches, OVERLAP, SEED))
cat("============================================================\n\n"); .flush()

truth <- compute_truth(cfg$n_truth, OVERLAP, EFFECT, seed = SEED)
cat(sprintf("True RD (no U): %.5f\n", truth$RD))
ref_dat <- generate_data(cfg$n_obs, OVERLAP, EFFECT, seed = SEED)
{
  ps <- predict(glm(reformulate(COVARIATES, "treatment"), ref_dat,
                    family = binomial()), type = "response")
  cat(sprintf("Reference PS range: [%.4f, %.4f]  frac<.10=%.1f%%  trt prev=%.3f\n\n",
              min(ps), max(ps), 100 * mean(ps < 0.10), mean(ref_dat$treatment)))
}
.flush()

# Batch lock seeds are spaced widely so the per-component rep seeds
# (lock$seed + rep_i [+ scenario/slope offsets]) never overlap across batches;
# this keeps the batches independent so the between-batch SD is a valid Monte
# Carlo standard error.
batch_seed_of <- function(b) SEED + b * 1000000L
batches <- vector("list", cfg$n_batches)
for (b in seq_len(cfg$n_batches)) {
  cat(sprintf("--- Batch %d/%d (lock seed = %d) ---\n", b, cfg$n_batches, batch_seed_of(b))); .flush()
  t0 <- Sys.time()
  batches[[b]] <- run_one_batch(ref_dat, batch_seed = batch_seed_of(b), reps = cfg$reps,
                                verbose = (MODE == "smoke"))
  cat(sprintf("    done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); .flush()
}

cand_order <- c("aggressive", "middle", "robust")
base_mat  <- sapply(batches, function(z) baseline_rmse(z$plas)[cand_order])
worst_mat <- sapply(batches, function(z) worst_rmse(z$dq)[cand_order])
rownames(base_mat) <- rownames(worst_mat) <- cand_order

# Per-cell (scenario x level x candidate) RMSE pooled across batches, with MC SE,
# for the degradation-gradient figure with Monte Carlo bands.
dq_cell_agg <- local({
  ml <- lapply(seq_along(batches), function(b) {
    m <- batches[[b]]$dq$metrics
    m <- m[m$scenario != "none", c("scenario", "level", "candidate", "rmse")]
    m$batch <- b; m
  })
  d <- do.call(rbind, ml)
  key <- with(d, paste(scenario, level, candidate, sep = "|"))
  ag <- do.call(rbind, lapply(split(d, key), function(g) {
    data.frame(scenario = g$scenario[1], level = g$level[1],
               candidate = g$candidate[1],
               rmse = mean(g$rmse),
               rmse_mcse = if (nrow(g) > 1) sd(g$rmse) / sqrt(nrow(g)) else NA_real_,
               n_batches = nrow(g), stringsAsFactors = FALSE)
  }))
  rownames(ag) <- NULL; ag
})
# Baseline RMSE per candidate pooled across batches (for the threshold lines).
baseline_pooled <- data.frame(candidate = cand_order,
  baseline_rmse = rowMeans(base_mat),
  baseline_mcse = if (cfg$n_batches > 1) apply(base_mat, 1, sd) / sqrt(cfg$n_batches) else NA_real_,
  stringsAsFactors = FALSE)

win_tbl <- do.call(rbind, lapply(seq_along(batches), function(b) {
  z <- batches[[b]]
  wmin <- suppressMessages(select_tmle_candidate(z$plas, rule = "min_rmse"))
  wmm  <- suppressMessages(select_tmle_candidate(z$plas, rule = "min_max_rmse",
                                                 dq_results = z$dq))
  data.frame(batch = b, win_min_rmse = wmin$candidate_id,
             win_min_max_rmse = wmm$candidate_id, stringsAsFactors = FALSE)
}))

B <- cfg$n_batches
agg <- data.frame(
  candidate = cand_order,
  baseline_rmse = rowMeans(base_mat),
  baseline_mcse = if (B > 1) apply(base_mat, 1, sd) / sqrt(B) else NA_real_,
  worst_rmse = rowMeans(worst_mat),
  worst_mcse = if (B > 1) apply(worst_mat, 1, sd) / sqrt(B) else NA_real_,
  stringsAsFactors = FALSE)
agg$baseline_threshold <- LOCK_THRESHOLDS$max_rmse_ratio * agg$baseline_rmse

cat("\n================  RESULTS  ================\n")
cat("\nBaseline RMSE per batch:\n");  print(round(base_mat, 5))
cat("\nWorst-case RMSE per batch:\n"); print(round(worst_mat, 5))
cat("\nAggregated (mean across batches, MC SE = sd/sqrt(B)):\n")
print(agg, row.names = FALSE)
cat("\nWinners per batch:\n"); print(win_tbl, row.names = FALSE)
cat("\nWorst threat per candidate (batch 1):\n")
for (cid in cand_order) {
  wt <- worst_threat(batches[[1]]$dq, cid)
  cat(sprintf("  %-11s %s [%s]  rmse=%.5f\n", cid, wt$scenario, wt$level, wt$rmse))
}

win_min_final <- names(sort(table(win_tbl$win_min_rmse), decreasing = TRUE))[1]
win_mm_final  <- names(sort(table(win_tbl$win_min_max_rmse), decreasing = TRUE))[1]
cat(sprintf("\nMODAL WINNERS:  min_rmse = %s   |   min_max_rmse = %s   (%s)\n",
            win_min_final, win_mm_final,
            if (win_min_final != win_mm_final) "DIVERGE" else "agree"))

# Separation between the two rule winners on worst-case RMSE.
sep <- unname(agg$worst_rmse[agg$candidate == win_min_final] -
              agg$worst_rmse[agg$candidate == win_mm_final])
comb_se <- if (B > 1) sqrt(agg$worst_mcse[agg$candidate == win_min_final]^2 +
                           agg$worst_mcse[agg$candidate == win_mm_final]^2) else NA_real_
cat(sprintf("Worst-case RMSE (min_rmse winner '%s' - min_max winner '%s') = %.5f",
            win_min_final, win_mm_final, sep))
if (B > 1) cat(sprintf("   combined MC SE = %.5f   (ratio = %.1f)", comb_se, sep / comb_se))
cat("\n")

deg1 <- summarize_dq_degradation(batches[[1]]$dq)

out <- list(
  mode = MODE,
  config = c(cfg, list(seed = SEED, overlap = OVERLAP, effect = EFFECT,
             effect_sizes = EFFECT_SIZES, pos_slopes = POS_SLOPES, u_ors = U_ORS,
             miss_fracs = MISS_FRACS, lock_thresholds = LOCK_THRESHOLDS)),
  truth = truth, candidates = make_candidates(),
  batch_seeds = vapply(seq_len(cfg$n_batches), batch_seed_of, numeric(1)),
  base_mat = base_mat, worst_mat = worst_mat, agg = agg, win_tbl = win_tbl,
  win_min_final = win_min_final, win_mm_final = win_mm_final,
  separation = sep, combined_mcse = comb_se,
  plas_metrics = batches[[1]]$plas$metrics, dq_metrics = batches[[1]]$dq$metrics,
  dq_cell_agg = dq_cell_agg, baseline_pooled = baseline_pooled,
  degradation = deg1,
  saved_at = Sys.time())
saveRDS(out, file.path(res_dir, sprintf("candidate_divergence_%s.rds", MODE)))
cat(sprintf("\nSaved results/candidate_divergence_%s.rds\nDone.\n", MODE))
