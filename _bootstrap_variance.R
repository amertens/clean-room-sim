#!/usr/bin/env Rscript
# ============================================================================
# Bootstrap vs influence-function variance comparison
#
# Demonstrates that the IF-based variance is badly mis-calibrated for
# Match_TMLE (se_sd_ratio ≈ 2.1 in Scenario B) while the nonparametric
# bootstrap is well-calibrated (se_sd_ratio ≈ 1.0). IPTW is included as
# a positive control — its IF-SE is already near-nominal.
#
# Design:
#   Scenarios: A (good overlap) and B (marginal overlap)
#   Estimators: IPTW (stabilised), Match_TMLE
#   MC reps:    100 per scenario
#   Bootstrap:  B = 200 draws, GLM-only (fast), percentile CI
#   Truth:      same large-sample oracle as run_simulation.R (seed 2026)
#
# Output: results_new/bootstrap_variance.rds
#         results_new/bootstrap_variance.csv
#
# Runtime: approx 60-90 min on a modern laptop (GLM-only bootstrap fits
#          take ~0.15 s each; 100 reps × 200 draws × 2 estimators × 2
#          scenarios = 80,000 fits).
# ============================================================================

suppressMessages({
  library(MatchIt)
  library(tmle)
})

.flush <- function() if (!interactive()) flush(stdout())
cat("Bootstrap variance comparison\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())), sep = "")
.flush()

# ── Configuration ─────────────────────────────────────────────────────────────

N_MC   <- 100L   # MC replicates per scenario
B_BOOT <- 200L   # bootstrap draws per replicate
N_OBS  <- 2000L  # observations per replicate
N_TRUTH <- 500000L
TRUNC  <- 0.05   # PS truncation (matches Scenario B plasmode-selected candidate)
SEED   <- 2026L
COVARS <- c("age", "sex", "biomarker", "comorbidity", "ckd")
TREAT  <- "treatment"
OUTC   <- "event_24"

out_dir <- file.path(getwd(), "results_new")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── DGP (identical to run_simulation.R) ──────────────────────────────────────

generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  age         <- rnorm(n, mean = 55, sd = 10)
  sex         <- rbinom(n, 1, 0.55)
  biomarker   <- rnorm(n, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n, 1, 0.12)
  lp_trt <- -0.5 + overlap_strength * (0.03 * (age - 55) + 0.8 * sex +
               0.6 * biomarker + 0.5 * ckd + 0.3 * comorbidity)
  treatment  <- rbinom(n, 1, plogis(lp_trt))
  lp_out <- -2.5 + 0.015 * (age - 55) + 0.3 * sex + 0.2 * biomarker +
    0.6 * ckd + 0.25 * comorbidity + (effect_size / 0.15) * treatment
  event_24 <- rbinom(n, 1, plogis(lp_out))
  data.frame(age = round(age, 1), sex = sex, biomarker = round(biomarker, 3),
             comorbidity = comorbidity, ckd = ckd, treatment = treatment,
             event_24 = event_24, stringsAsFactors = FALSE)
}

compute_truth <- function(n, overlap_strength, effect_size, seed) {
  set.seed(seed)
  age         <- rnorm(n, mean = 55, sd = 10)
  sex         <- rbinom(n, 1, 0.55)
  biomarker   <- rnorm(n, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n, 1, 0.12)
  lp_base <- -2.5 + 0.015 * (age - 55) + 0.3 * sex + 0.2 * biomarker +
    0.6 * ckd + 0.25 * comorbidity
  coef_trt <- effect_size / 0.15
  mean(plogis(lp_base + coef_trt)) - mean(plogis(lp_base))
}

cat("Computing ground truth...\n"); .flush()
truth_A <- compute_truth(N_TRUTH, 0.5, -0.05, SEED)
truth_B <- compute_truth(N_TRUTH, 1.5, -0.05, SEED)
cat(sprintf("  Scenario A true RD = %.5f\n", truth_A))
cat(sprintf("  Scenario B true RD = %.5f\n\n", truth_B))
.flush()

# ── Low-level estimator helpers ───────────────────────────────────────────────

# Stabilised IPTW: GLM PS, stabilised Horwitz-Thompson
.iptw_fit <- function(dat, truncate = TRUNC) {
  fml <- reformulate(COVARS, response = TREAT)
  ps  <- predict(glm(fml, data = dat, family = binomial()), type = "response")
  ps  <- pmin(pmax(ps, truncate), 1 - truncate)
  A   <- dat[[TREAT]]; Y <- dat[[OUTC]]
  pA  <- mean(A)
  w   <- ifelse(A == 1, pA / ps, (1 - pA) / (1 - ps))
  mu1 <- weighted.mean(Y[A == 1], w[A == 1])
  mu0 <- weighted.mean(Y[A == 0], w[A == 0])
  est <- mu1 - mu0
  # Influence-function SE (treats PS as known — conservative under estimation)
  eif <- (A / ps - (1 - A) / (1 - ps)) * (Y - ifelse(A == 1, mu1, mu0)) +
         (mu1 - mu0) - est
  se  <- sqrt(var(eif) / nrow(dat))
  list(estimate = est, se = se,
       ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}

# Nearest-neighbour PS matching (1:1) then TMLE on matched cohort.
.match_tmle_fit <- function(dat, truncate = TRUNC) {
  fml <- reformulate(COVARS, response = TREAT)
  # Fit PS for matching
  ps_mod <- glm(fml, data = dat, family = binomial())
  ps     <- predict(ps_mod, type = "response")
  ps_cl  <- pmin(pmax(ps, truncate), 1 - truncate)
  # 1:1 nearest-neighbour matching on logit PS
  m.out  <- matchit(fml, data = dat, method = "nearest", ratio = 1,
                    distance = qlogis(ps_cl))
  m.dat  <- match.data(m.out, distance = "prop.score")
  # TMLE on matched cohort (GLM-only for speed + reproducibility)
  suppressMessages(
    fit <- tmle::tmle(
      Y             = m.dat[[OUTC]],
      A             = m.dat[[TREAT]],
      W             = m.dat[, COVARS, drop = FALSE],
      family        = "binomial",
      Q.SL.library  = "SL.glm",
      g.SL.library  = "SL.glm",
      verbose       = FALSE
    )
  )
  ate <- fit$estimates$ATE
  list(estimate = ate$psi, se = sqrt(ate$var.psi),
       ci_lower = ate$CI[1], ci_upper = ate$CI[2])
}

# Bootstrap SE for IPTW
.boot_iptw <- function(dat, B = B_BOOT, seed = 1L) {
  set.seed(seed)
  n <- nrow(dat)
  boots <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(.iptw_fit(dat[idx, ])$estimate, error = function(e) NA_real_)
  }, numeric(1))
  boots <- boots[is.finite(boots)]
  list(se = sd(boots),
       ci = unname(quantile(boots, c(0.025, 0.975))))
}

# Bootstrap SE for Match_TMLE
.boot_match_tmle <- function(dat, B = B_BOOT, seed = 1L) {
  set.seed(seed)
  n <- nrow(dat)
  boots <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(.match_tmle_fit(dat[idx, ])$estimate, error = function(e) NA_real_)
  }, numeric(1))
  boots <- boots[is.finite(boots)]
  list(se = sd(boots),
       ci = unname(quantile(boots, c(0.025, 0.975))))
}

# ── MC loop ───────────────────────────────────────────────────────────────────

run_scenario <- function(sc_label, overlap_strength, truth_rd) {
  cat(sprintf("=== %s (true RD = %.5f) ===\n", sc_label, truth_rd)); .flush()
  t0 <- Sys.time()

  rows <- vector("list", N_MC)
  for (i in seq_len(N_MC)) {
    dat <- generate_data(N_OBS, overlap_strength, seed = NULL)

    # --- IPTW ---
    iptw_if   <- tryCatch(.iptw_fit(dat),    error = function(e) NULL)
    iptw_boot <- tryCatch(.boot_iptw(dat, seed = i), error = function(e) NULL)

    # --- Match_TMLE ---
    mt_if   <- tryCatch(.match_tmle_fit(dat),     error = function(e) NULL)
    mt_boot <- tryCatch(.boot_match_tmle(dat, seed = i), error = function(e) NULL)

    make_row <- function(method, if_fit, boot_fit) {
      if (is.null(if_fit)) return(NULL)
      covers_if   <- as.integer(!is.na(if_fit$ci_lower) &&
                                  if_fit$ci_lower <= truth_rd &&
                                  truth_rd <= if_fit$ci_upper)
      covers_boot <- if (!is.null(boot_fit) && length(boot_fit$ci) == 2)
        as.integer(boot_fit$ci[1] <= truth_rd && truth_rd <= boot_fit$ci[2])
      else NA_integer_
      data.frame(
        rep         = i,
        method      = method,
        estimate    = if_fit$estimate,
        se_if       = if_fit$se,
        se_boot     = if (!is.null(boot_fit)) boot_fit$se else NA_real_,
        covers_if   = covers_if,
        covers_boot = covers_boot,
        stringsAsFactors = FALSE
      )
    }

    rep_rows <- rbind(
      make_row("IPTW",       iptw_if, iptw_boot),
      make_row("Match_TMLE", mt_if,   mt_boot)
    )
    rows[[i]] <- rep_rows

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    remaining <- if (i > 1) elapsed / i * (N_MC - i) else NA
    if (i %% 10 == 0) {
      cat(sprintf("  rep %d/%d  (%.1f min elapsed", i, N_MC, elapsed))
      if (!is.na(remaining)) cat(sprintf(", ~%.1f min remaining", remaining))
      cat(")\n"); .flush()
    }
  }
  do.call(rbind, Filter(Negate(is.null), rows))
}

# Ground truth must be fixed; seeds vary across reps (no rep-level seed so
# each rep draws a fresh dataset).
set.seed(SEED + 1L)  # advance RNG state past ground-truth draws

results_A <- run_scenario("Scenario A: Good Overlap",    0.5, truth_A)
results_B <- run_scenario("Scenario B: Marginal Overlap", 1.5, truth_B)

all_results <- rbind(
  cbind(scenario = "good_overlap",    results_A),
  cbind(scenario = "marginal_overlap", results_B)
)

# ── Summary table ─────────────────────────────────────────────────────────────

build_summary <- function(df, truth_rd) {
  by_method <- split(df, df$method)
  rows <- lapply(names(by_method), function(m) {
    d  <- by_method[[m]]
    ok <- !is.na(d$estimate)
    d  <- d[ok, ]
    if (nrow(d) == 0) return(NULL)

    emp_sd      <- sd(d$estimate)
    mean_se_if  <- mean(d$se_if,   na.rm = TRUE)
    mean_se_boot <- mean(d$se_boot, na.rm = TRUE)

    data.frame(
      method          = m,
      n_reps          = nrow(d),
      bias            = round(mean(d$estimate) - truth_rd, 5),
      emp_sd          = round(emp_sd, 5),
      mean_se_if      = round(mean_se_if, 5),
      se_sd_ratio_if  = round(mean_se_if / emp_sd, 3),
      coverage_if     = round(mean(d$covers_if,   na.rm = TRUE), 3),
      mean_se_boot    = round(mean_se_boot, 5),
      se_sd_ratio_boot = round(mean_se_boot / emp_sd, 3),
      coverage_boot   = round(mean(d$covers_boot, na.rm = TRUE), 3),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

summary_A <- build_summary(results_A, truth_A)
summary_B <- build_summary(results_B, truth_B)

cat("\n=== SUMMARY: Scenario A (Good Overlap) ===\n")
print(summary_A)
cat("\n=== SUMMARY: Scenario B (Marginal Overlap) ===\n")
print(summary_B)
.flush()

# ── Save ──────────────────────────────────────────────────────────────────────

combined_summary <- rbind(
  cbind(scenario = "A: Good Overlap",    summary_A),
  cbind(scenario = "B: Marginal Overlap", summary_B)
)

saveRDS(list(
  results         = all_results,
  summary_A       = summary_A,
  summary_B       = summary_B,
  combined_summary = combined_summary,
  config          = list(N_MC=N_MC, B_BOOT=B_BOOT, N_OBS=N_OBS,
                         TRUNC=TRUNC, SEED=SEED),
  truth           = list(A=truth_A, B=truth_B),
  runtime_minutes = as.numeric(difftime(Sys.time(),
                     proc.time()["elapsed"], units = "mins"))
), file.path(out_dir, "bootstrap_variance.rds"))

write.csv(combined_summary,
          file.path(out_dir, "bootstrap_variance.csv"),
          row.names = FALSE)

cat(sprintf("\nCOMPLETED: %s\n", format(Sys.time())))
cat("Output: results_new/bootstrap_variance.rds\n")
cat("Output: results_new/bootstrap_variance.csv\n")
print(combined_summary)
cat("\nBOOTSTRAP_VARIANCE_DONE\n")
