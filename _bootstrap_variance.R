#!/usr/bin/env Rscript
# ============================================================================
# Bootstrap vs influence-function variance comparison
#
# Demonstrates that the IF-based variance is badly mis-calibrated for
# Match_TMLE (se_sd_ratio ≈ 2.1 in Scenario B) while the nonparametric
# bootstrap is well-calibrated (se_sd_ratio ≈ 1.0). IPTW is included as
# a positive control — its IF-SE is already near-nominal on this DGP.
#
# Design:
#   Scenarios: A (good overlap) and B (marginal overlap)
#   Estimators: IPTW (stabilised), Match_TMLE
#   MC reps:    100 per scenario
#   Bootstrap:  B = 200 draws, GLM-only, percentile CI
#
# Implementation: self-contained GLM + manual TMLE (no tmle:: package calls
# in the bootstrap loop; avoids subprocess issues on Windows). The manual
# TMLE targeting step closely replicates the tmle-package influence-function
# variance, allowing a fair IF vs bootstrap comparison.
#
# Output: results_new/bootstrap_variance.rds
#         results_new/bootstrap_variance.csv
# ============================================================================

suppressMessages(library(MatchIt))

.flush <- function() if (!interactive()) flush(stdout())
cat("Bootstrap variance comparison\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))
.flush()

# ── Configuration ─────────────────────────────────────────────────────────────

N_MC    <- 100L
B_BOOT  <- 200L
N_OBS   <- 2000L
N_TRUTH <- 500000L
TRUNC   <- 0.05
SEED    <- 2026L
COVARS  <- c("age", "sex", "biomarker", "comorbidity", "ckd")
TREAT   <- "treatment"
OUTC    <- "event_24"

out_dir <- file.path(getwd(), "results_new")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── DGP ───────────────────────────────────────────────────────────────────────

generate_data <- function(n, overlap_strength = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  age         <- rnorm(n, 55, 10)
  sex         <- rbinom(n, 1, 0.55)
  biomarker   <- rnorm(n)
  comorbidity <- sample(0:2, n, TRUE, c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n, 1, 0.12)
  lp_trt <- -0.5 + overlap_strength * (0.03*(age-55) + 0.8*sex +
               0.6*biomarker + 0.5*ckd + 0.3*comorbidity)
  treatment  <- rbinom(n, 1, plogis(lp_trt))
  lp_out <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*biomarker +
    0.6*ckd + 0.25*comorbidity + (-0.05/0.15)*treatment
  event_24 <- rbinom(n, 1, plogis(lp_out))
  data.frame(age=round(age,1), sex=sex, biomarker=round(biomarker,3),
             comorbidity=comorbidity, ckd=ckd, treatment=treatment,
             event_24=event_24, stringsAsFactors=FALSE)
}

compute_truth <- function(n, overlap_strength) {
  set.seed(SEED)
  age <- rnorm(n,55,10); sex <- rbinom(n,1,.55)
  bio <- rnorm(n); com <- sample(0:2,n,T,c(.5,.3,.2)); ckd <- rbinom(n,1,.12)
  lp  <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*bio + 0.6*ckd + 0.25*com
  mean(plogis(lp + (-0.05/0.15))) - mean(plogis(lp))
}

cat("Computing ground truth...\n"); .flush()
truth_A <- compute_truth(N_TRUTH, 0.5)
truth_B <- compute_truth(N_TRUTH, 1.5)
cat(sprintf("  Scenario A true RD = %.5f\n", truth_A))
cat(sprintf("  Scenario B true RD = %.5f\n\n", truth_B))
.flush()

# ── Lightweight TMLE helpers (no tmle:: package) ───────────────────────────────

# Manual TMLE: Q-fit, clever covariate, targeting step, ATE + IF SE.
# Returns list(estimate, se, ci_lower, ci_upper).
.tmle_manual <- function(Y, A, W, truncate = TRUNC) {
  Wdf  <- as.data.frame(W)
  n    <- length(Y)
  covs <- names(Wdf)

  # g-model (treatment mechanism)
  g_df  <- cbind(A = A, Wdf)
  g_mod <- glm(A ~ ., data = g_df, family = binomial())
  g     <- pmin(pmax(predict(g_mod, type = "response"), truncate), 1 - truncate)

  # Q-model (outcome mechanism)
  q_df  <- cbind(Y = Y, A = A, Wdf)
  q_mod <- glm(Y ~ ., data = q_df, family = binomial())
  Q1    <- predict(q_mod, newdata = cbind(A = 1, Wdf), type = "response")
  Q0    <- predict(q_mod, newdata = cbind(A = 0, Wdf), type = "response")
  Q_AW  <- A * Q1 + (1 - A) * Q0   # initial prediction at observed A

  # Clever covariate and targeting step
  H  <- A / g - (1 - A) / (1 - g)
  Q_AW_cl <- pmin(pmax(Q_AW, 1e-6), 1 - 1e-6)
  eps <- tryCatch(
    coef(glm(Y ~ H + offset(qlogis(Q_AW_cl)), family = binomial()))["H"],
    error = function(e) 0
  )

  # Targeted predictions
  Q1_star <- plogis(qlogis(pmin(pmax(Q1, 1e-6), 1-1e-6)) + eps / g)
  Q0_star <- plogis(qlogis(pmin(pmax(Q0, 1e-6), 1-1e-6)) - eps / (1 - g))
  ate     <- mean(Q1_star) - mean(Q0_star)

  # Influence-function SE
  eif <- (A / g - (1 - A) / (1 - g)) * (Y - Q_AW_cl) +
         (Q1_star - Q0_star) - ate
  se  <- sqrt(var(eif) / n)

  list(estimate = ate, se = se,
       ci_lower = ate - 1.96 * se,
       ci_upper = ate + 1.96 * se)
}

# ── Estimator fits ────────────────────────────────────────────────────────────

# IPTW: stabilised HT, IF-SE.
.iptw_fit <- function(dat) {
  A <- dat[[TREAT]]; Y <- dat[[OUTC]]
  W <- dat[, COVARS, drop = FALSE]
  g_mod <- glm(reformulate(COVARS, TREAT), data = dat, family = binomial())
  ps <- pmin(pmax(predict(g_mod, type = "response"), TRUNC), 1 - TRUNC)
  pA <- mean(A)
  w  <- ifelse(A == 1, pA / ps, (1 - pA) / (1 - ps))
  mu1 <- weighted.mean(Y[A == 1], w[A == 1])
  mu0 <- weighted.mean(Y[A == 0], w[A == 0])
  est <- mu1 - mu0
  eif <- (A / ps - (1 - A) / (1 - ps)) * (Y - ifelse(A == 1, mu1, mu0)) +
         (mu1 - mu0) - est
  se  <- sqrt(var(eif) / nrow(dat))
  list(estimate = est, se = se,
       ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}

# Match + TMLE: 1:1 NN matching on logit-PS, then manual TMLE on matched data.
.match_tmle_fit <- function(dat) {
  fml <- reformulate(COVARS, TREAT)
  g_mod <- glm(fml, data = dat, family = binomial())
  ps    <- pmin(pmax(predict(g_mod, type = "response"), TRUNC), 1 - TRUNC)
  suppressWarnings(
    m.out <- matchit(fml, data = dat, method = "nearest", ratio = 1,
                     distance = qlogis(ps))
  )
  m.dat <- match.data(m.out)
  .tmle_manual(m.dat[[OUTC]], m.dat[[TREAT]],
               m.dat[, COVARS, drop = FALSE])
}

# Bootstrap for either estimator (estimator_fn takes a data frame, returns list).
.bootstrap_se <- function(dat, estimator_fn, B = B_BOOT, seed = 1L) {
  set.seed(seed)
  n <- nrow(dat)
  boots <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(estimator_fn(dat[idx, , drop = FALSE])$estimate,
             error = function(e) NA_real_)
  }, numeric(1))
  boots <- boots[is.finite(boots)]
  if (length(boots) < 10L) return(list(se = NA_real_, ci = c(NA, NA)))
  list(se = sd(boots),
       ci = unname(quantile(boots, c(0.025, 0.975))))
}

# ── MC loop ───────────────────────────────────────────────────────────────────

run_scenario <- function(sc_label, overlap_strength, truth_rd) {
  cat(sprintf("=== %s (true RD = %.5f) ===\n", sc_label, truth_rd)); .flush()
  t0 <- Sys.time()
  rows <- vector("list", N_MC)

  for (i in seq_len(N_MC)) {
    dat <- generate_data(N_OBS, overlap_strength)

    iptw_if   <- tryCatch(.iptw_fit(dat),       error = function(e) NULL)
    iptw_boot <- tryCatch(.bootstrap_se(dat, .iptw_fit, seed = i),
                          error = function(e) NULL)

    mt_if   <- tryCatch(.match_tmle_fit(dat),   error = function(e) NULL)
    mt_boot <- tryCatch(.bootstrap_se(dat, .match_tmle_fit, seed = i),
                        error = function(e) NULL)

    make_row <- function(method, if_fit, boot_fit) {
      if (is.null(if_fit)) return(NULL)
      covers_if   <- as.integer(!is.na(if_fit$ci_lower) &&
                                  if_fit$ci_lower <= truth_rd &&
                                  truth_rd <= if_fit$ci_upper)
      covers_boot <- if (!is.null(boot_fit) && all(is.finite(boot_fit$ci)))
        as.integer(boot_fit$ci[1] <= truth_rd && truth_rd <= boot_fit$ci[2])
      else NA_integer_
      data.frame(rep = i, method = method,
                 estimate = if_fit$estimate,
                 se_if = if_fit$se,
                 se_boot = if (!is.null(boot_fit)) boot_fit$se else NA_real_,
                 covers_if = covers_if, covers_boot = covers_boot,
                 stringsAsFactors = FALSE)
    }

    rows[[i]] <- rbind(make_row("IPTW",       iptw_if, iptw_boot),
                       make_row("Match_TMLE", mt_if,   mt_boot))

    if (i %% 10 == 0) {
      el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      rem <- if (i > 1) el / i * (N_MC - i) else NA
      cat(sprintf("  rep %d/%d  (%.1f min", i, N_MC, el))
      if (!is.na(rem)) cat(sprintf(", ~%.0f min remaining", rem))
      cat(")\n"); .flush()
    }
  }
  do.call(rbind, Filter(Negate(is.null), rows))
}

set.seed(SEED + 99L)

results_A <- run_scenario("Scenario A: Good Overlap",     0.5, truth_A)
results_B <- run_scenario("Scenario B: Marginal Overlap",  1.5, truth_B)

# ── Summary ───────────────────────────────────────────────────────────────────

build_summary <- function(df, truth_rd) {
  do.call(rbind, lapply(unique(df$method), function(m) {
    d  <- df[df$method == m & !is.na(df$estimate), ]
    if (nrow(d) == 0) return(NULL)
    emp_sd    <- sd(d$estimate)
    mse_if    <- mean(d$se_if,   na.rm = TRUE)
    mse_boot  <- mean(d$se_boot, na.rm = TRUE)
    data.frame(
      method           = m,
      n_reps           = nrow(d),
      bias             = round(mean(d$estimate) - truth_rd, 5),
      emp_sd           = round(emp_sd,   5),
      mean_se_if       = round(mse_if,   5),
      se_sd_ratio_if   = round(mse_if  / emp_sd, 3),
      coverage_if      = round(mean(d$covers_if,   na.rm = TRUE), 3),
      mean_se_boot     = round(mse_boot, 5),
      se_sd_ratio_boot = round(mse_boot / emp_sd, 3),
      coverage_boot    = round(mean(d$covers_boot, na.rm = TRUE), 3),
      stringsAsFactors = FALSE
    )
  }))
}

summary_A <- build_summary(results_A, truth_A)
summary_B <- build_summary(results_B, truth_B)

cat("\n=== SUMMARY: Scenario A ===\n"); print(summary_A)
cat("\n=== SUMMARY: Scenario B ===\n"); print(summary_B)
.flush()

combined <- rbind(cbind(scenario = "A: Good Overlap",    summary_A),
                  cbind(scenario = "B: Marginal Overlap", summary_B))

saveRDS(list(results_A = results_A, results_B = results_B,
             summary_A = summary_A, summary_B = summary_B,
             combined  = combined,
             config    = list(N_MC=N_MC, B_BOOT=B_BOOT, N_OBS=N_OBS,
                              TRUNC=TRUNC, SEED=SEED),
             truth     = list(A=truth_A, B=truth_B)),
        file.path(out_dir, "bootstrap_variance.rds"))

write.csv(combined, file.path(out_dir, "bootstrap_variance.csv"),
          row.names = FALSE)

cat(sprintf("\nCOMPLETED: %s\n", format(Sys.time())))
cat("Output: results_new/bootstrap_variance.rds\n")
cat("Output: results_new/bootstrap_variance.csv\n")
cat("\nBOOTSTRAP_VARIANCE_DONE\n")
