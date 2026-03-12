# scripts/sim_scenario1_bad_overlap_tmle_ok.R
# (Create this file exactly, per instructions. Full standalone script.)

suppressPackageStartupMessages({
  requireNamespace("tmle", quietly = TRUE) || stop("Package 'tmle' is required.")
  requireNamespace("SuperLearner", quietly = TRUE) || stop("Package 'SuperLearner' is required.")
})

dir.create("outputs", showWarnings = FALSE)
out_dir <- file.path("outputs", "sim_scenario1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

expit <- function(x) 1 / (1 + exp(-x))

# -----------------------------
# Stage 1: Prespec parameters
# -----------------------------
params <- list(
  master_seed = 20260304L,
  N = 1500L,
  reps = 50L,
  truth_N = 50000L,         # increase later (e.g., 200000)
  strength_overlap = 2.6,   # increase to worsen overlap (e.g., 3.0–4.0)
  eps_clip = 0.005,         # enforce practical positivity (no deterministic A)
  g_trunc = c(0.01, 0.99),  # truncation bounds for diagnostics/IPTW
  effect_A = log(1.15),     # log-OR scale for outcome model (modifiable)
  target_event_rate = 0.10, # rough marginal target; achieved via alpha_Q tuning
  overlap_flag_rules = list(
    frac_extreme_thresh = 0.30,
    extreme_cut = c(0.05, 0.95),
    ess_frac_thresh = 0.20
  ),
  tmle = list(
    Q_SL = c("SL.glm", "SL.gam", "SL.glmnet", "SL.mean"),
    g_SL = c("SL.glm", "SL.glmnet", "SL.mean")
  ),
  iptw_bootstrap_B = 100L   # keep small for speed; can increase later
)

saveRDS(params, file.path(out_dir, "params_used.rds"))

# -----------------------------
# DGP: single replicate
# -----------------------------
generate_data <- function(N, seed, strength_overlap, eps_clip, effect_A, target_event_rate) {
  set.seed(seed)

  W1 <- rnorm(N, 0, 1)
  W2 <- rbinom(N, 1, 0.5)
  W3 <- rnorm(N, 0, 1)

  # Treatment: near-positivity (practical) via steep lp_g
  lp_g_raw <- strength_overlap * (3.0 * W1 + 1.5 * W2 - 1.0 * W3)
  # Center so P(A=1) ~ 0.5 in-sample
  alpha_g <- -mean(lp_g_raw)
  lp_g <- alpha_g + lp_g_raw
  g <- expit(lp_g)
  g <- pmin(pmax(g, eps_clip), 1 - eps_clip)
  A <- rbinom(N, 1, g)

  # Outcome: learnable, with mild nonlinearity capturable by GAM/GLMnet
  # Choose alpha_Q to hit approximate marginal event rate.
  # Use a quick one-step calibration based on mean risk under observed A.
  lp_Q_no_alpha <- effect_A * A + 1.0 * W1 + 0.5 * W2 - 0.5 * W3 + 0.3 * (W1^2)
  # Solve alpha_Q so mean(expit(alpha_Q + lp_Q_no_alpha)) ~ target_event_rate
  f <- function(alpha) mean(expit(alpha + lp_Q_no_alpha)) - target_event_rate
  alpha_Q <- uniroot(f, lower = -10, upper = 10)$root

  lp_Q <- alpha_Q + lp_Q_no_alpha
  pY <- expit(lp_Q)
  Y <- rbinom(N, 1, pY)

  data.frame(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y, g_true = g, pY_true = pY)
}

# -----------------------------
# Ground truth RD (Monte Carlo or exact risk)
# -----------------------------
compute_truth <- function(truth_N, seed, strength_overlap, eps_clip, effect_A, target_event_rate) {
  set.seed(seed)
  W1 <- rnorm(truth_N, 0, 1)
  W2 <- rbinom(truth_N, 1, 0.5)
  W3 <- rnorm(truth_N, 0, 1)

  # Define outcome model (same as DGP)
  # Calibrate alpha_Q on A=0.5 “mixture” to keep consistent with generator.
  A_mix <- rbinom(truth_N, 1, 0.5)
  lp_Q_no_alpha_mix <- effect_A * A_mix + 1.0 * W1 + 0.5 * W2 - 0.5 * W3 + 0.3 * (W1^2)
  f <- function(alpha) mean(expit(alpha + lp_Q_no_alpha_mix)) - target_event_rate
  alpha_Q <- uniroot(f, lower = -10, upper = 10)$root

  lp1 <- alpha_Q + effect_A * 1 + 1.0 * W1 + 0.5 * W2 - 0.5 * W3 + 0.3 * (W1^2)
  lp0 <- alpha_Q + effect_A * 0 + 1.0 * W1 + 0.5 * W2 - 0.5 * W3 + 0.3 * (W1^2)

  p1 <- expit(lp1)
  p0 <- expit(lp0)

  RD_true <- mean(p1) - mean(p0)
  list(RD_true = RD_true, alpha_Q = alpha_Q, truth_N = truth_N)
}

truth <- compute_truth(
  truth_N = params$truth_N,
  seed = params$master_seed + 999,
  strength_overlap = params$strength_overlap,
  eps_clip = params$eps_clip,
  effect_A = params$effect_A,
  target_event_rate = params$target_event_rate
)

saveRDS(truth, file.path(out_dir, "truth.rds"))

# -----------------------------
# Stage 2 diagnostics (outcome-blind)
# -----------------------------
fit_ps_glm <- function(df) {
  # intentionally simple, fast PS model for diagnostics
  mod <- glm(A ~ W1 + W2 + W3, family = binomial(), data = df)
  p <- as.numeric(stats::predict(mod, type = "response"))
  p
}

compute_overlap_diagnostics <- function(df, g_hat, g_trunc, rules) {
  g_hat_t <- pmin(pmax(g_hat, g_trunc[1]), g_trunc[2])
  A <- df$A
  w <- A / g_hat_t + (1 - A) / (1 - g_hat_t)

  frac_extreme <- mean(g_hat < rules$extreme_cut[1] | g_hat > rules$extreme_cut[2])
  ess <- (sum(w)^2) / sum(w^2)
  ess_frac <- ess / nrow(df)
  max_w <- max(w)

  overlap_flag <- (frac_extreme > rules$frac_extreme_thresh) || (ess_frac < rules$ess_frac_thresh)

  list(
    frac_extreme = frac_extreme,
    min_g = min(g_hat),
    max_g = max(g_hat),
    ess = ess,
    ess_frac = ess_frac,
    max_w = max_w,
    overlap_flag = overlap_flag
  )
}

# -----------------------------
# Estimators
# -----------------------------
tmle_rd <- function(df, Q_SL, g_SL) {
  W <- df[, c("W1", "W2", "W3")]
  fit <- tmle::tmle(
    Y = df$Y,
    A = df$A,
    W = W,
    family = "binomial",
    Q.SL.library = Q_SL,
    g.SL.library = g_SL
  )
  # For binary Y, tmle returns estimates for EY1 and EY0; RD is psi (ATE)
  est <- as.numeric(fit$estimates$ATE$psi)
  se  <- as.numeric(fit$estimates$ATE$se)
  ci  <- c(est - 1.96 * se, est + 1.96 * se)
  list(est = est, se = se, lcl = ci[1], ucl = ci[2])
}

crude_rd <- function(df) {
  A <- df$A
  Y <- df$Y
  p1 <- mean(Y[A == 1])
  p0 <- mean(Y[A == 0])
  n1 <- sum(A == 1)
  n0 <- sum(A == 0)
  est <- p1 - p0
  se <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  ci <- c(est - 1.96 * se, est + 1.96 * se)
  list(est = est, se = se, lcl = ci[1], ucl = ci[2])
}

iptw_rd <- function(df, g_hat, g_trunc, B = 0L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  A <- df$A
  Y <- df$Y

  g_hat_t <- pmin(pmax(g_hat, g_trunc[1]), g_trunc[2])
  w <- A / g_hat_t + (1 - A) / (1 - g_hat_t)

  mu1 <- sum(w * A * Y) / sum(w * A)
  mu0 <- sum(w * (1 - A) * Y) / sum(w * (1 - A))
  est <- mu1 - mu0

  if (B <= 0) {
    return(list(est = est, se = NA_real_, lcl = NA_real_, ucl = NA_real_))
  }

  # Fast bootstrap for CI (percentile)
  n <- nrow(df)
  boot_est <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    d2 <- df[idx, , drop = FALSE]
    g2 <- g_hat[idx]
    g2t <- pmin(pmax(g2, g_trunc[1]), g_trunc[2])
    A2 <- d2$A; Y2 <- d2$Y
    w2 <- A2 / g2t + (1 - A2) / (1 - g2t)
    mu1b <- sum(w2 * A2 * Y2) / sum(w2 * A2)
    mu0b <- sum(w2 * (1 - A2) * Y2) / sum(w2 * (1 - A2))
    boot_est[b] <- mu1b - mu0b
  }
  se <- stats::sd(boot_est)
  ci <- stats::quantile(boot_est, probs = c(0.025, 0.975), names = FALSE, type = 6)
  list(est = est, se = se, lcl = ci[1], ucl = ci[2])
}

# -----------------------------
# Stage 3: run replicates
# -----------------------------
set.seed(params$master_seed)

rep_results <- vector("list", params$reps)

for (r in seq_len(params$reps)) {
  seed_r <- params$master_seed + r
  df <- generate_data(
    N = params$N,
    seed = seed_r,
    strength_overlap = params$strength_overlap,
    eps_clip = params$eps_clip,
    effect_A = params$effect_A,
    target_event_rate = params$target_event_rate
  )

  # Outcome-blind PS diagnostics: fit g_hat from A,W only
  g_hat <- fit_ps_glm(df)

  diag <- compute_overlap_diagnostics(
    df = df,
    g_hat = g_hat,
    g_trunc = params$g_trunc,
    rules = params$overlap_flag_rules
  )

  # Estimation
  tm <- tmle_rd(df, Q_SL = params$tmle$Q_SL, g_SL = params$tmle$g_SL)
  cr <- crude_rd(df)
  ip <- iptw_rd(df, g_hat = g_hat, g_trunc = params$g_trunc,
                B = params$iptw_bootstrap_B, seed = seed_r + 10000)

  rep_results[[r]] <- data.frame(
    rep = r,
    N = nrow(df),
    event_rate = mean(df$Y),
    treat_rate = mean(df$A),
    frac_extreme = diag$frac_extreme,
    min_g = diag$min_g,
    max_g = diag$max_g,
    ess_frac = diag$ess_frac,
    max_w = diag$max_w,
    overlap_flag = diag$overlap_flag,
    tmle_est = tm$est, tmle_se = tm$se, tmle_lcl = tm$lcl, tmle_ucl = tm$ucl,
    crude_est = cr$est, crude_se = cr$se, crude_lcl = cr$lcl, crude_ucl = cr$ucl,
    iptw_est = ip$est, iptw_se = ip$se, iptw_lcl = ip$lcl, iptw_ucl = ip$ucl
  )
}

rep_df <- do.call(rbind, rep_results)
saveRDS(rep_df, file.path(out_dir, "replicate_results.rds"))

# -----------------------------
# Summaries: bias + coverage (overall and conditional on overlap_flag bad)
# -----------------------------
RD_true <- truth$RD_true

coverage <- function(est, lcl, ucl, truth) mean(lcl <= truth & truth <= ucl, na.rm = TRUE)

summarize_one <- function(df, prefix) {
  est <- df[[paste0(prefix, "_est")]]
  se  <- df[[paste0(prefix, "_se")]]
  lcl <- df[[paste0(prefix, "_lcl")]]
  ucl <- df[[paste0(prefix, "_ucl")]]

  data.frame(
    estimator = prefix,
    mean_est = mean(est, na.rm = TRUE),
    bias = mean(est, na.rm = TRUE) - RD_true,
    rmse = sqrt(mean((est - RD_true)^2, na.rm = TRUE)),
    emp_sd = stats::sd(est, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    coverage = coverage(est, lcl, ucl, RD_true)
  )
}

overall <- rbind(
  summarize_one(rep_df, "tmle"),
  summarize_one(rep_df, "crude"),
  summarize_one(rep_df, "iptw")
)

bad_only <- subset(rep_df, overlap_flag)
bad <- rbind(
  summarize_one(bad_only, "tmle"),
  summarize_one(bad_only, "crude"),
  summarize_one(bad_only, "iptw")
)
bad$subset <- "overlap_flag_bad"
overall$subset <- "overall"

diag_summary <- data.frame(
  RD_true = RD_true,
  reps = params$reps,
  N = params$N,
  prop_overlap_bad = mean(rep_df$overlap_flag),
  mean_frac_extreme = mean(rep_df$frac_extreme),
  mean_ess_frac = mean(rep_df$ess_frac),
  mean_max_w = mean(rep_df$max_w)
)

summary_df <- rbind(overall, bad)
summary_out <- merge(summary_df, diag_summary, by = NULL)
utils::write.csv(summary_out, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)

# session info
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

# console print
cat("\nScenario 1: Bad overlap but TMLE OK (tune strength_overlap to get desired behavior)\n")
print(diag_summary)
cat("\nSummary metrics:\n")
print(summary_df)

cat("\nWrote outputs to: ", out_dir, "\n", sep = "")
