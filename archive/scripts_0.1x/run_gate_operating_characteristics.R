#!/usr/bin/env Rscript
# ============================================================================
# Operating characteristics of the outcome-blind DQ-stress gate
# ============================================================================
#
# CALIBRATION experiment. For a grid of true unmeasured-confounder strengths
# s, compare
#   B_real(s)   = bias of the realised TMLE (adjusting for measured W only)
#                 under a TRUE unmeasured confounder of strength s acting on
#                 BOTH treatment and outcome (the ORACLE quantity).
#   B_screen(s) = bias the outcome-blind plasmode screen PREDICTS: fit the
#                 nuisance surfaces g0(W), Q0(W) on the analyst's cohort,
#                 inject a what-if confounder of strength s into BOTH the
#                 synthetic treatment and synthetic outcome, and score the
#                 same TMLE against the synthetic truth.
# Synthetic nuisances are fit two ways: a misspecified GLM and a flexible
# SuperLearner (glm + earth + ranger). Under a NONLINEAR treatment/outcome
# surface the GLM-based screen under-predicts the realised bias (the central
# plasmode pitfall; Schreck et al. 2024); the flexible screen tracks it.
#
# Decision rule: STOP if |bias| > max_abs_bias. Screen sensitivity,
# specificity and false-STOP rate vs the oracle follow from the per-cell
# STOP verdicts. Parallelised over a PSOCK cluster.
#
# Env overrides for a quick check: GATE_REPS, GATE_CORES.
# Output: results/gate_oc.rds.
# ============================================================================

suppressMessages(library(parallel))

config <- list(
  n            = 2000L,
  reps         = as.integer(Sys.getenv("GATE_REPS", "200")),
  s_grid       = c(1.0, 1.5, 2.0, 2.5, 3.0),
  max_abs_bias = 0.02,
  U_prev       = 0.3,
  delta_logit  = 0.6,
  n_cores      = as.integer(Sys.getenv("GATE_CORES",
                   as.character(max(1L, min(20L, parallel::detectCores() - 2L))))),
  results_dir  = "results",
  seed         = 2026L
)
if (!dir.exists(config$results_dir)) dir.create(config$results_dir, recursive = TRUE)

# ---- functions (defined on master, exported to workers) --------------------
clamp <- function(p) pmin(pmax(p, 1e-3), 1 - 1e-3)

g_logit <- function(W1, W2, nl) {
  if (nl) -0.3 + 0.7 * W1 - 0.5 * W1^2 + 0.6 * W1 * W2 + 0.4 * W2
  else    -0.3 + 0.5 * W1 + 0.4 * W2
}

q_logit <- function(W1, W2, nl) {
  if (nl) -1 + 0.9 * W1 - 0.7 * W1^2 + 0.8 * W1 * W2 + 0.6 * W2
  else    -1 + 0.6 * W1 + 0.5 * W2
}

gen_real <- function(n, s, nl, U_prev, delta) {
  W1 <- rnorm(n); W2 <- rbinom(n, 1, 0.5); U <- rbinom(n, 1, U_prev)
  A  <- rbinom(n, 1, plogis(g_logit(W1, W2, nl) + log(s) * U))
  Y  <- rbinom(n, 1, plogis(q_logit(W1, W2, nl) + delta * A + log(s) * U))
  data.frame(Y = Y, A = A, W1 = W1, W2 = W2)
}

# Realised estimator the analyst runs: TMLE with a (limited) GLM nuisance
# library, adjusting for the measured covariates only.
tmle_W <- function(d) {
  fit <- tmle::tmle(Y = d$Y, A = d$A, W = d[, c("W1", "W2")], family = "binomial",
                    Q.SL.library = "SL.glm", g.SL.library = "SL.glm")
  fit$estimates$ATE$psi
}

fit_nuis <- function(d, method) {
  if (method == "glm") {
    q0 <- predict(glm(Y ~ W1 + W2, d, family = binomial()), type = "response")
    g0 <- predict(glm(A ~ W1 + W2, d, family = binomial()), type = "response")
  } else {
    lib <- c("SL.glm", "SL.earth", "SL.ranger"); cv <- list(V = 5L)
    q0 <- SuperLearner::SuperLearner(Y = d$Y, X = d[, c("W1", "W2")],
            family = binomial(), SL.library = lib, cvControl = cv)$SL.predict[, 1]
    g0 <- SuperLearner::SuperLearner(Y = d$A, X = d[, c("W1", "W2")],
            family = binomial(), SL.library = lib, cvControl = cv)$SL.predict[, 1]
  }
  list(q0 = clamp(q0), g0 = clamp(g0))
}

# The plasmode DQ screen: inject a what-if confounder of strength s into BOTH
# the synthetic treatment and outcome, then score the analyst's TMLE.
screen_bias <- function(d, s, method, U_prev, delta) {
  nu <- fit_nuis(d, method)
  n  <- nrow(d); U <- rbinom(n, 1, U_prev)
  A_sim <- rbinom(n, 1, plogis(qlogis(nu$g0) + log(s) * U))
  Y_sim <- rbinom(n, 1, plogis(qlogis(nu$q0) + delta * A_sim + log(s) * U))
  synth_rd <- mean(plogis(qlogis(nu$q0) + delta + log(s) * U)) -
              mean(plogis(qlogis(nu$q0) + log(s) * U))
  tmle_W(data.frame(Y = Y_sim, A = A_sim, W1 = d$W1, W2 = d$W2)) - synth_rd
}

one_rep <- function(i, s, nl, trd, U_prev, delta, n) {
  d <- gen_real(n, s, nl, U_prev, delta)
  tryCatch(c(
    bias_real       = tmle_W(d) - trd,
    bias_screen_glm = screen_bias(d, s, "glm", U_prev, delta),
    bias_screen_sl  = screen_bias(d, s, "sl",  U_prev, delta)),
    error = function(e) c(bias_real = NA, bias_screen_glm = NA, bias_screen_sl = NA))
}

# Large-sample true marginal RD under the full (W,U) law.
true_rd <- function(s, nl, U_prev, delta, n = 2e5) {
  W1 <- rnorm(n); W2 <- rbinom(n, 1, 0.5); U <- rbinom(n, 1, U_prev)
  base <- q_logit(W1, W2, nl) + log(s) * U
  mean(plogis(base + delta)) - mean(plogis(base))
}

cat(sprintf("Gate operating characteristics: n=%d reps=%d cores=%d thr=%.3f\n\n",
            config$n, config$reps, config$n_cores, config$max_abs_bias))

cl <- makeCluster(config$n_cores)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, suppressMessages({
  library(tmle); library(SuperLearner); library(ranger); library(earth)
})))
clusterExport(cl, c("clamp", "g_logit", "q_logit", "gen_real", "tmle_W",
                    "fit_nuis", "screen_bias", "one_rep"))
clusterSetRNGStream(cl, config$seed)

rows <- list()
for (nl in c(FALSE, TRUE)) {
  surf <- if (nl) "nonlinear" else "linear"
  set.seed(config$seed)
  trd_cache <- sapply(config$s_grid, true_rd, nl = nl,
                      U_prev = config$U_prev, delta = config$delta_logit)
  for (k in seq_along(config$s_grid)) {
    s <- config$s_grid[k]; trd <- trd_cache[k]
    mat <- do.call(rbind, parLapply(cl, seq_len(config$reps), one_rep,
      s = s, nl = nl, trd = trd, U_prev = config$U_prev,
      delta = config$delta_logit, n = config$n))
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    mn <- colMeans(mat); se <- apply(mat, 2, sd) / sqrt(nrow(mat))
    rows[[length(rows) + 1]] <- data.frame(
      surface = surf, true_s = s, true_rd = round(trd, 4), reps_ok = nrow(mat),
      bias_real       = round(mn[["bias_real"]], 4),
      bias_screen_glm = round(mn[["bias_screen_glm"]], 4),
      bias_screen_sl  = round(mn[["bias_screen_sl"]], 4),
      se_real         = round(se[["bias_real"]], 4),
      se_screen_glm   = round(se[["bias_screen_glm"]], 4),
      se_screen_sl    = round(se[["bias_screen_sl"]], 4),
      stop_oracle     = as.integer(abs(mn[["bias_real"]])       > config$max_abs_bias),
      stop_screen_glm = as.integer(abs(mn[["bias_screen_glm"]]) > config$max_abs_bias),
      stop_screen_sl  = as.integer(abs(mn[["bias_screen_sl"]])  > config$max_abs_bias),
      row.names = NULL, stringsAsFactors = FALSE)
    cat(sprintf("  [%s] s=%.1f  real=%+.4f  glm=%+.4f  sl=%+.4f\n",
                surf, s, mn[["bias_real"]], mn[["bias_screen_glm"]], mn[["bias_screen_sl"]]))
  }
}
res <- do.call(rbind, rows)

oc <- function(scr, ora) c(
  sensitivity = if (sum(ora == 1) > 0) mean(scr[ora == 1] == 1) else NA,
  specificity = if (sum(ora == 0) > 0) mean(scr[ora == 0] == 0) else NA,
  false_stop  = if (sum(ora == 0) > 0) mean(scr[ora == 0] == 1) else NA)

cat("\n=== Results ===\n"); print(res, row.names = FALSE)
cat("\nScreen vs oracle (GLM nuisance):\n");    print(round(oc(res$stop_screen_glm, res$stop_oracle), 3))
cat("Screen vs oracle (flexible nuisance):\n"); print(round(oc(res$stop_screen_sl,  res$stop_oracle), 3))

saveRDS(list(config = config, results = res), file.path(config$results_dir, "gate_oc.rds"))
cat("\nSaved results/gate_oc.rds\n")
