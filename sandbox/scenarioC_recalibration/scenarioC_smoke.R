#!/usr/bin/env Rscript
# ============================================================================
# Scenario C recalibration smoke test
# ----------------------------------------------------------------------------
# Goal: find an unmeasured-confounding parameterization in which DETECTION and
# FRAGILITY align, i.e. at a *plausible* true confounder strength s*:
#   (1) the realised TMLE (adjusting for measured W only) is materially biased
#       -> |bias_real| > max_abs_bias  => oracle STOP  (detection is real);
#   (2) the realised Wald CI under-covers the true marginal RD (< 0.95);
#   (3) the outcome-blind plasmode screen also flags it
#       -> |bias_screen| > max_abs_bias => screen STOP (fragility flags it).
# The current manuscript Scenario C sits at bias_real ~ 0.007 (< 0.02): the
# oracle says GO, so the workflow's STOP is a false positive. We sweep
# U_prev x s to locate a cell that clears the threshold at a defensible
# E-value-equivalent strength.
#
# Levers: s   = OR of U on BOTH A and Y (log(s)*U on each linear predictor)
#         U_prev = prevalence of the binary confounder
#         delta  = treatment effect on the outcome logit
# E-value (RR-scale) for a confounder OR: EV = OR + sqrt(OR*(OR-1)).
#   s=1.5 -> EV 2.37 ; s=2.0 -> EV 3.41 ; s=2.5 -> EV 4.44 ; s=3.0 -> EV 5.45
#
# Env overrides: SC_REPS (default 40), SC_CORES, SC_NL (0/1 linear/nonlinear).
# Output: sandbox/scenarioC_recalibration/scenarioC_smoke.rds + printed table.
# ============================================================================

suppressMessages({ library(parallel) })

cfg <- list(
  n            = 2000L,
  reps         = as.integer(Sys.getenv("SC_REPS", "40")),
  s_grid       = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5),
  U_prev_grid  = c(0.3, 0.5),
  delta_logit  = 0.6,
  max_abs_bias = 0.02,
  nl           = as.integer(Sys.getenv("SC_NL", "0")) == 1L,
  n_cores      = as.integer(Sys.getenv("SC_CORES",
                   as.character(max(1L, min(20L, parallel::detectCores() - 2L))))),
  out_dir      = "sandbox/scenarioC_recalibration",
  seed         = 2026L
)
if (!dir.exists(cfg$out_dir)) dir.create(cfg$out_dir, recursive = TRUE)

# ---- DGP (shared with run_gate_operating_characteristics.R) ----------------
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

# Realised estimator: TMLE with a GLM nuisance library adjusting for measured
# W only. Returns point estimate AND the Wald CI so we can score coverage.
tmle_W <- function(d) {
  fit <- tmle::tmle(Y = d$Y, A = d$A, W = d[, c("W1", "W2")], family = "binomial",
                    Q.SL.library = "SL.glm", g.SL.library = "SL.glm")
  a <- fit$estimates$ATE
  c(psi = a$psi, lo = a$CI[1], hi = a$CI[2])
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

# Outcome-blind screen: inject a what-if confounder of strength s into BOTH the
# synthetic treatment and outcome, score the analyst's TMLE vs the synthetic
# truth. Returns the screen bias.
screen_bias <- function(d, s, method, U_prev, delta) {
  nu <- fit_nuis(d, method)
  n  <- nrow(d); U <- rbinom(n, 1, U_prev)
  A_sim <- rbinom(n, 1, plogis(qlogis(nu$g0) + log(s) * U))
  Y_sim <- rbinom(n, 1, plogis(qlogis(nu$q0) + delta * A_sim + log(s) * U))
  synth_rd <- mean(plogis(qlogis(nu$q0) + delta + log(s) * U)) -
              mean(plogis(qlogis(nu$q0) + log(s) * U))
  unname(tmle_W(data.frame(Y = Y_sim, A = A_sim, W1 = d$W1, W2 = d$W2))["psi"]) - synth_rd
}

one_rep <- function(i, s, nl, trd, U_prev, delta, n) {
  d <- gen_real(n, s, nl, U_prev, delta)
  tryCatch({
    est <- tmle_W(d)
    c(bias_real      = unname(est["psi"]) - trd,
      cover_real     = as.numeric(est["lo"] <= trd & trd <= est["hi"]),
      bias_screen_sl = screen_bias(d, s, "sl", U_prev, delta))
  }, error = function(e) c(bias_real = NA, cover_real = NA, bias_screen_sl = NA))
}

true_rd <- function(s, nl, U_prev, delta, n = 2e5) {
  W1 <- rnorm(n); W2 <- rbinom(n, 1, 0.5); U <- rbinom(n, 1, U_prev)
  base <- q_logit(W1, W2, nl) + log(s) * U
  mean(plogis(base + delta)) - mean(plogis(base))
}

e_value <- function(or) round(or + sqrt(or * (or - 1)), 2)  # RR-scale E-value

cat(sprintf("Scenario C recalibration: n=%d reps=%d cores=%d surface=%s thr=%.3f\n\n",
            cfg$n, cfg$reps, cfg$n_cores, if (cfg$nl) "nonlinear" else "linear",
            cfg$max_abs_bias))

cl <- makeCluster(cfg$n_cores)
on.exit(stopCluster(cl), add = TRUE)
invisible(clusterEvalQ(cl, suppressMessages({
  library(tmle); library(SuperLearner); library(ranger); library(earth)
})))
clusterExport(cl, c("clamp", "g_logit", "q_logit", "gen_real", "tmle_W",
                    "fit_nuis", "screen_bias", "one_rep"))
clusterSetRNGStream(cl, cfg$seed)

rows <- list()
for (U_prev in cfg$U_prev_grid) {
  set.seed(cfg$seed)
  trd_cache <- sapply(cfg$s_grid, true_rd, nl = cfg$nl, U_prev = U_prev,
                      delta = cfg$delta_logit)
  for (k in seq_along(cfg$s_grid)) {
    s <- cfg$s_grid[k]; trd <- trd_cache[k]
    mat <- do.call(rbind, parLapply(cl, seq_len(cfg$reps), one_rep,
      s = s, nl = cfg$nl, trd = trd, U_prev = U_prev,
      delta = cfg$delta_logit, n = cfg$n))
    mat <- mat[stats::complete.cases(mat), , drop = FALSE]
    mn  <- colMeans(mat)
    rows[[length(rows) + 1]] <- data.frame(
      U_prev = U_prev, s = s, e_value = e_value(s), true_rd = round(trd, 4),
      reps_ok = nrow(mat),
      bias_real       = round(mn[["bias_real"]], 4),
      cover_real      = round(mn[["cover_real"]], 3),
      bias_screen_sl  = round(mn[["bias_screen_sl"]], 4),
      stop_oracle     = as.integer(abs(mn[["bias_real"]])      > cfg$max_abs_bias),
      stop_screen_sl  = as.integer(abs(mn[["bias_screen_sl"]]) > cfg$max_abs_bias),
      aligned = as.integer(abs(mn[["bias_real"]]) > cfg$max_abs_bias &
                           abs(mn[["bias_screen_sl"]]) > cfg$max_abs_bias),
      row.names = NULL, stringsAsFactors = FALSE)
    cat(sprintf("  U_prev=%.1f s=%.1f (EV %.2f)  bias_real=%+.4f cover=%.2f  screen_sl=%+.4f  oracle=%d screen=%d\n",
                U_prev, s, e_value(s), mn[["bias_real"]], mn[["cover_real"]],
                mn[["bias_screen_sl"]],
                as.integer(abs(mn[["bias_real"]]) > cfg$max_abs_bias),
                as.integer(abs(mn[["bias_screen_sl"]]) > cfg$max_abs_bias)))
  }
}
res <- do.call(rbind, rows)

cat("\n=== Scenario C recalibration grid ===\n"); print(res, row.names = FALSE)

aligned <- res[res$aligned == 1L, ]
cat("\n=== Cells where detection AND fragility align (oracle STOP & screen STOP) ===\n")
if (nrow(aligned)) {
  # Prefer the SMALLEST confounder strength that clears the bar: the most
  # defensible "true Scenario C" is the weakest confounder that still biases
  # the realised analysis past threshold.
  aligned <- aligned[order(aligned$U_prev, aligned$s), ]
  print(aligned, row.names = FALSE)
  best <- aligned[1, ]
  cat(sprintf("\nRECOMMENDED true Scenario C: U_prev=%.1f, s=%.1f (E-value %.2f)\n",
              best$U_prev, best$s, best$e_value))
  cat(sprintf("  realised bias %.4f (> %.2f, oracle STOP), coverage %.2f (< 0.95)\n",
              best$bias_real, cfg$max_abs_bias, best$cover_real))
  cat(sprintf("  screen bias  %.4f (> %.2f, fragility STOP) -> detection & fragility aligned\n",
              best$bias_screen_sl, cfg$max_abs_bias))
} else {
  cat("  none in this grid; widen s_grid or raise U_prev.\n")
}

saveRDS(list(cfg = cfg, results = res), file.path(cfg$out_dir, "scenarioC_smoke.rds"))
cat(sprintf("\nSaved %s/scenarioC_smoke.rds\n", cfg$out_dir))
