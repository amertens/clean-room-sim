#!/usr/bin/env Rscript
# ============================================================================
# Plasmode synthetic-outcome fidelity diagnostic for the rescueCo case study.
#
# The DQ-stress plasmode reuses the empirical (W, A) structure and generates
# synthetic outcomes from a covariate-only model Q0(W) = E[Y | W]. This script
# quantifies how faithfully that generator reproduces the REAL outcome model,
# using CROSS-VALIDATED (out-of-sample) metrics so the numbers are honest
# rather than optimistic resubstitution:
#   (1) discrimination of Q0(W) on held-out data (10-fold cross-validated AUC);
#   (2) calibration of Q0(W) on held-out data (max abs predicted-minus-observed
#       deviation across risk deciles of the pooled out-of-fold predictions).
# The in-sample AUC is reported too, only to show the resubstitution optimism.
#
# Notes on what this does and does NOT establish:
#  - The marginal outcome prevalence is reproduced EXACTLY by construction (a
#    logistic GLM with an intercept forces mean(fitted) = mean(Y)), so a
#    prevalence match is a consistency check, not fidelity evidence, and is not
#    reported as such.
#  - Q0 mirrors the package generator (glm, predictions clamped to
#    [1e-3, 1-1e-3] as in R/plasmode_dq.R). It still cannot verify the synthetic
#    joint (A, Y) | W: the synthetic outcome carries a fixed risk difference and
#    no effect modification, the irreducible residual risk.
#
# Output: rescueCo/results/plasmode_fidelity.csv
# ============================================================================

rc  <- "rescueCo/results"
out <- file.path(rc, "plasmode_fidelity.csv")
lk  <- readRDS(file.path(rc, "stage1_lock_unmasked.rds"))

d  <- lk$data
W  <- lk$covariates
cc <- stats::complete.cases(d[, c(lk$outcome, lk$treatment, W)])
dd <- d[cc, , drop = FALSE]
Yc <- dd[[lk$outcome]]
n  <- nrow(dd)

clamp <- function(p) pmin(pmax(p, 1e-3), 1 - 1e-3)          # mirror plasmode_dq.R
fml   <- stats::reformulate(W, response = lk$outcome)        # backtick-safe

# Mann-Whitney AUC; NA (not NaN) when a class is empty.
auc <- function(y, s) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(s)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Max |predicted - observed| across risk bins, robust to tied quantile breaks.
calib_dev <- function(p, y, nbin = 10L) {
  br <- unique(stats::quantile(p, seq(0, 1, length.out = nbin + 1L)))
  if (length(br) < 3L) return(NA_real_)          # too few distinct predictions
  b  <- cut(p, breaks = br, include.lowest = TRUE, labels = FALSE)
  devs <- vapply(sort(unique(b)), function(k)
    mean(p[b == k]) - mean(y[b == k]), numeric(1))
  max(abs(devs))
}

# In-sample Q0 (resubstitution) for reference only.
q0_in <- clamp(stats::predict(
  stats::glm(fml, data = dd, family = stats::binomial()), type = "response"))

# 10-fold cross-validated out-of-fold predictions.
set.seed(2026)
K    <- 10L
fold <- sample(rep_len(seq_len(K), n))
p_cv <- numeric(n)
for (k in seq_len(K)) {
  tr <- fold != k; te <- fold == k
  fit <- stats::glm(fml, data = dd[tr, , drop = FALSE], family = stats::binomial())
  p_cv[te] <- clamp(stats::predict(fit, newdata = dd[te, , drop = FALSE],
                                   type = "response"))
}

res <- data.frame(
  metric = c("n_complete_case", "real_prevalence", "Q0_auc_insample",
             "Q0_auc_cv", "calibration_max_abs_dev_cv"),
  value  = c(n, round(mean(Yc), 4), round(auc(Yc, q0_in), 3),
             round(auc(Yc, p_cv), 3), round(calib_dev(p_cv, Yc), 3)),
  stringsAsFactors = FALSE)
write.csv(res, out, row.names = FALSE)
cat("Wrote", out, "\n"); print(res)
