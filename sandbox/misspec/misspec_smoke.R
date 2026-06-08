# ============================================================
# Reviewer item #1 (SMOKE): misspecified DGP
# ------------------------------------------------------------
# Shows that under a nonlinear + effect-modified data-generating process,
# main-effects-GLM estimators are biased / undercover for the marginal RD,
# while a flexible SuperLearner TMLE (gam + interactions) recovers it.
# This is the scenario the current linear-logistic DGP cannot exhibit.
#
# Usage: Rscript sandbox/misspec/misspec_smoke.R [reps] [n]
#   defaults: reps = 25, n = 1500
# ============================================================
suppressWarnings(suppressMessages({ library(SuperLearner); library(tmle) }))

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 25L
n    <- if (length(args) >= 2) as.integer(args[[2]]) else 1500L
tau  <- 0.5          # treatment coef on logit (base arm); effect modified by sex
set.seed(20260607)

# --- DGP: same-sign nonlinear confounding so the crude contrast is clearly
#          biased. Nonlinearity is a quadratic in age and a sex x biomarker
#          interaction, present in BOTH the PS and the outcome (so a main-
#          effects GLM is misspecified for both nuisances and double
#          robustness cannot rescue it). a^2 kept moderate in the PS to
#          avoid manufacturing a positivity violation. ------------------
# PS nonlinearity is a symmetric sex x biomarker interaction (a main-effects
# GLM cannot capture it, so double robustness still fails), but it does not
# pile the propensity score at 0/1 the way an always-positive a^2 would.
.lp_g <- function(a, sex, bio, ckd, com)
  -0.1 + 0.5 * a + 1.0 * sex * bio + 0.5 * ckd + 0.2 * com
.lp_y <- function(A, a, sex, bio, ckd, com)
  -0.6 + 0.4 * a + 0.7 * a^2 + 1.0 * sex * bio + 0.6 * ckd + 0.3 * com +
   tau * A * (1 + 0.4 * sex)

gen <- function(n) {
  age <- rnorm(n, 55, 10); a <- (age - 55) / 10
  sex <- rbinom(n, 1, 0.5); bio <- rnorm(n, 0, 1)
  ckd <- rbinom(n, 1, 0.2); com <- sample(0:2, n, TRUE, prob = c(0.5, 0.3, 0.2))
  ps  <- plogis(.lp_g(a, sex, bio, ckd, com))
  A   <- rbinom(n, 1, ps)
  Y   <- rbinom(n, 1, plogis(.lp_y(A, a, sex, bio, ckd, com)))
  list(W = data.frame(age = age, sex = sex, bio = bio, ckd = ckd, com = com),
       A = A, Y = Y, ps = ps)
}

# --- True marginal RD (large sample) + PS-range diagnostic -----------
big <- 3e5
age <- rnorm(big, 55, 10); a <- (age - 55) / 10
sex <- rbinom(big, 1, 0.5); bio <- rnorm(big, 0, 1)
ckd <- rbinom(big, 1, 0.2); com <- sample(0:2, big, TRUE, prob = c(0.5, 0.3, 0.2))
truth <- mean(plogis(.lp_y(1, a, sex, bio, ckd, com))) -
         mean(plogis(.lp_y(0, a, sex, bio, ckd, com)))
psd   <- plogis(.lp_g(a, sex, bio, ckd, com))
cat(sprintf("True marginal RD = %.5f\n", truth))
cat(sprintf("PS range [%.4f, %.4f]  frac<.05=%.1f%%  frac>.95=%.1f%%\n",
            min(psd), max(psd), 100 * mean(psd < .05), 100 * mean(psd > .95)))
rm(age, a, sex, bio, ckd, com, psd)

lib_glm <- "SL.glm"
lib_sl  <- c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.mean")

res <- vector("list", reps)
for (i in seq_len(reps)) {
  d <- gen(n)
  # Crude
  crude <- mean(d$Y[d$A == 1]) - mean(d$Y[d$A == 0])
  # GLM-only TMLE (misspecified nuisances)
  fg <- tryCatch(tmle(d$Y, d$A, d$W, family = "binomial",
                      Q.SL.library = lib_glm, g.SL.library = lib_glm),
                 error = function(e) NULL)
  # Flexible SL TMLE
  fs <- tryCatch(tmle(d$Y, d$A, d$W, family = "binomial",
                      Q.SL.library = lib_sl, g.SL.library = lib_sl),
                 error = function(e) NULL)
  grab <- function(f) if (is.null(f)) c(NA, NA, NA) else
    c(f$estimates$ATE$psi, f$estimates$ATE$CI[1], f$estimates$ATE$CI[2])
  res[[i]] <- data.frame(
    rep = i, crude = crude,
    glm = grab(fg)[1], glm_lo = grab(fg)[2], glm_hi = grab(fg)[3],
    sl  = grab(fs)[1], sl_lo  = grab(fs)[2], sl_hi  = grab(fs)[3])
  cat(sprintf("  rep %2d/%d  crude=%.4f  glmTMLE=%.4f  slTMLE=%.4f\n",
              i, reps, crude, res[[i]]$glm, res[[i]]$sl)); flush(stdout())
}
R <- do.call(rbind, res)

summ <- function(est, lo, hi, label) {
  ok <- is.finite(est)
  data.frame(estimator = label, n_ok = sum(ok),
    bias = mean(est[ok]) - truth,
    emp_sd = sd(est[ok]),
    coverage = mean(lo[ok] <= truth & truth <= hi[ok]),
    ci_width = mean(hi[ok] - lo[ok]))
}
tab <- rbind(
  data.frame(estimator = "Crude", n_ok = sum(is.finite(R$crude)),
             bias = mean(R$crude) - truth, emp_sd = sd(R$crude),
             coverage = NA, ci_width = NA),
  summ(R$glm, R$glm_lo, R$glm_hi, "TMLE (GLM nuisances)"),
  summ(R$sl,  R$sl_lo,  R$sl_hi,  "TMLE (SL: gam+interaction)"))
cat("\n=== Misspecified-DGP smoke summary (truth =", round(truth, 5), ") ===\n")
print(tab, row.names = FALSE, digits = 4)
cat("\nMISSPEC_SMOKE_DONE\n")
