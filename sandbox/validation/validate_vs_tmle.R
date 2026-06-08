# ============================================================
# Tier-2 validation: cleanTMLE estimators vs the reference `tmle` package
# and the analytic truth, on identical data.
# ------------------------------------------------------------
# Produces a small table showing that cleanTMLE's point estimate and SE agree
# with tmle::tmle (same nuisance libraries) and recover the analytic marginal
# risk difference. Establishes that the package's estimator machinery is a
# correct implementation, separate from its governance/DQ contributions.
#
# Usage: Rscript sandbox/validation/validate_vs_tmle.R [n]   (default n = 20000)
# Output: sandbox/validation/validation_vs_tmle.csv
# ============================================================
suppressWarnings(suppressMessages({ library(SuperLearner); library(tmle); library(cleanTMLE) }))

args <- commandArgs(trailingOnly = TRUE)
n    <- if (length(args) >= 1) as.integer(args[[1]]) else 20000L
set.seed(2026)
COVARS <- c("age", "sex", "biomarker", "comorbidity", "ckd")

# Linear-logistic DGP (matches run_simulation.R Scenario A), analytic truth.
gen <- function(n) {
  age<-rnorm(n,55,10); sex<-rbinom(n,1,.55); biomarker<-rnorm(n)
  comorbidity<-sample(0:2,n,T,c(.5,.3,.2)); ckd<-rbinom(n,1,.12)
  lp_trt<- -0.5 + 0.5*(0.03*(age-55)+0.8*sex+0.6*biomarker+0.5*ckd+0.3*comorbidity)
  treatment<-rbinom(n,1,plogis(lp_trt))
  lp_out<- -2.5 + 0.015*(age-55)+0.3*sex+0.2*biomarker+0.6*ckd+0.25*comorbidity+(-0.05/0.15)*treatment
  data.frame(age=age,sex=sex,biomarker=biomarker,comorbidity=comorbidity,ckd=ckd,
             treatment=treatment, event_24=rbinom(n,1,plogis(lp_out)))
}
truth_fn <- function(n=2e6){ set.seed(99)
  age<-rnorm(n,55,10); sex<-rbinom(n,1,.55); bio<-rnorm(n); com<-sample(0:2,n,T,c(.5,.3,.2)); ckd<-rbinom(n,1,.12)
  lp<- -2.5+0.015*(age-55)+0.3*sex+0.2*bio+0.6*ckd+0.25*com
  mean(plogis(lp+(-0.05/0.15)))-mean(plogis(lp)) }
truth <- truth_fn()
dat <- gen(n); W <- dat[, COVARS]; A <- dat$treatment; Y <- dat$event_24
lib <- c("SL.glm")

# Reference: tmle::tmle
ref <- tmle(Y, A, W, family = "binomial", Q.SL.library = lib, g.SL.library = lib)
ref_rd <- ref$estimates$ATE$psi; ref_se <- sqrt(ref$estimates$ATE$var.psi)

# cleanTMLE cross-fitted TMLE (same GLM library)
ct <- cleanTMLE::estimate_tmle_risk_point(
  data = dat, treatment = "treatment", outcome = "event_24",
  covariates = COVARS, sl_library = lib, truncate = 0.01, n_folds = 1L)
ct_rd <- ct$estimates$ATE$estimate; ct_se <- ct$estimates$ATE$se

tab <- data.frame(
  source     = c("Analytic truth", "tmle::tmle (GLM)", "cleanTMLE TMLE (GLM)"),
  RD         = round(c(truth, ref_rd, ct_rd), 5),
  SE         = round(c(NA, ref_se, ct_se), 5),
  diff_vs_tmle = round(c(NA, 0, ct_rd - ref_rd), 5))
cat(sprintf("n = %d   analytic truth RD = %.5f\n\n", n, truth))
print(tab, row.names = FALSE)
dir.create("sandbox/validation", showWarnings = FALSE, recursive = TRUE)
write.csv(tab, "sandbox/validation/validation_vs_tmle.csv", row.names = FALSE)
cat("\nVALIDATION_DONE\n")
