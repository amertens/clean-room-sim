# ============================================================
# Proportional-odds (PO) assumption test for the ordinal GOSE outcome
# ============================================================
# Produces rescueCo/results/gose_ordinal_po.csv.
#
# Reconciliation note (2026-07-28): an earlier clean-room artifact reported a
# borderline PO test (LR = 11.905, p = 0.064) computed with a MASS::polr-based
# omnibus PO test on the transfers-only cohort. The main pipeline reports a
# clearly non-rejecting nominal-effects LR test (p = 0.51). To reconcile, this
# script STANDARDISES on ordinal::nominal_test() — the likelihood-ratio test of
# proportionality used by the coauthors (analysis/11_meeting_followup_analyses.R)
# — and runs it on the ALIGNED primary cohort (interfacility transfers AND
# prior-care-elsewhere excluded). The common odds ratio is still reported from
# MASS::polr for interpretability. See reconciliation_2026-07-28.md (Discrepancy 2).
# ============================================================

suppressPackageStartupMessages({
  library(MASS)
  library(ordinal)
})

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")

cfg <- load_cr_config()
cr_log("=== Proportional-odds assumption test (ordinal GOSE) ===")

dat <- load_stage_output("stage1_cohort.rds")

# Analytic subset: observed GOSE and treatment on the aligned primary cohort.
pc <- dat[!is.na(dat$gose_score) & !is.na(dat$A), ]
cr_log(sprintf("PO test cohort: n = %d (observed GOSE + A on aligned primary cohort)", nrow(pc)))

# Explicit 1..8 levels keep nominal_test stable when some categories are absent.
pc$gose_f <- factor(pc$gose_score, levels = 1:8, ordered = TRUE)

# --- Common odds ratio from the proportional-odds model (MASS::polr) ---
m_polr    <- polr(gose_f ~ A, data = pc, Hess = TRUE)
beta_A    <- coef(m_polr)["A"]
se_A      <- sqrt(diag(vcov(m_polr)))["A"]
common_or <- exp(beta_A)
or_lo     <- exp(beta_A - 1.96 * se_A)
or_hi     <- exp(beta_A + 1.96 * se_A)
p_or      <- 2 * pnorm(abs(beta_A / se_A), lower.tail = FALSE)

# --- PO assumption test: ordinal::nominal_test (standardised LR test) ---
clm_fit <- clm(gose_f ~ A, data = pc)
nt      <- nominal_test(clm_fit)
lr_stat <- nt["A", "LRT"]
lr_df   <- nt["A", "Df"]
lr_p    <- nt["A", "Pr(>Chi)"]
cr_log(sprintf("nominal_test on A: LR = %.3f, df = %d, p = %.4f", lr_stat, lr_df, lr_p))

result <- data.frame(
  estimator  = "Unadjusted PO (MASS::polr); PO test = ordinal::nominal_test",
  common_OR  = round(common_or, 3),
  ci_lo      = round(or_lo, 3),
  ci_hi      = round(or_hi, 3),
  p_value    = round(p_or, 4),
  PO_test_LR = round(lr_stat, 3),
  PO_test_df = lr_df,
  PO_test_p  = round(lr_p, 4),
  n          = nrow(pc),
  stringsAsFactors = FALSE
)

write.csv(result, file.path(cfg$paths$results, "gose_ordinal_po.csv"), row.names = FALSE)
cr_log("Saved gose_ordinal_po.csv")
print(result)

verdict <- if (lr_p < 0.05) "PO assumption REJECTED (p < 0.05)" else
  "PO assumption NOT rejected (p >= 0.05)"
cr_log(verdict)
