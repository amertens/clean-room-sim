# Verify cleanTMLE modular TMLE reproduces the rescueCo local TMLE estimate
# on the same lock, for the full-cohort full-covariates spec.
suppressPackageStartupMessages({
  library(cleanTMLE)
})

lock <- load_lock("rescueCo/results/stage1_lock_unmasked.rds")
ps_fit <- readRDS("rescueCo/results/stage2_ct_ps_fit.rds")

# cleanTMLE modular TMLE
g <- fit_tmle_treatment_mechanism(lock, ps_fit)
Q <- fit_tmle_outcome_mechanism(lock, g, override_clean_room = TRUE)
upd <- run_tmle_targeting_step(g, Q)
est <- extract_tmle_estimate(upd)

ate <- est$estimates$ATE
cat("cleanTMLE modular TMLE on rescueCo lock:\n")
cat(sprintf("  RD = %.4f [%.4f, %.4f]  SE = %.4f  p = %.4f\n",
            ate$estimate, ate$ci_lower, ate$ci_upper, ate$se, ate$p_value))

# Compare to the local Full-cohort TMLE row from binary_outcome_comparison.csv
lc <- read.csv("rescueCo/results/binary_outcome_comparison.csv",
                stringsAsFactors = FALSE)
fc <- lc[grepl("Full-cohort TMLE", lc$method), ]
cat("\nrescueCo local Full-cohort TMLE (from binary_outcome_comparison.csv):\n")
cat(sprintf("  RD = %.4f [%.4f, %.4f]  SE = %.4f\n",
            fc$estimate, fc$ci_lower, fc$ci_upper, fc$se))

cat(sprintf("\nDifference in RD: %.5f\n",
            ate$estimate - fc$estimate))
