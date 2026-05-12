# Compute the true marginal causal RD = E[Y(1)] - E[Y(0)] under each DGP,
# accounting for the pmax(., 0.01) clamp on p_Y1 = p_Y0 + true_RD.
setwd(here::here("plasmode_selection_paper"))

K <- 5e5L
true_RD_target <- -0.05

# Replicate the outcome-probability surfaces from dgps.R without
# the rbinom() step so we get the *potential-outcome* probabilities.
linpred <- function(dgp, age, sex, biomarker, comorbidity, bmi_z, smoke) {
  switch(dgp,
    linear =
      -2.0 + 0.03*(age-60) + 0.20*sex + 0.40*biomarker +
        0.50*comorbidity + 0.10*bmi_z + 0.40*smoke,
    nonlinear_smooth =
      -2.0 + 0.03*(age-60) + 0.20*sex +
        0.40*biomarker + 0.40*biomarker^2 +
        0.50*comorbidity + 0.10*bmi_z + 0.20*sin(bmi_z) + 0.40*smoke,
    interactions =
      -2.0 + 0.03*(age-60) + 0.20*sex + 0.40*biomarker +
        0.50*comorbidity + 0.10*bmi_z + 0.40*smoke +
        0.40*biomarker*comorbidity + 0.30*bmi_z*smoke,
    sparse =
      -1.2 + 0.6*biomarker,
    high_dim_noise =
      -2.0 + 0.03*(age-60) + 0.20*sex + 0.40*biomarker +
        0.50*comorbidity + 0.10*bmi_z + 0.40*smoke
  )
}

set.seed(7L)
age         <- rnorm(K, 60, 10)
sex         <- rbinom(K, 1, 0.5)
biomarker   <- rnorm(K, 0, 1)
comorbidity <- rbinom(K, 1, 0.3)
bmi_z       <- rnorm(K, 0, 1)
smoke       <- rbinom(K, 1, 0.25)

for (dgp in c("linear","nonlinear_smooth","interactions","sparse",
              "high_dim_noise")) {
  lp <- linpred(dgp, age, sex, biomarker, comorbidity, bmi_z, smoke)
  p_Y0 <- plogis(lp)
  p_Y1 <- pmin(pmax(p_Y0 + true_RD_target, 0.01), 0.99)
  true_RD <- mean(p_Y1) - mean(p_Y0)
  frac_clamped_low  <- mean(p_Y0 + true_RD_target < 0.01)
  frac_clamped_high <- mean(p_Y0 + true_RD_target > 0.99)
  cat(sprintf("%-18s  E[p_Y0]=%.3f  E[p_Y1]=%.3f  true_RD=%+.5f  ",
              dgp, mean(p_Y0), mean(p_Y1), true_RD))
  cat(sprintf("low-clamp=%.1f%%  high-clamp=%.1f%%\n",
              100*frac_clamped_low, 100*frac_clamped_high))
}
