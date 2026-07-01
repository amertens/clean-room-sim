#!/usr/bin/env Rscript
# Generate and save the four main simulation figures from simulation_results.rds
# Outputs: reports/sim_figures/*.png

library(ggplot2)
outdir <- file.path("reports", "sim_figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

res <- readRDS("results_new/simulation_results.rds")

# ── Method colour palette ────────────────────────────────────────────────────
pal <- c(
  Crude       = "#D73027",
  IPTW        = "#FC8D59",
  "PS Match"  = "#FEE090",
  TMLE        = "#74ADD1",
  TMLE_CF     = "#313695",
  Match_TMLE  = "#4DAC26"
)

# ── 1. Boxplots: sampling distribution ──────────────────────────────────────
raw <- do.call(rbind, lapply(names(res$results), function(sc) {
  df <- as.data.frame(res$results[[sc]])
  df$scenario <- res$scenarios[[sc]]$label
  df$truth    <- res$truths[[sc]]$RD
  df
}))
raw$scenario <- factor(raw$scenario, levels = sapply(res$scenarios, `[[`, "label"))

p_box <- ggplot(raw, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.4) +
  geom_hline(aes(yintercept = truth), linetype = "dashed", colour = "red", linewidth = 0.7) +
  facet_wrap(~ scenario, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = pal, na.value = "grey70") +
  labs(x = NULL, y = "Risk-difference estimate",
       title = "Sampling distribution by estimator and scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "fig1_boxplots.png"), p_box,
       width = 7, height = 11, dpi = 150)
cat("Saved fig1_boxplots.png\n")

# ── 2. Forest plot ───────────────────────────────────────────────────────────
max_abs_bias <- 0.02
fp <- do.call(rbind, lapply(names(res$summaries), function(sc) {
  df <- res$summaries[[sc]]
  df$scenario <- res$scenarios[[sc]]$label
  df$truth    <- res$truths[[sc]]$RD
  df
}))
fp$lo <- fp$mean_est - 1.96 * fp$mc_se_bias  # mc_se of the MEAN, not coverage SE
fp$hi <- fp$mean_est + 1.96 * fp$mc_se_bias
# use mc_se from summary column
if ("mc_se" %in% names(fp)) {
  fp$lo <- fp$mean_est - 1.96 * fp$mc_se
  fp$hi <- fp$mean_est + 1.96 * fp$mc_se
} else {
  # fall back: emp_sd / sqrt(n_ok)
  fp$lo <- fp$mean_est - 1.96 * fp$emp_sd / sqrt(fp$n_ok)
  fp$hi <- fp$mean_est + 1.96 * fp$emp_sd / sqrt(fp$n_ok)
}
fp$scenario <- factor(fp$scenario, levels = sapply(res$scenarios, `[[`, "label"))

p_forest <- ggplot(fp, aes(x = mean_est, y = method, colour = method)) +
  geom_vline(aes(xintercept = truth), linetype = "dashed", colour = "grey40") +
  geom_vline(aes(xintercept = truth - max_abs_bias),
             linetype = "dashed", colour = "red", alpha = 0.5) +
  geom_vline(aes(xintercept = truth + max_abs_bias),
             linetype = "dashed", colour = "red", alpha = 0.5) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25) +
  facet_wrap(~ scenario, scales = "free_x", ncol = 1) +
  scale_colour_manual(values = pal, na.value = "grey70") +
  labs(x = "Mean RD estimate (Monte Carlo 95% interval)", y = NULL,
       title = "Forest plot: mean estimates vs decision thresholds (red dashed: ±0.02)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "fig2_forest.png"), p_forest,
       width = 8, height = 10, dpi = 150)
cat("Saved fig2_forest.png\n")

# ── 3. Coverage bar ──────────────────────────────────────────────────────────
plot_data <- do.call(rbind, lapply(names(res$summaries), function(sc) {
  df <- res$summaries[[sc]]
  df$scenario <- res$scenarios[[sc]]$label
  df
}))
plot_data$scenario <- factor(plot_data$scenario, levels = sapply(res$scenarios, `[[`, "label"))

p_cov <- ggplot(plot_data, aes(x = method, y = coverage, fill = scenario)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "red") +
  coord_cartesian(ylim = c(0.0, 1.0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "95% CI Coverage",
       title = "Coverage by estimator and scenario",
       fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "fig3_coverage.png"), p_cov,
       width = 8, height = 5, dpi = 150)
cat("Saved fig3_coverage.png\n")

# ── 4. SE calibration ────────────────────────────────────────────────────────
p_se <- ggplot(plot_data, aes(x = method, y = se_sd_ratio,
                               colour = scenario, shape = scenario)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  scale_y_continuous(limits = c(0.7, 1.4)) +
  labs(x = NULL, y = "Mean SE / Empirical SD",
       title = "Standard-error calibration (1.0 = well-calibrated)",
       colour = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "fig4_se_cal.png"), p_se,
       width = 8, height = 5, dpi = 150)
cat("Saved fig4_se_cal.png\n")

# ── 5. Bootstrap variance comparison ─────────────────────────────────────────
bv <- read.csv("results_new/bootstrap_variance.csv")
bv_long <- rbind(
  transform(bv, se_type = "IF-based",   se_sd_ratio = se_sd_ratio_if,   coverage = coverage_if),
  transform(bv, se_type = "Bootstrap",  se_sd_ratio = se_sd_ratio_boot, coverage = coverage_boot)
)
bv_long$scenario <- gsub(".*: ", "", bv_long$scenario)
bv_long$scenario <- factor(bv_long$scenario,
  levels = c("Very Good Overlap", "Good Overlap", "Marginal Overlap"))

p_boot_cov <- ggplot(bv_long, aes(x = method, y = coverage,
                                   fill = se_type)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "red") +
  facet_wrap(~ scenario, ncol = 3) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "95% CI Coverage", fill = "SE type",
       title = "Bootstrap vs IF-based CI coverage across overlap regimes") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(outdir, "fig5_boot_coverage.png"), p_boot_cov,
       width = 9, height = 5, dpi = 150)
cat("Saved fig5_boot_coverage.png\n")

cat("\nAll figures saved to", outdir, "\n")
