# Figures for the outcome-blind data-quality (DQ) stress simulation.
# Reads results_new/dq_stress_{good_overlap,marginal_overlap,unmeasured_conf}.rds
# Writes PNGs to reports/dq_figures/.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
})

outdir <- "reports/dq_figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

scen_files <- c(
  "Good overlap"            = "results_new/dq_stress_good_overlap.rds",
  "Marginal overlap"        = "results_new/dq_stress_marginal_overlap.rds",
  "Unmeasured confounding"  = "results_new/dq_stress_unmeasured_conf.rds"
)

threat_labels <- c(
  none            = "Baseline (none)",
  cov_miss        = "Covariate missingness",
  trt_misclass    = "Treatment misclass.",
  out_misclass    = "Outcome misclass.",
  unmeasured_U    = "Unmeasured confounder",
  near_positivity = "Near-positivity"
)

# Combine all three design scenarios into one tidy frame.
dat <- bind_rows(lapply(names(scen_files), function(nm) {
  m <- readRDS(scen_files[[nm]])$metrics
  m$design   <- nm
  m
}))

dat <- dat %>%
  mutate(
    design      = factor(design, levels = names(scen_files)),
    threat      = factor(threat_labels[scenario], levels = threat_labels),
    candidate   = factor(candidate,
                         levels = c("aggressive", "middle", "robust"),
                         labels = c("aggressive (trunc=0.001)",
                                    "middle (trunc=0.025)",
                                    "robust (trunc=0.20)"))
  )

# A within-threat ordinal "severity" rank so gradients plot left-to-right.
sev_order <- dat %>%
  distinct(scenario, level) %>%
  group_by(scenario) %>%
  arrange(level, .by_group = TRUE) %>%
  mutate(severity = row_number()) %>%
  ungroup()
dat <- left_join(dat, sev_order, by = c("scenario", "level"))

cand_cols <- c("aggressive (trunc=0.001)" = "#D55E00",
               "middle (trunc=0.025)"     = "#0072B2",
               "robust (trunc=0.20)"      = "#009E73")

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave2 <- function(file, plot, w, h) {
  ggsave(file.path(outdir, file), plot, width = w, height = h, dpi = 150, bg = "white")
  message("wrote ", file)
}

## ── Figure 1: near-positivity gradient, the headline divergence ──────────────
np <- dat %>% filter(scenario == "near_positivity")
np_base <- dat %>% filter(scenario == "none") %>%
  group_by(design) %>% summarise(rmse = mean(rmse), .groups = "drop")

f1 <- ggplot(np, aes(level, rmse, colour = candidate, group = candidate)) +
  geom_hline(data = np_base, aes(yintercept = rmse),
             linetype = "dashed", colour = "grey50", linewidth = .4) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.4) +
  facet_wrap(~design, scales = "free_y") +
  scale_colour_manual(values = cand_cols, name = NULL) +
  labs(title = "Candidate RMSE under near-positivity stress",
       subtitle = "Dashed line = baseline RMSE (no DQ threat). Light truncation runs away as overlap thins.",
       x = "Near-positivity slope multiplier", y = "RMSE") +
  base_theme
ggsave2("fig1_near_positivity_gradient.png", f1, 10, 4.2)

## ── Figure 2: minimax vs min-RMSE selection ──────────────────────────────────
# baseline RMSE (min-RMSE rule) vs worst-case RMSE across all DQ threats (minimax).
sel <- dat %>%
  group_by(design, candidate) %>%
  summarise(baseline   = rmse[scenario == "none"][1],
            worst_case = max(rmse),
            .groups = "drop") %>%
  pivot_longer(c(baseline, worst_case),
               names_to = "rule", values_to = "rmse") %>%
  mutate(rule = recode(rule,
                       baseline   = "Baseline RMSE (min-RMSE rule)",
                       worst_case = "Worst-case RMSE (minimax rule)"))

f2 <- ggplot(sel, aes(candidate, rmse, fill = rule)) +
  geom_col(position = position_dodge(width = .7), width = .65) +
  facet_wrap(~design, scales = "free_y") +
  scale_fill_manual(values = c("Baseline RMSE (min-RMSE rule)" = "#9ECAE1",
                               "Worst-case RMSE (minimax rule)" = "#08519C"),
                    name = NULL) +
  scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
  labs(title = "Min-RMSE vs minimax selection give different winners",
       subtitle = "Lowest baseline RMSE need not be the most robust across data-quality threats.",
       x = NULL, y = "RMSE") +
  base_theme +
  theme(axis.text.x = element_text(size = 9))
ggsave2("fig2_minimax_vs_minrmse.png", f2, 11, 4.6)

## ── Figure 3: RMSE heatmap across every threat level ─────────────────────────
heat <- dat %>%
  mutate(cell = paste(scenario, level)) %>%
  arrange(threat, severity) %>%
  mutate(cell = factor(cell, levels = unique(cell)))

f3 <- ggplot(heat, aes(candidate, cell, fill = rmse)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = formatC(rmse, format = "f", digits = 3)),
            size = 2.6, colour = "grey15") +
  facet_wrap(~design) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "RMSE") +
  scale_x_discrete(labels = c("aggr.", "middle", "robust")) +
  labs(title = "RMSE by candidate across the full DQ threat grid",
       x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.y = element_text(size = 7.5),
        legend.position = "right")
ggsave2("fig3_rmse_heatmap.png", f3, 11, 7)

## ── Figure 4: coverage degradation ───────────────────────────────────────────
f4 <- ggplot(dat, aes(severity, coverage, colour = candidate, group = candidate)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = .8) +
  geom_point(size = 1.8) +
  facet_grid(design ~ threat, scales = "free_x") +
  scale_colour_manual(values = cand_cols, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "95% CI coverage under escalating data-quality threats",
       subtitle = "Dotted line = nominal 0.95. Severity increases left to right within each threat family.",
       x = "Threat severity (ordinal within family)", y = "CI coverage") +
  base_theme +
  theme(strip.text.x = element_text(size = 8),
        panel.spacing.x = unit(.4, "lines"))
ggsave2("fig4_coverage_degradation.png", f4, 12, 6.5)

## ── Figure 5: bias-variance under the unmeasured-confounder sweep ─────────────
uc <- dat %>% filter(scenario == "unmeasured_U")
f5 <- ggplot(uc, aes(severity, bias, colour = candidate, group = candidate)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = .9) +
  geom_point(size = 2) +
  facet_wrap(~design, scales = "free_y") +
  scale_colour_manual(values = cand_cols, name = NULL) +
  scale_x_continuous(breaks = 1:5,
                     labels = c("OR2", "OR3", "OR4", "OR6", "OR8")) +
  labs(title = "Bias under the unmeasured-confounder sweep",
       subtitle = "Confounding strength rises as the U->treatment and U->outcome odds ratios grow.",
       x = "Unmeasured confounder strength (matched OR on treatment and outcome)",
       y = "Bias (estimate - truth)") +
  base_theme
ggsave2("fig5_unmeasured_bias.png", f5, 10, 4.2)

message("\nAll figures written to ", normalizePath(outdir))
