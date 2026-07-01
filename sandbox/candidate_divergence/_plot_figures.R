library(ggplot2)
library(pkgload)
pkgload::load_all("cleanTMLE", quiet = TRUE)

out    <- readRDS("sandbox/candidate_divergence/results/candidate_divergence_full.rds")
figdir <- "sandbox/candidate_divergence/figures"
dir.create(figdir, showWarnings = FALSE)

cand_levels <- c("aggressive", "middle", "robust")
cand_labels <- c("aggressive\n(trunc 0.001)", "middle\n(trunc 0.025)", "robust\n(trunc 0.20)")
cand_cols_id <- c(aggressive = "#D55E00", middle = "#009E73", robust = "#0072B2")
# Duplicate with display-label names so both factor and ID versions work
cand_cols  <- setNames(cand_cols_id, cand_labels)
cand_fills <- cand_cols
thr_ratio   <- out$config$lock_thresholds$max_rmse_ratio   # 2.0

cell <- out$dq_cell_agg
base_tbl <- out$baseline_pooled
base_rmse <- setNames(base_tbl$baseline_rmse, base_tbl$candidate)
cell$candidate <- factor(cell$candidate, levels = cand_levels, labels = cand_labels)
cell$baseline  <- base_rmse[as.character(cell$candidate)]

# fix: map labels back to IDs for baseline lookup
id_to_label <- setNames(cand_labels, cand_levels)
label_to_id <- setNames(cand_levels, cand_labels)
cell$baseline <- base_rmse[label_to_id[as.character(cell$candidate)]]
cell$ratio     <- cell$rmse / cell$baseline
cell$ratio_lo  <- pmax(0, (cell$rmse - cell$rmse_mcse) / cell$baseline)
cell$ratio_hi  <-     (cell$rmse + cell$rmse_mcse) / cell$baseline

# Per-candidate threshold lines (each candidate's own baseline × 2)
thresh_df <- data.frame(
  candidate = factor(cand_labels, levels = cand_labels),
  yint = thr_ratio   # all 2.0 × their own baselines
)

# ── Pull out separate threat subsets ──────────────────────────────────────────
pos  <- cell[cell$scenario == "near_positivity", ]
uu   <- cell[cell$scenario == "unmeasured_U", ]
miss <- cell[cell$scenario == "cov_miss", ]
pos$x_num <- as.numeric(sub("slope_x", "", pos$level))
uu$x_num  <- as.numeric(sub("OR_trt([0-9.]+)_out.*", "\\1", uu$level, perl=TRUE))
miss$x_num <- as.numeric(miss$level)

# Degradation table — all three threats, absolute RMSE
all_deg <- cell[cell$scenario != "none", ]
all_deg$x_num <- NA
all_deg$x_num[all_deg$scenario == "near_positivity"] <-
  as.numeric(sub("slope_x", "", all_deg$level[all_deg$scenario == "near_positivity"]))
all_deg$x_num[all_deg$scenario == "unmeasured_U"] <-
  as.numeric(sub("OR_trt([0-9.]+)_out.*", "\\1",
                 all_deg$level[all_deg$scenario == "unmeasured_U"], perl=TRUE))
all_deg$x_num[all_deg$scenario == "cov_miss"] <-
  as.numeric(all_deg$level[all_deg$scenario == "cov_miss"])

# Tidy agg for bar chart (baseline + worst case)
agg <- out$agg
agg$candidate <- factor(agg$candidate, levels = cand_levels, labels = cand_labels)

theme_pub <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: Baseline vs Worst-case RMSE — divergence bar chart
# ═══════════════════════════════════════════════════════════════════════════════
long_bar <- rbind(
  data.frame(candidate = agg$candidate, rmse = agg$baseline_rmse,
             se = agg$baseline_mcse, type = "Baseline RMSE"),
  data.frame(candidate = agg$candidate, rmse = agg$worst_rmse,
             se = agg$worst_mcse,    type = "Worst-case RMSE")
)
long_bar$type <- factor(long_bar$type, levels = c("Baseline RMSE", "Worst-case RMSE"))

# Add rule-winner annotations
winner_df <- data.frame(
  type      = factor(c("Baseline RMSE", "Worst-case RMSE"),
                     levels = c("Baseline RMSE", "Worst-case RMSE")),
  candidate = factor(c(id_to_label[out$win_min_final], id_to_label[out$win_mm_final]),
                     levels = cand_labels),
  label     = c("min_rmse\nwinner", "min_max_rmse\nwinner"),
  rmse      = c(agg$baseline_rmse[as.character(agg$candidate) == id_to_label[out$win_min_final]],
                agg$worst_rmse[as.character(agg$candidate)   == id_to_label[out$win_mm_final]])
)

p1 <- ggplot(long_bar, aes(candidate, rmse, fill = candidate)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_errorbar(aes(ymin = rmse - 1.96 * se, ymax = rmse + 1.96 * se),
                width = 0.2, colour = "grey30") +
  geom_text(data = winner_df, aes(label = label, y = rmse + 0.005),
            size = 3, colour = "grey20", hjust = 0.5, vjust = 0) +
  scale_fill_manual(values = cand_cols) +
  facet_wrap(~type, scales = "free_y") +
  labs(title = "Rule divergence: min_rmse and min_max_rmse select different candidates",
       subtitle = sprintf("min_rmse → %s    |    min_max_rmse → %s    |    separation %.4f (MC SE %.4f)",
                          id_to_label[out$win_min_final],
                          id_to_label[out$win_mm_final],
                          out$separation, out$combined_mcse),
       x = NULL, y = "RMSE", fill = NULL) +
  theme_pub + theme(legend.position = "none")

ggsave(file.path(figdir, "fig1_divergence_bar.png"), p1, width = 9, height = 5, dpi = 150)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: Near-positivity degradation — absolute RMSE with ribbon
# ═══════════════════════════════════════════════════════════════════════════════
base_lines_pos <- data.frame(
  candidate = factor(cand_labels, levels = cand_labels),
  base = base_rmse[cand_levels],
  threshold = base_rmse[cand_levels] * thr_ratio
)

p2 <- ggplot(pos, aes(x_num, rmse, colour = candidate, fill = candidate)) +
  geom_hline(data = base_lines_pos, aes(yintercept = threshold),
             linetype = "dashed", colour = "grey30") +
  geom_ribbon(aes(ymin = rmse - rmse_mcse, ymax = rmse + rmse_mcse),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  geom_hline(data = base_lines_pos, aes(yintercept = base),
             linetype = "dotted", colour = "grey60") +
  scale_colour_manual(values = cand_cols) +
  scale_fill_manual(values = cand_fills) +
  facet_wrap(~candidate, ncol = 3, scales = "fixed") +
  labs(title = "Near-positivity threat: RMSE by PS-slope amplification factor",
       subtitle = "Dashed = 2× baseline threshold   Dotted = baseline RMSE   Ribbon = ±1 MC SE",
       x = "PS slope amplification (severity)",
       y = "RMSE", colour = NULL, fill = NULL) +
  theme_pub + theme(legend.position = "none")

ggsave(file.path(figdir, "fig2_positivity_absolute.png"), p2, width = 9, height = 4, dpi = 150)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: Unmeasured-confounding degradation — absolute RMSE
# ═══════════════════════════════════════════════════════════════════════════════
base_lines_uu <- base_lines_pos   # same threshold structure

p3 <- ggplot(uu, aes(x_num, rmse, colour = candidate, fill = candidate)) +
  geom_hline(data = base_lines_uu, aes(yintercept = threshold),
             linetype = "dashed", colour = "grey30") +
  geom_ribbon(aes(ymin = rmse - rmse_mcse, ymax = rmse + rmse_mcse),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = cand_cols) +
  scale_fill_manual(values = cand_fills) +
  facet_wrap(~candidate, ncol = 3, scales = "fixed") +
  labs(title = "Unmeasured-confounding threat: RMSE by confounder strength (OR on A and Y)",
       subtitle = "Dashed = 2× baseline threshold   Ribbon = ±1 MC SE",
       x = "Unmeasured-confounder OR",
       y = "RMSE", colour = NULL, fill = NULL) +
  theme_pub + theme(legend.position = "none")

ggsave(file.path(figdir, "fig3_confounding_absolute.png"), p3, width = 9, height = 4, dpi = 150)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: RMSE ratio across ALL threats — candidates overlaid
# ═══════════════════════════════════════════════════════════════════════════════
# Create clean x-axis labels per scenario
pos$x_label  <- sprintf("×%.0f", pos$x_num)
uu$x_label   <- sprintf("OR %.0f", uu$x_num)
miss$x_label <- sprintf("%.0f%%", miss$x_num * 100)

pos$scen_lab  <- "Near-positivity"
uu$scen_lab   <- "Unmeasured confounding"
miss$scen_lab <- "Covariate missingness"

all_ratio <- do.call(rbind, list(
  data.frame(pos[, c("candidate","rmse","rmse_mcse","baseline","ratio",
                     "ratio_lo","ratio_hi","x_num","x_label","scen_lab")]),
  data.frame(uu[,  c("candidate","rmse","rmse_mcse","baseline","ratio",
                     "ratio_lo","ratio_hi","x_num","x_label","scen_lab")]),
  data.frame(miss[,c("candidate","rmse","rmse_mcse","baseline","ratio",
                     "ratio_lo","ratio_hi","x_num","x_label","scen_lab")])
))
all_ratio$scen_lab <- factor(all_ratio$scen_lab,
  levels = c("Near-positivity","Unmeasured confounding","Covariate missingness"))
all_ratio$x_ord <- interaction(all_ratio$scen_lab, all_ratio$x_num, drop = TRUE)

p4 <- ggplot(all_ratio, aes(x_label, ratio, colour = candidate, group = candidate)) +
  geom_hline(yintercept = thr_ratio, linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  annotate("text", x = -Inf, y = thr_ratio, label = " 2× threshold",
           hjust = 0, vjust = -0.4, size = 3, colour = "grey30") +
  geom_ribbon(aes(ymin = ratio_lo, ymax = ratio_hi, fill = candidate),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = cand_cols, labels = cand_labels) +
  scale_fill_manual(values = cand_fills,   labels = cand_labels) +
  facet_wrap(~scen_lab, scales = "free_x", nrow = 1) +
  labs(title = "RMSE ratio to each candidate's own baseline — all three threats",
       subtitle = "Values above 2.0 (dashed) exceed the prespecified locked threshold",
       x = "Threat severity",
       y = "RMSE / baseline RMSE",
       colour = "Candidate", fill = "Candidate") +
  theme_pub +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

ggsave(file.path(figdir, "fig4_ratio_all_threats.png"), p4, width = 11, height = 4.5, dpi = 150)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: Coverage degradation across all threats
# ═══════════════════════════════════════════════════════════════════════════════
dq_cov <- out$dq_metrics[out$dq_metrics$scenario != "none", ]
dq_cov$candidate <- factor(dq_cov$candidate, levels = cand_levels, labels = cand_labels)
dq_cov$x_num <- NA
dq_cov$x_num[dq_cov$scenario == "near_positivity"] <-
  as.numeric(sub("slope_x","", dq_cov$level[dq_cov$scenario == "near_positivity"]))
dq_cov$x_num[dq_cov$scenario == "unmeasured_U"] <-
  as.numeric(sub("OR_trt([0-9.]+)_out.*","\\1",
                 dq_cov$level[dq_cov$scenario == "unmeasured_U"], perl=TRUE))
dq_cov$x_num[dq_cov$scenario == "cov_miss"] <- as.numeric(dq_cov$level[dq_cov$scenario == "cov_miss"])
dq_cov$scen_lab <- recode_scen <- c(near_positivity = "Near-positivity",
                                    unmeasured_U    = "Unmeasured confounding",
                                    cov_miss        = "Covariate missingness")[dq_cov$scenario]
dq_cov$scen_lab <- factor(dq_cov$scen_lab,
  levels = c("Near-positivity","Unmeasured confounding","Covariate missingness"))

baseline_cov <- out$plas_metrics[, c("candidate","coverage")]
baseline_cov$candidate <- factor(baseline_cov$candidate, levels = cand_levels, labels = cand_labels)

min_cov_thresh <- out$config$lock_thresholds$min_coverage   # 0.88

p5 <- ggplot(dq_cov, aes(x_num, coverage, colour = candidate, group = candidate)) +
  geom_hline(yintercept = min_cov_thresh, linetype = "dashed",
             colour = "firebrick", linewidth = 0.7) +
  annotate("text", x = -Inf, y = min_cov_thresh,
           label = sprintf(" Locked floor (%.0f%%)", min_cov_thresh * 100),
           hjust = 0, vjust = -0.4, size = 3, colour = "firebrick") +
  geom_hline(yintercept = 0.95, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = cand_cols, labels = cand_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  facet_wrap(~scen_lab, scales = "free_x", nrow = 1) +
  labs(title = "95% CI coverage degradation across all three threats",
       subtitle = sprintf("Red dashed = locked minimum coverage (%.0f%%)   Dotted = nominal 95%%",
                          min_cov_thresh * 100),
       x = "Threat severity", y = "Coverage", colour = "Candidate") +
  theme_pub +
  guides(colour = guide_legend(nrow = 1))

ggsave(file.path(figdir, "fig5_coverage.png"), p5, width = 11, height = 4.5, dpi = 150)

# ═══════════════════════════════════════════════════════════════════════════════
# Figure 6: Per-batch stability — worst-case RMSE across 5 independent batches
# ═══════════════════════════════════════════════════════════════════════════════
worst_long <- do.call(rbind, lapply(seq_len(ncol(out$worst_mat)), function(b) {
  data.frame(
    batch     = b,
    candidate = factor(rownames(out$worst_mat), levels = cand_levels, labels = cand_labels),
    worst_rmse = out$worst_mat[, b]
  )
}))
# colour the min-max winner per batch
win_mm_label <- id_to_label[out$win_mm_final]
win_min_label <- id_to_label[out$win_min_final]

p6 <- ggplot(worst_long, aes(factor(batch), worst_rmse, colour = candidate, group = candidate)) +
  geom_line(linewidth = 0.8, alpha = 0.6) +
  geom_point(aes(size = (candidate == win_mm_label)), shape = 19) +
  scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2), guide = "none") +
  scale_colour_manual(values = cand_cols, labels = cand_labels) +
  labs(title = "Worst-case RMSE is stable across 5 independent batches",
       subtitle = sprintf("Large dots = min_max_rmse winner (%s)   rule changes decision in all 5 batches",
                          win_mm_label),
       x = "Batch", y = "Worst-case RMSE", colour = "Candidate") +
  theme_pub +
  guides(colour = guide_legend(nrow = 1))

ggsave(file.path(figdir, "fig6_batch_stability.png"), p6, width = 8, height = 4.5, dpi = 150)

cat("All 6 figures saved to", figdir, "\n")
