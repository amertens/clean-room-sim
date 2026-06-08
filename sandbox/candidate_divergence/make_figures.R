#!/usr/bin/env Rscript
# Build the deliverable figures from candidate_divergence_<mode>.rds:
#   figures/degradation_gradient.png  - RMSE-ratio-to-baseline vs threat
#       severity, one line per candidate, horizontal line at the locked
#       max_rmse_ratio threshold. Left panel: near-positivity (aggressive
#       crosses, robust stays under). Right panel: unmeasured confounding
#       (robust climbs highest, middle lowest -> minimax picks middle).
#   figures/selection_table.png       - baseline RMSE, worst-case RMSE, worst
#       threat, and the two rule winners, each with Monte Carlo SE.
suppressWarnings(suppressMessages({
  library(ggplot2)
  has_grid <- requireNamespace("gridExtra", quietly = TRUE)
}))

.this_dir <- tryCatch({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
res_dir <- file.path(.this_dir, "results"); fig_dir <- file.path(.this_dir, "figures")
mode <- { a <- commandArgs(TRUE); if (length(a) >= 1L) a[[1]] else "full" }

out <- readRDS(file.path(res_dir, sprintf("candidate_divergence_%s.rds", mode)))
cand_levels <- c("aggressive", "middle", "robust")
cand_cols   <- c(aggressive = "#D55E00", middle = "#009E73", robust = "#0072B2")
thr <- out$config$lock_thresholds$max_rmse_ratio

cell <- out$dq_cell_agg
base <- setNames(out$baseline_pooled$baseline_rmse, out$baseline_pooled$candidate)
cell$baseline <- base[cell$candidate]
cell$ratio      <- cell$rmse / cell$baseline
cell$ratio_lo   <- (cell$rmse - ifelse(is.na(cell$rmse_mcse), 0, cell$rmse_mcse)) / cell$baseline
cell$ratio_hi   <- (cell$rmse + ifelse(is.na(cell$rmse_mcse), 0, cell$rmse_mcse)) / cell$baseline
cell$candidate  <- factor(cell$candidate, levels = cand_levels)

pos <- cell[cell$scenario == "near_positivity", ]
pos$x <- as.numeric(sub("slope_x", "", pos$level))
uu  <- cell[cell$scenario == "unmeasured_U", ]
uu$x <- as.numeric(sub("OR_trt([0-9.]+)_out.*", "\\1", uu$level))

base_panel <- function(df, xlab, title) {
  ggplot(df, aes(x, ratio, colour = candidate, fill = candidate)) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = "grey30") +
    annotate("text", x = min(df$x), y = thr, vjust = -0.6, hjust = 0,
             label = sprintf("locked threshold (%.1fx baseline)", thr),
             size = 3, colour = "grey30") +
    geom_ribbon(aes(ymin = ratio_lo, ymax = ratio_hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    scale_colour_manual(values = cand_cols) + scale_fill_manual(values = cand_cols) +
    labs(x = xlab, y = "RMSE ratio to baseline", title = title, colour = NULL, fill = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

pA <- base_panel(pos, "Near-positivity severity (PS slope amplification)",
                 "A. Near-positivity threat")
pB <- base_panel(uu, "Unmeasured-confounding strength (OR on A and Y)",
                 "B. Unmeasured-confounding threat")

if (has_grid) {
  g <- gridExtra::arrangeGrob(pA, pB, ncol = 2)
  ggsave(file.path(fig_dir, "degradation_gradient.png"), g,
         width = 11, height = 5, dpi = 150)
} else {
  ggsave(file.path(fig_dir, "degradation_gradient_A.png"), pA, width = 6, height = 5, dpi = 150)
  ggsave(file.path(fig_dir, "degradation_gradient_B.png"), pB, width = 6, height = 5, dpi = 150)
}

# ── Selection table ─────────────────────────────────────────────────────────
agg <- out$agg
worst_threat <- sapply(cand_levels, function(cid) {
  s <- cell[cell$candidate == cid, ]; w <- s[which.max(s$rmse), ]
  sprintf("%s [%s]", w$scenario, w$level)
})
fmt <- function(m, se) ifelse(is.na(se), sprintf("%.4f", m),
                              sprintf("%.4f (%.4f)", m, se))
tab <- data.frame(
  Candidate     = agg$candidate,
  `PS truncation` = c("0.001", "0.025", "0.20")[match(agg$candidate, cand_levels)],
  `Baseline RMSE (MC SE)` = fmt(agg$baseline_rmse, agg$baseline_mcse),
  `Worst-case RMSE (MC SE)` = fmt(agg$worst_rmse, agg$worst_mcse),
  `Worst threat` = worst_threat[agg$candidate],
  `min_rmse`     = ifelse(agg$candidate == out$win_min_final, "<-- selected", ""),
  `min_max_rmse` = ifelse(agg$candidate == out$win_mm_final, "<-- selected", ""),
  check.names = FALSE, stringsAsFactors = FALSE)

cap <- sprintf(paste0("min_rmse selects '%s' (lowest baseline RMSE); ",
  "min_max_rmse selects '%s' (lowest worst-case RMSE). Worst-case separation ",
  "between the two = %.4f (combined MC SE %.4f)."),
  out$win_min_final, out$win_mm_final, out$separation, out$combined_mcse)

if (has_grid) {
  tt <- gridExtra::ttheme_default(base_size = 9, core = list(fg_params = list(hjust = 0, x = 0.02)))
  tg <- gridExtra::tableGrob(tab, rows = NULL, theme = tt)
  capg <- grid::textGrob(paste(strwrap(cap, width = 110), collapse = "\n"),
                         gp = grid::gpar(fontsize = 8), just = "left", x = 0.01, hjust = 0)
  g2 <- gridExtra::arrangeGrob(tg, capg, ncol = 1, heights = c(4, 1))
  ggsave(file.path(fig_dir, "selection_table.png"), g2, width = 11, height = 3.2, dpi = 150)
}
write.csv(tab, file.path(res_dir, "selection_table.csv"), row.names = FALSE)
cat("Figures written to", fig_dir, "\n")
cat(cap, "\n")
print(tab, row.names = FALSE)
