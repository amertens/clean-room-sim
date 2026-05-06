# ============================================================
# Clean-Room Workflow: Plotting Functions
# ============================================================

#' PS overlap density plot
plot_ps_overlap <- function(ps, A, file_path = NULL) {
  require(ggplot2)
  df <- data.frame(
    ps = ps[!is.na(ps) & !is.na(A)],
    group = factor(ifelse(A[!is.na(ps) & !is.na(A)] == 1,
                          "Rescue.Co (treated)", "Non-Rescue.Co (control)"))
  )

  p <- ggplot(df, aes(x = ps, fill = group)) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Propensity Score Overlap",
      x = "Propensity Score",
      y = "Density",
      fill = "Group"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("steelblue", "coral"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 5)
    message("Saved: ", file_path)
  }
  p
}

#' Love plot (SMD before and after matching)
plot_love <- function(love_data, file_path = NULL) {
  require(ggplot2)

  df <- data.frame(
    covariate = rep(love_data$covariate, 2),
    smd = c(abs(love_data$smd_pre), abs(love_data$smd_post)),
    stage = rep(c("Before matching", "After matching"), each = nrow(love_data))
  )
  df$stage <- factor(df$stage, levels = c("Before matching", "After matching"))

  # Sort by pre-matching SMD
  ord <- order(abs(love_data$smd_pre))
  df$covariate <- factor(df$covariate, levels = love_data$covariate[ord])

  # Limit to top 30 covariates for readability
  top_covs <- love_data$covariate[order(-abs(love_data$smd_pre))][1:min(30, nrow(love_data))]
  df <- df[df$covariate %in% top_covs, ]

  p <- ggplot(df, aes(x = smd, y = covariate, color = stage, shape = stage)) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray50") +
    labs(
      title = "Covariate Balance: Love Plot",
      x = "Absolute Standardized Mean Difference",
      y = NULL,
      color = NULL, shape = NULL
    ) +
    theme_minimal() +
    scale_color_manual(values = c("Before matching" = "coral",
                                   "After matching" = "steelblue"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 10, height = 8)
    message("Saved: ", file_path)
  }
  p
}

#' Weight distribution plot
plot_weight_distribution <- function(weights, A, file_path = NULL) {
  require(ggplot2)
  df <- data.frame(
    weight = weights,
    group = factor(ifelse(A == 1, "Treated", "Control"))
  )

  p <- ggplot(df, aes(x = weight, fill = group)) +
    geom_histogram(alpha = 0.6, bins = 50, position = "identity") +
    labs(
      title = "IPTW Weight Distribution",
      x = "Weight",
      y = "Count",
      fill = "Group"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("steelblue", "coral"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 5)
    message("Saved: ", file_path)
  }
  p
}

#' Survival curves plot
plot_survival_curves <- function(surv_tmle_results, file_path = NULL) {
  require(ggplot2)

  if (is.null(surv_tmle_results) || is.null(surv_tmle_results$estimates)) {
    warning("No survival TMLE results to plot")
    return(invisible(NULL))
  }

  est <- surv_tmle_results$estimates
  df <- data.frame(
    time = rep(est$time, 2),
    risk = c(est$risk_1, est$risk_0),
    group = rep(c("Rescue.Co (A=1)", "Non-Rescue.Co (A=0)"), each = nrow(est))
  )

  p <- ggplot(df, aes(x = time, y = risk, color = group)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    labs(
      title = "Counterfactual Cumulative Risk Curves (Survival TMLE)",
      x = "Time (days)",
      y = "Cumulative Risk",
      color = "Group"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("steelblue", "coral"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 5)
    message("Saved: ", file_path)
  }
  p
}

#' Simulation performance comparison plot
plot_simulation_performance <- function(sim_metrics, file_path = NULL) {
  require(ggplot2)

  df <- sim_metrics[, c("estimator", "bias", "rmse", "coverage")]
  df_long <- reshape(df, varying = c("bias", "rmse", "coverage"),
                     v.names = "value", timevar = "metric",
                     times = c("bias", "rmse", "coverage"),
                     direction = "long")

  p <- ggplot(df_long, aes(x = estimator, y = value, fill = estimator)) +
    geom_col() +
    facet_wrap(~ metric, scales = "free_y") +
    labs(title = "Simulation Performance Comparison",
         x = NULL, y = "Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 10, height = 5)
    message("Saved: ", file_path)
  }
  p
}

#' Forest plot comparing estimators
plot_estimator_comparison <- function(results_df, file_path = NULL) {
  require(ggplot2)

  # Expect columns: method, estimate, ci_lower, ci_upper
  p <- ggplot(results_df, aes(x = estimate, y = method)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Estimator Comparison: Treatment Effect Estimates",
      x = "Risk Difference",
      y = NULL
    ) +
    theme_minimal()

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 5)
    message("Saved: ", file_path)
  }
  p
}

# ============================================================
# Extended Diagnostics
# ============================================================

#' Three-way love plot: unweighted, matched, and IPTW-weighted SMDs
plot_weighted_love <- function(smd_unweighted, smd_matched, smd_weighted,
                               file_path = NULL) {
  require(ggplot2)

  df <- data.frame(
    covariate = rep(smd_unweighted$covariate, 3),
    smd       = c(abs(smd_unweighted$smd), abs(smd_matched$smd), abs(smd_weighted$smd)),
    method    = rep(c("Unweighted", "PS matched", "IPTW weighted"),
                    each = nrow(smd_unweighted)),
    stringsAsFactors = FALSE
  )
  df$method <- factor(df$method, levels = c("Unweighted", "PS matched", "IPTW weighted"))

  ord <- order(abs(smd_unweighted$smd))
  df$covariate <- factor(df$covariate, levels = smd_unweighted$covariate[ord])

  top_covs <- smd_unweighted$covariate[order(-abs(smd_unweighted$smd))][
    1:min(30, nrow(smd_unweighted))
  ]
  df <- df[df$covariate %in% top_covs, ]

  p <- ggplot(df, aes(x = smd, y = covariate, color = method, shape = method)) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray50") +
    labs(
      title = "Covariate Balance: Unweighted vs. Matched vs. IPTW",
      x     = "Absolute Standardized Mean Difference",
      y     = NULL, color = NULL, shape = NULL
    ) +
    theme_minimal() +
    scale_color_manual(values = c("Unweighted" = "coral",
                                   "PS matched" = "steelblue",
                                   "IPTW weighted" = "forestgreen"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 10, height = 8)
    message("Saved: ", file_path)
  }
  p
}

#' Influence curve diagnostic plot
plot_influence_curve <- function(ic_values, A, file_path = NULL) {
  require(ggplot2)

  df <- data.frame(
    ic    = ic_values,
    group = factor(ifelse(A == 1, "Treated", "Control"))
  )

  p <- ggplot(df, aes(x = ic, fill = group)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Influence Curve (IC) Distribution by Treatment Group",
      x     = "Influence Curve Value",
      y     = "Density",
      fill  = "Group"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("steelblue", "coral"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 5)
    message("Saved: ", file_path)
  }
  p
}

#' Clever covariate (H_n) diagnostic plot
plot_clever_covariate <- function(H_values, A, file_path = NULL) {
  require(ggplot2)

  df <- data.frame(
    H     = H_values,
    group = factor(ifelse(A == 1, "Treated", "Control"))
  )

  p <- ggplot(df, aes(x = group, y = H, fill = group)) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "Clever Covariate (H_n) by Treatment Group",
      x     = NULL,
      y     = "H_n = A/g(1|W) - (1-A)/g(0|W)",
      fill  = "Group"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("steelblue", "coral"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 7, height = 5)
    message("Saved: ", file_path)
  }
  p
}

#' Survival risk curves with confidence bands
plot_survival_curves_ci <- function(surv_tmle_results, file_path = NULL) {
  require(ggplot2)

  if (is.null(surv_tmle_results) || is.null(surv_tmle_results$estimates)) {
    warning("No survival TMLE results to plot")
    return(invisible(NULL))
  }

  est <- surv_tmle_results$estimates
  has_se <- "rd_se" %in% names(est)

  df <- data.frame(
    time  = rep(est$time, 2),
    risk  = c(est$risk_1, est$risk_0),
    group = rep(c("Rescue.Co (A=1)", "Non-Rescue.Co (A=0)"), each = nrow(est))
  )

  if (has_se) {
    df$lower <- df$risk - 1.96 * rep(est$rd_se, 2)
    df$upper <- df$risk + 1.96 * rep(est$rd_se, 2)
    df$lower <- pmax(df$lower, 0)
    df$upper <- pmin(df$upper, 1)
  }

  p <- ggplot(df, aes(x = time, y = risk, color = group)) +
    geom_line(size = 1.2) +
    geom_point(size = 3)

  if (has_se) {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, fill = group),
                         alpha = 0.15, color = NA)
  }

  p <- p +
    labs(
      title = "Counterfactual Cumulative Risk Curves (Survival TMLE)",
      x     = "Time (days)",
      y     = "Cumulative Risk",
      color = "Group", fill = "Group"
    ) +
    theme_minimal() +
    scale_color_manual(values = c("steelblue", "coral")) +
    scale_fill_manual(values = c("steelblue", "coral"))

  if (!is.null(file_path)) {
    ggsave(file_path, p, width = 8, height = 5)
    message("Saved: ", file_path)
  }
  p
}
