# Run the manuscript's key analysis chunks against the smoke results
# to validate the R code paths (headline table, selection cross-tab,
# FIORD sensitivity, forest plot). Catches type/column bugs that the
# placeholder render would silently skip.
options(warn = 1L)
setwd(here::here("plasmode_selection_paper"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

sim <- readRDS("results/sim_smoke_all.rds")
true_RD <- -0.05
sim$bias <- sim$estimate - true_RD
sim$covers <- (sim$ci_lower <= true_RD) & (sim$ci_upper >= true_RD)

cat("=== HEADLINE TABLE ===\n")
summary_tbl <- sim |>
  group_by(dgp, workflow) |>
  summarise(
    n_valid  = sum(!is.na(estimate)),
    bias     = mean(bias, na.rm = TRUE),
    abs_bias = mean(abs(bias), na.rm = TRUE),
    emp_sd   = sd(estimate, na.rm = TRUE),
    mean_se  = mean(se, na.rm = TRUE),
    se_sd    = mean(se, na.rm = TRUE) / sd(estimate, na.rm = TRUE),
    rmse     = sqrt(mean(bias^2, na.rm = TRUE)),
    coverage = mean(covers, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(dgp, workflow)
print(as.data.frame(summary_tbl), digits = 3)

cat("\n=== SELECTION CROSS-TAB (min_rmse, complete columns) ===\n")
all_cands <- c("param_t01", "paramplus_t01",
               "rich_t01", "rich_t05",
               "rich_screener_t01", "rich_screener_t05")
pick <- sim |>
  filter(workflow == "plasmode_min_rmse") |>
  mutate(
    fallback = grepl("_FALLBACK$", selected_candidate),
    sel      = factor(sub("_FALLBACK$", "", selected_candidate),
                      levels = all_cands)
  ) |>
  group_by(dgp, sel) |>
  summarise(n = n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = sel, values_from = n,
                     values_fill = 0L, names_expand = TRUE)
print(as.data.frame(pick))

cat("\n=== FIORD SENSITIVITY ===\n")
fiord_compare <- sim |>
  filter(workflow %in% c("plasmode_min_rmse",
                          "plasmode_fiord_two_stage")) |>
  group_by(dgp, workflow) |>
  summarise(
    abs_bias = mean(abs(bias), na.rm = TRUE),
    rmse     = sqrt(mean(bias^2, na.rm = TRUE)),
    coverage = mean(covers, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(fiord_compare), digits = 3)

cat("\n=== Forest plot data ok? ===\n")
fp <- sim |>
  filter(workflow %in% c("fixed_parametric","fixed_rich",
                          "plasmode_min_rmse"))
cat("Forest-plot rows:", nrow(fp),
    "  workflows:", paste(unique(fp$workflow), collapse = ", "), "\n")

cat("\nAll manuscript chunks executed without error.\n")
