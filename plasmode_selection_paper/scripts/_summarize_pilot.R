# Build a plain-text summary of the pilot results so the headline
# numbers can be committed alongside the gitignored sim_pilot.rds.
options(warn = 1L)
setwd(here::here("plasmode_selection_paper"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})
sim <- readRDS("results/sim_pilot.rds")
true_RD <- -0.05
sim$bias   <- sim$estimate - true_RD
sim$covers <- (sim$ci_lower <= true_RD) & (sim$ci_upper >= true_RD)
n_reps_per_cell <- length(unique(sim$rep))

con <- file("results/sim_pilot_summary.txt", open = "wt")
on.exit(close(con))
say <- function(...) cat(..., "\n", sep = "", file = con)

say("Plasmode-selection paper -- pilot summary")
say("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
say("Source: results/sim_pilot.rds (", nrow(sim), " rows; ",
    n_reps_per_cell, " outer reps x 4 workflows x ",
    length(unique(sim$dgp)), " DGPs)")
say("Pilot config: n=500, inner_reps=10, true marginal RD=-0.05")
say("Monte Carlo SE of coverage at n_reps=", n_reps_per_cell,
    ": ~", sprintf("%.3f", sqrt(0.5 * 0.5 / n_reps_per_cell)),
    " (pilot is NOT inferentially defensible)")
say("")
say("==== Headline operating characteristics ====")
sumtbl <- sim |>
  group_by(dgp, workflow) |>
  summarise(
    bias     = mean(bias, na.rm = TRUE),
    rmse     = sqrt(mean(bias^2, na.rm = TRUE)),
    emp_sd   = sd(estimate, na.rm = TRUE),
    mean_se  = mean(se, na.rm = TRUE),
    se_sd    = mean(se, na.rm = TRUE) / sd(estimate, na.rm = TRUE),
    coverage = mean(covers, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(dgp, workflow)
capture.output(print(as.data.frame(sumtbl), digits = 3),
               file = con, append = TRUE)
say("")

say("==== Candidate-selection cross-tab (min_rmse) ====")
all_cands <- c("param_t01", "paramplus_t01",
               "rich_t01", "rich_t05",
               "rich_screener_t01", "rich_screener_t05")
pick <- sim |>
  filter(workflow == "plasmode_min_rmse") |>
  mutate(sel = factor(sub("_FALLBACK$", "", selected_candidate),
                      levels = all_cands)) |>
  group_by(dgp, sel) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = sel, values_from = n,
              values_fill = 0L, names_expand = TRUE)
capture.output(print(as.data.frame(pick)),
               file = con, append = TRUE)
say("")

say("==== Candidate-selection cross-tab (fiord_two_stage) ====")
pick2 <- sim |>
  filter(workflow == "plasmode_fiord_two_stage") |>
  mutate(sel = factor(sub("_FALLBACK$", "", selected_candidate),
                      levels = all_cands)) |>
  group_by(dgp, sel) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = sel, values_from = n,
              values_fill = 0L, names_expand = TRUE)
capture.output(print(as.data.frame(pick2)),
               file = con, append = TRUE)
say("")

cat("Wrote results/sim_pilot_summary.txt\n")
