# Validate the R chunk bodies edited in the manuscript against the real files.
setwd(file.path(dirname(normalizePath(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))), "."))

# --- setup chunk helpers (M2/M3) ---
.bv_def <- function(boot_var, scenario_pat, method, col) {
  r <- boot_var[grepl(scenario_pat, boot_var$scenario, ignore.case=TRUE) &
                boot_var$method == method, col, drop=TRUE]
  if (length(r)==0 || all(is.na(r))) NA else r[1]
}
.bv_cov_mcse <- function(boot_var, scenario_pat, method, cov_col) {
  p <- .bv_def(boot_var, scenario_pat, method, cov_col)
  n <- .bv_def(boot_var, scenario_pat, method, "n_reps")
  if (is.na(p) || is.na(n)) NA else sqrt(p * (1 - p) / n)
}
boot_path <- file.path("..", "results_new", "bootstrap_variance.csv")
stopifnot(file.exists(boot_path))
boot_var <- read.csv(boot_path)
mt_A_cov_if   <- .bv_def(boot_var, "^A", "Match_TMLE", "coverage_if")
mt_A_cov_boot <- .bv_def(boot_var, "^A", "Match_TMLE", "coverage_boot")
mt_A_cov_if_se   <- .bv_cov_mcse(boot_var, "^A", "Match_TMLE", "coverage_if")
mt_A_cov_boot_se <- .bv_cov_mcse(boot_var, "^A", "Match_TMLE", "coverage_boot")
mt_C_cov_if_se   <- .bv_cov_mcse(boot_var, "Very Good", "Match_TMLE", "coverage_if")
iptw_B_cov_if_se <- .bv_cov_mcse(boot_var, "Marginal", "IPTW", "coverage_if")
cat(sprintf("[M2] mt_A_cov_if=%.2f (MCSE %.3f) boot=%.2f (MCSE %.3f)  mt_C_if_se=%.3f  iptw_B_if_se=%.3f\n",
            mt_A_cov_if, mt_A_cov_if_se, mt_A_cov_boot, mt_A_cov_boot_se,
            mt_C_cov_if_se, iptw_B_cov_if_se))

# bootstrap table disp build (M2)
cov_if_mcse   <- sqrt(boot_var$coverage_if   * (1 - boot_var$coverage_if)   / boot_var$n_reps)
cov_boot_mcse <- sqrt(boot_var$coverage_boot * (1 - boot_var$coverage_boot) / boot_var$n_reps)
disp <- data.frame(
  Scenario = boot_var$scenario, Method = boot_var$method,
  `Coverage (IF)`   = sprintf("%.2f (%.2f)", boot_var$coverage_if, cov_if_mcse),
  `Coverage (boot)` = sprintf("%.2f (%.2f)", boot_var$coverage_boot, cov_boot_mcse),
  check.names = FALSE)
cat("[M2] bootstrap table coverage cells:\n"); print(disp)

# --- divergence chunks (M1) ---
div_path <- file.path("..", "sandbox", "candidate_divergence", "results", "candidate_divergence_full.rds")
div_fig  <- file.path("..", "sandbox", "candidate_divergence", "figures", "degradation_gradient.png")
div_tab  <- file.path("..", "sandbox", "candidate_divergence", "results", "selection_table.csv")
stopifnot(file.exists(div_path), file.exists(div_fig), file.exists(div_tab))
div   <- readRDS(div_path)
div_b <- setNames(div$agg$baseline_rmse, div$agg$candidate)
div_w <- setNames(div$agg$worst_rmse,    div$agg$candidate)
cat(sprintf("[M1] baseline robust=%.4f middle=%.4f aggr=%.4f\n",
            div_b["robust"], div_b["middle"], div_b["aggressive"]))
cat(sprintf("[M1] worst aggr=%.4f robust=%.4f middle=%.4f  sep=%.4f mcse=%.4f win_min=%s win_mm=%s\n",
            div_w["aggressive"], div_w["robust"], div_w["middle"],
            div$separation, div$combined_mcse, div$win_min_final, div$win_mm_final))
sel <- read.csv(div_tab, check.names = FALSE)
cat("[M1] selection table:\n"); print(sel)
cat("\nALL CHUNKS VALID\n")
