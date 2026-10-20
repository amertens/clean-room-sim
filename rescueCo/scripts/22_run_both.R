# Runs the design stage for both runs in sequence (one R job): first four
# months, then the final cohort. 100 repetitions (extended to 300 where a cell
# could change the decision), 6 parallel workers. Spec section 6, call 12.
rs <- file.path(R.home("bin"), "Rscript")
for (run in c("first4m", "final")) {
  st <- system2(rs, c("rescueCo/scripts/22_design_stage.R", run, "run", "100", "6"))
  if (st != 0) stop("Run ", run, " failed with status ", st)
}
