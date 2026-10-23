# Runs the design stage for both runs in sequence (one R job): the final
# cohort first (its first lock is the expected confirmatory comparison), then
# the first four months. Each lock is checkpointed, so a restart skips finished
# locks. 100 repetitions (extended to 300 where a cell
# could change the decision), 4 parallel workers (6 until 22 October, when a
# worker ran out of memory on the larger comparisons; the worker count does not
# change any result, since each repetition has its own random stream). Spec
# section 6, call 12.
rs <- file.path(R.home("bin"), "Rscript")
for (run in c("final", "first4m")) {
  st <- system2(rs, c("rescueCo/scripts/22_design_stage.R", run, "run", "100", "4"))
  if (st != 0) stop("Run ", run, " failed with status ", st)
}
