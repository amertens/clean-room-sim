SIM <- "C:/Users/andre/OneDrive/Documents/clean-room-sim/simulation"
.libPaths(c(file.path(SIM, "lib"), .libPaths()))
options(sim.dir = SIM, width = 200)
for (f in c("dgp.R", "map.R", "key.R", "run.R", "summarize.R")) source(file.path(SIM, "R", f))
cal <- readRDS(file.path(SIM, "pilot", "calibration.rds"))
scen <- cal$scenarios
# Answer-key sanity check on 40 datasets (idx 1 to 40), not cached.
for (cfg in list(c("moderate_inside_effect", 500), c("severe_outside_effect", 2000),
                 c("good_inside_effect", 500))) {
  s <- scen[scen$id == cfg[1], ]
  n <- as.integer(cfg[2])
  r <- do.call(rbind, lapply(1:40, function(i) run_key(s, n, i)))
  ak <- answer_key(r, truth_table(s))
  cat("\n", cfg[1], "n =", n, "\n")
  print(format(ak$candidates[c("estimand", "candidate", "M", "failed", "truth", "bias",
                               "bias_mcse", "coverage", "halfwidth", "pass")], digits = 3),
        row.names = FALSE)
}
# Where does one map repetition spend its time? (n = 500)
s <- scen[scen$id == "moderate_inside_effect", ]
d <- draw_data(s, 500, 1)
plan <- make_plan(plan_seed(s, 500, 1), 5L)
lock <- cleanTMLE::create_analysis_lock(d[c("id", "a", COVARIATES)], "a", COVARIATES, plan)
des <- cleanTMLE::assess_design(lock)
X <- cleanTMLE:::.design_matrix(lock$data, COVARIATES)
surf <- cleanTMLE:::.surfaces(plan, X, des$g_world, plan$seed)
pf <- tempfile()
Rprof(pf, interval = 0.01)
t0 <- Sys.time()
for (r in 1:5) invisible(cleanTMLE:::.one_rep(r, X, des$g_world, surf, plan))
Rprof(NULL)
cat("\nseconds per repetition:", as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5, "\n")
sp <- summaryRprof(pf)$by.total
print(head(sp[grep("tmle|glm|target|fit_Q|fit_g|estimate_row|implaus|bound|ato", rownames(sp)), ], 15))
