# Pilot for the estimability simulation (BUILD-BRIEF.md, "Pilot script", and
# PLAN.md "Revisions after pilot steps 1 and 2").
# Run from the clean-room-sim folder:
#   Rscript simulation/pilot.R 1          # calibration and DGP diagnostics
#   Rscript simulation/pilot.R 2          # smoke runs (as run before the revisions: commit 82ad787)
#   Rscript simulation/pilot.R par        # parallel stability test (12 short map runs)
#   Rscript simulation/pilot.R 3 <W>      # contrast pilot at n = 500 on W workers
#   Rscript simulation/pilot.R 3a         # contrast-pilot analysis (pilot_check)

SIM <- "C:/Users/andre/OneDrive/Documents/clean-room-sim/simulation"
.libPaths(c(file.path(SIM, "lib"), .libPaths()))
options(sim.dir = SIM, width = 220)
for (f in c("dgp.R", "map.R", "key.R", "run.R", "summarize.R", "pilot_check.R"))
  source(file.path(SIM, "R", f))
stopifnot(utils::packageVersion("cleanTMLE") == "0.4.1",
          dirname(find.package("cleanTMLE")) == normalizePath(file.path(SIM, "lib"), "/"))

args <- commandArgs(trailingOnly = TRUE)
step <- args[1]
PILOT <- file.path(SIM, "pilot")
dir.create(PILOT, showWarnings = FALSE)
secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

if (step == "1") {
  t0 <- Sys.time()
  Xp <- pop_covariates()
  grid <- seq(0.25, 4, by = 0.25)
  cal <- do.call(rbind, lapply(c("inside", "near", "outside"), calibrate,
                               gamma_grid = grid, Xp = Xp))
  cat("\nCalibration (share of population with true pi outside 0.05 to 0.95; b0 per model):\n")
  wide <- reshape(cal[c("model", "gamma", "share_outside", "b0")], idvar = c("gamma"),
                  timevar = "model", direction = "wide")
  wide <- wide[, c("gamma", "share_outside.inside", "b0.inside", "b0.near", "b0.outside")]
  names(wide)[2] <- "share_outside"
  print(format(wide, digits = 4), row.names = FALSE)
  levels <- overlap_levels(OVERLAP_GAMMAS, Xp)
  cat("\nOverlap levels:\n"); print(levels, row.names = FALSE)
  scen <- scenario_table(levels, Xp)
  diag <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) dgp_diagnostics(scen[i, ], Xp)))
  cat("\nScenarios and DGP diagnostics (Table S2):\n")
  print(format(diag[, -1], digits = 3), row.names = FALSE)
  saveRDS(list(calibration = cal, levels = levels, scenarios = scen, diagnostics = diag,
               beta = BETA, pop_n = POP_N, pop_seed = POP_SEED,
               cleanTMLE = as.character(utils::packageVersion("cleanTMLE")),
               created = format(Sys.time(), tz = "UTC", usetz = TRUE)),
          file.path(PILOT, "calibration.rds"))
  utils::write.csv(cal, file.path(PILOT, "calibration_grid.csv"), row.names = FALSE)
  utils::write.csv(diag, file.path(PILOT, "dgp_diagnostics.csv"), row.names = FALSE)
  cat(sprintf("\nStep 1 took %.1f s.\n", secs(t0)))
}

cal <- if (file.exists(file.path(PILOT, "calibration.rds"))) readRDS(file.path(PILOT, "calibration.rds"))
levels <- cal$levels
scen <- cal$scenarios

if (step == "2") {
  # Two datasets per map variant and key, at one level and n each.
  smoke <- data.frame(level = c("out16", "out41"), n = c(500L, 2000L), stringsAsFactors = FALSE)
  for (j in seq_len(nrow(smoke))) {
    lv <- levels[levels$id == smoke$level[j], ]
    for (st in c("map_generate", "map_sample")) {
      t0 <- Sys.time()
      print(run_stage(st, lv, smoke$n[j], 1:2, reps = 200L))
      cat(sprintf("%s %s n=%d: %.1f s\n", st, lv$id, smoke$n[j], secs(t0)))
    }
    print(run_stage("key", scen[scen$overlap == lv$id, ], smoke$n[j], 1:2))
  }
}

if (step == "par") {
  # Stability and speed-up of parallel workers: 12 short map runs (reps = 20)
  # at n = 500, first sequentially, then on 6 workers, each into its own
  # scratch cache so neither reads the other's results; then 2 runs at
  # n = 2,000 sequentially for the scaling with n.
  W <- if (length(args) > 1) as.integer(args[2]) else 6L
  lv <- levels[levels$id %in% c("out04", "out16", "out26", "out41"), ]
  base <- file.path(SIM, "cache", "_partest", format(Sys.time(), "%Y%m%d_%H%M%S"))
  t0 <- Sys.time()
  seq_rows <- run_stage("map_generate", lv, 500L, 9001:9003, reps = 20L,
                        cache_dir = file.path(base, "seq"))
  t_seq <- secs(t0)
  cat(sprintf("Sequential: 12 runs in %.1f s\n", t_seq))
  t0 <- Sys.time()
  par_rows <- run_stage("map_generate", lv, 500L, 9001:9003, reps = 20L, workers = W,
                        cache_dir = file.path(base, sprintf("par%d", W)))
  t_par <- secs(t0)
  cat(sprintf("Parallel (%d workers): 12 runs in %.1f s; speed-up %.2f\n", W, t_par,
              t_seq / t_par))
  same <- isTRUE(all.equal(seq_rows[setdiff(names(seq_rows), c("sec_assess", "sec_total"))],
                           par_rows[setdiff(names(par_rows), c("sec_assess", "sec_total"))],
                           check.attributes = FALSE))
  cat("Parallel rows equal sequential rows:", same, "\n")
  t0 <- Sys.time()
  big <- run_stage("map_generate", lv[lv$id == "out16", ], 2000L, 9001:9002, reps = 20L,
                   cache_dir = file.path(base, "seq"))
  t_big <- secs(t0)
  per_rep <- c(n500 = t_seq / (12 * 20), n2000 = t_big / (2 * 20))
  cat(sprintf("Seconds per repetition (generate, one truncation): n=500 %.2f, n=2000 %.2f\n",
              per_rep[1], per_rep[2]))
  saveRDS(list(workers = W, t_seq = t_seq, t_par = t_par, speedup = t_seq / t_par,
               identical = same, per_rep = per_rep, seq_rows = seq_rows, par_rows = par_rows),
          file.path(PILOT, sprintf("parallel_test_w%d.rds", W)))
}

if (step == "3") {
  # Contrast pilot: n = 500; every overlap level; R = 100 map runs per level for
  # both variants; M = 200 key datasets per scenario (the map datasets are the
  # first 100 of them).
  W <- if (length(args) > 1) as.integer(args[2]) else 1L
  for (st in c("key", "map_generate", "map_sample")) {
    t0 <- Sys.time()
    units <- if (st == "key") scen else levels
    idx <- if (st == "key") 1:200 else 1:100
    invisible(run_stage(st, units, 500L, idx, reps = 200L, workers = W))
    cat(sprintf("%s: %d cells in %.1f s (%d workers)\n", st, nrow(units) * length(idx),
                secs(t0), W))
  }
}

if (step == "3a") {
  key <- collect("key"); key <- key[key$n == 500 & key$idx <= 200, ]
  maps <- rbind(collect("map_generate"), collect("map_sample"))
  maps <- maps[maps$n == 500 & maps$idx <= 100 & maps$reps == 200, ]
  truths <- truth_table(scen)
  ak <- answer_key(key, truths)
  rt <- add_diagnostic_calls(run_table(maps, ak$estimable, scen))
  agree <- agreement_table(rt)
  con <- table1_contrasts(rt, agree, k = 5)
  c3 <- claim3_rows(key, maps, truths, scen)
  c3s <- claim3_summary(c3, k = 5)
  c4 <- claim4_rows(key, maps, scen)
  c4s <- claim4_summary(c4)
  saveRDS(list(answer_key = ak, run_table = rt, thresholds = attr(rt, "thresholds"),
               agreement = agree, contrasts = con, claim3 = c3, claim3_summary = c3s,
               claim4 = c4, claim4_summary = c4s),
          file.path(PILOT, "step3.rds"))
  utils::write.csv(ak$candidates, file.path(PILOT, "step3_answer_key.csv"), row.names = FALSE)
  utils::write.csv(agree, file.path(PILOT, "step3_agreement.csv"), row.names = FALSE)
  utils::write.csv(con, file.path(PILOT, "step3_contrasts.csv"), row.names = FALSE)
  utils::write.csv(c3s, file.path(PILOT, "step3_claim3.csv"), row.names = FALSE)
  utils::write.csv(c4s, file.path(PILOT, "step3_claim4.csv"), row.names = FALSE)
  cat("\nAnswer key (estimable) by scenario and estimand:\n")
  print(stats::xtabs(estimable ~ scenario + estimand, ak$estimable))
  cat("\nMap runtimes (s):\n")
  print(stats::aggregate(sec_total ~ variant + overlap, unique(maps[c("variant", "overlap", "idx", "sec_total")]), mean))
  cat("\nThresholds:\n"); print(attr(rt, "thresholds"))
  cat("\nAgreement:\n"); print(format(agree, digits = 3), row.names = FALSE)
  cat("\nTable 1 contrasts:\n")
  print(format(con[c("estimand", "relation", "measure", "comparison", "other", "a_value",
                     "b_value", "gap", "sd_d", "pilot_z", "n_sim_safe")], digits = 3),
        row.names = FALSE)
  cat("\nClaim 3:\n"); print(format(c3s[c("overlap", "effect", "measure", "prespecified",
                                           "smallest_p", "datasets", "gap", "sd_d", "pilot_z",
                                           "n_sim_safe")], digits = 3), row.names = FALSE)
  cat("\nClaim 4:\n"); print(format(c4s, digits = 3), row.names = FALSE)
}
