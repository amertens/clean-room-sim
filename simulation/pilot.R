# Pilot for the estimability simulation (BUILD-BRIEF.md, "Pilot script").
# Run from the clean-room-sim folder:
#   Rscript simulation/pilot.R 1   # calibration and DGP diagnostics
#   Rscript simulation/pilot.R 2   # smoke runs and runtime projection
# Step 3 (the contrast pilot) is not run here until its cost is agreed.

SIM <- "C:/Users/andre/OneDrive/Documents/clean-room-sim/simulation"
.libPaths(c(file.path(SIM, "lib"), .libPaths()))
options(sim.dir = SIM, width = 200)
for (f in c("dgp.R", "map.R", "key.R", "run.R", "summarize.R")) source(file.path(SIM, "R", f))
stopifnot(utils::packageVersion("cleanTMLE") == "0.4.1",
          dirname(find.package("cleanTMLE")) == normalizePath(file.path(SIM, "lib"), "/"))

step <- as.integer(commandArgs(trailingOnly = TRUE)[1])
PILOT <- file.path(SIM, "pilot")
dir.create(PILOT, showWarnings = FALSE)

if (step == 1L) {
  t0 <- Sys.time()
  Xp <- pop_covariates()
  grid <- seq(0.25, 4, by = 0.25)
  cal <- do.call(rbind, lapply(c("inside", "near", "outside"), calibrate,
                               gamma_grid = grid, Xp = Xp))
  gammas <- pick_gamma(cal)
  cat("\nCalibration (share of population with true pi outside 0.05 to 0.95; b0 per model):\n")
  wide <- reshape(cal[c("model", "gamma", "share_outside", "b0")], idvar = c("gamma"),
                  timevar = "model", direction = "wide")
  wide <- wide[, c("gamma", "share_outside.inside", "b0.inside", "b0.near", "b0.outside")]
  names(wide)[2] <- "share_outside"
  print(format(wide, digits = 4), row.names = FALSE)
  cat("\nChosen gamma:\n"); print(gammas)
  if (anyNA(gammas)) stop("An overlap level cannot be reached on this grid.")
  scen <- scenario_table(gammas, Xp)
  diag <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) dgp_diagnostics(scen[i, ], Xp)))
  cat("\nScenarios and DGP diagnostics (Table S2):\n")
  print(format(diag[, -1], digits = 3), row.names = FALSE)
  saveRDS(list(calibration = cal, gammas = gammas, scenarios = scen, diagnostics = diag,
               beta = BETA, pop_n = POP_N, pop_seed = POP_SEED,
               cleanTMLE = as.character(utils::packageVersion("cleanTMLE")),
               created = format(Sys.time(), tz = "UTC", usetz = TRUE)),
          file.path(PILOT, "calibration.rds"))
  utils::write.csv(cal, file.path(PILOT, "calibration_grid.csv"), row.names = FALSE)
  utils::write.csv(diag, file.path(PILOT, "dgp_diagnostics.csv"), row.names = FALSE)
  cat(sprintf("\nStep 1 took %.1f s.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

if (step == 2L) {
  cal <- readRDS(file.path(PILOT, "calibration.rds"))
  scen <- cal$scenarios
  smoke <- data.frame(id = c("moderate_inside_effect", "severe_outside_effect"),
                      n = c(500L, 2000L), stringsAsFactors = FALSE)
  REPS <- 200L
  out <- list()
  for (j in seq_len(nrow(smoke))) {
    s <- scen[scen$id == smoke$id[j], ]
    for (st in STAGES) {
      t0 <- Sys.time()
      r <- run_stage(st, s, smoke$n[j], 1:2, reps = REPS)
      cat(sprintf("%s %s n=%d: %.1f s for 2 datasets\n", st, s$id, smoke$n[j],
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      out[[st]] <- rbind(out[[st]], r)
    }
  }
  # A cache hit returns the stored rows without recomputing.
  t0 <- Sys.time()
  again <- run_stage("key", scen[scen$id == smoke$id[1], ], smoke$n[1], 1L)
  hit_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  same <- isTRUE(all.equal(again, out$key[out$key$scenario == smoke$id[1] &
                                            out$key$n == smoke$n[1] & out$key$idx == 1L, ],
                           check.attributes = FALSE))
  cat(sprintf("Cache hit: %.2f s, identical rows: %s\n", hit_sec, same))

  maps <- rbind(out$map_generate, out$map_sample)
  cat("\nMap verdicts and diagnostics:\n")
  print(format(maps[, c("scenario", "n", "idx", "variant", "estimand", "status", "selected",
                        "ci_halfwidth", "unresolved", "share_outside", "max_weight",
                        "ess_frac", "pref_share", "sec_total")], digits = 3), row.names = FALSE)
  cat("\nAnswer-key fits:\n")
  print(format(out$key[, c("scenario", "n", "idx", "estimand", "candidate", "estimate", "se",
                           "ci_lower", "ci_upper", "failed", "p_value", "min_p", "sec_total")],
               digits = 3), row.names = FALSE)

  # Summaries run end to end (two datasets: the numbers mean nothing yet).
  truths <- truth_table(scen[scen$id %in% smoke$id, ])
  ak <- answer_key(out$key, truths)
  ag <- agreement(maps, ak$estimable, scen)
  cat("\nAnswer key on 2 datasets (code check only):\n"); print(ak$estimable)
  cat("\nAgreement on 2 datasets (code check only):\n"); print(ag)

  # Runtime per call (one dataset), by n, and projected hours.
  per_call <- function(rows, stage) {
    u <- unique(rows[c("n", "idx", "sec_total")])
    stats::aggregate(sec_total ~ n, u, mean)
  }
  tg <- per_call(out$map_generate); ts <- per_call(out$map_sample); tk <- per_call(out$key)
  tim <- data.frame(n = tg$n, map_generate = tg$sec_total, map_sample = ts$sec_total,
                    key = tk$sec_total)
  cat("\nSeconds per call (mean of 2 datasets):\n"); print(tim)
  n_scen <- nrow(scen)
  R <- 100; M <- 500
  per_n_map <- n_scen * R * (tim$map_generate + tim$map_sample) / 3600
  per_n_key <- n_scen * M * tim$key / 3600
  proj <- data.frame(n = tim$n, map_hours = per_n_map, key_hours = per_n_key,
                     total_hours = per_n_map + per_n_key)
  cat(sprintf("\nProjected sequential hours, full grid (%d scenarios x 2 n), R = %d map runs per variant, M = %d key datasets:\n",
              n_scen, R, M))
  print(proj); cat(sprintf("Total: %.1f h\n", sum(proj$total_hours)))
  # Map runs are outcome-free: if X and A were seeded by (overlap, n, idx), one
  # map run would serve all six model x effect scenarios of an overlap level.
  shared <- 3 * R * sum(tim$map_generate + tim$map_sample) / 3600
  cat(sprintf("Map hours if map datasets were shared across models and effects: %.1f h\n", shared))
  contrast <- data.frame(n = tim$n, map100_hours = 100 * (tim$map_generate + tim$map_sample) / 3600,
                         key200_hours = 200 * tim$key / 3600)
  cat("\nContrast pilot: 100 map runs (both variants) and 200 key datasets, per scenario:\n")
  print(contrast)
  saveRDS(list(rows = out, maps = maps, timing = tim, projection = proj,
               shared_map_hours = shared, contrast = contrast, cache_hit_sec = hit_sec,
               cache_identical = same, answer_key = ak, agreement = ag),
          file.path(PILOT, "smoke.rds"))
}
