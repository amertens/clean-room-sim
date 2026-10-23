# ============================================================
# Export the Rescue.Co re-execution's design-stage numbers for the registry
# ============================================================
# Reads the outcome-free design stage (design data, assessments, simulation
# and negative-control objects) and writes aggregate tables to
# writing/results-export/, which writing/registry.json points to. No patient
# identifier and no primary outcome is read or written. Rerun from the
# project root after each design-stage run:
#   Rscript writing/export_results.R
# ============================================================

DESIGN <- "rescueCo/data/design"
RES <- "rescueCo/results/reexec_050"
OUT <- "writing/results-export"
dir.create(OUT, showWarnings = FALSE)
RUNS <- c("final", "first4m")

# The population gate (rescueCo/scripts/23_registration_bundle.R), applied here
# to every estimand below the ATE rather than only to the primary estimand.
GATE_KEEP <- 0.80
GATE_SHIFT <- 0.25
GATE_COVS <- c("age", "eiss_ais", "head_injury_severity",
               "intent_Intentional..assault.homicide..", "inj_place_Private.House..Home")
ATO_ALLOWED <- c("L3", "S1", "S4")

gate_rows <- list()
for (run in RUNS) {
  dd <- file.path(RES, run, "design")
  locks <- sub("_populations\\.csv$", "", list.files(dd, "_populations\\.csv$"))
  for (lk in locks) {
    pop <- read.csv(file.path(dd, paste0(lk, "_populations.csv")))
    cnt <- read.csv(file.path(dd, paste0(lk, "_trimmed_counts.csv")))
    share <- function(g) cnt$kept[cnt$group == g] / (cnt$kept[cnt$group == g] + cnt$removed[cnt$group == g])
    named <- pop$covariate %in% GATE_COVS
    for (e in c("trimmed_ATE", "ATT", "ATO")) {
      sh <- abs(pop[[paste0(e, "_shift")]][named])
      kt <- if (e == "trimmed_ATE") share("treated") else 1
      kc <- if (e == "trimmed_ATE") share("control") else 1
      # Shifts are not shown when trimming removes fewer than five patients; the
      # gate then counts the shift as zero, as script 23 does.
      ms <- if (all(is.na(sh))) NA_real_ else max(sh, 0, na.rm = TRUE)
      why <- c(if (kt < GATE_KEEP || kc < GATE_KEEP) "keeps under 80% of an arm",
               if (isTRUE(ms > GATE_SHIFT)) "shift over 0.25 SD on a named covariate",
               if (e == "ATO" && !lk %in% ATO_ALLOWED) "overlap-weighted effect not acceptable here")
      gate_rows[[length(gate_rows) + 1]] <- data.frame(
        run = run, lock = lk, estimand = e, named_covariates_found = sum(named),
        kept_treated = kt, kept_comparator = kc, max_abs_shift = ms,
        passes = !length(why), reason = if (length(why)) paste(why, collapse = "; ") else "passes")
    }
  }
}
gate <- do.call(rbind, gate_rows)
write.csv(gate, file.path(OUT, "population_gate.csv"), row.names = FALSE)

# Numerators for the percentages: patients outside the band, and patients
# lost to follow-up by arm (shares in the source are exact to 4 decimals, which
# recovers each count exactly for arms under 10,000 patients).
ov_rows <- list()
for (run in RUNS) {
  dd <- file.path(RES, run, "design")
  for (f in list.files(dd, "_overlap\\.csv$")) {
    ov <- read.csv(file.path(dd, f))
    k <- round(ov$n * ov$share_outside_band)
    ov_rows[[length(ov_rows) + 1]] <- data.frame(run = run, lock = sub("_overlap\\.csv$", "", f),
      outside_treated = k[ov$arm == "treated"], outside_comparator = k[ov$arm == "control"], outside = sum(k))
  }
}
write.csv(do.call(rbind, ov_rows), file.path(OUT, "overlap_counts.csv"), row.names = FALSE)
loss_rows <- list()
for (f in list.files(file.path(RES, "final"), "^loss_by_arm_.*\\.csv$", full.names = TRUE)) {
  l <- read.csv(f)
  stopifnot(all(l$patients < 10000))
  loss_rows[[length(loss_rows) + 1]] <- data.frame(run = "final", lock = l$lock, arm = l$arm, patients = l$patients,
    six_month_not_recorded = round(l$patients * l$six_month_not_recorded),
    disposition_not_recorded = round(l$patients * l$disposition_not_recorded))
}
write.csv(do.call(rbind, loss_rows), file.path(OUT, "loss_counts.csv"), row.names = FALSE)

# Cohort counts from the design data (registry-wide and per run; no outcome).
flow_rows <- list()
for (run in RUNS) {
  dd <- readRDS(file.path(DESIGN, paste0("design_data_", run, ".rds")))
  D <- dd$data
  rc <- D$arm %in% "rescueco"
  add <- function(item, n) flow_rows[[length(flow_rows) + 1]] <<- data.frame(run = run, item = item, n = n)
  add("patients", nrow(D))
  for (a in c("rescueco", "other_operator", "other_means")) add(paste0("arm_", a), sum(D$arm %in% a))
  add("arm_unknown", sum(is.na(D$arm)))
  add("ambulance_operator_unknown", sum(D$operator_unknown == 1, na.rm = TRUE))
  add("rescueco_no_dispatch_record", sum(rc & D$crossref_tier %in% 5L))
  add("rescueco_probable_match", sum(rc & D$crossref_tier %in% 3L))
  dir_l3 <- D$arm %in% c("rescueco", "other_operator") & D$direct_arrival == 1 & !(D$crossref_tier %in% 5L)
  add("rescueco_probable_match_direct", sum(rc & dir_l3 & D$crossref_tier %in% 3L))
  if ("six_month_recorded" %in% names(D))
    add("rescueco_probable_match_direct_followed", sum(rc & dir_l3 & D$crossref_tier %in% 3L & D$six_month_recorded == 1))
  add("transfer_or_prior_care", sum(D$direct_arrival == 0, na.rm = TRUE))
  add("walked_in", sum(D$walked_in == 1, na.rm = TRUE))
  add("covariate_columns", length(dd$covariates))
  add("arrival_state_columns", length(dd$arrival_state))
}
write.csv(do.call(rbind, flow_rows), file.path(OUT, "cohort_counts.csv"), row.names = FALSE)

# The outcome-blind simulation, one row per lock, estimand and candidate,
# summarized over the simulated outcome models (worst case across them).
sim_rows <- list(); verdict_rows <- list(); ev_rows <- list(); nc_rows <- list()
for (run in RUNS) {
  ld <- file.path(DESIGN, "locks", run)
  for (f in list.files(ld, "_simulation\\.rds$")) {
    lk <- sub("_simulation\\.rds$", "", f)
    s <- readRDS(file.path(ld, f))
    m <- s$metrics
    for (key in unique(paste(m$estimand, m$candidate))) {
      x <- m[paste(m$estimand, m$candidate) == key, ]
      sim_rows[[length(sim_rows) + 1]] <- data.frame(
        run = run, lock = lk, estimand = x$estimand[1], candidate = x$candidate[1],
        library = x$library[1], truncation = x$truncation[1], models = nrow(x),
        reps_min = min(x$reps), reps_max = max(x$reps), failed_max = max(x$fail_rate),
        max_abs_bias = max(abs(x$bias)), min_coverage = min(x$coverage), max_coverage = max(x$coverage),
        min_se_mean = min(x$se_mean), max_se_mean = max(x$se_mean),
        min_emp_se = min(x$emp_se), max_emp_se = max(x$emp_se),
        min_halfwidth = min(x$ci_halfwidth), max_halfwidth = max(x$ci_halfwidth),
        max_rmse = max(x$rmse), max_bias_mcse = max(x$bias_mcse),
        min_coverage_lower = min(x$coverage_lower), max_coverage_upper = max(x$coverage_upper),
        unresolved = sum(x$unresolved), pass_bias_coverage = all(x$pass_bc), pass_halfwidth = all(x$pass_hw),
        pass = all(x$pass), any_unresolved = any(x$unresolved))
    }
    v <- s$verdict
    verdict_rows[[length(verdict_rows) + 1]] <- data.frame(run = run, lock = lk, v, row.names = NULL)
    if (is.data.frame(s$expected_events))
      ev_rows[[length(ev_rows) + 1]] <- data.frame(run = run, lock = lk, s$expected_events, row.names = NULL)
  }
  for (f in list.files(ld, "_nc\\.rds$")) {
    lk <- sub("_nc\\.rds$", "", f)
    n <- readRDS(file.path(ld, f))
    if (is.data.frame(n$table))
      nc_rows[[length(nc_rows) + 1]] <- data.frame(run = run, lock = lk, n$table, row.names = NULL)
    verdict_nc <- n$verdict
    nc_rows_v <- data.frame(run = run, lock = lk, estimand = names(verdict_nc), verdict = unname(verdict_nc),
                            null_band_lower = n$criteria$null_band[1], null_band_upper = n$criteria$null_band[2])
    write.csv(nc_rows_v, file.path(OUT, paste0("nc_verdict_", run, "_", lk, ".csv")), row.names = FALSE)
  }
}
if (length(sim_rows)) {
  sim <- do.call(rbind, sim_rows)
  write.csv(sim, file.path(OUT, "simulation_summary.csv"), row.names = FALSE)
  # Ranges by learner library (over truncation levels), by estimand and over all estimands.
  span <- function(x, by) do.call(rbind, lapply(split(x, x[by], drop = TRUE), function(y) data.frame(
    y[1, by, drop = FALSE], worst_bias_min = min(y$max_abs_bias), worst_bias_max = max(y$max_abs_bias),
    coverage_min = min(y$min_coverage), coverage_max = max(y$max_coverage),
    se_mean_min = min(y$min_se_mean), se_mean_max = max(y$max_se_mean),
    emp_se_min = min(y$min_emp_se), emp_se_max = max(y$max_emp_se),
    halfwidth_min = min(y$min_halfwidth), halfwidth_max = max(y$max_halfwidth),
    bias_mcse_max = max(y$max_bias_mcse), coverage_lower_min = min(y$min_coverage_lower),
    coverage_upper_max = max(y$max_coverage_upper), cells = sum(y$models), unresolved_cells = sum(y$unresolved),
    pass_bias_coverage = all(y$pass_bias_coverage), pass_halfwidth = all(y$pass_halfwidth), row.names = NULL)))
  write.csv(span(sim, c("run", "lock", "library", "estimand")),
            file.path(OUT, "simulation_by_library_estimand.csv"), row.names = FALSE)
  write.csv(span(sim, c("run", "lock", "library")), file.path(OUT, "simulation_by_library.csv"), row.names = FALSE)
}
if (length(verdict_rows)) write.csv(do.call(rbind, verdict_rows), file.path(OUT, "simulation_verdicts.csv"), row.names = FALSE)
if (length(ev_rows)) {
  ev <- do.call(rbind, ev_rows)
  write.csv(ev, file.path(OUT, "expected_events.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(split(ev, ev[c("run", "lock")], drop = TRUE), function(y) data.frame(
    run = y$run[1], lock = y$lock[1], treated_min = min(y$treated), treated_max = max(y$treated),
    control_min = min(y$control), control_max = max(y$control)))),
    file.path(OUT, "expected_events_range.csv"), row.names = FALSE)
}
if (length(nc_rows)) {
  nc <- do.call(rbind, nc_rows)
  write.csv(nc, file.path(OUT, "negative_controls.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(split(nc, nc[c("run", "lock", "rung", "estimand")], drop = TRUE), function(y) data.frame(
    run = y$run[1], lock = y$lock[1], rung = y$rung[1], estimand = y$estimand[1], controls = nrow(y),
    estimated = sum(y$status == "estimated"), inestimable = sum(y$status != "estimated"),
    inside = sum(y$reading %in% "inside"), outside = sum(y$reading %in% "outside"),
    inconclusive = sum(y$reading %in% "inconclusive")))),
    file.path(OUT, "negative_controls_summary.csv"), row.names = FALSE)
}
cat("Wrote", paste(list.files(OUT), collapse = ", "), "\n")
