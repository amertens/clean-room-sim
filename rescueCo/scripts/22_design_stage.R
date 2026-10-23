# ============================================================
# 22. Design stage on cleanTMLE 0.5.0: plans, locks, dossiers
# ============================================================
# For each lock of one run: the plan, the variable dictionary, the lock, the
# design assessment, the outcome-free simulation, the negative-control ladder
# and the rendered dossier. Checkpointed per lock: a lock whose dossier exists
# is skipped. Reads design data only; no outcome exists in this process.
#
# Usage:
#   Rscript rescueCo/scripts/22_design_stage.R <run> pilot      # time 2 reps on L3 and L1
#   Rscript rescueCo/scripts/22_design_stage.R <run> run <reps> <workers>  # every lock of the run
#   Rscript rescueCo/scripts/22_design_stage.R <run> run <reps> <workers> <lock>  # one lock
#   Rscript rescueCo/scripts/22_design_stage.R <run> verify <reps>   # rebuild finished locks, compare hashes
#   run = "first4m" or "final"
#
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/reexec_guards.R")   # rx_log(): rescueCo/logs/pipeline.log is not opened (ledger, section 5)
# cleanTMLE 0.5.0 from the project library (rescueCo/lib, built from the
# cleanTMLE submodule), so the installed 0.3.0 other projects use is untouched.
.libPaths(c("rescueCo/lib", .libPaths()))
suppressMessages(library(cleanTMLE))
stopifnot(utils::packageVersion("cleanTMLE") >= "0.5.0")
Sys.setenv(QUARTO_PATH = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe")
cfg <- load_cr_config()
args <- commandArgs(trailingOnly = TRUE)
RUN <- args[1]; MODE <- args[2]
REPS <- if (length(args) >= 3) as.integer(args[3]) else 100L
WORKERS <- if (length(args) >= 4) as.integer(args[4]) else 1L
ONLY <- if (length(args) >= 5) args[5] else NULL   # one lock only (the supervisor sizes workers per lock)
stopifnot(RUN %in% c("first4m", "final"), MODE %in% c("pilot", "run", "check", "verify"))
rx_log(paste("=== 22 design stage, run", RUN, "mode", MODE,
             "cleanTMLE", as.character(utils::packageVersion("cleanTMLE")), "==="))

dd <- readRDS(file.path(cfg$data$design_dir, paste0("design_data_", RUN, ".rds")))
D <- dd$data
OUT <- file.path(cfg$data$design_dir, "locks", RUN)       # patient-level objects
RES <- file.path(cfg$paths$results, "reexec_050", RUN)    # aggregate tables and dossiers
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Locks (spec section 4)
# ------------------------------------------------------------
both_ops <- c("rescueco", "other_operator")
# Compute decisions of 2026-10-23 (spec section 10), made before any of these
# locks was built and without any outcome: `outcome_sl = FALSE` drops the super
# learner outcome candidate (the super learner still fits the propensity score),
# and `extend = FALSE` sets max_reps to reps. Both apply to the comparisons that
# cannot carry the confirmatory claim; L5, the one with moderate overlap, keeps
# both candidates. L3, S1 and S4 carry neither flag, so their plans, and L3's
# lock hashes, are unchanged.
LOCKS <- list(
  L1 = list(label = "Rescue.Co vs all other transport, all arrivals (protocol)",
            treated = "rescueco", control = c("other_operator", "other_means"), direct = FALSE,
            outcome_sl = FALSE, extend = FALSE),
  L2 = list(label = "Rescue.Co vs all other transport, direct arrivals",
            treated = "rescueco", control = c("other_operator", "other_means"), direct = TRUE,
            outcome_sl = FALSE, extend = FALSE),
  L3 = list(label = "Rescue.Co vs other ambulance operators, direct arrivals",
            treated = "rescueco", control = "other_operator", direct = TRUE),
  L4 = list(label = "Rescue.Co vs arrival by other means, direct arrivals",
            treated = "rescueco", control = "other_means", direct = TRUE,
            outcome_sl = FALSE, extend = FALSE),
  L5 = list(label = "Any ambulance vs arrival by other means, direct arrivals",
            treated = both_ops, control = "other_means", direct = TRUE, extend = FALSE),
  S1 = list(label = "As L3, arrival state added as a stand-in",
            treated = "rescueco", control = "other_operator", direct = TRUE, arrival = TRUE),
  S2 = list(label = "As L4, walk-ins excluded",
            treated = "rescueco", control = "other_means", direct = TRUE, no_walkin = TRUE,
            outcome_sl = FALSE, extend = FALSE),
  S3 = list(label = "As L5, walk-ins excluded",
            treated = both_ops, control = "other_means", direct = TRUE, no_walkin = TRUE,
            outcome_sl = FALSE, extend = FALSE),
  S4 = list(label = "As L3, probable dispatch matches (tier 3) excluded",
            treated = "rescueco", control = "other_operator", direct = TRUE, no_tier3 = TRUE))
# Run order: the expected confirmatory comparison first, so results arrive in a
# useful order (compute call 12, spec section 6).
LOCKS <- LOCKS[c("L3", "L4", "L5", "L2", "L1", "S1", "S4", "S2", "S3")]
# The first four months: the confirmatory comparison (done) and the protocol's
# comparison only (spec section 10).
if (RUN == "first4m") LOCKS <- LOCKS[c("L3", "L1")]
COMPLETE_CASE <- RUN == "final"

# ------------------------------------------------------------
# Plan pieces shared by every lock (spec section 5)
# ------------------------------------------------------------
TIME_ZERO <- list(
  event = "transport decision (dispatch, or departure from the scene)",
  measured_as = "time of injury, because arrivals by other means have no dispatch time",
  follow_up_from = "arrival at hospital",
  response = paste("The six-month horizon is counted from arrival; the time from injury",
                   "to arrival is short relative to six months."))
ELIG_ALIVE <- list(
  description = "Arrived alive at one of the three participating hospitals and was enrolled in the registry",
  assessed = "after",
  response = paste("The target-trial estimand is the effect in injured patients for whom a",
                   "transport choice was made; the registry identifies it only if transport does",
                   "not change survival to arrival. A two-sided bias analysis for deaths before",
                   "arrival and a restriction to patients who would almost surely have arrived",
                   "alive are registered before unblinding (D6)."))
ELIG_ARM <- list(
  description = paste("Transport arm known: operator recorded for ambulance arrivals, and Rescue.Co",
                      "transport confirmed against dispatch records"),
  assessed = "after",
  response = paste("This classifies the treatment rather than selecting on what it changes; the",
                   "ambulance arrivals with no recorded operator and the Rescue.Co patients with",
                   "no dispatch record are excluded."))
ELIG_DIRECT <- list(
  description = "Arrived directly from the scene: not an interfacility transfer, and no care sought elsewhere before arrival",
  assessed = "after",
  response = paste("Route, and care sought elsewhere, may follow the transport decision; the",
                   "protocol comparison over all arrivals (L1) is reported without this",
                   "restriction, with the share of each arm it removes."))
ELIG_WALKIN <- list(description = "Did not walk in", assessed = "at",
                    response = NULL)
ELIG_CC <- list(
  description = "Six-month status recorded (complete case)",
  assessed = "after",
  response = paste("Complete case is the primary handling of missing outcomes and the estimator",
                   "the design stage evaluates; the share lost by arm is reported, and a",
                   "longitudinal TMLE is a sensitivity analysis."))
DOMAINS <- c(nc_household_urban = "household and wealth",
             nc_cooking_wood_charcoal = "household and wealth",
             nc_hypertension = "health and access to care",
             nc_diabetes_insulin = "health and access to care",
             nc_hiv_art = "health and access to care")

make_plan <- function(spec, covs, data, reps, max_reps = 3L * reps) {
  el <- list(arrived_alive = ELIG_ALIVE, arm_known = ELIG_ARM)
  if (isTRUE(spec$direct)) el$direct_arrival <- ELIG_DIRECT
  if (isTRUE(spec$no_walkin)) el$no_walk_in <- ELIG_WALKIN
  if (COMPLETE_CASE) el$complete_case <- ELIG_CC
  prev <- colMeans(data[, names(DOMAINS)])
  band <- max(0.02, 0.25 * stats::median(prev))
  ex <- list()
  cook <- grep("^cooking_fuel_", covs, value = TRUE)
  if (length(cook)) ex$nc_cooking_wood_charcoal <- cook
  chron <- grep("^chronic_illness_|^medical_history_2$", covs, value = TRUE)
  if (length(chron)) for (k in c("nc_hypertension", "nc_diabetes_insulin", "nc_hiv_art"))
    ex[[k]] <- chron
  analysis_plan(
    outcome = "severe_disability", estimands = c("ATE", "trimmed_ATE", "ATT", "ATO"),
    time_zero = TIME_ZERO, eligibility = el,
    ps_library = "default",
    candidates = list(library = if (isFALSE(spec$outcome_sl)) list(glm = "glm") else
                        list(glm = "glm", sl = "default"), truncation = c(0.01, 0.05)),
    tolerance = list(bias = 0.0275, coverage = 0.90, ci_halfwidth = 0.055),
    surfaces = list(baseline_risk = 0.22),
    trim_band = c(0.05, 0.95),
    negative_controls = DOMAINS, nc_exclude = if (length(ex)) ex else NULL,
    nc_criteria = list(null_band = c(-band, band), rule = "interval", min_per_domain = 1L),
    roles = list(analyst = "Andrew Mertens", data_steward = "Andrew Mertens (steward scripts 20 and 24)",
                 reviewers = "Andrew Mertens (the demonstration has no outside approver)"),
    K = 5L, V = 5L, reps = reps, max_reps = max(reps, max_reps), seed = 20261020L)
}

lock_rows <- function(spec) {
  keep <- D$arm %in% c(spec$treated, spec$control) & !(D$crossref_tier %in% 5L)
  if (isTRUE(spec$direct)) keep <- keep & D$direct_arrival == 1
  if (isTRUE(spec$no_walkin)) keep <- keep & D$walked_in == 0
  if (isTRUE(spec$no_tier3)) keep <- keep & !(D$crossref_tier %in% 3L)
  keep
}

build_lock <- function(name, spec, reps, max_reps = 3L * reps) {
  if (isFALSE(spec$extend)) max_reps <- reps
  keep <- lock_rows(spec)
  x <- D[keep, , drop = FALSE]
  x$A <- as.integer(x$arm %in% spec$treated)
  if (COMPLETE_CASE) {
    loss <- do.call(rbind, lapply(c(1L, 0L), function(a) {
      s <- x[x$A == a, ]
      data.frame(lock = name, arm = if (a == 1) "treated" else "comparator", patients = nrow(s),
                 disposition_not_recorded = round(mean(s$disposition_recorded == 0), 4),
                 six_month_not_recorded = round(mean(s$six_month_recorded == 0), 4))
    }))
    utils::write.csv(loss, file.path(RES, paste0("loss_by_arm_", name, ".csv")), row.names = FALSE)
    x <- x[x$six_month_recorded == 1, , drop = FALSE]
  }
  covs <- c(dd$covariates, if (isTRUE(spec$arrival)) dd$arrival_state)
  obs <- if (COMPLETE_CASE) c("disposition_recorded", "six_month_recorded") else character()
  x <- x[, c("patient_id", "A", covs, dd$negative_controls, obs)]
  dic <- dd$dictionary[dd$dictionary$column %in% c(covs, dd$negative_controls, obs), ]
  dic$role[dic$column %in% dd$arrival_state] <- "stand_in"
  dic <- rbind(dic,
    data.frame(column = "patient_id", variable = "patient_id", source = "registry",
               follow_up = FALSE, recorded = "after", fixed = "after", role = "id",
               reason = "", sensitivity = ""),
    data.frame(column = "A", variable = "transport", follow_up = FALSE,
               source = "registry transport section, confirmed against dispatch records",
               recorded = "after", fixed = "at", role = "treatment", reason = "", sensitivity = ""))
  plan <- make_plan(spec, covs, x, reps, max_reps)
  create_analysis_lock(x, treatment = "A", covariates = covs, plan = plan,
                       id = "patient_id", dictionary = dic)
}

# ------------------------------------------------------------
# Check: build every lock (the dictionary and plan checks run) and print each
# lock's size and identification flags; no assessment or simulation.
# ------------------------------------------------------------
if (MODE == "check") {
  for (nm in names(LOCKS)) {
    lk <- build_lock(nm, LOCKS[[nm]], reps = 2L, max_reps = 2L)
    idt <- cleanTMLE:::.identification_table(lk)
    rx_log(sprintf("CHECK %s: n = %d (%d treated, %d comparator); %d covariates; flags: %s", nm,
                   nrow(lk$data), sum(lk$data$A), sum(lk$data$A == 0), length(lk$covariates),
                   paste(paste(idt$component, idt$item, sep = ": ")[idt$flag == "flag"], collapse = " | ")))
  }
  quit(save = "no")
}

# ------------------------------------------------------------
# Verify: rebuild every lock of this run that already has a design report and
# compare its hash with the lock on disk, so that no edit to this script can
# change a lock that has been simulated. Run before every relaunch.
# ------------------------------------------------------------
if (MODE == "verify") {
  bad <- character()
  for (nm in names(LOCKS)) {
    if (!file.exists(file.path(OUT, paste0(nm, "_dossier.rds")))) next
    old <- load_lock(file.path(OUT, paste0(nm, "_lock.rds")))$lock_hash
    new <- build_lock(nm, LOCKS[[nm]], REPS)$lock_hash
    rx_log(sprintf("VERIFY %s %s: %s", RUN, nm, if (identical(old, new)) "hash reproduced" else "HASH DIFFERS"))
    if (!identical(old, new)) bad <- c(bad, nm)
  }
  if (length(bad)) stop("Rebuilt locks differ from those on disk: ", paste(bad, collapse = ", "), call. = FALSE)
  quit(save = "no")
}

# ------------------------------------------------------------
# Pilot: time two repetitions on L3 and L1
# ------------------------------------------------------------
if (MODE == "pilot") {
  for (nm in intersect(c("L3", "L1"), names(LOCKS))) {
    lk <- build_lock(nm, LOCKS[[nm]], reps = 2L, max_reps = 2L)
    t0 <- Sys.time()
    ds <- assess_design(lk)
    t1 <- Sys.time()
    sm <- simulate_design(lk, ds)
    t2 <- Sys.time()
    rx_log(sprintf("PILOT %s: n = %d (%d treated); assess %.1f s; simulate (2 reps) %.1f s; per rep %.1f s",
                   nm, nrow(lk$data), sum(lk$data$A), as.numeric(t1 - t0, units = "secs"),
                   as.numeric(t2 - t1, units = "secs"), as.numeric(t2 - t1, units = "secs") / 2))
  }
  quit(save = "no")
}

# ------------------------------------------------------------
# Run: every lock, checkpointed
# ------------------------------------------------------------
if (WORKERS > 1L) {
  future::plan(future::multisession, workers = WORKERS)
  rx_log(paste("Repetitions run on", WORKERS, "parallel workers"))
}
# ------------------------------------------------------------
# ------------------------------------------------------------
if (!is.null(ONLY)) stopifnot(ONLY %in% names(LOCKS))
for (nm in names(LOCKS)) {
  if (!is.null(ONLY) && nm != ONLY) next
  done <- file.path(OUT, paste0(nm, "_dossier.rds"))
  if (file.exists(done)) { rx_log(paste(nm, "already done; skipped")); next }
  spec <- LOCKS[[nm]]
  t0 <- Sys.time()
  lk <- build_lock(nm, spec, REPS)
  save_lock(lk, file.path(OUT, paste0(nm, "_lock.rds")))
  ds <- assess_design(lk); saveRDS(ds, file.path(OUT, paste0(nm, "_design.rds")))
  sm <- simulate_design(lk, ds); saveRDS(sm, file.path(OUT, paste0(nm, "_simulation.rds")))
  nc <- negative_control_ladder(lk, ds, sm); saveRDS(nc, file.path(OUT, paste0(nm, "_nc.rds")))
  dos <- design_report(lk, ds, sm, nc, file = file.path(RES, paste0(nm, "_dossier.html")))
  saveRDS(dos, done)
  utils::write.csv(dos$identification, file.path(RES, paste0(nm, "_identification.csv")), row.names = FALSE)
  rx_log(sprintf("%s (%s): n = %d; primary %s; NC %s; lock %s; %.1f min", nm, spec$label,
                 nrow(lk$data), dos$decision$primary, dos$decision$nc_verdict,
                 substr(lk$lock_hash, 1, 12), as.numeric(Sys.time() - t0, units = "mins")))
}
rx_log(paste("Run", RUN, "complete."))
