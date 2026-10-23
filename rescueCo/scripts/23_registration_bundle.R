# ============================================================
# 23. Registration bundle: population gate, confirmatory comparison, hashes
# ============================================================
# Reads the final run's design reports (no outcome exists yet) and applies, by
# the rules written in the plan (spec section 5): the population gate to each
# lock's primary estimand, and the outcome-free rule that picks the comparison
# carrying the confirmatory claim. Runs twice:
#   pick   once L1 to L5 all have design reports. Records the confirmatory
#          comparison for script 22b, and saves and hashes the cohort of the
#          longitudinal TMLE (the comparison's rows before the complete-case
#          restriction) before any outcome is joined.
#   final  once every sensitivity lock and the 180-day lock has a design
#          report, or 7 days after the last of L1 to L5 has one (decision 7):
#          a lock unfinished then, or finished later, is reported as not run.
#          Writes the tables, the bundle with every hash, and the approvals
#          template that Andrew completes (decision 3).
# Copies of the bundle and the hashes go to rescueCo/registration/, which is
# committed and pushed before script 24 builds the outcome store (decision 4).
#
# Usage: Rscript rescueCo/scripts/23_registration_bundle.R pick|final
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md and
# the decisions of 2026-10-23.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/reexec_guards.R")
.libPaths(c("rescueCo/lib", .libPaths()))
suppressMessages(library(cleanTMLE))
cfg <- load_cr_config()
MODE <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(MODE %in% c("pick", "final"))
rx_log(paste("=== 23 registration bundle,", MODE, "==="))

LOCKD <- file.path(cfg$data$design_dir, "locks", "final")
RES <- file.path(cfg$paths$results, "reexec_050")
ORDER <- paste0("L", 1:5)
SENS <- paste0("S", 1:4)
GATE_KEEP <- 0.80
GATE_SHIFT <- 0.25
GATE_COVS <- c("age", "eiss_ais", "head_injury_severity",
               "intent_Intentional..assault.homicide..", "inj_place_Private.House..Home")
ATO_ALLOWED <- c("L3", "S1", "S4")
CUTOFF_DAYS <- 7
dir.create(REG_DIR, showWarnings = FALSE)

gate <- function(nm, d) {
  e <- d$decision$primary
  if (is.na(e)) return(data.frame(lock = nm, estimand = NA, kept_treated = NA, kept_control = NA,
                                  max_shift = NA, passes = FALSE, why = "no feasible estimand"))
  if (e == "ATE") return(data.frame(lock = nm, estimand = e, kept_treated = 1, kept_control = 1,
                                    max_shift = 0, passes = TRUE, why = "whole population"))
  pop <- d$design$populations
  sh <- abs(pop[[paste0(e, "_shift")]][pop$covariate %in% GATE_COVS])
  cnt <- d$design$trimmed_profile$counts
  kt <- if (e == "trimmed_ATE") cnt$kept[cnt$group == "treated"] /
    (cnt$kept[cnt$group == "treated"] + cnt$removed[cnt$group == "treated"]) else 1
  kc <- if (e == "trimmed_ATE") cnt$kept[cnt$group == "control"] /
    (cnt$kept[cnt$group == "control"] + cnt$removed[cnt$group == "control"]) else 1
  why <- c(if (kt < GATE_KEEP || kc < GATE_KEEP) "keeps under 80% of an arm",
           if (max(sh, 0, na.rm = TRUE) > GATE_SHIFT) "shift over 0.25 SD on a named covariate",
           if (e == "ATO" && !sub("_death180$", "", nm) %in% ATO_ALLOWED)   # the 180-day lock is its comparison's
             "overlap-weighted effect not declared acceptable here")
  data.frame(lock = nm, estimand = e, kept_treated = round(kt, 3), kept_control = round(kc, 3),
             max_shift = round(max(sh, 0, na.rm = TRUE), 3), passes = !length(why),
             why = if (length(why)) paste(why, collapse = "; ") else "passes")
}

dossier_path <- function(nm) file.path(LOCKD, paste0(nm, "_dossier.rds"))
created <- function(nm) as.POSIXct(readRDS(dossier_path(nm))$created,
                                   format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
lock_outcome <- function(nm) load_lock(file.path(LOCKD, paste0(nm, "_lock.rds")))$plan$outcome
summary_row <- function(nm, d, g) data.frame(
  lock = nm, lock_hash = d$lock_hash, dossier_hash = d$dossier_hash, outcome = lock_outcome(nm),
  overlap = d$design$overlap_grade, share_outside_band = round(d$design$share_outside_band, 3),
  primary = d$decision$primary, nc_verdict = d$decision$nc_verdict,
  status = paste(names(d$decision$status), unlist(d$decision$status), sep = ": ", collapse = "; "),
  gate = g$why, stringsAsFactors = FALSE)

miss <- ORDER[!file.exists(dossier_path(ORDER))]
if (length(miss))
  stop("Script 23 needs a design report for every primary comparison; missing: ",
       paste(miss, collapse = ", "), ".", call. = FALSE)

rows <- list(); gates <- list()
for (nm in ORDER) {
  d <- readRDS(dossier_path(nm))
  gates[[nm]] <- gate(nm, d)
  rows[[nm]] <- summary_row(nm, d, gates[[nm]])
}
tab <- do.call(rbind, rows); gtab <- do.call(rbind, gates)

# Confirmatory comparison: the first of L1 to L5, in declared order, with good or
# moderate overlap, a primary estimand that passes the gate, and no STOP.
ok <- with(tab, overlap %in% c("good", "moderate") & lock %in% gtab$lock[gtab$passes] &
             !(nc_verdict %in% "STOP"))
pick <- tab$lock[ok][1]
rx_log(paste("Confirmatory comparison by the rule in the plan:", if (is.na(pick)) "none qualifies" else pick))

# ------------------------------------------------------------
# pick: record the comparison, and freeze the longitudinal TMLE's cohort
# ------------------------------------------------------------
if (MODE == "pick") {
  utils::write.csv(data.frame(lock = pick, primary = if (is.na(pick)) NA else tab$primary[tab$lock == pick],
                              decided = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
                   file.path(REG_DIR, "confirmatory.csv"), row.names = FALSE)
  if (is.na(pick)) {
    rx_log("No comparison qualifies: no 180-day lock and no longitudinal TMLE cohort.")
    quit(save = "no")
  }
  dd <- readRDS(file.path(cfg$data$design_dir, "design_data_final.rds"))
  lk <- load_lock(file.path(LOCKD, paste0(pick, "_lock.rds")))
  D <- dd$data[rx_primary_rows(dd$data, pick), , drop = FALSE]
  D$A <- as.integer(D$arm %in% RX_PRIMARY[[pick]]$treated)
  co <- D[, c("patient_id", "A", lk$covariates, "disposition_recorded", "six_month_recorded")]
  cc <- co[co$six_month_recorded == 1, , drop = FALSE]
  m <- match(lk$data[[lk$id]], cc$patient_id)
  if (nrow(cc) != nrow(lk$data) || anyNA(m) || !identical(cc$A[m], as.integer(lk$data[[lk$treatment]])))
    stop("The cohort rebuilt for the longitudinal TMLE does not match ", pick,
         "'s complete-case rows; check RX_PRIMARY against script 22.", call. = FALSE)
  cf <- file.path(cfg$data$design_dir, paste0("lmtp_cohort_", pick, ".rds"))
  saveRDS(co, cf)
  utils::write.csv(data.frame(lock = pick, file = cf, n = nrow(co), n_complete_case = nrow(cc),
                              md5 = rx_md5(cf)),
                   file.path(REG_DIR, "lmtp_cohort.csv"), row.names = FALSE)
  rx_log(sprintf("Longitudinal TMLE cohort for %s: %d patients (%d complete case); md5 %s",
                 pick, nrow(co), nrow(cc), rx_md5(cf)))
  quit(save = "no")
}

# ------------------------------------------------------------
# final: the sensitivity and 180-day locks, the cut-off, the bundle
# ------------------------------------------------------------
recorded <- utils::read.csv(file.path(REG_DIR, "confirmatory.csv"), stringsAsFactors = FALSE)$lock[1]
if (!identical(is.na(recorded), is.na(pick)) || (!is.na(pick) && recorded != pick))
  stop("The rule now picks ", pick, " but pick mode recorded ", recorded, ".", call. = FALSE)
SECONDARY <- if (is.na(pick)) character() else paste0(pick, "_death180")
OPT <- c(SENS, SECONDARY)
cutoff <- max(do.call(c, lapply(ORDER, created))) + CUTOFF_DAYS * 86400
state <- vapply(OPT, function(nm) {
  if (!file.exists(dossier_path(nm))) "unfinished at the cut-off"
  else if (created(nm) > cutoff) "finished after the cut-off" else "done"
}, character(1))
if (any(!file.exists(dossier_path(OPT))) && Sys.time() < cutoff)
  stop("Unfinished: ", paste(OPT[!file.exists(dossier_path(OPT))], collapse = ", "),
       ". Script 23 (final) waits until they finish or until the cut-off, ",
       format(cutoff, "%Y-%m-%d %H:%M UTC"), ".", call. = FALSE)

for (nm in OPT[state == "done"]) {
  d <- readRDS(dossier_path(nm))
  gates[[nm]] <- gate(nm, d)
  rows[[nm]] <- summary_row(nm, d, gates[[nm]])
}
tab <- do.call(rbind, rows); gtab <- do.call(rbind, gates)
role <- function(nm) if (identical(nm, pick)) "confirmatory" else if (nm %in% ORDER) "descriptive" else
  if (nm %in% SENS) "sensitivity" else "secondary"
tab$role <- vapply(tab$lock, role, character(1))
not_run <- OPT[state != "done"]
if (length(not_run))
  rx_log(paste("Reported as not run (decision 7):", paste(not_run, state[state != "done"], collapse = "; ")))

ap <- rbind(
  data.frame(lock = tab$lock, role = tab$role, lock_hash = tab$lock_hash, dossier_hash = tab$dossier_hash,
             outcome = tab$outcome, primary = tab$primary, nc_verdict = tab$nc_verdict,
             approved_by = "", approved_on = "", note = ""),
  if (length(not_run))
    data.frame(lock = not_run, role = "not run", lock_hash = NA, dossier_hash = NA, outcome = NA,
               primary = NA, nc_verdict = NA, approved_by = "", approved_on = "",
               note = state[state != "done"]))

note <- c("# Rescue.Co re-execution: registration bundle",
          "", sprintf("Built %s from the final run's design reports, before any outcome was joined.",
                      format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
          "", sprintf("Confirmatory comparison by the rule in the plan: **%s**.",
                      if (is.na(pick)) "none qualifies" else pick),
          sprintf("Cut-off for the sensitivity and 180-day locks: %s (%d days after the last primary design report).",
                  format(cutoff, "%Y-%m-%d %H:%M UTC"), CUTOFF_DAYS),
          if (length(not_run)) sprintf("Not run: %s.", paste(not_run, state[state != "done"], sep = ", ", collapse = "; ")),
          "", "| Lock | Role | Outcome | Lock hash | Design-report hash | Overlap | Primary | Negative controls | Gate |",
          "|---|---|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |", tab$lock, tab$role, tab$outcome,
                  tab$lock_hash, tab$dossier_hash, tab$overlap, tab$primary, tab$nc_verdict, tab$gate))
for (dir in c(RES, REG_DIR)) {
  utils::write.csv(tab, file.path(dir, "registration_locks.csv"), row.names = FALSE)
  utils::write.csv(gtab, file.path(dir, "population_gate.csv"), row.names = FALSE)
  writeLines(note, file.path(dir, "registration_bundle.md"))
}
utils::write.csv(ap, file.path(REG_DIR, "approvals_final_template.csv"), row.names = FALSE)
rx_log("Wrote the bundle and approvals_final_template.csv; complete the template as approvals_final.csv, commit and push.")
