# ============================================================
# 23. Registration bundle: population gate, confirmatory comparison, hashes
# ============================================================
# Reads the final run's dossiers (no outcome exists yet) and applies, by the
# rules written in the plan (spec section 5): the population gate to each
# lock's primary estimand, and the outcome-free rule that picks the comparison
# carrying the confirmatory claim. Writes the tables and a registration note
# with every lock hash, for the OSF registration (D9).
#
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
.libPaths(c("rescueCo/lib", .libPaths()))
suppressMessages(library(cleanTMLE))
cfg <- load_cr_config()
cr_log("=== 23 registration bundle ===")

LOCKD <- file.path(cfg$data$design_dir, "locks", "final")
RES <- file.path(cfg$paths$results, "reexec_050")
ORDER <- paste0("L", 1:5)
ALL <- c(ORDER, paste0("S", 1:4))
GATE_KEEP <- 0.80
GATE_SHIFT <- 0.25
GATE_COVS <- c("age", "eiss_ais", "head_injury_severity",
               "intent_Intentional..assault.homicide..", "inj_place_Private.House..Home")
ATO_ALLOWED <- c("L3", "S1", "S4")

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
           if (e == "ATO" && !nm %in% ATO_ALLOWED) "overlap-weighted effect not declared acceptable here")
  data.frame(lock = nm, estimand = e, kept_treated = round(kt, 3), kept_control = round(kc, 3),
             max_shift = round(max(sh, 0, na.rm = TRUE), 3), passes = !length(why),
             why = if (length(why)) paste(why, collapse = "; ") else "passes")
}

rows <- list(); gates <- list()
for (nm in ALL) {
  f <- file.path(LOCKD, paste0(nm, "_dossier.rds"))
  if (!file.exists(f)) { cr_log(paste(nm, "has no dossier yet")); next }
  d <- readRDS(f)
  g <- gate(nm, d); gates[[nm]] <- g
  rows[[nm]] <- data.frame(
    lock = nm, lock_hash = d$lock_hash, dossier_hash = d$dossier_hash,
    overlap = d$design$overlap_grade, share_outside_band = round(d$design$share_outside_band, 3),
    primary = d$decision$primary, nc_verdict = d$decision$nc_verdict,
    status = paste(names(d$decision$status), unlist(d$decision$status), sep = ": ", collapse = "; "),
    gate = g$why, stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows); gtab <- do.call(rbind, gates)
utils::write.csv(tab, file.path(RES, "registration_locks.csv"), row.names = FALSE)
utils::write.csv(gtab, file.path(RES, "population_gate.csv"), row.names = FALSE)

# Confirmatory comparison: the first of L1 to L5, in declared order, with good or
# moderate overlap, a primary estimand that passes the gate, and no STOP.
ok <- with(tab[tab$lock %in% ORDER, ], overlap %in% c("good", "moderate") &
             lock %in% gtab$lock[gtab$passes] & !(nc_verdict %in% "STOP"))
pick <- tab$lock[tab$lock %in% ORDER][ok][1]
cr_log(paste("Confirmatory comparison by the registered rule:",
             if (is.na(pick)) "none qualifies" else pick))

note <- c("# Rescue.Co re-execution: registration bundle",
          "", sprintf("Built %s from the final run's dossiers, before any outcome was joined.",
                      format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
          "", sprintf("Confirmatory comparison by the registered rule: **%s**.",
                      if (is.na(pick)) "none qualifies" else pick),
          "", "| Lock | Lock hash | Dossier hash | Overlap | Primary | Negative controls | Gate |",
          "|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %s | %s | %s | %s | %s |", tab$lock, tab$lock_hash, tab$dossier_hash,
                  tab$overlap, tab$primary, tab$nc_verdict, tab$gate))
writeLines(note, file.path(RES, "registration_bundle.md"))
cr_log("Wrote registration_locks.csv, population_gate.csv, registration_bundle.md")
