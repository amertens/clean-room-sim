# ============================================================
# 24. Data steward: the outcome store (only after approval)
# ============================================================
# Run by the data steward only after the dossiers are approved and the lock
# hashes are registered (D9). Refuses to run without the approval record
# rescueCo/data/outcomes/APPROVAL.txt, which names the approver, the date, and
# the registration link. Scores the outcomes as the unblinded analysis does
# (R/steward.R) and writes the outcome store. The log carries registry-wide,
# arm-blind counts of each scoring rule, never an outcome by arm.
#
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/steward.R")
cfg <- load_cr_config()
approval <- file.path(cfg$data$outcome_dir, "APPROVAL.txt")
if (!file.exists(approval))
  stop("No approval record at ", approval, ". The outcome store is built only after ",
       "the dossiers are approved and the lock hashes registered (D9).", call. = FALSE)
cr_log("=== 24 steward: outcome store ===")
cr_log(paste("Approval record:", paste(readLines(approval, warn = FALSE), collapse = " | ")))

reg <- haven::read_dta(cfg$data$trauma_registry_dta,
                       col_select = tidyselect::all_of(c("patient_id", "disposition", "arrival_date")))
fu <- steward_read_followup(cfg$data$sixmo_followup_dta)
fu <- steward_link_followup(fu, as.character(reg$patient_id),
                            cfg$data$followup_typo_audit, cfg$data$followup_recon_actions)
sc <- steward_score(reg, fu)
rm(fu)

# Arm-blind checks: the indicators the design stage used must match the store.
ind <- readRDS(file.path(cfg$data$design_dir, "observation_indicators.rds"))
chk <- steward_indicators(sc)
stopifnot(identical(chk$six_month_recorded, ind$six_month_recorded[match(chk$patient_id, ind$patient_id)]))
cr_log(paste("Registry-wide scoring counts: GOSE overridden by a hospital death",
             sum(sc$gose_death_reconciled), "; death date implausible", sum(sc$death_date_implausible),
             "; primary outcome known", sum(!is.na(sc$severe_disability)),
             "; death within 180 days known", sum(!is.na(sc$death_by_6mo))))

dir.create(cfg$data$outcome_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sc[, c("patient_id", "severe_disability", "death_by_6mo", "death_in_hospital", "gose")],
        file.path(cfg$data$outcome_dir, "outcome_store.rds"))
cr_log("Wrote outcome_store.rds.")
