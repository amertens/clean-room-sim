# ============================================================
# 24. Data steward: the outcome store (only after approval)
# ============================================================
# Run by the data steward only after Andrew has approved every design report
# of the final run that is to be unblinded: rescueCo/registration/
# approvals_final.csv lists each lock with its lock and design-report hashes,
# the approver and the date, and it must be committed and on a remote branch
# (decisions 3 and 4 of 2026-10-23; the demonstration has no outside approver
# or registration, spec 2026-10-18, "Governance in the demonstration"). Scores
# the outcomes as the unblinded analysis does (R/steward.R) and writes the
# outcome store. The log carries registry-wide, arm-blind counts of each
# scoring rule, never an outcome by arm. Never run before approval, not even
# as a rehearsal; the steward functions are tested on a synthetic follow-up
# file instead.
#
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/steward.R")
source("rescueCo/R/reexec_guards.R")
cfg <- load_cr_config()
ap <- rx_read_approvals(file.path(cfg$data$design_dir, "locks", "final"))
if (file.exists(file.path(cfg$data$outcome_dir, "outcome_store.rds")))
  stop("The outcome store exists already; script 24 runs once.", call. = FALSE)
rx_log(paste("=== 24 steward: outcome store; approved:", paste(ap$lock, collapse = ", "), "==="))

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
rx_log(paste("Registry-wide scoring counts: GOSE overridden by a hospital death",
             sum(sc$gose_death_reconciled), "; death date implausible", sum(sc$death_date_implausible),
             "; primary outcome known", sum(!is.na(sc$severe_disability)),
             "; death within 180 days known", sum(!is.na(sc$death_by_6mo))))

dir.create(cfg$data$outcome_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sc[, c("patient_id", "severe_disability", "death_by_6mo", "death_in_hospital", "gose")],
        file.path(cfg$data$outcome_dir, "outcome_store.rds"))
rx_log("Wrote outcome_store.rds.")
