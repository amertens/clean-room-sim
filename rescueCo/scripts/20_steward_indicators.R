# ============================================================
# 20. Data steward: the two outcome-observation indicators (D3)
# ============================================================
# Run by the data steward before the design stage. Links the six-month
# follow-up to the registry, scores the outcomes in memory, and writes ONLY two
# indicators per patient: disposition recorded, and six-month status recorded
# (the primary outcome is known). No outcome value is written, printed or
# logged; the log carries registry-wide totals of the two indicators, which D3
# approves for the design stage (the blinding ledger records what they reveal).
#
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/steward.R")
cfg <- load_cr_config()
cr_log("=== 20 steward: outcome-observation indicators ===")

reg <- haven::read_dta(cfg$data$trauma_registry_dta,
                       col_select = tidyselect::all_of(c("patient_id", "disposition",
                                                         "arrival_date")))
fu <- steward_read_followup(cfg$data$sixmo_followup_dta)
fu <- steward_link_followup(fu, as.character(reg$patient_id),
                            cfg$data$followup_typo_audit, cfg$data$followup_recon_actions)
n_linked <- sum(as.character(reg$patient_id) %in% fu$patients_id_gose)
cr_log(paste("Registry patients with a linked follow-up record:", n_linked,
             "(the unblinded analysis links 7,784)"))
if (n_linked != 7784L)
  stop("Linkage differs from the unblinded analysis (", n_linked, " vs 7,784).", call. = FALSE)

ind <- steward_indicators(steward_score(reg, fu))
rm(fu)
cr_log(paste("Disposition recorded:", sum(ind$disposition_recorded), "of", nrow(ind),
             "; six-month status recorded:", sum(ind$six_month_recorded), "of", nrow(ind)))

dir.create(cfg$data$design_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(ind, file.path(cfg$data$design_dir, "observation_indicators.rds"))
cr_log("Wrote observation_indicators.rds (id and two indicators only).")
