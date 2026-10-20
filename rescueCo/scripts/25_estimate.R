# ============================================================
# 25. After approval: unblind each lock and estimate
# ============================================================
# For each lock of the final run: unblind() with the approved dossier and the
# steward's outcome store, then estimate_effect() with the estimand and the
# candidate the dossier selected (complete case, the estimator the design
# stage evaluated). Then the longitudinal TMLE sensitivity analysis for the
# confirmatory comparison, with lmtp, which the design stage did not evaluate.
#
# Usage: Rscript rescueCo/scripts/25_estimate.R "<approver name>"
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
.libPaths(c("rescueCo/lib", .libPaths()))
suppressMessages(library(cleanTMLE))
cfg <- load_cr_config()
approved_by <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(approved_by) || !nzchar(approved_by))
  stop("Give the approver's name as the first argument.", call. = FALSE)
cr_log("=== 25 estimation after approval ===")

LOCKD <- file.path(cfg$data$design_dir, "locks", "final")
RES <- file.path(cfg$paths$results, "reexec_050")
store <- readRDS(file.path(cfg$data$outcome_dir, "outcome_store.rds"))

fits <- list()
for (f in list.files(LOCKD, pattern = "_lock\\.rds$", full.names = TRUE)) {
  nm <- sub("_lock\\.rds$", "", basename(f))
  lk <- load_lock(f)
  d <- readRDS(file.path(LOCKD, paste0(nm, "_dossier.rds")))
  oc <- store[store$patient_id %in% lk$data$patient_id, c("patient_id", "severe_disability")]
  ub <- unblind(lk, d, oc, approved_by = approved_by)
  saveRDS(ub, file.path(LOCKD, paste0(nm, "_unblinded.rds")))
  fit <- estimate_effect(ub)
  saveRDS(fit, file.path(LOCKD, paste0(nm, "_fit.rds")))
  export_design_log(ub, fit, file = file.path(RES, paste0(nm, "_design_log.csv")), format = "muntner")
  fits[[nm]] <- as.data.frame(fit)
  cr_log(paste(nm, "estimated:", fit$estimand))
}
utils::write.csv(do.call(rbind, fits), file.path(RES, "estimates_complete_case.csv"), row.names = FALSE)

# Longitudinal TMLE sensitivity analysis (lmtp), confirmatory comparison only.
# Nodes as in the unblinded analysis: disposition recorded, death at
# disposition (sets the outcome), six-month status recorded, outcome. Run with
# the design stage's covariates fixed before transport; arrival state enters
# the later regressions as a node after transport is NOT available in the
# package, so here it is omitted (decision 4). Labelled as not evaluated in the
# design stage.
# To be completed when the confirmatory comparison is known (script 23).
cr_log("Complete-case estimates written; longitudinal TMLE sensitivity analysis pending the confirmatory comparison.")
