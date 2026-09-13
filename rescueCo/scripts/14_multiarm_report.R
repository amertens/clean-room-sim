# ============================================================
# Multi-arm Stage 5: sensitivity, artifacts, reconciliation, report
# ============================================================
# Computes E-values from the TMLE-adjusted arm risks and the additive
# bias-to-null distances on the ladder estimates, regenerates the
# manuscript artifacts from the multiarm outputs (never edited by hand),
# rebuilds the reconciliation table, and renders the estimability case
# study. Run after scripts 10-13.
# ============================================================

suppressMessages(library(cleanTMLE))
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                 "bootstrap.R"))
source("rescueCo/R/utils.R")
cr_log("=== Multi-arm Stage 5: sensitivity, artifacts, report ===")

OUT <- "rescueCo/results/multiarm"
ART <- "rescueCo/results/manuscript_artifacts"
est <- utils::read.csv(file.path(OUT, "ladder_estimates.csv"),
                       stringsAsFactors = FALSE)
locks <- readRDS(file.path(OUT, "multiarm_locks.rds"))

# ── E-values from adjusted arm risks, and the bias that moves the CI to
#    null (the Diaz and van der Laan additive reading), per binomial row ──
sens <- do.call(rbind, lapply(seq_len(nrow(est)), function(i) {
  r <- est[i, ]
  if (!identical(r$family, "binomial") || !is.finite(r$risk_treated) ||
      !is.finite(r$risk_control) || r$risk_control <= 0)
    return(NULL)
  rr <- r$risk_treated / r$risk_control
  # RR confidence bound on the same scale, delta-method on the log-RR from
  # the RD interval mapped through the control risk (approximation recorded
  # as such; the exact RR interval lives in the tmle fit objects).
  rr_lo <- max((r$risk_treated - 1.96 * r$se) / r$risk_control, 1e-6)
  rr_hi <- (r$risk_treated + 1.96 * r$se) / r$risk_control
  ev <- tryCatch(compute_evalue(rr, ci_bound = if (rr < 1) rr_hi else rr_lo),
                 error = function(e) list(evalue_point = NA_real_,
                                          evalue_ci = NA_real_))
  # Additive bias-to-null: |estimate| to move the point to 0, and the
  # smaller CI distance to 0 (the bound that crosses first).
  d_point <- abs(r$estimate)
  d_ci <- if (r$ci_lower > 0) r$ci_lower else if (r$ci_upper < 0)
    -r$ci_upper else 0
  data.frame(contrast = r$contrast, outcome = r$outcome,
             estimand = r$estimand, rr = round(rr, 4),
             evalue_point = round(as.numeric(ev$evalue_point %||% ev[[1]]), 3),
             evalue_ci = round(as.numeric(ev$evalue_ci %||% NA_real_), 3),
             bias_to_null_point = round(d_point, 5),
             bias_to_null_ci = round(abs(d_ci), 5),
             floor_declared = 0.0035,
             stringsAsFactors = FALSE)
}))
if (!is.null(sens)) {
  utils::write.csv(sens, file.path(OUT, "sensitivity_evalues.csv"),
                   row.names = FALSE)
  cr_log(paste("Wrote sensitivity_evalues.csv (", nrow(sens), "rows )"))
}

# ── Manuscript artifacts, regenerated from the pipeline outputs ─────────────
dir.create(ART, showWarnings = FALSE, recursive = TRUE)

# Effects table: the ladder estimates.
utils::write.csv(est[, c("contrast", "outcome", "estimand", "is_primary",
                         "estimate", "se", "ci_lower", "ci_upper", "p_value",
                         "n", "support_verdict", "implausible")],
                 file.path(ART, "table_effects_multiarm.csv"),
                 row.names = FALSE)

# Support and gate path.
file.copy(file.path(OUT, "support_by_contrast.csv"),
          file.path(ART, "table_support_by_contrast.csv"), overwrite = TRUE)
file.copy(file.path(OUT, "estimand_feasibility.csv"),
          file.path(ART, "table_estimand_feasibility.csv"), overwrite = TRUE)
if (file.exists(file.path(OUT, "nc_ladder.csv")))
  file.copy(file.path(OUT, "nc_ladder.csv"),
            file.path(ART, "table_nc_ladder.csv"), overwrite = TRUE)
if (file.exists(file.path(OUT, "ladder_design_log.csv")))
  file.copy(file.path(OUT, "ladder_design_log.csv"),
            file.path(ART, "decision_summary_multiarm.csv"),
            overwrite = TRUE)
for (cn in names(locks)) {
  f <- file.path(OUT, sprintf("fig_support_%s.png", cn))
  if (file.exists(f))
    file.copy(f, file.path(ART, sprintf("fig_support_%s.png", cn)),
              overwrite = TRUE)
}

# Case-study metadata, computed.
meta <- list(
  case_study = "Rescue.Co multi-arm estimability",
  cleanTMLE_version = as.character(utils::packageVersion("cleanTMLE")),
  data_cut = "2026-03-27",
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  locks = lapply(names(locks), function(cn) list(
    contrast = cn, label = locks[[cn]]$contrast$label,
    n = nrow(locks[[cn]]$data),
    n_treated = sum(locks[[cn]]$data[[locks[[cn]]$treatment]]),
    lock_hash = locks[[cn]]$lock_hash)),
  primary_outcome = "death_by_6mo",
  ladder = locks[[1]]$estimand_ladder[c("primary", "fallbacks", "trigger")])
if (requireNamespace("yaml", quietly = TRUE)) {
  writeLines(yaml::as.yaml(meta),
             file.path(ART, "case_study_metadata_multiarm.yml"))
} else {
  saveRDS(meta, file.path(ART, "case_study_metadata_multiarm.rds"))
}

# ── Reconciliation table, regenerated ────────────────────────────────────────
source("rescueCo/scripts/build_reconciliation_table.R")

# ── Render the case study ────────────────────────────────────────────────────
qy <- Sys.which("quarto")
if (nzchar(qy)) {
  cr_log("Rendering rescueco_estimability_case_study.qmd")
  status <- system2(qy, c("render",
                          "rescueCo/reports/rescueco_estimability_case_study.qmd"),
                    stdout = TRUE, stderr = TRUE)
  writeLines(status, file.path(OUT, "render_log.txt"))
} else {
  cr_log("quarto not on PATH; render the report manually")
}
cr_log("Multi-arm Stage 5 complete.")
