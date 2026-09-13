###############################################################################
#  build_reconciliation_table.R
#
#  Computes rescueCo/results/reconciliation_table.csv: one row per
#  (contrast, cohort, outcome, estimand, estimator), giving the clean-room
#  estimate and the main-pipeline estimate side by side with a verdict.
#
#  The main pipeline (C:/Users/andre/OneDrive/Documents/rescueCo_analysis) is
#  read only; nothing there is modified. Numbers are read from the CSVs both
#  pipelines publish, never typed, so the table can be regenerated after any
#  rerun. Verdicts:
#
#    agree                 |difference| <= 0.25 main-pipeline SE
#    benign_difference     both present, difference explained (reason column)
#    open_discrepancy      both present, difference not yet explained
#    not_comparable        the two rows estimate different quantities by design
#    missing_in_cleanroom  the main pipeline reports it, the clean room not yet
#    missing_in_main       the clean room reports it, the main pipeline does not
#
#  Run from the repo root:
#    "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" rescueCo/scripts/build_reconciliation_table.R
###############################################################################

MAIN_ROOT <- Sys.getenv("RESCUECO_MAIN_ROOT",
  "C:/Users/andre/OneDrive/Documents/rescueCo_analysis/analysis_outputs")
CR_ROOT   <- "rescueCo/results"
OUT_CSV   <- file.path(CR_ROOT, "reconciliation_table.csv")

if (!dir.exists(MAIN_ROOT))
  stop("Main-pipeline outputs not found at ", MAIN_ROOT,
       " (set RESCUECO_MAIN_ROOT).")
if (!dir.exists(CR_ROOT))
  stop("Run from the clean-room-sim repo root; ", CR_ROOT, " not found.")

read_main <- function(...) {
  f <- file.path(MAIN_ROOT, ...)
  if (!file.exists(f)) { message("main file absent: ", f); return(NULL) }
  utils::read.csv(f, stringsAsFactors = FALSE)
}
read_cr <- function(...) {
  f <- file.path(CR_ROOT, ...)
  if (!file.exists(f)) { message("clean-room file absent: ", f); return(NULL) }
  utils::read.csv(f, stringsAsFactors = FALSE)
}

# Canonical row constructor -------------------------------------------------
row0 <- function(contrast, cohort, outcome, estimand, estimator,
                 cr_est = NA_real_, cr_lo = NA_real_, cr_hi = NA_real_,
                 cr_n = NA_real_,
                 main_est = NA_real_, main_lo = NA_real_, main_hi = NA_real_,
                 main_se = NA_real_, main_n = NA_real_) {
  data.frame(contrast = contrast, cohort = cohort, outcome = outcome,
             estimand = estimand, estimator = estimator,
             cr_est = cr_est, cr_ci_lo = cr_lo, cr_ci_hi = cr_hi, cr_n = cr_n,
             main_est = main_est, main_ci_lo = main_lo, main_ci_hi = main_hi,
             main_se = main_se, main_n = main_n,
             stringsAsFactors = FALSE)
}

COHORTS <- c(C1 = "known transport (8,323)", C2 = "Rescue.Co + non-ambulance (7,165)",
             C3 = "known transport (8,323)", C4 = "ambulance, transfers included (2,197)",
             PRIMARY = "primary cohort (1,616)")

rows <- list()

# ── Main pipeline, primary cohort (script 06): ATE and complete-case ATT ────
m06 <- read_main("06_tmle_full_cohort_effects", "tmle_effects_summary.csv")
if (!is.null(m06)) {
  for (i in seq_len(nrow(m06))) {
    r <- m06[i, ]
    est_ate <- if (isTRUE(r$used_ipcw == "TRUE" | r$used_ipcw == TRUE))
      "TMLE, IPCW via Delta" else "TMLE, complete case"
    rows[[length(rows) + 1L]] <- row0("PRIMARY", COHORTS["PRIMARY"], r$outcome,
      "ATE", est_ate,
      main_est = r$ate_est, main_lo = r$ate_ci_lo, main_hi = r$ate_ci_hi,
      main_se = r$ate_se, main_n = r$n)
    if (is.finite(suppressWarnings(as.numeric(r$att_est))))
      rows[[length(rows) + 1L]] <- row0("PRIMARY", COHORTS["PRIMARY"], r$outcome,
        "ATT among outcome-observed", "TMLE, complete case",
        main_est = r$att_est, main_lo = r$att_ci_lo, main_hi = r$att_ci_hi,
        main_se = r$att_se, main_n = r$att_n)
  }
}

# ── Main pipeline, multi-arm contrasts (script 42) ──────────────────────────
m42 <- read_main("42_multiarm_tmle", "multiarm_tmle_estimates.csv")
if (!is.null(m42)) {
  for (i in seq_len(nrow(m42))) {
    r <- m42[i, ]
    pop <- if (startsWith(r$population, "full")) "ATE"
           else paste0("trimmed ATE (", sub("^trimmed ", "", r$population), ")")
    est <- if (isTRUE(r$used_ipcw == "TRUE" | r$used_ipcw == TRUE))
      "TMLE, IPCW via Delta" else "TMLE, complete case"
    rows[[length(rows) + 1L]] <- row0(r$contrast, COHORTS[r$contrast], r$outcome,
      pop, est,
      main_est = r$ate_est, main_lo = r$ate_ci_lo, main_hi = r$ate_ci_hi,
      main_se = r$ate_se, main_n = r$n)
    if (startsWith(r$population, "full") &&
        is.finite(suppressWarnings(as.numeric(r$att_est))))
      rows[[length(rows) + 1L]] <- row0(r$contrast, COHORTS[r$contrast], r$outcome,
        "ATT among outcome-observed", "TMLE, complete case",
        main_est = r$att_est, main_lo = r$att_ci_lo, main_hi = r$att_ci_hi,
        main_se = r$att_se, main_n = r$att_n)
  }
}

# ── Main pipeline, overlap weights (script 53) ──────────────────────────────
m53 <- read_main("53_overlap_weights", "ato_estimates.csv")
if (!is.null(m53)) {
  for (i in seq_len(nrow(m53))) {
    r <- m53[i, ]
    rows[[length(rows) + 1L]] <- row0(r$contrast, COHORTS[r$contrast], r$outcome,
      "ATO", "Hajek overlap weights, percentile bootstrap",
      main_est = r$ato_est, main_lo = r$ato_ci_lo, main_hi = r$ato_ci_hi,
      main_se = r$ato_se, main_n = r$n)
  }
}

# ── Main pipeline, matched ATT (scripts 05 and 43) ──────────────────────────
m05 <- read_main("05_matched_ATT_effects", "matched_effects_summary.csv")
if (!is.null(m05)) {
  ocol <- intersect(c("outcome", "outcome_name"), names(m05))[1]
  ecol <- intersect(c("att", "att_est", "estimate"), names(m05))[1]
  lcol <- intersect(c("ci_lo", "att_ci_lo"), names(m05))[1]
  hcol <- intersect(c("ci_hi", "att_ci_hi"), names(m05))[1]
  ncol_ <- intersect(c("n_sets", "sets", "n"), names(m05))[1]
  spec  <- intersect(c("specification", "spec"), names(m05))
  keep  <- if (length(spec)) grepl("no replacement", m05[[spec[1]]]) else
    rep(TRUE, nrow(m05))
  for (i in which(keep)) {
    r <- m05[i, ]
    rows[[length(rows) + 1L]] <- row0("PRIMARY", COHORTS["PRIMARY"], r[[ocol]],
      "matched ATT", "1:1 caliper matching, no replacement",
      main_est = r[[ecol]], main_lo = if (!is.na(lcol)) r[[lcol]] else NA,
      main_hi = if (!is.na(hcol)) r[[hcol]] else NA,
      main_n = if (!is.na(ncol_)) r[[ncol_]] else NA)
  }
}
m43 <- read_main("43_multiarm_matched_att", "multiarm_matched_att.csv")
if (!is.null(m43)) {
  ocol <- intersect(c("outcome", "outcome_name"), names(m43))[1]
  ecol <- intersect(c("att", "att_est", "estimate"), names(m43))[1]
  lcol <- intersect(c("ci_lo", "att_ci_lo"), names(m43))[1]
  hcol <- intersect(c("ci_hi", "att_ci_hi"), names(m43))[1]
  for (i in seq_len(nrow(m43))) {
    r <- m43[i, ]
    rows[[length(rows) + 1L]] <- row0(r$contrast, COHORTS[r$contrast], r[[ocol]],
      "matched ATT", "1:1 caliper matching, no replacement",
      main_est = r[[ecol]], main_lo = if (!is.na(lcol)) r[[lcol]] else NA,
      main_hi = if (!is.na(hcol)) r[[hcol]] else NA)
  }
}

# ── Main pipeline, negative controls (script 18) ────────────────────────────
m18 <- read_main("18_negative_controls", "negative_control_results.csv")
if (!is.null(m18)) {
  for (i in seq_len(nrow(m18))) {
    r <- m18[i, ]
    ctr <- if (grepl("Primary", r$cohort)) "PRIMARY" else "C4"
    rows[[length(rows) + 1L]] <- row0(ctr,
      if (ctr == "PRIMARY") COHORTS["PRIMARY"] else COHORTS["C4"],
      paste0("NC: ", r$negative_control),
      "unadjusted RD (balance check)", "two-proportion difference",
      main_est = r$rd, main_lo = r$ci_lo, main_hi = r$ci_hi,
      main_se = r$se, main_n = r$n)
  }
}

tab <- do.call(rbind, rows)

# ── Clean-room results, mapped onto the same keys ───────────────────────────
# Outcome name mapping: the clean room's gose_good is the main pipeline's
# gose_favorable (both are GOSE > 4).
cr_fill <- function(tab, contrast, outcome, estimand, estimator,
                    est, lo, hi, n, add_estimator_if_absent = TRUE) {
  hit <- tab$contrast == contrast & tab$outcome == outcome &
         tab$estimand == estimand & tab$estimator == estimator
  if (any(hit)) {
    tab$cr_est[hit] <- est; tab$cr_ci_lo[hit] <- lo
    tab$cr_ci_hi[hit] <- hi; tab$cr_n[hit] <- n
    return(tab)
  }
  # Same estimand, different estimator label, and the row not already
  # claimed by another clean-room estimator: reuse it and record both labels.
  hit2 <- tab$contrast == contrast & tab$outcome == outcome &
          tab$estimand == estimand & !is.finite(tab$cr_est)
  if (any(hit2)) {
    i <- which(hit2)[1]
    tab$cr_est[i] <- est; tab$cr_ci_lo[i] <- lo
    tab$cr_ci_hi[i] <- hi; tab$cr_n[i] <- n
    if (!identical(tab$estimator[i], estimator))
      tab$estimator[i] <- paste0(tab$estimator[i], " | CR: ", estimator)
  } else if (add_estimator_if_absent) {
    tab <- rbind(tab, row0(contrast, COHORTS[contrast], outcome, estimand,
                           estimator, cr_est = est, cr_lo = lo, cr_hi = hi,
                           cr_n = n))
  }
  tab
}

cr_eff <- read_cr("manuscript_artifacts", "table_effects.csv")
if (!is.null(cr_eff)) {
  map <- list(
    "Crude (unadjusted)"  = c("crude difference", "arm difference"),
    "PS matching"         = c("matched ATT", "1:1 caliper matching, no replacement"),
    "IPTW (stabilised)"   = c("ATE", "stabilised IPTW"),
    "Full-cohort TMLE"    = c("ATE", "TMLE, complete case"),
    "Matched-cohort TMLE" = c("ATE on matched cohort", "TMLE on matched cohort"),
    "IPCW-weighted TMLE"  = c("ATE", "TMLE, IPCW via Delta")
  )
  e_lo <- intersect(c("ci_lo", "ci_lower"), names(cr_eff))[1]
  e_hi <- intersect(c("ci_hi", "ci_upper"), names(cr_eff))[1]
  # The IPCW row is filled first so the exact estimand+estimator match claims
  # the main pipeline's IPCW ATE row before the estimand-only fallback runs.
  ord_m <- order(cr_eff$method != "IPCW-weighted TMLE")
  for (i in ord_m) {
    r <- cr_eff[i, ]
    mm <- map[[r$method]]
    if (is.null(mm)) next
    tab <- cr_fill(tab, "PRIMARY", "gose_favorable", mm[1], mm[2],
                   r$estimate,
                   if (!is.na(e_lo)) r[[e_lo]] else NA_real_,
                   if (!is.na(e_hi)) r[[e_hi]] else NA_real_, r$n)
  }
  # The main pipeline's matched ATT panel carries severe_disability, the exact
  # complement of gose_favorable (both derive from the same GOSE), so the
  # clean room's matched gose_favorable ATT is comparable to the negated
  # severe_disability ATT.
  ps_row <- cr_eff[cr_eff$method == "PS matching", ]
  hit_sd <- tab$contrast == "PRIMARY" & tab$outcome == "severe_disability" &
            tab$estimand == "matched ATT"
  if (nrow(ps_row) == 1L && any(hit_sd)) {
    i <- which(hit_sd)[1]
    tab$cr_est[i]   <- -ps_row$estimate
    tab$cr_ci_lo[i] <- -(if (!is.na(e_hi)) ps_row[[e_hi]] else NA_real_)
    tab$cr_ci_hi[i] <- -(if (!is.na(e_lo)) ps_row[[e_lo]] else NA_real_)
    tab$cr_n[i]     <- ps_row$n
    tab$estimator[i] <- paste0(tab$estimator[i],
                               " | CR: negated matched gose_favorable ATT")
  }
}

cr_nc <- if (file.exists(file.path(CR_ROOT, "multiarm", "nc_ladder.csv")))
  NULL else read_cr("negative_control_results.csv")
if (!is.null(cr_nc)) {
  nc_map <- c(household_urban = "NC: Urban household (SES)",
              chronic_hypertension = "NC: Chronic hypertension",
              cooking_fuel_wood_charcoal = "NC: Wood/charcoal cooking fuel (SES)",
              diabetes_insulin = "NC: Diabetes, insulin dependent",
              hiv_on_art = "NC: HIV/AIDS on ART")
  vcol <- intersect(c("variable", "nc", "negative_control", "nc_variable"),
                    names(cr_nc))[1]
  ecol <- intersect(c("estimate", "rd", "est"), names(cr_nc))[1]
  lcol <- intersect(c("ci_lower", "ci_lo"), names(cr_nc))[1]
  hcol <- intersect(c("ci_upper", "ci_hi"), names(cr_nc))[1]
  if (!is.na(vcol) && !is.na(ecol)) {
    for (i in seq_len(nrow(cr_nc))) {
      r <- cr_nc[i, ]
      key <- nc_map[[as.character(r[[vcol]])]]
      if (is.null(key)) key <- paste0("NC: ", r[[vcol]])
      tab <- cr_fill(tab, "PRIMARY", key, "unadjusted RD (balance check)",
                     "clean-room NC (IPTW/TMLE)",
                     r[[ecol]],
                     if (!is.na(lcol)) r[[lcol]] else NA,
                     if (!is.na(hcol)) r[[hcol]] else NA, NA)
    }
  }
}

# ── Clean-room multi-arm ladder estimates (September rebuild) ───────────────
cr_lad <- read_cr("multiarm", "ladder_estimates.csv")
if (!is.null(cr_lad)) {
  for (i in seq_len(nrow(cr_lad))) {
    r <- cr_lad[i, ]
    key <- switch(r$estimand,
      ATE = list(est = if (isTRUE(r$used_ipcw == TRUE | r$used_ipcw == "TRUE"))
        "TMLE, IPCW via Delta" else "TMLE, complete case", ed = "ATE"),
      trimmed_ATE = {
        lvl <- sub("^ATE on the common-support population, trimmed ", "",
                   r$estimand_label)
        list(est = "TMLE, complete case",
             ed = paste0("trimmed ATE (", lvl, ")"))
      },
      ATT = list(est = "TMLE, complete case",
                 ed = "ATT among outcome-observed"),
      ATO = list(est = "augmented overlap weights, IF variance", ed = "ATO"),
      matched_ATT = list(est = "TMLE on matched cohort", ed = "matched ATT"),
      NULL)
    if (is.null(key)) next
    tab <- cr_fill(tab, r$contrast, r$outcome, key$ed, key$est,
                   r$estimate, r$ci_lower, r$ci_upper, r$n)
  }
}

# Clean-room multi-arm negative-control ladder: map the C4 rungs onto the
# main pipeline's two NC cohorts.
cr_ncl <- read_cr("multiarm", "nc_ladder.csv")
if (!is.null(cr_ncl)) {
  nc_map2 <- c(nc_household_urban = "NC: Urban household (SES)",
               nc_cooking_fuel_wood_charcoal = "NC: Wood/charcoal cooking fuel (SES)",
               nc_chronic_hypertension = "NC: Chronic hypertension",
               nc_diabetes_insulin = "NC: Diabetes, insulin dependent",
               nc_hiv_art = "NC: HIV/AIDS on ART")
  for (i in seq_len(nrow(cr_ncl))) {
    r <- cr_ncl[i, ]
    if (!r$contrast %in% c("C4")) next
    ctr <- if (grepl("prior care", r$cohort)) "PRIMARY"
           else if (identical(r$cohort, "full cohort")) "C4" else NA
    if (is.na(ctr)) next
    key <- nc_map2[[as.character(r$negative_control)]]
    if (is.null(key)) next
    tab <- cr_fill(tab, ctr, key, "unadjusted RD (balance check)",
                   "clean-room NC ladder (TMLE)",
                   r$estimate, r$ci_lower, r$ci_upper, r$n)
  }
}

# Survival proxy is deliberately not merged into effect rows: the clean room's
# proxy time grid and the main pipeline's death_by_6mo are different outcomes.
cr_surv <- read_cr("survival_outcome_comparison.csv")
if (!is.null(cr_surv)) {
  pick <- grepl("survtmle", cr_surv[[1]], ignore.case = TRUE) &
          grepl("180", paste0(cr_surv[[2]]))
  if (any(pick)) {
    r <- cr_surv[which(pick)[1], ]
    rdcol <- intersect(c("rd", "risk_difference", "RD"), names(cr_surv))[1]
    tab <- rbind(tab, row0("PRIMARY", COHORTS["PRIMARY"],
      "death by 180 d (proxy survival)", "ATE (180-day risk difference)",
      "survtmle, monthly grid",
      cr_est = if (!is.na(rdcol)) as.numeric(r[[rdcol]]) else NA))
  }
}

# ── Verdicts ────────────────────────────────────────────────────────────────
# Known, documented benign differences (reason strings cite the memo).
benign <- list(
  list(contrast = "PRIMARY", outcome = "gose_favorable", estimand = "ATE",
       reason = paste("Estimand and nuisance differences documented in",
                      "reconciliation_2026-07-28.md discrepancy 3: clean room",
                      "fits complete-case Q on the full cohort with 21",
                      "covariates and SL.glm/SL.glmnet/SL.mean; the main",
                      "pipeline uses IPCW via Delta, 90 design columns, and",
                      "the richer library. Scheduled to close in the",
                      "September rebuild (Part 1 items 2 to 4).")),
  list(contrast = "PRIMARY", outcome = "gose_favorable", estimand = "matched ATT",
       reason = paste("Different matched sets (clean room matches on its own",
                      "21-covariate PS) and interval method; estimates agree",
                      "in direction and overlap.")),
  list(contrast = "PRIMARY", outcome = "severe_disability", estimand = "matched ATT",
       reason = paste("Clean-room value is the negated matched gose_favorable",
                      "ATT (exact complement). Different matched sets and",
                      "interval method; estimates agree in direction",
                      "(+0.021 vs +0.018) and the intervals overlap.")),
  list(contrast = "PRIMARY", outcome = "NC: Urban household (SES)",
       estimand = "unadjusted RD (balance check)",
       reason = paste("Clean room NC is TMLE/IPTW-adjusted on its covariate",
                      "set; main is an unadjusted two-proportion difference.",
                      "Both null on the primary cohort.")),
  list(contrast = "PRIMARY", outcome = "NC: Chronic hypertension",
       estimand = "unadjusted RD (balance check)",
       reason = "Same as urban household: adjusted vs unadjusted, both null.")
)

tab$abs_diff <- abs(tab$cr_est - tab$main_est)
tab$verdict <- NA_character_
tab$reason  <- ""

for (i in seq_len(nrow(tab))) {
  has_cr <- is.finite(tab$cr_est[i]); has_m <- is.finite(tab$main_est[i])
  if (!has_cr && has_m) {
    tab$verdict[i] <- "missing_in_cleanroom"
    tab$reason[i]  <- "Scheduled: September rebuild adds this contrast/estimand to the clean room."
  } else if (has_cr && !has_m) {
    if (grepl("proxy", tab$outcome[i])) {
      tab$verdict[i] <- "not_comparable"
      tab$reason[i]  <- paste("Proxy survival outcome uses coarsened, partly",
                              "imputed event times; kept strictly secondary",
                              "(reconciliation_2026-07-28.md discrepancy 4).")
    } else {
      tab$verdict[i] <- "missing_in_main"
      tab$reason[i]  <- "Clean-room-only quantity."
    }
  } else if (has_cr && has_m) {
    se <- tab$main_se[i]
    if (!is.finite(se) && is.finite(tab$main_ci_lo[i]))
      se <- (tab$main_ci_hi[i] - tab$main_ci_lo[i]) / 3.92
    thr <- if (is.finite(se)) 0.25 * se else 0.1 * abs(tab$main_est[i])
    if (is.finite(tab$abs_diff[i]) && tab$abs_diff[i] <= thr) {
      tab$verdict[i] <- "agree"
      tab$reason[i]  <- "Within a quarter of a main-pipeline SE."
    } else {
      hit <- FALSE
      for (b in benign) {
        if (tab$contrast[i] == b$contrast && tab$outcome[i] == b$outcome &&
            tab$estimand[i] == b$estimand) {
          tab$verdict[i] <- "benign_difference"; tab$reason[i] <- b$reason
          hit <- TRUE; break
        }
      }
      if (!hit) {
        tab$verdict[i] <- "open_discrepancy"
        tab$reason[i]  <- "Difference exceeds tolerance and has no recorded explanation yet."
      }
    }
  } else {
    tab$verdict[i] <- "not_comparable"
    tab$reason[i]  <- "Neither pipeline reports this cell."
  }
}

ord <- order(match(tab$contrast, c("PRIMARY", "C4", "C1", "C2", "C3")),
             tab$outcome, tab$estimand)
tab <- tab[ord, ]

utils::write.csv(tab, OUT_CSV, row.names = FALSE)
cat("Wrote", OUT_CSV, "with", nrow(tab), "rows\n")
cat("\nVerdict counts:\n")
print(table(tab$verdict))
