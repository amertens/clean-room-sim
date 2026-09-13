# ============================================================
# Multi-arm Stage 1: cohort, outcomes, wide W, and five locks
# ============================================================
# Builds the known-transport cohort (all modes, expected n = 8,323), the ten
# study outcomes, the wide design matrix, and five analysis locks: the four
# multi-arm contrasts (C1 Rescue.Co vs everyone else, the protocol
# comparator; C2 Rescue.Co vs non-ambulance; C3 any ambulance vs
# non-ambulance; C4 Rescue.Co vs other ambulance, transfers included) and
# the primary cohort (transfers and prior care excluded).
#
# Construction rules follow the main pipeline
# (../rescueCo_analysis/analysis/02_define_exposure_outcomes_covariates.R and
# 00_config.R), reimplemented here, not copied: factor expansion with rare
# levels under 2 percent pooled (Missing included), checkbox items under
# 1 percent pooled into <family>_other_rare, missingness indicators above
# 1 percent, near-constant binaries excluded, all prevalence filters
# computed on the primary cohort, no hosp_code (an outcome proxy), and the
# derived timing covariates (months since start, hour sine and cosine,
# strike windows, weekend).
#
# OUTCOME BLINDING IS PHYSICAL HERE. The ten outcome columns are written to
# multiarm_outcomes.rds and are ABSENT from every lock: each lock carries a
# masked placeholder for its primary outcome (death_by_6mo) and nothing
# else outcome-shaped. The design stage (script 11) and the support
# simulation (script 12) run on the locks; only the estimation stage
# (script 13) joins the outcomes back by patient id.
#
# NO treatment-outcome association is computed in this script.
# ============================================================

library(yaml)
library(haven)
suppressMessages(library(cleanTMLE))

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                 "bootstrap.R"))
source("rescueCo/R/utils.R")
cfg <- load_cr_config()
cr_log("=== Multi-arm Stage 1: cohort, outcomes, W, locks ===")

OUT <- "rescueCo/results/multiarm"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

as_num <- function(x) suppressWarnings(as.numeric(x))

# ------------------------------------------------------------
# A. Load raw data
# ------------------------------------------------------------
reg <- read_dta(cfg$data$trauma_registry_dta)
fu  <- read_dta(cfg$data$sixmo_followup_dta)
cr_log(paste("Registry:", nrow(reg), "x", ncol(reg),
             "; follow-up:", nrow(fu), "x", ncol(fu)))

# ------------------------------------------------------------
# B. Exposure and transport mode (main analysis/02 sections A-B)
# ------------------------------------------------------------
reg$A <- as.integer(as_num(reg$ems_ambulance))
mode_num <- as_num(reg$hosp_transport)             # 1 = Ambulance
idx_fill <- is.na(reg$A) & !is.na(mode_num) & mode_num == 1
reg$A[idx_fill] <- 0L
cr_log(paste("Filled A = 0 for", sum(idx_fill),
             "ambulance arrivals with a blank ems_ambulance flag"))
rid <- trimws(as.character(reg$rescue_id))
rid[rid %in% c("", "NA")] <- NA_character_
idx_rid <- is.na(reg$A) & !is.na(rid)
if (any(idx_rid)) {
  reg$A[idx_rid] <- 1L
  cr_log(paste("Filled A = 1 for", sum(idx_rid),
               "patients with a Rescue.Co ID but no checkbox"))
}
reg$rescueco_id_present <- !is.na(rid)

mode_present <- !is.na(reg$hosp_transport) &
  nzchar(trimws(as.character(reg$hosp_transport)))
reg$transport3 <- ifelse(!is.na(reg$A) & reg$A == 1, "Rescue.Co ambulance",
                  ifelse(!is.na(reg$A) & reg$A == 0, "Non-Rescue.Co ambulance",
                  ifelse(mode_present, "Non-ambulance", NA_character_)))
tt <- table(reg$transport3, useNA = "ifany")
cr_log(paste("transport3:", paste(names(tt), tt, sep = "=", collapse = ", ")))

# Reference counts from the main pipeline (run_log of analysis/42):
# Non-ambulance 6126, Non-Rescue.Co ambulance 1158, Rescue.Co 1039.
ref_counts <- c(`Non-ambulance` = 6126, `Non-Rescue.Co ambulance` = 1158,
                `Rescue.Co ambulance` = 1039)
for (nmm in names(ref_counts)) {
  got <- sum(reg$transport3 == nmm, na.rm = TRUE)
  if (got != ref_counts[[nmm]])
    cr_log(paste("NOTE: transport3 level", nmm, "has", got,
                 "patients here vs", ref_counts[[nmm]],
                 "in the main pipeline; record in the reconciliation."))
}

# Transfers and prior care (main analysis/02: inter_facility_transfer code 99
# treated as not a transfer; prior_care == 1 is Yes).
xfer_num <- as_num(reg$inter_facility_transfer)
reg$is_transfer <- as.integer(!is.na(xfer_num) & xfer_num == 1)
pc_num <- as_num(reg$prior_care)
reg$sought_prior_care <- as.integer(!is.na(pc_num) & pc_num == 1)
reg$pre_hospital_care <- as.integer(reg$is_transfer == 1 |
                                      reg$sought_prior_care == 1)
reg$in_primary_cohort <- !is.na(reg$A) & reg$pre_hospital_care == 0
cr_log(paste("Primary cohort:", sum(reg$in_primary_cohort), "(",
             sum(reg$in_primary_cohort & reg$A == 1), "Rescue.Co /",
             sum(reg$in_primary_cohort & reg$A == 0), "other ambulance )"))

# ------------------------------------------------------------
# C. Datetimes and the three time outcomes (main analysis/02 section E)
# ------------------------------------------------------------
mk_dt <- function(date_col, time_col) {
  d <- suppressWarnings(as.Date(as.character(date_col)))
  t <- trimws(as.character(time_col))
  t[!grepl("^\\d{1,2}:\\d{2}", t)] <- NA
  out <- as.POSIXct(paste(d, t), format = "%Y-%m-%d %H:%M", tz = "UTC")
  out[is.na(d) | is.na(t)] <- NA
  out
}
reg$injury_datetime  <- mk_dt(reg$inj_date, reg$inj_time)
reg$arrival_datetime <- mk_dt(reg$arrival_date, reg$arrival_time)

MAX_MIN <- 30 * 24 * 60   # 30 days, the main pipeline's plausibility cap
guard_minutes <- function(x, label) {
  bad <- !is.na(x) & (x < 0 | x > MAX_MIN)
  if (any(bad))
    cr_log(paste(label, "discarded (out of [0, 30d]):", sum(bad),
                 "( Rescue.Co", sum(bad & reg$A %in% 1),
                 ", other amb", sum(bad & reg$A %in% 0),
                 ", non-amb", sum(bad & is.na(reg$A)), ")"))
  x[bad] <- NA_real_
  x
}
reg$time_to_arrival_minutes <- guard_minutes(as.numeric(difftime(
  reg$arrival_datetime, reg$injury_datetime, units = "mins")),
  "time_to_arrival")

# Definitive interventions: first timestamp across the definitive-care set
# (main analysis/00_config.R DEFINITIVE_TX_MAP; analgesic excluded as
# near-universal).
def_map <- list(
  c("treat_blooddate", "treat_bloodtime"),
  c("treat_tranexaciddate", "treat_tranexacidtime"),
  c("treat_antibioticdate", "treat_antibiotictime"),
  c("treat_splitdate", "treat_splittime_2"),
  c("treat_reductdate", "treat_reducttime"),
  c("treat_debridedate", "treat_debridetime"),
  c("treat_anticoagulationdate", "treat_anticoagulationtime"),
  c("treat_antitetanusdate", "treat_antitetanustime_2"),
  c("treat_colloiddate", "treat_colloidtime"),
  c("treat_crystalloiddate", "treat_crystalloidtime"),
  c("treat_ppidate", "treat_ppitime"),
  c("treat_nasogastricdate", "treat_nasogastrictime_2"),
  c("treat_catheterdate", "treat_cathetertime_2"))
def_epochs <- sapply(def_map, function(p) {
  if (!all(p %in% names(reg))) return(rep(NA_real_, nrow(reg)))
  as.numeric(mk_dt(reg[[p[1]]], reg[[p[2]]]))
})
reg$first_def_tx_epoch <- apply(def_epochs, 1, function(r)
  if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE))
first_dt <- as.POSIXct(reg$first_def_tx_epoch, origin = "1970-01-01",
                       tz = "UTC")
reg$time_to_def_intervention_minutes <- guard_minutes(as.numeric(difftime(
  first_dt, reg$injury_datetime, units = "mins")), "time_to_def_intervention")
reg$time_arrival_to_intervention_minutes <- guard_minutes(as.numeric(difftime(
  first_dt, reg$arrival_datetime, units = "mins")),
  "time_arrival_to_intervention")

# ------------------------------------------------------------
# D. Follow-up merge and GOSE (walk-down as in scripts/01, reconciliation
#    rules as in main analysis/02 section D)
# ------------------------------------------------------------
fu_id <- cfg$merge$followup_id
fu <- fu[!duplicated(fu[[fu_id]]), ]
fu_keep <- intersect(c(fu_id, "consciousness", "independence_home",
                       "need_assistance", "shop", "travel_locally", "work",
                       "majorchange_work", "social_restriction",
                       "change_personality", "change_personality_new",
                       "dailylife", "dailylife_new", "patient_alive", "dod"),
                     names(fu))
dat <- merge(reg, fu[, fu_keep, drop = FALSE],
             by.x = cfg$merge$registry_id, by.y = fu_id, all.x = TRUE)
cr_log(paste("Merged follow-up:", sum(!is.na(dat$patient_alive)),
             "patients with 6-month vital status"))

disp_num <- as_num(dat$disposition)
# death_in_hospital: disposition 5 (Died) ONLY. Code 8, discharged home to
# die, is a prognosis at discharge, not a death (main analysis/02: of 8 such
# patients registry-wide, 5 were confirmed alive at 6 months).
dat$death_in_hospital <- ifelse(!is.na(disp_num) & disp_num == 5, 1L,
                                ifelse(!is.na(disp_num), 0L, NA_integer_))

# GOSE walk-down (same item logic as scripts/01_stage1_build_cohort.R).
dat$gose <- NA_integer_
dead_fu <- !is.na(dat$patient_alive) & as_num(dat$patient_alive) == 0
dat$gose[dead_fu] <- 1L
alive <- !is.na(dat$patient_alive) & as_num(dat$patient_alive) == 1
no_consc <- alive & !is.na(dat$consciousness) & as_num(dat$consciousness) == 0
dat$gose[no_consc] <- 2L
consc <- alive & !is.na(dat$consciousness) & as_num(dat$consciousness) == 1
needs_help <- consc & !is.na(dat$independence_home) &
  as_num(dat$independence_home) == 1
dat$gose[needs_help & !is.na(dat$need_assistance) &
           as_num(dat$need_assistance) == 1] <- 3L
dat$gose[needs_help & (is.na(dat$need_assistance) |
                         as_num(dat$need_assistance) == 0)] <- 4L
indep <- consc & (is.na(dat$independence_home) |
                    as_num(dat$independence_home) == 0)
cant <- indep & ((!is.na(dat$shop) & as_num(dat$shop) == 0) |
                 (!is.na(dat$travel_locally) &
                    as_num(dat$travel_locally) == 0) |
                 (!is.na(dat$work) & as_num(dat$work) == 0))
dat$gose[cant & is.na(dat$gose)] <- 5L
reduced <- indep & !cant &
  ((!is.na(dat$majorchange_work) & as_num(dat$majorchange_work) == 1) |
   (!is.na(dat$social_restriction) & as_num(dat$social_restriction) >= 1))
dat$gose[reduced & is.na(dat$gose)] <- 6L
pers <- indep & !cant & !reduced &
  ((!is.na(dat$change_personality) & as_num(dat$change_personality) == 1 &
    !is.na(dat$change_personality_new) &
    as_num(dat$change_personality_new) == 1) |
   (!is.na(dat$dailylife) & as_num(dat$dailylife) == 1 &
    !is.na(dat$dailylife_new) & as_num(dat$dailylife_new) == 1))
dat$gose[pers & is.na(dat$gose)] <- 7L
dat$gose[consc & is.na(dat$gose)] <- 8L

# Reconciliation rules (main analysis/02):
# 1. Back-fill: a hospital death with no follow-up GOSE is GOSE 1.
bf <- is.na(dat$gose) & dat$death_in_hospital %in% 1L
dat$gose[bf] <- 1L
# 2. Contradiction override: a hospital death with a non-1 follow-up GOSE is
#    GOSE 1 (flagged).
cx <- dat$death_in_hospital %in% 1L & !is.na(dat$gose) & dat$gose != 1L
dat$gose_death_reconciled <- as.integer(cx)
dat$gose[cx] <- 1L
if (any(cx)) cr_log(paste("GOSE contradiction override on", sum(cx), "rows"))
# 3. Late-death release: a death dated more than 180 days after arrival, and
#    not an in-hospital death, is not a 6-month death; its GOSE 1 is
#    released. A death dated before arrival is impossible: date discarded,
#    death flag retained.
dod_d <- suppressWarnings(as.Date(as.character(dat$dod)))
arr_d <- as.Date(dat$arrival_datetime)
days_to_death <- as.numeric(dod_d - arr_d)
late <- !is.na(days_to_death) & days_to_death > 180 &
  !(dat$death_in_hospital %in% 1L)
dat$death_date_implausible <- as.integer(!is.na(days_to_death) &
                                           days_to_death < 0)
if (any(late)) {
  cr_log(paste("Late-death release on", sum(late),
               "rows (death after 180 days)"))
  dat$gose[late & dat$gose %in% 1L] <- NA_integer_
}

# death_by_6mo: union of follow-up vital status, GOSE 1, hospital death,
# minus the late-death releases.
dat$death_by_6mo <- ifelse(
  (dead_fu & !late) | dat$gose %in% 1L | dat$death_in_hospital %in% 1L, 1L,
  ifelse(!is.na(dat$patient_alive) | !is.na(dat$gose), 0L, NA_integer_))
dat$severe_disability <- ifelse(!is.na(dat$gose),
                                as.integer(dat$gose <= 4L), NA_integer_)
dat$gose_favorable <- ifelse(!is.na(dat$gose),
                             as.integer(dat$gose > 4L), NA_integer_)

# Transfers out (main analysis/02: disposition 7; higher care needs the
# disposition_transfer reason checkbox 1).
hc <- as_num(dat$disposition_transfer___1)
dat$transferred_out <- ifelse(!is.na(disp_num) & disp_num == 7, 1L,
                              ifelse(!is.na(disp_num), 0L, NA_integer_))
dat$transferred_out_higher_care <- ifelse(
  !is.na(disp_num) & disp_num == 7 & !is.na(hc) & hc == 1, 1L,
  ifelse(!is.na(disp_num), 0L, NA_integer_))

OUTCOME_VARS <- c("death_by_6mo", "death_in_hospital", "severe_disability",
                  "gose_favorable", "gose", "transferred_out",
                  "transferred_out_higher_care", "time_to_arrival_minutes",
                  "time_to_def_intervention_minutes",
                  "time_arrival_to_intervention_minutes")
for (v in OUTCOME_VARS)
  cr_log(sprintf("  outcome %-38s non-NA %5d  mean %.4f", v,
                 sum(!is.na(dat[[v]])), mean(dat[[v]], na.rm = TRUE)))

# ------------------------------------------------------------
# E. Negative controls (registered as cohort columns BEFORE any W filter,
#    so no variance screen can drop them; main analysis/18)
# ------------------------------------------------------------
lab_chr <- function(x) tolower(trimws(as.character(haven::as_factor(x))))
dat$nc_household_urban <- as.integer(grepl("urban", lab_chr(dat$household_area)))
dat$nc_cooking_fuel_wood_charcoal <- as.integer(
  grepl("wood|charcoal|firewood", lab_chr(dat$cooking_fuel)))
dat$nc_chronic_hypertension <- as.integer(as_num(dat$chronic_illness___2) == 1)
dat$nc_diabetes_insulin <- as.integer(as_num(dat$chronic_illness___0) == 1)
dat$nc_hiv_art <- as.integer(as_num(dat$chronic_illness___3) == 1)
NC_VARS <- c("nc_household_urban", "nc_cooking_fuel_wood_charcoal",
             "nc_chronic_hypertension", "nc_diabetes_insulin", "nc_hiv_art")

# ------------------------------------------------------------
# F. Derived timing covariates (main analysis/02 section F)
# ------------------------------------------------------------
index_date <- as.Date(dat$injury_datetime)
index_date[is.na(index_date)] <- as.Date(dat$arrival_datetime)[is.na(index_date)]
dat$months_since_start <- as.numeric(index_date -
                                       min(index_date, na.rm = TRUE)) / 30.44
ih <- suppressWarnings(as.integer(format(dat$injury_datetime, "%H")))
ih[is.na(ih)] <- suppressWarnings(
  as.integer(format(dat$arrival_datetime, "%H")))[is.na(ih)]
dat$hour_sin <- sin(2 * pi * ih / 24)
dat$hour_cos <- cos(2 * pi * ih / 24)
wd <- as.POSIXlt(index_date)$wday
dat$arrival_weekend <- as.integer(wd %in% c(0, 6))

# Strike windows (main analysis/00_config.R STRIKE_WINDOWS).
strikes <- list(c("2023-09-13", "2023-10-17"), c("2024-03-14", "2024-05-08"),
                c("2025-02-27", "2025-04-08"), c("2025-07-07", "2025-08-18"),
                c("2025-12-19", "2026-01-03"), c("2026-01-17", "2026-03-27"))
dat$strike_period <- 0L
for (w in strikes)
  dat$strike_period[!is.na(index_date) & index_date >= as.Date(w[1]) &
                      index_date <= as.Date(w[2])] <- 1L

# Blood pressure split with plausibility screens (main analysis/02:
# sbp in [50, 260], dbp in [20, 200], pulse pressure >= 10, both forced
# missing together so one indicator is emitted).
bp <- gsub("-", "/", as.character(dat$vitals_1_bp))
parts <- strsplit(bp, "/", fixed = TRUE)
dat$sbp <- as_num(vapply(parts, function(p) if (length(p) >= 1) p[1] else NA_character_,
                         character(1)))
dat$dbp <- as_num(vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_,
                         character(1)))
bad_bp <- (!is.na(dat$sbp) & (dat$sbp < 50 | dat$sbp > 260)) |
  (!is.na(dat$dbp) & (dat$dbp < 20 | dat$dbp > 200)) |
  (!is.na(dat$sbp) & !is.na(dat$dbp) & (dat$sbp - dat$dbp) < 10)
dat$sbp[bad_bp] <- NA_real_
dat$dbp[bad_bp] <- NA_real_
dat$dbp[is.na(dat$sbp)] <- NA_real_
dat$sbp[is.na(dat$dbp)] <- NA_real_

# n_regions_injured: AIS regions above zero; NA when none recorded.
ais_cols <- c("general_ais", "face_ais", "ais_head", "torso_ais", "abd_ais",
              "extremities_ais")
ais <- sapply(ais_cols, function(v) as_num(dat[[v]]))
dat$n_regions_injured <- apply(ais, 1, function(r)
  if (all(is.na(r))) NA_real_ else sum(r > 0, na.rm = TRUE))

# Head-injury severity (main analysis/02: mech_inj 9/10/11 are mild,
# moderate, severe head injury; one ordinal keeps the rare severe levels in
# the design where the prevalence filter would drop their indicators).
mi <- function(k) ifelse(is.na(as_num(dat[[paste0("mech_inj___", k)]])), 0,
                         as_num(dat[[paste0("mech_inj___", k)]]))
dat$head_injury_severity <- pmax(1 * (mi(9) == 1), 2 * (mi(10) == 1),
                                 3 * (mi(11) == 1))

# Sub-county normalisation (main analysis/02 .norm_sc).
norm_sc <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v <- gsub("[^A-Z ]", "", v)
  v <- gsub("[[:space:]]+", " ", v)
  v <- sub(" COUNTY$", "", v)
  v <- sub("^KAMKUNJI$", "KAMUKUNJI", v)
  v[v %in% c("", "NA", "UNKNOWN", "NOT APPLICABLE")] <- NA
  v
}
dat$inj_sub_county_norm <- norm_sc(dat$inj_sub_county)

# ------------------------------------------------------------
# G. Wide design matrix (main analysis/02 sections G-H rules)
# ------------------------------------------------------------
ref <- which(dat$in_primary_cohort)   # all filters anchored on the 1,616
MISS_IND_MIN <- 0.01
RARE_LEVEL_PREV <- 0.02
RARE_BINARY_PREV <- 0.01
MIN_DESIGN_PREV <- 0.02

CONTINUOUS <- c("age", "inj_distance", "eiss_ais", "n_regions_injured",
                "vitals_1_hr", "eyes_gcs", "verbal_gcs", "motor_gcs",
                "sbp", "dbp", "months_since_start", "hour_sin", "hour_cos",
                "head_injury_severity")
BINARY <- c("strike_period", "arrival_weekend")
NOMINAL <- c("sex", "marital_status", "education", "television", "cable_tv",
             "mixer", "cooking_fuel", "occupation", "home_ownership",
             "inj_place", "inj_activity", "alcohol_use", "inj_alcohol_use",
             "inj_sub_county_norm", "tobacco_use", "intent",
             "penetrating_wound", "pulse")
CHECKBOX <- c("mech_inj", "medical_history", "chronic_illness",
              "prior_surgery", "dispo_payment")
# household_area, cellphone, signs_of_life: near-constant, excluded up front
# (main analysis/00_config.R NEAR_CONSTANT_EXCLUDED). hosp_code: outcome
# proxy, never in the adjustment set (main README).

enc <- list()
# Continuous: missingness indicator above 1 percent (in ref), then overall
# median imputation (not within arm: within-arm imputation writes treatment
# information into W).
for (v in CONTINUOUS) {
  xv <- as_num(dat[[v]])
  miss <- is.na(xv)
  if (mean(miss[ref]) >= MISS_IND_MIN)
    enc[[paste0(v, "_missing")]] <- as.integer(miss)
  med <- stats::median(xv[ref], na.rm = TRUE)
  if (is.na(med)) med <- stats::median(xv, na.rm = TRUE)
  xv[miss] <- med
  enc[[v]] <- xv
}
for (v in BINARY) enc[[v]] <- as.integer(dat[[v]])

# Nominal: label, blanks to Missing, rare levels (under 2 percent in ref,
# Missing included) pooled to Other, levels ordered by ref frequency so the
# modal level is the omitted reference.
for (v in NOMINAL) {
  f <- as.character(haven::as_factor(dat[[v]]))
  f[is.na(f) | !nzchar(trimws(f))] <- "Missing"
  tab_ref <- table(f[ref])
  rare <- names(tab_ref)[tab_ref / length(ref) < RARE_LEVEL_PREV]
  f[f %in% rare] <- "Other"
  lev <- names(sort(table(f[ref]), decreasing = TRUE))
  lev <- c(lev, setdiff(unique(f), lev))
  fac <- factor(f, levels = lev)
  if (nlevels(fac) < 2L) next
  mm <- stats::model.matrix(~ fac - 1)[, -1, drop = FALSE]
  colnames(mm) <- paste0(v, "_", make.names(sub("^fac", "", colnames(mm))))
  for (j in seq_len(ncol(mm))) enc[[colnames(mm)[j]]] <- as.numeric(mm[, j])
}

# Checkbox families: keep items, pool the under-1-percent items per family.
for (fam in CHECKBOX) {
  cols <- grep(paste0("^", fam, "___[0-9]+$"), names(dat), value = TRUE)
  if (!length(cols)) next
  fam_rare <- rep(0L, nrow(dat))
  n_rare <- 0L
  for (cc in cols) {
    xv <- as.integer(as_num(dat[[cc]]) == 1)
    xv[is.na(xv)] <- 0L
    prev <- mean(xv[ref])
    nm2 <- sub("___", "_", cc)
    if (prev < RARE_BINARY_PREV) {
      fam_rare <- pmax(fam_rare, xv)
      n_rare <- n_rare + 1L
    } else enc[[nm2]] <- xv
  }
  if (n_rare >= 2L) enc[[paste0(fam, "_other_rare")]] <- fam_rare
}

W_full <- as.data.frame(enc, check.names = TRUE)

# Filters, all on the ref cohort: zero variance, exact duplicates,
# near-constant binaries (< 2 percent or > 98 percent prevalence).
sd_ref <- vapply(W_full, function(x) stats::sd(x[ref]), numeric(1))
drop_zv <- names(W_full)[!is.finite(sd_ref) | sd_ref < 1e-10]
W_full <- W_full[, setdiff(names(W_full), drop_zv), drop = FALSE]
dup <- duplicated(t(as.matrix(W_full[ref, ])))
drop_dup <- names(W_full)[dup]
W_full <- W_full[, setdiff(names(W_full), drop_dup), drop = FALSE]
is_bin <- vapply(W_full, function(x) all(x %in% c(0, 1)), logical(1))
prev <- colMeans(W_full[ref, , drop = FALSE])
drop_nzv <- names(W_full)[is_bin & (prev < MIN_DESIGN_PREV |
                                      prev > 1 - MIN_DESIGN_PREV)]
W <- W_full[, setdiff(names(W_full), drop_nzv), drop = FALSE]
cr_log(paste("Design matrix:", ncol(W), "columns after filters (dropped",
             length(drop_zv), "zero-variance,", length(drop_dup),
             "duplicate,", length(drop_nzv), "near-constant)"))
utils::write.csv(
  data.frame(column = names(W),
             prevalence_ref = round(colMeans(W[ref, , drop = FALSE]), 4)),
  file.path(OUT, "design_columns.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(dropped = c(drop_zv, drop_dup, drop_nzv),
             reason = c(rep("zero_variance", length(drop_zv)),
                        rep("duplicate", length(drop_dup)),
                        rep("near_constant", length(drop_nzv)))),
  file.path(OUT, "design_columns_dropped.csv"), row.names = FALSE)

# ------------------------------------------------------------
# H. Assemble the analysis frame, the outcome store, and the locks
# ------------------------------------------------------------
keep_transport <- !is.na(dat$transport3)
core <- data.frame(patient_id = dat[[cfg$merge$registry_id]],
                   transport3 = dat$transport3,
                   is_transfer = dat$is_transfer,
                   sought_prior_care = dat$sought_prior_care,
                   pre_hospital_care = dat$pre_hospital_care,
                   in_primary_cohort = dat$in_primary_cohort,
                   rescueco_id_present = dat$rescueco_id_present,
                   stringsAsFactors = FALSE)
core <- cbind(core, dat[, NC_VARS], W)
core <- core[keep_transport, , drop = FALSE]
outcomes <- cbind(data.frame(patient_id = dat[[cfg$merge$registry_id]]),
                  dat[, OUTCOME_VARS])
outcomes <- outcomes[keep_transport, , drop = FALSE]
cr_log(paste("Known-transport cohort:", nrow(core)))

saveRDS(outcomes, file.path(OUT, "multiarm_outcomes.rds"))
saveRDS(core, file.path(OUT, "multiarm_cohort.rds"))

# The locks: a masked placeholder for the primary outcome. The real outcome
# columns live only in multiarm_outcomes.rds until Stage 4.
core$death_by_6mo <- NA_real_
covariate_cols <- names(W)

contrasts <- list(
  C1 = list(treated = "Rescue.Co ambulance",
            control = c("Non-Rescue.Co ambulance", "Non-ambulance"),
            label = "Rescue.Co vs everyone else (protocol comparator)"),
  C2 = list(treated = "Rescue.Co ambulance", control = "Non-ambulance",
            label = "Rescue.Co vs non-ambulance"),
  C3 = list(treated = c("Rescue.Co ambulance", "Non-Rescue.Co ambulance"),
            control = "Non-ambulance",
            label = "Any ambulance vs non-ambulance"),
  C4 = list(treated = "Rescue.Co ambulance",
            control = "Non-Rescue.Co ambulance",
            label = "Rescue.Co vs other ambulance (transfers included)"))

locks <- create_contrast_locks(
  core, "transport3", contrasts,
  outcome = "death_by_6mo", covariates = covariate_cols,
  sl_library = c("SL.glm", "SL.mean"),   # per-fit libraries set downstream
  seed = as.integer(cfg$seed %||% 42L),
  negative_controls = NC_VARS)

# PRIMARY: the 1,616 cohort as its own lock, same W, same controls.
prim <- core[core$in_primary_cohort, , drop = FALSE]
prim$.A_contrast <- as.integer(prim$transport3 == "Rescue.Co ambulance")
lock_p <- create_analysis_lock(
  prim, ".A_contrast", "death_by_6mo", covariate_cols,
  sl_library = c("SL.glm", "SL.mean"), seed = as.integer(cfg$seed %||% 42L))
lock_p$contrast <- list(name = "PRIMARY",
                        label = "Rescue.Co vs other ambulance, transfers and prior care excluded",
                        treated = "Rescue.Co ambulance",
                        control = "Non-Rescue.Co ambulance",
                        source_column = "transport3")
for (nc in NC_VARS) lock_p <- define_negative_control(lock_p, nc)
locks$PRIMARY <- lock_p
class(locks) <- c("contrast_locks", "list")

# Ladder, primary outcome decision, and pending amendments on every lock.
for (cn in names(locks)) {
  lk <- locks[[cn]]
  lk$.outcome_masked <- TRUE
  lk <- declare_estimand_ladder(lk, primary = "ATE",
                                fallbacks = c("trimmed_ATE", "ATT", "ATO"),
                                trigger = "SEVERE",
                                bias_to_null_floor = 0.0035)
  lk <- cleanTMLE:::.log_design_decision(lk, "outcome",
    "Primary outcome death_by_6mo (about 100 events in the primary cohort); severe_disability replaces gose_good so the rare category is the event; IPCW via Delta is the primary missing-outcome handling, complete case the sensitivity.")
  lk <- cleanTMLE:::.log_design_decision(lk, "pending_amendment",
    "Rescue.Co is checking 76 treated patients in the primary cohort without a Rescue.Co ID (meeting 2026-09-10); any reclassification reruns from this stage.")
  lk <- cleanTMLE:::.log_design_decision(lk, "pending_amendment",
    "The study team recalls improving Rescue.Co identification partway through enrolment; part of the rising Rescue.Co share may be ascertainment. Date pending.")
  locks[[cn]] <- lk
}

saveRDS(locks, file.path(OUT, "multiarm_locks.rds"))
utils::write.csv(
  do.call(rbind, lapply(names(locks), function(cn) {
    lk <- locks[[cn]]
    data.frame(contrast = cn, label = lk$contrast$label,
               n = nrow(lk$data), n_treated = sum(lk$data[[lk$treatment]]),
               n_control = sum(lk$data[[lk$treatment]] == 0),
               n_covariates = length(lk$covariates))
  })), file.path(OUT, "lock_summary.csv"), row.names = FALSE)

cr_log("Multi-arm Stage 1 complete.")
for (cn in names(locks))
  cr_log(sprintf("  %-8s n = %5d (%4d treated / %4d control)", cn,
                 nrow(locks[[cn]]$data),
                 sum(locks[[cn]]$data[[locks[[cn]]$treatment]]),
                 sum(locks[[cn]]$data[[locks[[cn]]$treatment]] == 0)))
