# ============================================================
# 21. Design data for the re-execution on cleanTMLE 0.5.0
# ============================================================
# Builds, for one design-stage run, the outcome-free design data: arms,
# cohort flags, enrolment date, the covariates fixed before transport, the
# arrival-state block (sensitivity analyses only), the five negative controls,
# and the steward's two outcome-observation indicators; plus the variable
# dictionary and an aggregate column manifest.
#
# Usage: Rscript rescueCo/scripts/21_build_design_data.R <run>
#   run = "first4m" (patients enrolled 19 February to 18 June 2024, the window
#   of the protocol's own outcome-blind check) or "final" (all patients).
#
# This script never opens the follow-up file and derives no outcome. It reads
# registry columns from a whitelist; disposition, payment, hospital reached,
# distance, arrival times, dispatch and scene times, signs of life, scene care
# and operator name are not read. Encoding follows the unblinded analysis
# (rescueCo_analysis/analysis/02_define_exposure_outcomes_covariates.R).
#
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
cfg <- load_cr_config()
args <- commandArgs(trailingOnly = TRUE)
RUN <- if (length(args)) args[1] else "final"
stopifnot(RUN %in% c("first4m", "final"))
cr_log(paste("=== 21 design data, run", RUN, "==="))

ENROL_START <- as.Date("2024-02-19")
FIRST4M_END <- as.Date("2024-06-18")
as_num <- function(x) if (inherits(x, "haven_labelled")) as.numeric(unclass(x)) else
  suppressWarnings(as.numeric(as.character(x)))

# ------------------------------------------------------------
# A. Registry columns (whitelist)
# ------------------------------------------------------------
REG_SCALAR <- c("patient_id", "ems_ambulance", "rescue_id", "hosp_transport",
                "inter_facility_transfer", "prior_care", "inj_date", "inj_time",
                "arrival_date", "age", "sex", "marital_status", "education",
                "occupation", "home_ownership", "television", "cable_tv", "mixer",
                "cooking_fuel", "household_area", "alcohol_use", "tobacco_use",
                "inj_alcohol_use", "inj_sub_county", "inj_place", "inj_activity",
                "intent", "penetrating_wound", "eiss_ais", "general_ais", "face_ais",
                "ais_head", "torso_ais", "abd_ais", "extremities_ais",
                "vitals_1_bp", "vitals_1_hr", "eyes_gcs", "verbal_gcs", "motor_gcs",
                "pulse")
REG_FAMILIES <- c("mech_inj", "medical_history", "chronic_illness", "prior_surgery")
hdr <- names(haven::read_dta(cfg$data$trauma_registry_dta, n_max = 0))
fam_cols <- unlist(lapply(REG_FAMILIES, function(f)
  grep(paste0("^", f, "___[0-9]+$"), hdr, value = TRUE)))
stopifnot(all(REG_SCALAR %in% hdr))
reg <- haven::read_dta(cfg$data$trauma_registry_dta,
                       col_select = tidyselect::all_of(c(REG_SCALAR, fam_cols)))
reg$patient_id <- as.character(reg$patient_id)
cr_log(paste("Registry:", nrow(reg), "patients;", ncol(reg), "whitelisted columns"))

# ------------------------------------------------------------
# B. Arms (the unblinded analysis's rules; inventory row 6)
# ------------------------------------------------------------
ems <- as_num(reg$ems_ambulance)
rid <- trimws(as.character(reg$rescue_id)); rid[rid %in% c("", "NA")] <- NA
mode <- as_num(reg$hosp_transport)
rescueco <- (!is.na(ems) & ems == 1) | (is.na(ems) & !is.na(rid))
other_op <- !is.na(ems) & ems == 0 & !rescueco
ambulance <- !is.na(mode) & mode == 1
reg$operator_unknown <- as.integer(ambulance & !rescueco & !other_op)
other_means <- !is.na(mode) & !(mode %in% c(1, 99)) & !rescueco & !other_op
reg$arm <- ifelse(rescueco, "rescueco", ifelse(other_op, "other_operator",
                  ifelse(other_means, "other_means", NA_character_)))
reg$walked_in <- as.integer(other_means & mode == 0)
tab <- table(reg$arm, useNA = "ifany")
cr_log(paste("Arms:", paste(names(tab), tab, sep = "=", collapse = ", "),
             "; operator unknown:", sum(reg$operator_unknown)))
expect <- c(rescueco = 1039, other_operator = 1146, other_means = 6118)
for (k in names(expect)) if (sum(reg$arm == k, na.rm = TRUE) != expect[[k]])
  stop("Arm ", k, " differs from the unblinded analysis.", call. = FALSE)

# Rescue.co's cross-reference tier for Rescue.Co patients without an id
# (1, 2, 3 = matched or probable; 5 = no dispatch record). Read without names.
xr_path <- cfg$data$rescue_crossref_xlsx
xr_hdr <- names(readxl::read_excel(xr_path, sheet = "2. Cross-reference", n_max = 0))
xr <- readxl::read_excel(xr_path, sheet = "2. Cross-reference",
                         col_types = ifelse(xr_hdr %in% c("patient_id", "match_status"),
                                            "text", "skip"))
xr_tier <- as.integer(substr(xr$match_status, 1, 1))
reg$crossref_tier <- xr_tier[match(reg$patient_id, trimws(xr$patient_id))]
stopifnot(sum(reg$crossref_tier == 5, na.rm = TRUE) == 14)

# ------------------------------------------------------------
# C. Cohort flags and dates
# ------------------------------------------------------------
xf <- as_num(reg$inter_facility_transfer); pc <- as_num(reg$prior_care)
reg$transfer <- as.integer(!is.na(xf) & xf == 1)
reg$prior_care_elsewhere <- as.integer(!is.na(pc) & pc == 1)
reg$direct_arrival <- as.integer(reg$transfer == 0 & reg$prior_care_elsewhere == 0)
reg$enrolment_date <- as.Date(reg$arrival_date)
stopifnot(min(reg$enrolment_date, na.rm = TRUE) == ENROL_START)
in_run <- if (RUN == "first4m") !is.na(reg$enrolment_date) &
  reg$enrolment_date <= FIRST4M_END else rep(TRUE, nrow(reg))

# ------------------------------------------------------------
# D. Derived covariates fixed before transport
# ------------------------------------------------------------
idx_date <- as.Date(reg$inj_date)
idx_date[is.na(idx_date)] <- reg$enrolment_date[is.na(idx_date)]   # same day for nearly all
reg$months_since_start <- as.numeric(idx_date - ENROL_START) / 30.44
ih <- suppressWarnings(as.integer(sub("^(\\d{1,2}):.*$", "\\1", trimws(as.character(reg$inj_time)))))
ih[!is.na(ih) & (ih < 0 | ih > 23)] <- NA                         # no fallback to arrival hour
reg$hour_sin <- sin(2 * pi * ih / 24)
reg$hour_cos <- cos(2 * pi * ih / 24)
reg$arrival_weekend <- as.integer(as.POSIXlt(idx_date)$wday %in% c(0, 6))   # weekend of injury
strikes <- list(c("2023-09-13", "2023-10-17"), c("2024-03-14", "2024-05-08"),
                c("2025-02-27", "2025-04-08"), c("2025-07-07", "2025-08-18"),
                c("2025-12-19", "2026-01-03"), c("2026-01-17", "2026-03-27"))
reg$strike_period <- 0L
for (w in strikes) reg$strike_period[!is.na(idx_date) & idx_date >= as.Date(w[1]) &
                                       idx_date <= as.Date(w[2])] <- 1L
ais <- sapply(c("general_ais", "face_ais", "ais_head", "torso_ais", "abd_ais",
                "extremities_ais"), function(v) as_num(reg[[v]]))
reg$n_regions_injured <- apply(ais, 1, function(r) if (all(is.na(r))) NA_real_ else sum(r > 0, na.rm = TRUE))
mi <- function(k) { x <- as_num(reg[[paste0("mech_inj___", k)]]); ifelse(is.na(x), 0, x) }
reg$head_injury_severity <- pmax(1 * (mi(9) == 1), 2 * (mi(10) == 1), 3 * (mi(11) == 1))
norm_sc <- function(x) {
  v <- toupper(trimws(as.character(x))); v <- gsub("[^A-Z ]", "", v)
  v <- gsub("[[:space:]]+", " ", v); v <- sub(" COUNTY$", "", v)
  v <- sub("^KAMKUNJI$", "KAMUKUNJI", v)
  v[v %in% c("", "NA", "UNKNOWN", "NOT APPLICABLE")] <- NA
  v
}
reg$inj_sub_county <- norm_sc(reg$inj_sub_county)

# Arrival state (stand-ins for the sensitivity analysis only).
bp <- strsplit(gsub("-", "/", as.character(reg$vitals_1_bp)), "/", fixed = TRUE)
reg$sbp <- as_num(vapply(bp, function(p) if (length(p) >= 1) p[1] else NA_character_, ""))
reg$dbp <- as_num(vapply(bp, function(p) if (length(p) >= 2) p[2] else NA_character_, ""))
bad <- (!is.na(reg$sbp) & (reg$sbp < 50 | reg$sbp > 260)) |
  (!is.na(reg$dbp) & (reg$dbp < 20 | reg$dbp > 200)) |
  (!is.na(reg$sbp) & !is.na(reg$dbp) & reg$sbp - reg$dbp < 10)
reg$sbp[bad | is.na(reg$dbp)] <- NA; reg$dbp[bad | is.na(reg$sbp)] <- NA

# Negative controls (an unchecked box is absent).
lab <- function(x) tolower(trimws(as.character(haven::as_factor(x))))
chk <- function(v) { x <- as_num(reg[[v]]); as.integer(!is.na(x) & x == 1) }
reg$nc_household_urban <- as.integer(grepl("urban", lab(reg$household_area)))
reg$nc_cooking_wood_charcoal <- as.integer(grepl("wood|charcoal|firewood", lab(reg$cooking_fuel)))
reg$nc_hypertension <- chk("chronic_illness___2")
reg$nc_diabetes_insulin <- chk("chronic_illness___0")
reg$nc_hiv_art <- chk("chronic_illness___3")
NC <- c("nc_household_urban", "nc_cooking_wood_charcoal", "nc_hypertension",
        "nc_diabetes_insulin", "nc_hiv_art")

# ------------------------------------------------------------
# E. Encoding, anchored on this run's primary cohort (direct arrivals carried
#    by Rescue.Co or another operator; the 14 without a dispatch record kept
#    for anchoring, as in the unblinded analysis)
# ------------------------------------------------------------
ref <- which(in_run & reg$direct_arrival == 1 & reg$arm %in% c("rescueco", "other_operator"))
cr_log(paste("Encoding anchor (this run's primary cohort):", length(ref), "patients"))
MISS_DROP <- 0.40; MISS_IND <- 0.01; RARE_LEVEL <- 0.02; RARE_BOX <- 0.01; NEAR_CONST <- 0.02

encode <- function(continuous, nominal, families, binary = character()) {
  enc <- list(); src <- character()
  put <- function(nm, x, v) { enc[[nm]] <<- x; src[[nm]] <<- v }
  for (v in continuous) {
    x <- as_num(reg[[v]]); miss <- is.na(x)
    if (mean(miss[ref]) > MISS_DROP) next
    if (mean(miss[ref]) >= MISS_IND) put(paste0(v, "_missing"), as.integer(miss), v)
    med <- stats::median(x[ref], na.rm = TRUE)
    x[miss] <- med; put(v, x, v)
  }
  for (v in binary) put(v, as.integer(reg[[v]]), v)
  for (v in nominal) {
    f <- if (is.character(reg[[v]])) reg[[v]] else as.character(haven::as_factor(reg[[v]]))
    f[is.na(f) | !nzchar(trimws(f))] <- "Missing"
    tr <- table(f[ref]); f[f %in% names(tr)[tr / length(ref) < RARE_LEVEL]] <- "Other"
    lev <- names(sort(table(f[ref]), decreasing = TRUE)); lev <- c(lev, setdiff(unique(f), lev))
    fac <- factor(f, levels = lev)
    if (nlevels(fac) < 2L) next
    mm <- stats::model.matrix(~ fac - 1)[, -1, drop = FALSE]
    for (j in seq_len(ncol(mm)))
      put(paste0(v, "_", make.names(sub("^fac", "", colnames(mm)[j]))), as.numeric(mm[, j]), v)
  }
  for (fam in families) {
    rare <- rep(0L, nrow(reg)); n_rare <- 0L
    for (cc in grep(paste0("^", fam, "___[0-9]+$"), names(reg), value = TRUE)) {
      x <- as.integer(as_num(reg[[cc]]) == 1); x[is.na(x)] <- 0L
      if (mean(x[ref]) < RARE_BOX) { rare <- pmax(rare, x); n_rare <- n_rare + 1L }
      else put(sub("___", "_", cc), x, fam)
    }
    if (n_rare >= 2L) put(paste0(fam, "_other_rare"), rare, fam)
  }
  W <- as.data.frame(enc, check.names = FALSE)
  sdr <- vapply(W, function(x) stats::sd(x[ref]), 0)
  keep <- is.finite(sdr) & sdr > 1e-10
  W <- W[, keep, drop = FALSE]
  W <- W[, !duplicated(t(as.matrix(W[ref, , drop = FALSE]))), drop = FALSE]
  is_bin <- vapply(W, function(x) all(x %in% c(0, 1)), TRUE)
  pr <- colMeans(W[ref, , drop = FALSE])
  W <- W[, !(is_bin & (pr < NEAR_CONST | pr > 1 - NEAR_CONST)), drop = FALSE]
  list(W = W, variable = src[names(W)])
}

pre <- encode(
  continuous = c("age", "eiss_ais", "n_regions_injured", "head_injury_severity",
                 "months_since_start", "hour_sin", "hour_cos"),
  nominal = c("sex", "marital_status", "education", "occupation", "home_ownership",
              "television", "cable_tv", "mixer", "cooking_fuel", "inj_place",
              "inj_activity", "alcohol_use", "inj_alcohol_use", "inj_sub_county",
              "tobacco_use", "intent", "penetrating_wound"),
  families = c("mech_inj", "medical_history", "chronic_illness", "prior_surgery"),
  binary = c("strike_period", "arrival_weekend"))
arr <- encode(continuous = c("sbp", "dbp", "vitals_1_hr", "eyes_gcs", "verbal_gcs", "motor_gcs"),
              nominal = "pulse", families = character())
arr$W <- arr$W[, setdiff(names(arr$W), names(pre$W)), drop = FALSE]
cr_log(paste("Covariates fixed before transport:", ncol(pre$W), "columns from",
             length(unique(pre$variable)), "variables; arrival state:", ncol(arr$W), "columns"))

# ------------------------------------------------------------
# F. Assemble, join the steward's indicators, write
# ------------------------------------------------------------
ind <- readRDS(file.path(cfg$data$design_dir, "observation_indicators.rds"))
core <- data.frame(patient_id = reg$patient_id, arm = reg$arm,
                   operator_unknown = reg$operator_unknown, walked_in = reg$walked_in,
                   crossref_tier = reg$crossref_tier, transfer = reg$transfer,
                   prior_care_elsewhere = reg$prior_care_elsewhere,
                   direct_arrival = reg$direct_arrival, enrolment_date = reg$enrolment_date,
                   stringsAsFactors = FALSE)
core <- cbind(core, reg[, NC], pre$W, arr$W)
m <- match(core$patient_id, ind$patient_id)
stopifnot(!anyNA(m))
core$disposition_recorded <- ind$disposition_recorded[m]
core$six_month_recorded <- ind$six_month_recorded[m]
core <- core[in_run, , drop = FALSE]
rownames(core) <- NULL

# Dictionary rows for the covariates, stand-ins, controls and indicators.
# Timing: every registry field is recorded after time zero (the registry is
# completed at hospital); `fixed` separates the variables (spec section 3).
HOSP_CODED <- c("eiss_ais", "n_regions_injured", "penetrating_wound")
dict_row <- function(col, variable, source, fixed, role, reason = "", sensitivity = "",
                     follow_up = FALSE, recorded = "after")
  data.frame(column = col, variable = variable, source = source, follow_up = follow_up,
             recorded = recorded, fixed = fixed, role = role, reason = reason,
             sensitivity = sensitivity, stringsAsFactors = FALSE)
pre_reason <- function(v)
  if (v %in% HOSP_CODED) "Coded at hospital, but describes the injury itself." else
  if (v == "head_injury_severity") paste("Graded from the injury-form head-injury checkboxes; whether they",
                                         "were graded from the arrival Glasgow Coma Scale is an open question to co-authors.") else ""
pre_source <- function(v)
  if (v %in% c("inj_sub_county", "inj_place", "inj_activity", "intent", "inj_alcohol_use",
               "mech_inj", "head_injury_severity", "months_since_start", "hour_sin",
               "hour_cos", "arrival_weekend", "strike_period")) "registry injury form" else
  if (v %in% HOSP_CODED) "registry hospital assessment" else "registry enrolment form"
dict <- do.call(rbind, c(
  lapply(names(pre$W), function(cn) { v <- pre$variable[[cn]]
    dict_row(cn, v, pre_source(v), "before", "covariate", pre_reason(v)) }),
  lapply(names(arr$W), function(cn) dict_row(cn, arr$variable[[cn]],
    "registry arrival vital signs", "after", "stand_in",
    "Closest record of how ill the patient was at the scene, which decides whether an ambulance is called.",
    "The primary lock for the same comparison, which leaves arrival state out.")),
  lapply(NC, function(cn) dict_row(cn, cn, "registry enrolment form", "before", "negative_control")),
  list(dict_row("disposition_recorded", "disposition_recorded", "registry disposition (steward)",
                "after", "observation", "Approved design-stage input (D3, 2026-10-18); see the blinding ledger."),
       dict_row("six_month_recorded", "six_month_recorded", "six-month follow-up (steward)",
                "after", "observation", "Approved design-stage input (D3, 2026-10-18); see the blinding ledger.",
                follow_up = TRUE))))

OUTD <- cfg$data$design_dir
saveRDS(list(data = core, dictionary = dict, covariates = names(pre$W),
             arrival_state = names(arr$W), negative_controls = NC, run = RUN),
        file.path(OUTD, paste0("design_data_", RUN, ".rds")))
RES <- file.path(cfg$paths$results, "reexec_050")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(dict, file.path(RES, paste0("dictionary_", RUN, ".csv")), row.names = FALSE)
utils::write.csv(data.frame(column = c(names(pre$W), names(arr$W)),
                            block = c(rep("before transport", ncol(pre$W)), rep("arrival state", ncol(arr$W))),
                            mean_in_anchor = round(colMeans(cbind(pre$W, arr$W)[ref, , drop = FALSE]), 4)),
                 file.path(RES, paste0("design_columns_", RUN, ".csv")), row.names = FALSE)
cr_log(paste("Wrote design_data_", RUN, ".rds: ", nrow(core), " patients", sep = ""))
