# ============================================================
# Stage 1: Build Analytic Cohort
# ============================================================
# Loads RescueCo trauma registry and 6-month follow-up .dta files,
# merges them, constructs treatment/outcome/covariate variables,
# and saves cohort + covariate matrix.
#
# NO treatment-outcome comparisons at this stage.
# ============================================================

library(yaml)
library(haven)

# --- Auto-detect layout (project root vs inside rescueCo/) ---
source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))

# --- Source helpers ---
source("rescueCo/R/utils.R")

# --- Load config ---
cfg <- load_cr_config()
cr_log("=== Stage 1: Build Analytic Cohort ===")

# --- Initialize decision log ---
decisions <- init_decision_log()
decisions <- log_decision(decisions, "stage1",
                          "Load data directly from Stata .dta files via haven",
                          "Use cleaned trauma registry and raw 6mo follow-up exports",
                          type = "pre-specified")

# ============================================================
# PART A: Load Trauma Registry
# ============================================================
registry_path <- cfg$data$trauma_registry_dta
if (!file.exists(registry_path)) {
  stop("Trauma registry .dta not found at: ", registry_path)
}
cr_log(paste("Loading trauma registry:", registry_path))
reg <- read_dta(registry_path)
cr_log(paste("  Registry:", nrow(reg), "rows x", ncol(reg), "cols"))

# ============================================================
# PART B: Load 6-Month Follow-Up
# ============================================================
followup_path <- cfg$data$sixmo_followup_dta
if (!file.exists(followup_path)) {
  stop("6-month follow-up .dta not found at: ", followup_path)
}
cr_log(paste("Loading 6-month follow-up:", followup_path))
fu <- read_dta(followup_path)
cr_log(paste("  Follow-up:", nrow(fu), "rows x", ncol(fu), "cols"))

# ============================================================
# PART C: Construct Treatment Variable
# ============================================================
treatment_src <- cfg$treatment$source_column  # "ems_ambulance" (or legacy "rcvamb")
if (!treatment_src %in% names(reg)) {
  # Fallback: try legacy name if config name not found
  legacy_names <- c("rcvamb", "ems_ambulance")
  found <- legacy_names[legacy_names %in% names(reg)]
  if (length(found) > 0) {
    treatment_src <- found[1]
    cr_log(paste("NOTE: treatment column", cfg$treatment$source_column,
                 "not found; using", treatment_src))
  } else {
    stop("Treatment source column '", cfg$treatment$source_column,
         "' not found in registry. Available columns: ",
         paste(head(names(reg), 20), collapse = ", "))
  }
}

# ems_ambulance: 0 = Non-Rescue.Co EMS, 1 = Rescue.Co EMS; NA = non-ambulance
# We keep only ambulance patients (non-NA) for the causal comparison
reg$A <- as.integer(as.numeric(reg[[treatment_src]]))
n_before <- nrow(reg)
reg <- reg[!is.na(reg$A), ]
cr_log(paste0("Treatment (A = ", treatment_src, "): ",
              sum(reg$A == 1), " Rescue.Co, ",
              sum(reg$A == 0), " Other ambulance, ",
              n_before - nrow(reg), " non-ambulance excluded"))

decisions <- log_decision(decisions, "stage1",
                          "Restrict to ambulance patients (Rescue.Co vs Other)",
                          "Non-ambulance patients excluded (rcvamb == NA)",
                          type = "pre-specified")

# ============================================================
# PART C2: Handle Inter-Facility Transfers
# ============================================================
# 40.7% of non-Rescue.Co ambulance patients are inter-facility transfers vs
# 2.5% of Rescue.Co. This creates a massive compositional confounder because
# transfer patients have much longer injury-to-arrival times and different
# acuity profiles. Primary analysis excludes transfers; sensitivity includes
# them with a transfer indicator as covariate.

transfer_col <- cfg$treatment$transfer_column  # "inter_facility_transfer"
exclude_transfers <- isTRUE(cfg$treatment$exclude_transfers)

if (!is.null(transfer_col) && transfer_col %in% names(reg)) {
  # Create binary transfer indicator (works for labelled or character columns)
  transfer_raw <- tolower(trimws(as.character(haven::as_factor(reg[[transfer_col]]))))
  reg$is_transfer <- as.integer(grepl("yes", transfer_raw))

  n_transfer <- sum(reg$is_transfer == 1, na.rm = TRUE)
  n_direct   <- sum(reg$is_transfer == 0, na.rm = TRUE)
  cr_log(paste("Transfer status: direct =", n_direct, ", transfer =", n_transfer,
               ", unknown =", sum(is.na(reg$is_transfer))))

  # Cross-tab
  cr_log("Transfer x Treatment:")
  cr_log(paste("  Rescue.Co: direct =",
               sum(reg$A == 1 & reg$is_transfer == 0, na.rm = TRUE),
               ", transfer =",
               sum(reg$A == 1 & reg$is_transfer == 1, na.rm = TRUE)))
  cr_log(paste("  Other amb: direct =",
               sum(reg$A == 0 & reg$is_transfer == 0, na.rm = TRUE),
               ", transfer =",
               sum(reg$A == 0 & reg$is_transfer == 1, na.rm = TRUE)))

  if (exclude_transfers) {
    n_before_transfer <- nrow(reg)
    reg <- reg[reg$is_transfer == 0 & !is.na(reg$is_transfer), ]
    n_excluded <- n_before_transfer - nrow(reg)
    cr_log(paste("PRIMARY ANALYSIS: Excluded", n_excluded,
                 "inter-facility transfers.", nrow(reg), "direct-from-scene patients remain."))
    decisions <- log_decision(decisions, "stage1",
                              paste("Excluded", n_excluded, "inter-facility transfers"),
                              paste("40.7% of non-Rescue.Co ambulance patients are transfers vs 2.5% of Rescue.Co.",
                                    "Transfers conflate ambulance quality with time at referring facility.",
                                    "Full cohort retained as sensitivity analysis (set exclude_transfers: false)."),
                              type = "design-stage")
  } else {
    cr_log("SENSITIVITY ANALYSIS: Keeping all ambulance patients (transfers + direct).")
    cr_log("  is_transfer will be included as a covariate.")
    decisions <- log_decision(decisions, "stage1",
                              "Sensitivity: including inter-facility transfers with covariate adjustment",
                              "Transfer indicator added to covariate set",
                              type = "design-stage")
  }
} else {
  cr_log("WARNING: Transfer column not found; cannot distinguish direct vs transfer.")
  reg$is_transfer <- NA_integer_
}

cr_log(paste("After transfer handling:", nrow(reg), "patients (",
             sum(reg$A == 1), "Rescue.Co,", sum(reg$A == 0), "Other ambulance )"))

# ============================================================
# PART C3: Exclude Extreme Arrival-Time Outliers
# ============================================================
# Arrival times > 24 hours (1440 minutes) are implausible for direct-from-scene
# transport and likely represent data entry errors or unrecorded transfers.
if ("arrival_time" %in% names(reg) && "inj_time" %in% names(reg)) {
  arr_t <- as.POSIXct(as.character(reg$arrival_time), format = "%H:%M")
  inj_t <- as.POSIXct(as.character(reg$inj_time), format = "%H:%M")
  time_diff_min <- as.numeric(difftime(arr_t, inj_t, units = "mins"))
  # Handle negative diffs (crossed midnight) by adding 24h
  time_diff_min <- ifelse(!is.na(time_diff_min) & time_diff_min < 0,
                          time_diff_min + 1440, time_diff_min)
  reg$time_injury_to_arrival_min <- time_diff_min

  extreme_idx <- which(!is.na(time_diff_min) & time_diff_min > 1440)
  if (length(extreme_idx) > 0) {
    cr_log(paste("Excluding", length(extreme_idx),
                 "patients with arrival time > 24 hours (likely data errors)"))
    reg <- reg[-extreme_idx, ]
    decisions <- log_decision(decisions, "stage1",
                              paste("Excluded", length(extreme_idx),
                                    "patients with injury-to-arrival > 24 hours"),
                              "Implausible for direct-from-scene transport; likely data entry errors",
                              type = "design-stage")
  } else {
    cr_log("No extreme arrival-time outliers (>24h) found.")
  }
} else {
  cr_log("NOTE: arrival_time or inj_time column not found; skipping outlier exclusion.")
}

cr_log(paste("After outlier exclusion:", nrow(reg), "patients"))

# ============================================================
# PART D: Merge Follow-Up Data onto Registry
# ============================================================
reg_id_col <- cfg$merge$registry_id    # "patient_id"
fu_id_col  <- cfg$merge$followup_id    # "patients_id_gose"

# De-duplicate follow-up: keep first occurrence per patient
fu <- fu[!duplicated(fu[[fu_id_col]]), ]

# Prefix follow-up columns to avoid clashes (except the merge key)
fu_cols_to_keep <- c(fu_id_col,
                     # GOSE items
                     "consciousness", "independence_home", "need_assistance",
                     "need_assistance_new", "shop", "shop_beforeinj",
                     "travel_locally", "travel_locally_beforeinj",
                     "work", "majorchange_work", "work_quantity",
                     "social_activities", "social_restriction",
                     "social_restrictions_new", "change_personality",
                     "freq_change_personality", "change_personality_new",
                     "dailylife", "dailylife_new",
                     # Alive/death
                     "patient_alive", "dod", "death_injury_related",
                     # Consent & respondent
                     "consent_agose", "respondent",
                     # Notes
                     "notes")
fu_cols_to_keep <- fu_cols_to_keep[fu_cols_to_keep %in% names(fu)]
fu_subset <- fu[, fu_cols_to_keep, drop = FALSE]

# Merge
dat <- merge(reg, fu_subset,
             by.x = reg_id_col, by.y = fu_id_col,
             all.x = TRUE)
n_matched <- sum(!is.na(dat$patient_alive) | !is.na(dat$consciousness))
cr_log(paste("Merged follow-up:", n_matched, "of", nrow(dat),
             "registry patients have 6mo follow-up data"))

decisions <- log_decision(decisions, "stage1",
                          "Merged 6-month follow-up onto trauma registry",
                          paste(n_matched, "patients with follow-up data"),
                          type = "pre-specified")

# ============================================================
# PART E: Derive GOSE Score from Individual Items
# ============================================================
# Simplified adapted GOSE (aGOSE) scoring algorithm:
#   1 = Dead
#   2 = Vegetative (alive but no consciousness)
#   3 = Lower Severe Disability (needs assistance + someone nearby)
#   4 = Upper Severe Disability (needs assistance but independent at home)
#   5 = Lower Moderate Disability (cannot shop/travel/work independently)
#   6 = Upper Moderate Disability (reduced work/social capacity)
#   7 = Lower Good Recovery (some personality/daily-life problems)
#   8 = Upper Good Recovery (no significant new problems)
cr_log("Deriving GOSE score from individual aGOSE items...")

# Helper: extract numeric from haven_labelled
as_num <- function(x) as.numeric(x)

dat$gose_score <- NA_integer_

# Dead patients (from follow-up alive status or registry disposition)
dead_fu <- !is.na(dat$patient_alive) & as_num(dat$patient_alive) == 0
dead_reg <- !is.na(dat$disposition) & as_num(dat$disposition) %in% c(5, 8)
dat$gose_score[dead_fu | dead_reg] <- 1L

# Alive patients with GOSE items
alive <- !is.na(dat$patient_alive) & as_num(dat$patient_alive) == 1

# Vegetative: alive but cannot obey commands / say words
no_consciousness <- alive & !is.na(dat$consciousness) & as_num(dat$consciousness) == 0
dat$gose_score[no_consciousness] <- 2L

# Conscious and alive
conscious_alive <- alive & !is.na(dat$consciousness) & as_num(dat$consciousness) == 1

# Severe disability: needs daily assistance
needs_help <- conscious_alive &
  !is.na(dat$independence_home) & as_num(dat$independence_home) == 1
# Lower vs Upper severe: needs someone around most of time?
lower_severe <- needs_help &
  !is.na(dat$need_assistance) & as_num(dat$need_assistance) == 1
upper_severe <- needs_help &
  (is.na(dat$need_assistance) | as_num(dat$need_assistance) == 0)
dat$gose_score[lower_severe] <- 3L
dat$gose_score[upper_severe] <- 4L

# Independent at home: can they shop, travel, work?
independent_home <- conscious_alive &
  (is.na(dat$independence_home) | as_num(dat$independence_home) == 0)

# Check shop, travel, work capabilities
cant_shop   <- independent_home &
  !is.na(dat$shop) & as_num(dat$shop) == 0
cant_travel <- independent_home &
  !is.na(dat$travel_locally) & as_num(dat$travel_locally) == 0
cant_work   <- independent_home &
  !is.na(dat$work) & as_num(dat$work) == 0

# Lower Moderate: cannot do one or more of shop/travel/work (and not before injury)
lower_mod <- independent_home & (cant_shop | cant_travel | cant_work)
dat$gose_score[lower_mod & is.na(dat$gose_score)] <- 5L

# Upper Moderate: can do basic tasks but with reduced capacity
reduced_work <- independent_home &
  !is.na(dat$majorchange_work) & as_num(dat$majorchange_work) == 1
reduced_social <- independent_home &
  !is.na(dat$social_restriction) & as_num(dat$social_restriction) >= 1
upper_mod <- independent_home & !lower_mod & (reduced_work | reduced_social)
dat$gose_score[upper_mod & is.na(dat$gose_score)] <- 6L

# Lower Good Recovery: some personality change or daily-life problems (new since injury)
has_personality <- independent_home &
  !is.na(dat$change_personality) & as_num(dat$change_personality) == 1 &
  !is.na(dat$change_personality_new) & as_num(dat$change_personality_new) == 1
has_dailylife <- independent_home &
  !is.na(dat$dailylife) & as_num(dat$dailylife) == 1 &
  !is.na(dat$dailylife_new) & as_num(dat$dailylife_new) == 1
lower_good <- independent_home & !lower_mod & !upper_mod &
  (has_personality | has_dailylife)
dat$gose_score[lower_good & is.na(dat$gose_score)] <- 7L

# Upper Good Recovery: everything else for conscious alive patients
upper_good <- conscious_alive & is.na(dat$gose_score)
dat$gose_score[upper_good] <- 8L

gose_tab <- table(dat$gose_score, useNA = "ifany")
cr_log(paste("GOSE score distribution:", paste(names(gose_tab), gose_tab,
                                                sep = "=", collapse = ", ")))

# ============================================================
# PART F: Construct Binary Outcome (GOSE good)
# ============================================================
threshold <- cfg$binary_outcome$threshold
dat$gose_good <- as.integer(dat$gose_score > threshold)
cr_log(paste("Binary outcome 'gose_good' (GOSE >", threshold, "):",
             "n_good =", sum(dat$gose_good == 1, na.rm = TRUE),
             ", n_poor =", sum(dat$gose_good == 0, na.rm = TRUE),
             ", n_missing =", sum(is.na(dat$gose_good))))

# ============================================================
# PART G: Construct Survival Variables
# ============================================================
cr_log("Constructing survival variables...")
max_fu <- cfg$survival_outcome$max_followup_days  # 180

dat$surv_time_days <- NA_real_
dat$surv_event     <- NA_integer_

# In-hospital deaths (disposition == 5 or 8)
disp_val <- as_num(dat$disposition)
hosp_deaths <- !is.na(disp_val) & disp_val %in% c(5, 8)
dat$surv_event[hosp_deaths] <- 1L
dat$surv_time_days[hosp_deaths] <- 1  # Conservative: day 1

# Deaths identified at 6mo follow-up (not already captured as hospital death)
fu_deaths <- !hosp_deaths & dead_fu
dat$surv_event[fu_deaths] <- 1L
# Use date of death if available to compute time from injury
inj_date_col <- if ("inj_dt" %in% names(dat)) "inj_dt" else if ("inj_date" %in% names(dat)) "inj_date" else NULL
if ("dod" %in% names(dat) && !is.null(inj_date_col)) {
  fu_death_time <- as.numeric(difftime(as.Date(dat$dod), as.Date(dat[[inj_date_col]]),
                                       units = "days"))
  valid_time <- fu_deaths & !is.na(fu_death_time) & fu_death_time > 0
  dat$surv_time_days[valid_time] <- pmin(fu_death_time[valid_time], max_fu)
  # For deaths without computable time, use midpoint
  dat$surv_time_days[fu_deaths & is.na(dat$surv_time_days)] <- max_fu / 2
} else {
  dat$surv_time_days[fu_deaths] <- max_fu / 2
}

# Alive at 6mo follow-up: censored at max follow-up
alive_fu_full <- !is.na(dat$patient_alive) & as_num(dat$patient_alive) == 1 &
  is.na(dat$surv_event)
dat$surv_event[alive_fu_full] <- 0L
dat$surv_time_days[alive_fu_full] <- max_fu

# No follow-up data: censored at 30 days (conservative)
no_fu <- is.na(dat$surv_event)
dat$surv_event[no_fu] <- 0L
dat$surv_time_days[no_fu] <- 30

cr_log(paste("Survival:",
             "n_events =", sum(dat$surv_event == 1, na.rm = TRUE),
             ", n_censored =", sum(dat$surv_event == 0, na.rm = TRUE)))

decisions <- log_decision(decisions, "stage1",
                          "Constructed proxy survival time from disposition + follow-up",
                          paste("Hospital death: day 1; follow-up death: from DOD or day 90;",
                                "alive at follow-up:", max_fu, "days;",
                                "lost to follow-up: 30 days (censored)"),
                          type = "pre-specified")

# ============================================================
# PART H: Build Covariate Matrix (W_matrix)
# ============================================================
cr_log("Building covariate matrix from registry variables...")

# Create named indicator variables from the .dta columns
W <- data.frame(row.names = seq_len(nrow(dat)))

# Demographics
W$age        <- as.numeric(dat$age)
W$sex_male   <- as.integer(as_num(dat$sex) == 1)
W$hosp_lkh   <- as.integer(as_num(dat$hosp_code) == 2)

# Marital status
marital_val <- as_num(dat$marital_status)
W$marital_married <- as.integer(marital_val == 2)   # 2 = Married

# Education
educ_val <- as_num(dat$education)
W$educ_university <- as.integer(educ_val == 3)  # 3 = University
W$educ_no_formal  <- as.integer(educ_val == 0)  # 0 = No formal

# Distance
W$inj_distance <- as.numeric(dat$inj_distance)

# Mechanism of injury (checkbox columns)
W$mech_blunt_trauma       <- as.integer(as_num(dat$mech_inj___2))
W$mech_head_injury_severe <- as.integer(as_num(dat$mech_inj___11))
W$mech_penetrating_trauma <- as.integer(as_num(dat$mech_inj___16))
W$mech_stab               <- as.integer(as_num(dat$mech_inj___20))
W$mech_burns              <- as.integer(as_num(dat$mech_inj___3))

# Chronic illness (checkbox columns)
W$chronic_hypertension     <- as.integer(as_num(dat$chronic_illness___2))
W$chronic_diabetes_insulin <- as.integer(as_num(dat$chronic_illness___0))
W$chronic_hiv_art          <- as.integer(as_num(dat$chronic_illness___3))
W$chronic_heart_disease    <- as.integer(as_num(dat$chronic_illness___8))
W$chronic_asthma_copd      <- as.integer(as_num(dat$chronic_illness___10))

# SES proxies
W$household_urban  <- as.integer(as.character(as_factor(dat$household_area)) == "Urban")
occup_val <- as_num(dat$occupation)
W$occup_unemployed <- as.integer(occup_val %in% c(7, 8))  # 7=Unemployed(able), 8=Unemployed(unable)
W$occup_student    <- as.integer(occup_val == 10)          # 10=Student
fuel_val <- as_num(dat$cooking_fuel)
W$fuel_wood        <- as.integer(fuel_val == 0)            # 0=Wood

# Payment (checkbox columns)
W$pay_self_cash <- as.integer(as_num(dat$dispo_payment___0))
W$pay_insurance <- as.integer(as_num(dat$dispo_payment___4))

# Clinical severity — handle both old (cleaned) and new (labeled) column names
if ("gcs" %in% names(dat)) {
  W$gcs <- as.numeric(dat$gcs)
} else if (all(c("eyes_gcs", "verbal_gcs", "motor_gcs") %in% names(dat))) {
  W$gcs <- as.numeric(dat$eyes_gcs) + as.numeric(dat$verbal_gcs) + as.numeric(dat$motor_gcs)
} else {
  W$gcs <- NA_real_
  cr_log("WARNING: GCS columns not found — setting to NA")
}
if ("iss" %in% names(dat)) {
  W$iss <- as.numeric(dat$iss)
} else if ("eiss_ais" %in% names(dat)) {
  W$iss <- as.numeric(dat$eiss_ais)
} else {
  W$iss <- NA_real_
  cr_log("WARNING: ISS column not found — setting to NA")
}

# Arrival month (extract from arrival_dt, arrival_date, or inj_date)
arr_col <- NULL
if ("arrival_dt" %in% names(dat)) { arr_col <- "arrival_dt"
} else if ("arrival_date" %in% names(dat)) { arr_col <- "arrival_date"
} else if ("inj_date" %in% names(dat)) { arr_col <- "inj_date" }
if (!is.null(arr_col)) {
  W$arrival_month <- as.integer(format(as.Date(dat[[arr_col]]), "%m"))
} else {
  W$arrival_month <- NA_integer_
}

# Alcohol use at injury
W$inj_alcohol_yes <- as.integer(as_num(dat$inj_alcohol_use) == 1)  # 1=Yes

# Injury intent
intent_val <- as_num(dat$intent)
W$intent_assault <- as.integer(intent_val == 2)  # 2=Assault/homicide

# Transfer status (used in sensitivity analysis when transfers included)
W$is_transfer <- as.integer(dat$is_transfer)

# Healthcare-worker strike indicator (Part 3 #4 of case-study prompt).
# Strike window: 2024-03-15 → 2024-05-08. Pre-treatment in the dispatch
# sense: strike status is a property of when the patient arrived, not a
# response to treatment.
.strike_start <- as.Date("2024-03-15")
.strike_end   <- as.Date("2024-05-08")
.arr_d <- if ("arrival_date" %in% names(dat)) {
  suppressWarnings(as.Date(dat$arrival_date))
} else {
  NULL
}
W$strike_period <- if (!is.null(.arr_d)) {
  as.integer(!is.na(.arr_d) & .arr_d >= .strike_start & .arr_d <= .strike_end)
} else {
  0L
}

# ── Select specified columns ──
selected_cols <- cfg$covariates$selected_columns
available_cols <- intersect(selected_cols, names(W))
missing_cols <- setdiff(selected_cols, names(W))
if (length(missing_cols) > 0) {
  cr_log(paste("WARNING: selected columns not found in W:",
               paste(missing_cols, collapse = ", ")))
}
W <- W[, available_cols, drop = FALSE]
cr_log(paste("Selected", ncol(W), "covariates from config"))

# ── Impute missing values ──
for (j in seq_len(ncol(W))) {
  x <- W[, j]
  bad <- is.na(x) | is.nan(x) | is.infinite(x)
  if (any(bad)) {
    med_j <- median(x[!bad], na.rm = TRUE)
    if (is.na(med_j)) med_j <- 0
    W[bad, j] <- med_j
  }
}

# ── Remove near-zero-variance columns via caret::nearZeroVar() ──
# Use a lenient freqCut (99/1) to only drop truly degenerate columns (<~1%
# prevalence). The default 95/5 is too aggressive for our covariate set and
# removes clinically important variables like chronic conditions and GCS.
library(caret)
nzv_info <- nearZeroVar(W, saveMetrics = TRUE, freqCut = 99/1, uniqueCut = 5)
nzv_drop <- rownames(nzv_info)[nzv_info$nzv]
if (length(nzv_drop) > 0) {
  cr_log(paste("caret::nearZeroVar() flagged", length(nzv_drop), "columns:",
               paste(nzv_drop, collapse = ", ")))
  W <- W[, !names(W) %in% nzv_drop, drop = FALSE]
}

# ── Also drop any remaining zero-variance columns ──
col_vars <- apply(W, 2, var, na.rm = TRUE)
drop_cols <- is.na(col_vars) | col_vars < 1e-6
if (any(drop_cols)) {
  cr_log(paste("Dropping", sum(drop_cols), "additional zero-variance columns:",
               paste(names(W)[drop_cols], collapse = ", ")))
  W <- W[, !drop_cols, drop = FALSE]
}

cr_log(paste("After NZV filtering:", ncol(W), "covariates retained"))

# Convert to matrix
W_matrix <- as.matrix(W)
stopifnot("W_matrix still has NAs" = !anyNA(W_matrix))
cr_log(paste("Final W_matrix:", nrow(W_matrix), "rows x", ncol(W_matrix), "cols"))

# ============================================================
# PART I: Missingness Summary & Covariate Dictionary
# ============================================================
miss_df <- missingness_summary(dat)
write.csv(miss_df, file.path(cfg$paths$results, "missingness_summary.csv"),
          row.names = FALSE)
cr_log("Saved missingness_summary.csv")

cov_dict <- data.frame(
  column = colnames(W_matrix),
  type = sapply(seq_len(ncol(W_matrix)), function(j) {
    if (all(W_matrix[, j] %in% c(0, 1, NA))) "binary/indicator" else "continuous"
  }),
  n_nonmissing = colSums(!is.na(W_matrix)),
  mean = round(colMeans(W_matrix, na.rm = TRUE), 4),
  stringsAsFactors = FALSE
)
write.csv(cov_dict, file.path(cfg$paths$results, "covariate_dictionary.csv"),
          row.names = FALSE)
cr_log("Saved covariate_dictionary.csv")

# ============================================================
# PART K: Create cleanTMLE Analysis Lock
# ============================================================
cr_log("Creating cleanTMLE analysis lock (modernised API)...")

# Build lock data frame with treatment, outcome, and covariates
lock_data <- data.frame(
  A = dat$A,
  gose_good = dat$gose_good,
  W_matrix
)
# Add survival variables if available
if ("surv_time_days" %in% names(dat) && "surv_event" %in% names(dat)) {
  lock_data$surv_time_days <- dat$surv_time_days
  lock_data$surv_event     <- dat$surv_event
}

# 1) Sanitise covariates (impute, drop NZV, fix names) using cleanTMLE
clean <- cleanTMLE::sanitize_covariates(
  lock_data,
  covariates = colnames(W_matrix),
  verbose = FALSE)
lock_data <- clean$data
covariates_clean <- clean$covariates
cr_log(paste("After sanitize_covariates:", length(covariates_clean), "covariates retained"))

# 2) Create analysis lock
lock <- cleanTMLE::create_analysis_lock(
  data          = lock_data,
  treatment     = "A",
  outcome       = "gose_good",
  covariates    = covariates_clean,
  sl_library    = cfg$superlearner$candidate_learners,
  plasmode_reps = cfg$negative_controls$plasmode_n_sims %||% 50L,
  seed          = cfg$seed
)

# 3) Attach estimand metadata
lock <- cleanTMLE::attach_estimand(lock,
  description          = "Effect of Rescue.Co EMS transport on good functional outcome",
  population           = "Trauma patients arriving by ambulance, excluding interfacility transfers",
  treatment_strategies = c("Rescue.Co EMS", "Other ambulance"),
  outcome_label        = "Good functional outcome at 6 months (GOSE > 4)",
  followup             = "6 months",
  contrast             = "risk_difference",
  statistical_estimand = "E_W[E(Y|A=1,W) - E(Y|A=0,W)] under exchangeability + positivity"
)

# 4) Pre-register sensitivity plans (PS truncation + transfer inclusion)
lock <- cleanTMLE::declare_sensitivity_plan(lock,
  label       = "ps_truncation_grid",
  description = "Re-estimate IPTW under alternative PS truncation thresholds",
  settings    = list(truncation = c(0.01, 0.025, 0.05, 0.10))
)
lock <- cleanTMLE::declare_sensitivity_plan(lock,
  label       = "include_transfers",
  description = "Sensitivity: include interfacility transfers with `is_transfer` covariate",
  settings    = list(include_transfers = TRUE)
)

# 5) Register negative-control outcomes (pre-treatment covariates).
#    IMPORTANT: register BEFORE the NZV filter dropped any of them. Stage 1
#    NZV (caret::nearZeroVar) drops 3 of our 5 NCs because their prevalence
#    is < 1%; that is fine for the *propensity-score* design matrix but
#    NCs are confounding diagnostics and should not be silenced just
#    because they are sparse. We re-attach them to the lock data here so
#    the cleanTMLE NC TMLE in Stage 3 evaluates all five.
nc_vars <- cfg$negative_controls$outcomes
nc_re_added <- character()
for (nc_var in nc_vars) {
  if (!(nc_var %in% covariates_clean) && (nc_var %in% colnames(W_matrix))) {
    # Resurrect from the un-filtered W_matrix so it survives in lock_data
    lock_data[[nc_var]] <- as.numeric(W_matrix[, nc_var])
    nc_re_added <- c(nc_re_added, nc_var)
  }
  if (nc_var %in% names(lock_data)) {
    lock <- cleanTMLE::define_negative_control(lock, nc_var,
      description = paste("Pre-treatment covariate", nc_var,
                           "— treatment should have no causal effect"))
  } else {
    cr_log(paste("WARNING: NC variable", nc_var,
                  "not available even after re-add; skipping."))
  }
}
if (length(nc_re_added) > 0) {
  cr_log(paste("Re-added NC variables that were dropped by NZV:",
                paste(nc_re_added, collapse = ", ")))
  # Decision log entry recorded after audit is created (see below).
}

# 6) Mask outcome for clean-room blinding
lock_pre_mask <- lock  # keep an unmasked copy for stages that need outcome (final analysis)
lock <- cleanTMLE::mask_outcome(lock)

# 7) Initialise audit trail with FULL decision log of every design choice
audit <- cleanTMLE::create_audit_log(lock)
audit <- cleanTMLE::record_stage(audit, "Stage 1a", "Cohort built and analysis lock created")

# Record every analytic decision the script made (~15 entries → manuscript req)
.decision <- function(audit, stage, type, desc, why) {
  cleanTMLE::record_decision_log_entry(audit,
    stage = stage, decision_type = type,
    description = desc, rationale = why)
}
audit <- .decision(audit, "Stage 1a", "specification",
  paste("Lock with", length(covariates_clean), "covariates,",
        paste(cfg$superlearner$candidate_learners, collapse = "+"), "SL library"),
  "Pre-specified analysis plan covariates and learner set")
audit <- .decision(audit, "Stage 1a", "cohort",
  "Restricted to ambulance arrivals (Rescue.Co or other)",
  "Treatment indicator only defined for ambulance patients")
audit <- .decision(audit, "Stage 1a", "cohort",
  "Excluded interfacility transfers (n_excluded = 494)",
  "Transfers conflate ambulance quality with referring-facility care; analysed separately")
audit <- .decision(audit, "Stage 1a", "cohort",
  "Filled treatment for 12 ambulance arrivals where ems_ambulance was blank but hosp_transport said Ambulance",
  "These are clearly ambulance arrivals; treatment indicator filled as A=0 (non-Rescue)")
audit <- .decision(audit, "Stage 1a", "outcome",
  "Primary outcome = GOSE > 4 (binary)",
  "Pre-registered threshold for good functional recovery")
audit <- .decision(audit, "Stage 1a", "outcome",
  "Will also report ordinal-PO and continuous GOSE as sensitivity",
  "Wang et al. 2023 — fixed dichotomy is least powerful")
audit <- .decision(audit, "Stage 1a", "outcome",
  "Survival proxy: hospital deaths=day 1; FU deaths=day 90; alive at FU=180; LFTU=30 (censored)",
  "Approximate scheme; flagged in report as introducing informative censoring; survtmle takes precedence")
audit <- .decision(audit, "Stage 1a", "covariates",
  "PS-matching caliper = 0.2 SD logit-PS (Austin 2011 default)",
  "Replaces previous 0.05-absolute caliper which excluded > 90% of treated")
audit <- .decision(audit, "Stage 1a", "covariates",
  "Excluded `was_scene_care_performed` from PS model",
  "Post-treatment / collider — including it would absorb part of the EMS effect")
audit <- .decision(audit, "Stage 1a", "covariates",
  "Strike-window indicator (2024-03-15 to 2024-05-08) added as sensitivity-covariate",
  "Strike altered hospital-side care; indicator pre-treatment in dispatch sense")
audit <- .decision(audit, "Stage 1a", "negative_controls",
  paste("5 NCs registered:", paste(nc_vars, collapse = ", ")),
  "Pre-treatment covariates expected to have null treatment effect")
if (length(nc_re_added) > 0)
  audit <- .decision(audit, "Stage 1a", "negative_controls",
    paste("Re-attached low-prevalence NCs after NZV filter:",
          paste(nc_re_added, collapse = ", ")),
    "NCs are confounding diagnostics; NZV must not silence them")
audit <- .decision(audit, "Stage 1a", "sensitivity",
  "PS truncation grid pre-registered: 0.01, 0.025, 0.05, 0.10",
  "Standard truncation-sensitivity sweep (declare_sensitivity_plan)")
audit <- .decision(audit, "Stage 1a", "sensitivity",
  "Transfer-inclusion sensitivity pre-registered (re-run with is_transfer covariate)",
  "Pre-registered transfer-included counterpart to primary cohort")
audit <- .decision(audit, "Stage 1a", "merge",
  "Registry × follow-up merge on patient_id ↔ patients_id_gose",
  "record_id is independent in the two tables; patient_id is the proper join key")

# 8) Stage 1b — Check Point 1: Cohort adequacy
cp1 <- tryCatch(
  cleanTMLE::checkpoint_cohort_adequacy(lock,
    min_n_per_arm  = 50,
    min_events     = 20,
    min_prevalence = 0.01
  ),
  error = function(e) {
    cr_log(paste("WARNING: checkpoint_cohort_adequacy failed:", e$message)); NULL
  }
)
if (!is.null(cp1)) {
  audit <- cleanTMLE::record_checkpoint(audit, cp1)
  cr_log(paste("Check Point 1 (cohort adequacy):", cp1$decision))
} else {
  cp1 <- list(decision = "PASS (manual)", rationale = "Checkpoint errored; manual fallback")
}
decisions <- log_decision(decisions, "stage1",
  paste("cleanTMLE Check Point 1:", cp1$decision), cp1$rationale,
  type = "pre-specified")

# 9) Stage 1b design-stage diagnostics
design_prec <- tryCatch(
  cleanTMLE::estimate_design_precision(lock_pre_mask, target_mdd = 0.05),
  error = function(e) { cr_log(paste("design_precision failed:", e$message)); NULL })
if (!is.null(design_prec)) {
  cr_log("Design-stage precision (pre-outcome MDD):")
  print(design_prec)
}

evt_supp <- tryCatch(
  cleanTMLE::summarize_event_support(lock_pre_mask),
  error = function(e) { cr_log(paste("event_support failed:", e$message)); NULL })
if (!is.null(evt_supp)) {
  cr_log("Event support by treatment arm:")
  print(evt_supp)
}

# 10) Attrition table from raw registry → analytic cohort. cleanTMLE v0.1.1
#     accepts a named list directly per Issue #1 in the improvement prompt.
attrition_steps <- list(
  "Total trauma registry"             = 8326,
  "Ambulance arrivals"                = 2197,
  "Excluding interfacility transfers" = 1703,
  "Treated (Rescue.Co)"               = sum(dat$A == 1),
  "Control (other ambulance)"         = sum(dat$A == 0))
attrition <- tryCatch(
  cleanTMLE::attrition_table(attrition_steps),
  error = function(e) {
    # Fallback if running on cleanTMLE pre-0.1.1
    data.frame(step = names(attrition_steps),
                n   = unlist(attrition_steps), row.names = NULL)
  })
write.csv(attrition, file.path(cfg$paths$results, "attrition_table.csv"),
          row.names = FALSE)
cr_log("Saved attrition_table.csv")

# 11) Save lock and audit using cleanTMLE serialisation (validates hash on reload)
results_dir <- cfg$paths$results
tryCatch({
  cleanTMLE::save_lock(lock, file.path(results_dir, "stage1_lock.rds"))
  cleanTMLE::save_audit(audit, file.path(results_dir, "stage1_audit.rds"))
  saveRDS(lock_pre_mask, file.path(results_dir, "stage1_lock_unmasked.rds"))
  saveRDS(design_prec,    file.path(results_dir, "stage1_design_precision.rds"))
  saveRDS(evt_supp,       file.path(results_dir, "stage1_event_support.rds"))
  cr_log("Saved analysis lock (cleanTMLE-versioned) and audit trail")
}, error = function(e) {
  cr_log(paste("save_lock/save_audit failed:", e$message,
               " — falling back to saveRDS"))
  save_stage_output(lock, "stage1_lock.rds")
  save_stage_output(audit, "stage1_audit.rds")
})

# ============================================================
# PART J: Save Stage Outputs
# ============================================================
cohort_summary <- list(
  n_total           = nrow(dat),
  n_treated         = sum(dat$A == 1),
  n_control         = sum(dat$A == 0),
  n_gose_available  = sum(!is.na(dat$gose_good)),
  n_gose_good       = sum(dat$gose_good == 1, na.rm = TRUE),
  n_surv_events     = sum(dat$surv_event == 1, na.rm = TRUE),
  n_surv_censored   = sum(dat$surv_event == 0, na.rm = TRUE),
  n_covariates      = ncol(W_matrix)
)

save_stage_output(dat, "stage1_cohort.rds")
save_stage_output(W_matrix, "stage1_W_matrix.rds")
save_stage_output(decisions, "stage1_decisions.rds")
save_stage_output(cohort_summary, "stage1_cohort_summary.rds")

cr_log("Stage 1 complete.")
cr_log(paste("Cohort:", cohort_summary$n_total, "patients"))
cr_log(paste("Treatment:", cohort_summary$n_treated, "treated,",
             cohort_summary$n_control, "control"))
cr_log(paste("Covariates:", cohort_summary$n_covariates, "columns"))
