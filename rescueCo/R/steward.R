# Data-steward functions for the Rescue.Co re-execution (spec
# docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md).
#
# Only the steward scripts source this file: 20_steward_indicators.R, before
# the design stage, which writes two outcome-observation indicators and
# nothing else; and 24_steward_outcomes.R, after the dossier is approved,
# which writes the outcome store. The analyst's design build never sources it.
#
# The linkage and the scoring follow the unblinded analysis
# (rescueCo_analysis/analysis/01_ingest_merge_clean.R,
# 02_define_exposure_outcomes_covariates.R, shared_functions.R), reimplemented
# here. Nothing in this file prints or logs a value for an individual patient.

STEWARD_GOSE_SEVERE_CUT <- 4L      # GOSE 4 or below: severe disability or death
STEWARD_SHOP_TRAVEL_CODE <- 4L     # unable to shop or travel unaided: GOSE 4
STEWARD_HORIZON_DAYS <- 180        # six months, counted from arrival

.st_num <- function(x) {
  if (inherits(x, "haven_labelled")) return(as.numeric(unclass(x)))
  suppressWarnings(as.numeric(as.character(x)))
}

# Follow-up columns the steward reads: the id, the timestamp used to break ties,
# vital status, date of death, and the GOSE items. No name or phone column.
STEWARD_FU_COLUMNS <- c(
  "patients_id_gose", "patient_alive", "dod", "consciousness", "respondent",
  "independence_home", "need_assistance", "shop", "shop_beforeinj",
  "travel_locally", "travel_locally_beforeinj", "work", "majorchange_work",
  "social_activities", "social_restriction", "social_restrictions_new",
  "change_personality", "freq_change_personality", "change_personality_new",
  "dailylife", "dailylife_new")

# Read the follow-up file with a column whitelist (plus any survey timestamp).
steward_read_followup <- function(path) {
  nms <- names(haven::read_dta(path, n_max = 0))
  ts <- grep("^survey_timestamp", nms, value = TRUE)
  keep <- intersect(c(STEWARD_FU_COLUMNS, ts), nms)
  miss <- setdiff(STEWARD_FU_COLUMNS, nms)
  if (length(miss))
    stop("Follow-up file lacks: ", paste(miss, collapse = ", "), call. = FALSE)
  haven::read_dta(path, col_select = tidyselect::all_of(keep))
}

# Link follow-up records to registry ids as the unblinded analysis does:
# accepted typo corrections, then de-duplication (prefer a row with vital
# status, then the earliest timestamp), then the saved reconciliation actions.
# The two audit files hold study ids only.
steward_link_followup <- function(fu, reg_ids, typo_audit, recon_actions) {
  fid <- as.character(fu$patients_id_gose)
  ty <- utils::read.csv(typo_audit, stringsAsFactors = FALSE)[, c("bad_id", "resolved_to", "decision")]
  ty <- ty[grepl("^ACCEPT", ty$decision) & !is.na(ty$resolved_to), ]
  hit <- match(fid, ty$bad_id)
  fid[!is.na(hit)] <- ty$resolved_to[hit[!is.na(hit)]]
  fu$patients_id_gose <- fid

  ts <- grep("^survey_timestamp", names(fu), value = TRUE)[1]
  ord <- order(fu$patients_id_gose, is.na(.st_num(fu$patient_alive)),
               if (!is.na(ts)) as.POSIXct(as.character(fu[[ts]]), tz = "UTC") else seq_len(nrow(fu)))
  fu <- fu[ord, ]
  fu <- fu[!duplicated(fu$patients_id_gose), ]

  act <- utils::read.csv(recon_actions, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  lnk <- act[act$action == "LINKED (Isaac-confirmed)", ]
  for (r in seq_len(nrow(lnk))) {
    h <- which(fu$patients_id_gose == lnk$id[r])
    if (length(h) == 1L && lnk$resolved_to[r] %in% reg_ids)
      fu$patients_id_gose[h] <- lnk$resolved_to[r]
  }
  unlink_ids <- act$id[grepl("^UNLINKED", act$action)]
  fu <- fu[!fu$patients_id_gose %in% unlink_ids, ]
  fu <- fu[fu$patients_id_gose %in% reg_ids & !duplicated(fu$patients_id_gose), ]
  fu
}

# GOSE from the structured interview, as in the unblinded analysis's
# compute_gose_score(): "cannot say any words" is a vegetative state only when
# a proxy answered and no item shows independence; unable to shop or travel
# unaided (when able before the injury) is GOSE 4.
.st_gose_one <- function(alive, conscious, respondent, indep_home, need_assist,
                         shop_v, shop_pre, travel_v, travel_pre, work_v, work_major,
                         social_v, social_restr, social_new,
                         person_v, person_freq, person_new, daily_v, daily_new) {
  ok <- function(x, v) !is.na(x) && x == v
  if (is.na(alive)) {
    if (!is.na(conscious) || !is.na(indep_home) || !is.na(shop_v) ||
        !is.na(work_v) || !is.na(social_v)) alive <- 1 else return(NA_integer_)
  }
  if (alive == 0) return(1L)
  independent <- ok(indep_home, 0) || ok(shop_v, 1) || ok(travel_v, 1) ||
    ok(work_v, 1) || ok(social_v, 1)
  if (ok(conscious, 0) && ok(respondent, 2) && !independent) return(2L)
  if (ok(indep_home, 1)) return(if (ok(need_assist, 1)) 3L else 4L)
  if (all(is.na(c(indep_home, shop_v, travel_v, work_v, social_v, person_v, daily_v))))
    return(NA_integer_)
  if (ok(shop_v, 0) && ok(shop_pre, 1)) return(STEWARD_SHOP_TRAVEL_CODE)
  if (ok(travel_v, 0) && ok(travel_pre, 1)) return(STEWARD_SHOP_TRAVEL_CODE)
  if (ok(work_v, 0) && ok(work_major, 1)) return(5L)
  residual <- (ok(social_v, 0) && ok(social_new, 1)) ||
    (ok(person_v, 1) && ok(person_new, 1)) || (ok(daily_v, 1) && ok(daily_new, 1))
  if (residual) {
    if (ok(work_v, 0) || ok(social_restr, 3) || ok(person_freq, 3)) return(6L)
    return(7L)
  }
  if (ok(work_v, 0) && (is.na(work_major) || work_major == 0)) return(7L)
  8L
}

# The outcomes, one row per registry id. `reg` needs patient_id, disposition and
# arrival_date; `fu` is the linked follow-up from steward_link_followup().
steward_score <- function(reg, fu) {
  d <- data.frame(patient_id = as.character(reg$patient_id),
                  disp = .st_num(reg$disposition),
                  arrival_date = suppressWarnings(as.Date(as.character(reg$arrival_date))),
                  stringsAsFactors = FALSE)
  m <- match(d$patient_id, as.character(fu$patients_id_gose))
  g <- function(v) .st_num(fu[[v]])[m]
  d$followup_linked <- !is.na(m)
  d$disp[d$disp %in% 99] <- NA
  d$death_in_hospital <- ifelse(is.na(d$disp), NA_integer_, as.integer(d$disp == 5))

  alive <- g("patient_alive")
  d$gose <- mapply(.st_gose_one, alive, g("consciousness"), g("respondent"),
                   g("independence_home"), g("need_assistance"), g("shop"),
                   g("shop_beforeinj"), g("travel_locally"), g("travel_locally_beforeinj"),
                   g("work"), g("majorchange_work"), g("social_activities"),
                   g("social_restriction"), g("social_restrictions_new"),
                   g("change_personality"), g("freq_change_personality"),
                   g("change_personality_new"), g("dailylife"), g("dailylife_new"))
  died_h <- d$death_in_hospital %in% 1L
  d$gose[is.na(d$gose) & died_h] <- 1L                 # back-fill
  d$gose_death_reconciled <- as.integer(died_h & !is.na(d$gose) & d$gose != 1L)
  d$gose[died_h] <- 1L                                 # hospital record authoritative

  d$death_by_6mo <- ifelse(!is.na(alive) & alive == 0, 1L,
                           ifelse(!is.na(alive) & alive == 1, 0L, NA_integer_))
  d$death_by_6mo[d$gose %in% 1L] <- 1L
  d$death_by_6mo[died_h] <- 1L
  dod <- suppressWarnings(as.Date(as.character(fu$dod)))[m]
  days <- as.numeric(dod - d$arrival_date)
  is_d <- d$death_by_6mo %in% 1L
  late <- is_d & !is.na(days) & days > STEWARD_HORIZON_DAYS & !died_h
  d$death_by_6mo[late] <- 0L
  d$gose[late & d$gose %in% 1L] <- NA_integer_         # late-death release
  d$death_date_implausible <- as.integer(late | (is_d & !is.na(days) & days < 0))
  d$severe_disability <- ifelse(is.na(d$gose), NA_integer_,
                                as.integer(d$gose <= STEWARD_GOSE_SEVERE_CUT))
  d$disp <- NULL
  d
}

# The two outcome-observation indicators approved in D3, from the scored data:
# disposition recorded, and six-month status recorded (the primary outcome is
# known: a death at disposition, or a follow-up giving a GOSE score).
steward_indicators <- function(scored) {
  data.frame(patient_id = scored$patient_id,
             disposition_recorded = as.integer(!is.na(scored$death_in_hospital)),
             six_month_recorded = as.integer(!is.na(scored$severe_disability)),
             stringsAsFactors = FALSE)
}
