# ============================================================
# Multi-arm Stage 2: the design stage, outcome-blind
# ============================================================
# For every lock from script 10: the SuperLearner propensity fit (rwe_wide
# library, 10-fold CV) with the fitted ensemble weights recorded, the
# support assessment with its graded verdict, the estimand feasibility
# table, the profile of who falls outside the band, the negative-control
# ladder (full ambulance cohort, transfers excluded, transfers and prior
# care excluded, on the ambulance contrasts), the care-process collider
# check, and the assembled design report.
#
# No outcome column is present in any lock (script 10 stores outcomes
# separately), so nothing here can touch the treatment-outcome association.
#
# Checkpointed per lock: rerun resumes; set MULTIARM_FORCE=1 to refit.
# ============================================================

suppressMessages(library(cleanTMLE))
suppressMessages(library(SuperLearner))

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                 "bootstrap.R"))
source("rescueCo/R/utils.R")
cr_log("=== Multi-arm Stage 2: design (outcome-blind) ===")

OUT <- "rescueCo/results/multiarm"
locks <- readRDS(file.path(OUT, "multiarm_locks.rds"))
FORCE <- nzchar(Sys.getenv("MULTIARM_FORCE"))

# The design-stage nuisance library: the smooth wide-registry preset.
SL_LIB <- build_sl_library(role = "g", n_eff = 3000, preset = "rwe_wide")$library
cr_log(paste("PS library:", paste(vapply(SL_LIB, paste, character(1),
                                         collapse = "+"), collapse = ", ")))

# Interpretable variables for the who-is-unsupported profile (the script 54
# idea: standardised differences of removed vs kept on variables a reader
# can interpret).
PROFILE_VARS <- c("age", "sex_Female", "mech_inj_2", "inj_place_Private.House..Home",
                  "inj_distance", "eiss_ais", "head_injury_severity",
                  "n_regions_injured", "months_since_start",
                  "vitals_1_hr_missing", "sbp_missing", "pulse_Not.examined",
                  "intent_Intentional..assault.homicide..",
                  "is_transfer", "sought_prior_care")

support_rows <- list(); feas_rows <- list(); nd_rows <- list()
vr_rows <- list(); who_rows <- list(); wt_rows <- list()

for (cn in names(locks)) {
  lk <- locks[[cn]]
  ck <- file.path(OUT, paste0("design_", cn, ".rds"))
  if (file.exists(ck) && !FORCE) {
    cr_log(paste("[", cn, "] checkpoint exists; loading"))
    st <- readRDS(ck)
  } else {
    cr_log(paste("[", cn, "] fitting propensity (n =", nrow(lk$data), ")"))
    lk_fit <- lk
    lk_fit$sl_library <- SL_LIB
    # V matches the main pipeline's standalone propensity fits: 5 folds on
    # the large multi-arm cohorts (analysis/42), 10 on the primary cohort
    # (analysis/04).
    v_folds <- if (nrow(lk$data) > 4000) 5L else 10L
    t0 <- Sys.time()
    psf <- fit_ps_superlearner(lk_fit, truncate = 0.01, cv_folds = v_folds)
    cr_log(paste("[", cn, "] PS fit in",
                 round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
                 "min"))

    sup <- assess_support(psf, tree_search = TRUE, tree_min_n = 50L)
    fea <- estimand_feasibility(psf)
    who <- who_is_unsupported(psf, vars = PROFILE_VARS)

    # Negative-control ladder on the ambulance contrasts; single full-cohort
    # rung elsewhere (the restriction ladder is about the ambulance cohort).
    ncl <- NULL
    if (cn %in% c("C4", "PRIMARY")) {
      restr <- if (cn == "C4") list(
        `transfers excluded` = lk$data$is_transfer == 0,
        `transfers and prior care excluded` = lk$data$pre_hospital_care == 0)
      else list(`(already restricted)` = rep(TRUE, nrow(lk$data)))
      ncl <- tryCatch(run_negative_control_ladder(
        lk, restrictions = restr, method = "tmle", ps_method = "glm",
        verbose = FALSE), error = function(e) {
          cr_log(paste("[", cn, "] NC ladder failed:", conditionMessage(e)))
          NULL })
    } else {
      ncl <- tryCatch(run_negative_control_ladder(
        lk, restrictions = list(), method = "tmle", ps_method = "glm",
        verbose = FALSE), error = function(e) NULL)
    }

    # Care-process collider check, conditional on severity.
    chk <- tryCatch(check_process_indicators(
      lk, indicators = c("vitals_1_hr_missing", "sbp_missing",
                         "pulse_Not.examined"),
      condition_on = c("age", "eiss_ais", "n_regions_injured",
                       "head_injury_severity", "eyes_gcs", "verbal_gcs",
                       "motor_gcs")), error = function(e) NULL)

    rep_ <- design_report(lk, sup, fea, nc_ladder = ncl,
                          extra = list(collider_check = chk))
    st <- list(contrast = cn, label = lk$contrast$label, ps = psf,
               support = sup, feasibility = fea, who = who, nc_ladder = ncl,
               collider = chk, report = rep_,
               sl_weights = tryCatch(psf$sl_fit$coef, error = function(e) NULL))
    saveRDS(st, ck)
    cr_log(paste("[", cn, "] design checkpoint saved"))
  }

  s <- st$support$summary
  support_rows[[cn]] <- cbind(contrast = cn, label = st$label, s,
                              escalated = st$support$escalated,
                              stringsAsFactors = FALSE)
  feas_rows[[cn]] <- cbind(contrast = cn, st$feasibility$table,
                           stringsAsFactors = FALSE)
  if (!is.null(st$support$near_deterministic))
    nd_rows[[cn]] <- cbind(contrast = cn, st$support$near_deterministic,
                           stringsAsFactors = FALSE)
  if (!is.null(st$support$violation_regions))
    vr_rows[[cn]] <- cbind(contrast = cn, st$support$violation_regions,
                           stringsAsFactors = FALSE)
  if (!is.null(st$who))
    who_rows[[cn]] <- cbind(contrast = cn, st$who, stringsAsFactors = FALSE)
  if (!is.null(st$sl_weights))
    wt_rows[[cn]] <- data.frame(contrast = cn,
                                learner = names(st$sl_weights),
                                weight = round(as.numeric(st$sl_weights), 4),
                                stringsAsFactors = FALSE)

  # Mirrored support plot for the dossier.
  p <- plot(st$support)
  ggplot2::ggsave(file.path(OUT, paste0("fig_support_", cn, ".png")), p,
                  width = 7.5, height = 4.5, dpi = 150)
  cr_log(paste("[", cn, "]", st$support$verdict, "-",
               st$report$recommendation))
}

utils::write.csv(do.call(rbind, support_rows),
                 file.path(OUT, "support_by_contrast.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, feas_rows),
                 file.path(OUT, "estimand_feasibility.csv"),
                 row.names = FALSE)
if (length(nd_rows))
  utils::write.csv(do.call(rbind, nd_rows),
                   file.path(OUT, "near_deterministic_covariates.csv"),
                   row.names = FALSE)
if (length(vr_rows))
  utils::write.csv(do.call(rbind, vr_rows),
                   file.path(OUT, "violation_regions.csv"),
                   row.names = FALSE)
if (length(who_rows))
  utils::write.csv(do.call(rbind, who_rows),
                   file.path(OUT, "who_is_unsupported.csv"),
                   row.names = FALSE)
if (length(wt_rows))
  utils::write.csv(do.call(rbind, wt_rows),
                   file.path(OUT, "sl_propensity_weights.csv"),
                   row.names = FALSE)

# Negative-control ladders and collider checks, gathered.
ncl_rows <- list(); chk_rows <- list()
for (cn in names(locks)) {
  st <- readRDS(file.path(OUT, paste0("design_", cn, ".rds")))
  if (!is.null(st$nc_ladder))
    ncl_rows[[cn]] <- cbind(contrast = cn, st$nc_ladder$table,
                            stringsAsFactors = FALSE)
  if (!is.null(st$collider))
    chk_rows[[cn]] <- cbind(contrast = cn, st$collider,
                            stringsAsFactors = FALSE)
}
if (length(ncl_rows))
  utils::write.csv(do.call(rbind, ncl_rows),
                   file.path(OUT, "nc_ladder.csv"), row.names = FALSE)
if (length(chk_rows))
  utils::write.csv(do.call(rbind, chk_rows),
                   file.path(OUT, "collider_check.csv"), row.names = FALSE)

cr_log("Multi-arm Stage 2 complete.")
