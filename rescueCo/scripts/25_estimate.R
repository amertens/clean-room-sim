# ============================================================
# 25. After approval: unblind each lock and estimate
# ============================================================
# For each final lock in rescueCo/registration/approvals_final.csv: unblind()
# with its approved design report and the steward's outcome store, then
# estimate_effect() with the estimand and the candidate the design report
# selected (complete case, the estimator the design stage evaluated). A lock
# that unblind() refuses (no feasible estimand, or a STOP negative-control
# verdict) is reported as refused and never overridden (decision 1). A lock
# with a missing outcome on any of its rows is not estimated, since
# estimate_effect() would drop those rows and redraw the folds; the run then
# ends with an error so that the decision is made by hand. Then the
# longitudinal TMLE sensitivity analysis (lmtp) for the confirmatory
# comparison, which the design stage did not evaluate. No comparator is run
# (decision 6), and the candidate the simulation did not select is not
# estimated (decision 2).
#
# Usage:
#   Rscript rescueCo/scripts/25_estimate.R rehearse   # synthetic outcomes, scratch folder
#   Rscript rescueCo/scripts/25_estimate.R run        # after script 24
# Spec: docs/superpowers/specs/2026-10-20-rescueco-reexecution-0.5.0.md,
# section 7, and the decisions of 2026-10-23.
# ============================================================

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R", "bootstrap.R"))
source("rescueCo/R/utils.R")
source("rescueCo/R/reexec_guards.R")
.libPaths(c("rescueCo/lib", .libPaths()))
suppressMessages(library(cleanTMLE))
cfg <- load_cr_config()
MODE <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(MODE %in% c("rehearse", "run"))
REHEARSE <- MODE == "rehearse"
rx_check_fingerprint()
rx_log(paste("=== 25 estimation,", MODE, "==="))

LOCKD <- file.path(cfg$data$design_dir, "locks", "final")
if (REHEARSE) {
  OUT <- RES <- file.path("rescueCo/data/rehearsal", format(Sys.time(), "%Y%m%d_%H%M%S"))
} else {
  OUT <- "rescueCo/data/estimates/final"                 # patient-level objects
  RES <- file.path(cfg$paths$results, "estimates_050")   # aggregate tables, not reexec_050
  if (length(list.files(RES)))
    stop("Estimates exist already in ", RES, "; script 25 runs once.", call. = FALSE)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# The locks to unblind. A rehearsal takes them from script 23's template, under
# a made-up approver; the real run takes the approvals file, with every hash
# checked and the file committed and pushed.
if (REHEARSE) {
  ap <- utils::read.csv(file.path(REG_DIR, "approvals_final_template.csv"), colClasses = "character")
  ap <- ap[ap$role != "not run", , drop = FALSE]
  ap$approved_by <- "rehearsal (synthetic outcomes)"
} else {
  ap <- rx_read_approvals(LOCKD)
}
cohort_file <- file.path(REG_DIR, "lmtp_cohort.csv")
co_rec <- if (file.exists(cohort_file)) utils::read.csv(cohort_file, colClasses = "character") else NULL

if (REHEARSE) {
  # Synthetic outcomes for the plumbing only, independent of treatment and
  # covariates, so that each arm's risk is known and lmtp's convention can be
  # checked. No real outcome is read.
  ids <- unique(c(unlist(lapply(ap$lock, function(nm)
    load_lock(file.path(LOCKD, paste0(nm, "_lock.rds")))$data$patient_id)),
    if (!is.null(co_rec)) readRDS(co_rec$file[1])$patient_id))
  set.seed(20261023)
  dh <- stats::rbinom(length(ids), 1, 0.03)
  store <- data.frame(patient_id = ids, death_in_hospital = dh,
                      severe_disability = pmax(dh, stats::rbinom(length(ids), 1, 0.2)),
                      death_by_6mo = pmax(dh, stats::rbinom(length(ids), 1, 0.03)),
                      stringsAsFactors = FALSE)
} else {
  store <- readRDS(file.path(cfg$data$outcome_dir, "outcome_store.rds"))
}

# ------------------------------------------------------------
# Complete case, every approved lock
# ------------------------------------------------------------
res <- list(); missing_outcome <- character()
for (i in seq_len(nrow(ap))) {
  nm <- ap$lock[i]
  lk <- load_lock(file.path(LOCKD, paste0(nm, "_lock.rds")))
  d <- readRDS(file.path(LOCKD, paste0(nm, "_dossier.rds")))
  y <- lk$plan$outcome
  row <- data.frame(lock = nm, role = ap$role[i], outcome = y, n_lock = nrow(lk$data))
  oc <- store[store$patient_id %in% lk$data[[lk$id]], c("patient_id", y)]
  ub <- tryCatch(unblind(lk, d, oc, approved_by = ap$approved_by[i]), error = function(e) e)
  if (inherits(ub, "error")) {
    if (!startsWith(conditionMessage(ub), "unblind refused")) stop(ub)
    res[[nm]] <- cbind(row, status = conditionMessage(ub))
    rx_log(paste(nm, conditionMessage(ub)))
    next
  }
  saveRDS(ub, file.path(OUT, paste0(nm, "_unblinded.rds")))
  # unblind() counts design rows with no row in the store, not rows whose
  # outcome is missing, so the missing outcomes are counted here.
  k <- sum(is.na(ub$outcomes[[y]]))
  if (k > 0) {
    res[[nm]] <- cbind(row, status = sprintf("not estimated: outcome missing for %d of %d rows", k, nrow(lk$data)))
    missing_outcome <- c(missing_outcome, nm)
    rx_log(paste(nm, res[[nm]]$status))
    next
  }
  fit <- estimate_effect(ub)
  saveRDS(fit, file.path(OUT, paste0(nm, "_fit.rds")))
  export_design_log(ub, fit, file = file.path(RES, paste0(nm, "_design_log.csv")), format = "muntner")
  res[[nm]] <- cbind(row, status = "estimated", as.data.frame(fit))
  rx_log(paste(nm, "estimated:", fit$estimand))
}
est <- rx_bind_fill(res)
utils::write.csv(est, file.path(RES, "estimates_complete_case.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Longitudinal TMLE (lmtp), confirmatory comparison only
# ------------------------------------------------------------
# Every setting below is fixed before the outcome store exists (decision of
# 2026-10-23). Nodes in time order: disposition recorded (C.1), death at
# disposition (Y.1), six-month status recorded (C.2), the primary outcome
# (Y.2). As lmtp requires for a survival outcome, a death at disposition is
# carried forward to Y.2, and censoring is monotone: a patient whose
# disposition is not recorded is censored from the first node on. The cohort
# is the comparison's rows before the complete-case restriction, saved and
# hashed by script 23 (pick). lmtp 1.5.4 estimates the probability of no event
# by the last node, so each risk is one minus its estimate and the risk
# difference is minus the contrast; the rehearsal checks this convention.
#
# Learners (supervised-learning skill): one type per family, as in the design
# stage's library. Elastic net at three mixing values; MARS of degree 2;
# boosted trees over depth 1 to 5, learning rate 0.1 and 0.3, and 100 or 500
# trees, with no leaf minimum; a neural net with 5 units over four weight
# decays, one step below the design stage's grid, whose smallest decay the
# propensity score usually chose. Each setting enters separately, so its
# cross-validated risk can be checked for edges afterwards, and the ensemble,
# not lmtp's default discrete super learner, is used, as in the design stage.
# Folds: 5 for cross-fitting and 5 inside each super learner, as K and V in
# the plans. Density ratios are trimmed at their 0.999 quantile, lmtp's
# default.
SL.nnet.cr <- function(Y, X, newX, family, obsWeights, size = 5, decay = 0.01, ...) {
  X <- as.matrix(X); newX <- as.matrix(newX)
  mu <- colMeans(X)
  sdv <- apply(X, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  sc <- function(M) sweep(sweep(M, 2, mu), 2, sdv, "/")
  bin <- identical(family$family, "binomial")
  fit <- withr::with_seed(1L, nnet::nnet(x = sc(X), y = Y, size = size, decay = decay, entropy = bin,
                                          linout = !bin, weights = obsWeights, maxit = 500L,
                                          trace = FALSE, MaxNWts = 100000L))
  obj <- structure(list(object = fit, mu = mu, sd = sdv), class = "SL.nnet.cr")
  list(pred = as.numeric(stats::predict(fit, newdata = sc(newX), type = "raw")), fit = obj)
}
predict.SL.nnet.cr <- function(object, newdata, ...) {
  M <- sweep(sweep(as.matrix(newdata), 2, object$mu), 2, object$sd, "/")
  as.numeric(stats::predict(object$object, newdata = M, type = "raw"))
}
g_lrn <- SuperLearner::create.Learner("SL.glmnet", params = list(nfolds = 5),
                                      tune = list(alpha = c(0, 0.5, 1)), env = globalenv())
x_lrn <- SuperLearner::create.Learner("SL.xgboost", params = list(minobspernode = 0, nthread = 1),
                                      tune = list(max_depth = 1:5, shrinkage = c(0.1, 0.3), ntrees = c(100, 500)),
                                      env = globalenv())
e_lrn <- SuperLearner::create.Learner("SL.earth", params = list(degree = 2), env = globalenv())
n_lrn <- SuperLearner::create.Learner("SL.nnet.cr", tune = list(decay = c(1e-4, 1e-3, 1e-2, 1e-1)),
                                      env = globalenv())
SL_LIB <- c(g_lrn$names, e_lrn$names, x_lrn$names, n_lrn$names)
grid <- rbind(
  data.frame(learner = g_lrn$names, setting = paste("alpha", g_lrn$grid$alpha)),
  data.frame(learner = e_lrn$names, setting = "degree 2"),
  data.frame(learner = x_lrn$names, setting = paste("depth", x_lrn$grid$max_depth, "rate",
                                                    x_lrn$grid$shrinkage, "trees", x_lrn$grid$ntrees)),
  data.frame(learner = n_lrn$names, setting = paste("decay", n_lrn$grid$decay)))

conf <- ap[ap$role == "confirmatory", , drop = FALSE]
lm_why <- if (nrow(conf) != 1) "no confirmatory comparison" else
  if (!identical(est$status[est$lock == conf$lock], "estimated")) "the confirmatory comparison was not estimated" else
  if (!identical(conf$primary, "ATE")) "the confirmatory estimand is not the ATE, which lmtp's static interventions estimate here" else
  if (is.null(co_rec) || !conf$lock %in% co_rec$lock) "script 23 recorded no cohort for it" else NA

if (is.na(lm_why)) {
  cr <- co_rec[co_rec$lock == conf$lock, ]
  if (!identical(rx_md5(cr$file), cr$md5))
    stop("The longitudinal TMLE cohort changed after script 23 hashed it.", call. = FALSE)
  co <- readRDS(cr$file)
  covs <- setdiff(names(co), c("patient_id", "A", "disposition_recorded", "six_month_recorded"))
  s <- store[match(co$patient_id, store$patient_id), ]
  if (anyNA(s$patient_id))
    stop("Patients of the longitudinal TMLE cohort are missing from the outcome store.", call. = FALSE)
  C1 <- co$disposition_recorded
  Y1 <- ifelse(C1 == 1, s$death_in_hospital, NA)
  C2 <- C1 * co$six_month_recorded
  Y2 <- ifelse(C2 == 1, ifelse(Y1 %in% 1, 1L, s$severe_disability), NA)
  if (anyNA(Y1[C1 == 1]) || anyNA(Y2[C2 == 1]) || any(Y1 %in% 1 & C2 != 1))
    stop("The outcome store disagrees with the observation indicators the design stage used.", call. = FALSE)
  L <- data.frame(co[, c("A", covs)], C.1 = C1, Y.1 = Y1, C.2 = C2, Y.2 = Y2)
  rx_log(sprintf("lmtp cohort %s: %d patients; censored at the first node %d (of whom %d had six-month status recorded), at the second %d",
                 conf$lock, nrow(L), sum(C1 == 0), sum(C1 == 0 & co$six_month_recorded == 1),
                 sum(C1 == 1 & C2 == 0)))

  ctl <- lmtp::lmtp_control(.learners_outcome_folds = 5, .learners_trt_folds = 5, .trim = 0.999,
                            .discrete = FALSE, .return_full_fits = TRUE)
  run_lmtp <- function(shift) {
    set.seed(20261020)
    lmtp::lmtp_tmle(L, trt = "A", outcome = c("Y.1", "Y.2"), baseline = covs, cens = c("C.1", "C.2"),
                    shift = shift, mtp = FALSE, outcome_type = "survival",
                    learners_outcome = SL_LIB, learners_trt = SL_LIB, folds = 5, control = ctl)
  }
  f1 <- run_lmtp(lmtp::static_binary_on)
  f0 <- run_lmtp(lmtp::static_binary_off)
  saveRDS(list(on = f1, off = f0), file.path(OUT, paste0(conf$lock, "_lmtp.rds")))
  t1 <- generics::tidy(f1); t0 <- generics::tidy(f0)
  ct <- lmtp::lmtp_contrast(f1, ref = f0, type = "additive")$estimates
  need <- c("estimate", "std.error", "conf.low", "conf.high")
  if (!all(need %in% names(t1)) || !all(need %in% names(ct)))
    stop("lmtp's output columns are not as expected: ", paste(names(ct), collapse = ", "), ".", call. = FALSE)
  lm <- data.frame(
    lock = conf$lock, estimand = "ATE", n = nrow(L),
    method = "longitudinal TMLE (lmtp); not evaluated by the design stage",
    risk_treated = 1 - t1$estimate, risk_treated_low = 1 - t1$conf.high, risk_treated_high = 1 - t1$conf.low,
    risk_comparator = 1 - t0$estimate, risk_comparator_low = 1 - t0$conf.high, risk_comparator_high = 1 - t0$conf.low,
    risk_difference = -ct$estimate, std_error = ct$std.error, ci_low = -ct$conf.high, ci_high = -ct$conf.low)
  utils::write.csv(lm, file.path(RES, paste0("lmtp_", conf$lock, ".csv")), row.names = FALSE)

  # Each super learner's cross-validated risk by setting, for the edge check.
  sl_risks <- function(x, path) {
    if (inherits(x, "SuperLearner"))
      return(data.frame(fit = path, learner = names(x$cvRisk), cv_risk = unname(x$cvRisk)))
    if (is.list(x) && !is.data.frame(x) && length(x))
      return(do.call(rbind, Map(sl_risks, x, paste0(path, "/", seq_along(x)))))
    NULL
  }
  er <- rbind(sl_risks(f1, "on"), sl_risks(f0, "off"))
  utils::write.csv(grid, file.path(RES, "lmtp_learner_grid.csv"), row.names = FALSE)
  if (is.null(er)) rx_log("No super learner fits found in the lmtp objects; the edge check cannot run.") else
    utils::write.csv(er, file.path(RES, paste0("lmtp_cv_risk_", conf$lock, ".csv")), row.names = FALSE)

  if (REHEARSE) {
    # The synthetic outcomes are independent of treatment, so both risks
    # should be close to the synthetic risk; a value near one minus it means
    # lmtp's convention differs from the one assumed above.
    truth <- mean(s$severe_disability)
    if (any(abs(c(lm$risk_treated, lm$risk_comparator) - truth) > 0.05))
      stop(sprintf("lmtp risks %.3f and %.3f against a synthetic risk of %.3f; check the survival convention.",
                   lm$risk_treated, lm$risk_comparator, truth), call. = FALSE)
    rx_log(sprintf("Rehearsal: lmtp risks %.3f and %.3f against a synthetic risk of %.3f",
                   lm$risk_treated, lm$risk_comparator, truth))
  }
  rx_log(paste("Longitudinal TMLE written for", conf$lock))
} else {
  rx_log(paste("Longitudinal TMLE not run:", lm_why))
}

# Hashes of every output, committed and pushed before the unblinded analysis
# is opened for the comparison (report recommendation 22).
if (!REHEARSE) {
  files <- c(list.files(RES, full.names = TRUE), list.files(OUT, full.names = TRUE))
  utils::write.csv(data.frame(file = files, md5 = rx_md5(files)),
                   file.path(REG_DIR, "estimates_hashes.csv"), row.names = FALSE)
  rx_log("Wrote rescueCo/registration/estimates_hashes.csv; commit and push it before script 26.")
}
if (length(missing_outcome)) {
  rx_log(paste("Not estimated because of missing outcomes:", paste(missing_outcome, collapse = ", "),
               "; the decision is Andrew's."))
  quit(save = "no", status = 2)
}
rx_log("Script 25 finished.")
