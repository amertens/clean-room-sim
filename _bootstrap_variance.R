#!/usr/bin/env Rscript
# ============================================================================
# Bootstrap vs influence-function variance comparison
#
# Three scenarios that together tell a coherent story:
#   C (very good overlap, overlap_strength=0.25): Match_TMLE IF-SE is most
#     anti-conservative; bootstrap provides the largest gain. This is the
#     canonical Abadie-Imbens (2008) setting.
#   A (good overlap, overlap_strength=0.5):  moderate anti-conservatism;
#     bootstrap corrects modestly.
#   B (marginal overlap, overlap_strength=1.5): caliper self-selects a
#     well-overlapped sub-cohort; Match_TMLE IF-SE approximately calibrated;
#     IPTW IF-SE theoretically conservative (PS treated as known).
#
# Design:
#   Estimators: IPTW (stabilised), TMLE, TMLE_CF (cross-fitted), Match_TMLE
#   MC reps:    100 per scenario
#   Bootstrap:  B = 500 draws, percentile CI
#
# TMLE and TMLE_CF are included because they are the estimators that are
# anti-conservative under marginal overlap (influence-function SE/SD ~ 0.8);
# the study tests whether the nonparametric bootstrap repairs their coverage.
#
# Implementation: fully self-contained GLM + pure-R nearest-neighbour
# matching (same algorithm as cleanTMLE::run_match_workflow) + manual
# TMLE targeting step. No subprocess-prone packages inside any loop.
#
# Output: results_new/bootstrap_variance.rds
#         results_new/bootstrap_variance.csv
# ============================================================================

.flush <- function() if (!interactive()) flush(stdout())
cat("Bootstrap variance comparison\n")
cat(sprintf("Started: %s\n\n", format(Sys.time())))
.flush()

# ── Configuration ─────────────────────────────────────────────────────────────

N_MC    <- 100L
B_BOOT  <- 200L   # reduced from 500 for tractable runtime; 200 draws suffice
                  # for a percentile interval and a coverage estimate
N_OBS   <- 2000L
N_TRUTH <- 500000L
TRUNC   <- 0.05
SEED    <- 2026L
COVARS  <- c("age", "sex", "biomarker", "comorbidity", "ckd")
TREAT   <- "treatment"
OUTC    <- "event_24"

out_dir <- file.path(getwd(), "results_new")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── DGP ───────────────────────────────────────────────────────────────────────

generate_data <- function(n, overlap_strength = 0.5) {
  age         <- rnorm(n, 55, 10)
  sex         <- rbinom(n, 1, 0.55)
  biomarker   <- rnorm(n)
  comorbidity <- sample(0:2, n, TRUE, c(0.5, 0.3, 0.2))
  ckd         <- rbinom(n, 1, 0.12)
  lp_trt <- -0.5 + overlap_strength *
    (0.03*(age-55) + 0.8*sex + 0.6*biomarker + 0.5*ckd + 0.3*comorbidity)
  treatment  <- rbinom(n, 1, plogis(lp_trt))
  lp_out <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*biomarker +
    0.6*ckd + 0.25*comorbidity + (-0.05/0.15)*treatment
  event_24 <- rbinom(n, 1, plogis(lp_out))
  data.frame(age=round(age,1), sex=sex, biomarker=round(biomarker,3),
             comorbidity=comorbidity, ckd=ckd, treatment=treatment,
             event_24=event_24, stringsAsFactors=FALSE)
}

# True marginal RD (does not depend on overlap_strength — only outcome model)
compute_truth <- function(n) {
  set.seed(SEED)
  age <- rnorm(n,55,10); sex <- rbinom(n,1,.55)
  bio <- rnorm(n); com <- sample(0:2,n,T,c(.5,.3,.2)); ckd <- rbinom(n,1,.12)
  lp  <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*bio + 0.6*ckd + 0.25*com
  mean(plogis(lp + (-0.05/0.15))) - mean(plogis(lp))
}

cat("Computing ground truth...\n"); .flush()
truth_rd <- compute_truth(N_TRUTH)
cat(sprintf("  true RD = %.5f  (same for A and B by DGP)\n\n", truth_rd))
.flush()

# ── Pure-R 1:1 nearest-neighbour matching on logit-PS ────────────────────────
# Algorithm mirrors cleanTMLE::run_match_workflow exactly.

.match_nn <- function(A, logit_ps, caliper = NULL) {
  if (is.null(caliper)) caliper <- 0.2 * sd(logit_ps)
  treated_idx  <- which(A == 1)
  control_idx  <- which(A == 0)
  n_treated    <- length(treated_idx)
  used_ctrl    <- logical(length(control_idx))
  matched_ctrl <- integer(n_treated)

  for (i in seq_len(n_treated)) {
    dists <- abs(logit_ps[treated_idx[i]] - logit_ps[control_idx])
    dists[used_ctrl] <- Inf
    best <- which.min(dists)
    if (dists[best] <= caliper) {
      matched_ctrl[i] <- control_idx[best]
      used_ctrl[best] <- TRUE
    } else {
      matched_ctrl[i] <- NA_integer_
    }
  }
  valid <- !is.na(matched_ctrl)
  c(treated_idx[valid], matched_ctrl[valid])
}

# ── Manual TMLE (GLM only, no packages) ──────────────────────────────────────
# Returns list(estimate, se, ci_lower, ci_upper).

.tmle_glm <- function(Y, A, W, truncate = TRUNC) {
  Wdf <- as.data.frame(W)
  n   <- length(Y)

  g_mod <- glm(A ~ ., data = cbind(A=A, Wdf), family = binomial())
  g     <- pmin(pmax(predict(g_mod, type="response"), truncate), 1-truncate)

  q_mod <- glm(Y ~ ., data = cbind(Y=Y, A=A, Wdf), family = binomial())
  Q1    <- predict(q_mod, newdata = cbind(A=1, Wdf), type="response")
  Q0    <- predict(q_mod, newdata = cbind(A=0, Wdf), type="response")
  Q_AW  <- pmin(pmax(A*Q1 + (1-A)*Q0, 1e-6), 1-1e-6)

  H   <- A/g - (1-A)/(1-g)
  eps <- tryCatch(
    coef(glm(Y ~ H + offset(qlogis(Q_AW)), family=binomial()))["H"],
    error = function(e) 0
  )

  Q1s <- plogis(qlogis(pmin(pmax(Q1,1e-6),1-1e-6)) + eps/g)
  Q0s <- plogis(qlogis(pmin(pmax(Q0,1e-6),1-1e-6)) - eps/(1-g))
  ate <- mean(Q1s) - mean(Q0s)

  eif <- (A/g - (1-A)/(1-g)) * (Y - Q_AW) + (Q1s - Q0s) - ate
  se  <- sqrt(var(eif) / n)
  list(estimate=ate, se=se, ci_lower=ate-1.96*se, ci_upper=ate+1.96*se)
}

# ── Estimator fits ────────────────────────────────────────────────────────────

.iptw_fit <- function(dat) {
  A <- dat[[TREAT]]; Y <- dat[[OUTC]]; W <- dat[,COVARS,drop=FALSE]
  g_mod <- glm(reformulate(COVARS, TREAT), data=dat, family=binomial())
  ps <- pmin(pmax(predict(g_mod, type="response"), TRUNC), 1-TRUNC)
  pA <- mean(A)
  w  <- ifelse(A==1, pA/ps, (1-pA)/(1-ps))
  mu1 <- weighted.mean(Y[A==1], w[A==1])
  mu0 <- weighted.mean(Y[A==0], w[A==0])
  est <- mu1 - mu0
  eif <- (A/ps-(1-A)/(1-ps))*(Y-ifelse(A==1,mu1,mu0)) + (mu1-mu0) - est
  se  <- sqrt(var(eif)/nrow(dat))
  list(estimate=est, se=se, ci_lower=est-1.96*se, ci_upper=est+1.96*se)
}

.match_tmle_fit <- function(dat) {
  A <- dat[[TREAT]]; Y <- dat[[OUTC]]; W <- dat[,COVARS,drop=FALSE]
  g_mod <- glm(reformulate(COVARS, TREAT), data=dat, family=binomial())
  ps    <- pmin(pmax(predict(g_mod, type="response"), TRUNC), 1-TRUNC)
  idx   <- .match_nn(A, qlogis(ps))
  if (length(idx) < 10) stop("Too few matched pairs")
  m     <- dat[idx, ]
  .tmle_glm(m[[OUTC]], m[[TREAT]], m[,COVARS,drop=FALSE])
}

# Plain TMLE (GLM nuisances, full sample) — wraps the manual TMLE.
.tmle_fit <- function(dat) {
  .tmle_glm(dat[[OUTC]], dat[[TREAT]], dat[, COVARS, drop=FALSE])
}

# Cross-fitted TMLE (2-fold, GLM nuisances, pooled targeting). Added so the
# variance study covers the estimators that are anti-conservative under
# marginal overlap, not only IPTW and Match_TMLE.
.tmle_cf_fit <- function(dat, V = 2L) {
  Y <- dat[[OUTC]]; A <- dat[[TREAT]]; W <- dat[, COVARS, drop=FALSE]; n <- length(Y)
  folds <- sample(rep(seq_len(V), length.out = n))
  g <- Q1 <- Q0 <- numeric(n)
  for (v in seq_len(V)) {
    tr <- folds != v; te <- folds == v
    gm <- glm(A ~ ., data = cbind(A=A, W)[tr, , drop=FALSE], family=binomial())
    g[te]  <- predict(gm, newdata = W[te, , drop=FALSE], type="response")
    qm <- glm(Y ~ ., data = cbind(Y=Y, A=A, W)[tr, , drop=FALSE], family=binomial())
    Q1[te] <- predict(qm, newdata = cbind(A=1, W[te, , drop=FALSE]), type="response")
    Q0[te] <- predict(qm, newdata = cbind(A=0, W[te, , drop=FALSE]), type="response")
  }
  g    <- pmin(pmax(g, TRUNC), 1-TRUNC)
  Q_AW <- pmin(pmax(A*Q1 + (1-A)*Q0, 1e-6), 1-1e-6)
  H    <- A/g - (1-A)/(1-g)
  eps  <- tryCatch(coef(glm(Y ~ H + offset(qlogis(Q_AW)), family=binomial()))["H"],
                   error = function(e) 0)
  Q1s <- plogis(qlogis(pmin(pmax(Q1,1e-6),1-1e-6)) + eps/g)
  Q0s <- plogis(qlogis(pmin(pmax(Q0,1e-6),1-1e-6)) - eps/(1-g))
  ate <- mean(Q1s) - mean(Q0s)
  eif <- H*(Y - Q_AW) + (Q1s - Q0s) - ate
  se  <- sqrt(var(eif)/n)
  list(estimate=ate, se=se, ci_lower=ate-1.96*se, ci_upper=ate+1.96*se)
}

# Generic bootstrap SE wrapper
.boot_se <- function(dat, fn, B=B_BOOT, seed=1L) {
  set.seed(seed)
  n <- nrow(dat)
  boots <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n, n, replace=TRUE)
    tryCatch(fn(dat[idx,,drop=FALSE])$estimate, error=function(e) NA_real_)
  }, numeric(1))
  boots <- boots[is.finite(boots)]
  if (length(boots) < 10L) return(list(se=NA_real_, ci=c(NA_real_, NA_real_)))
  list(se=sd(boots), ci=unname(quantile(boots, c(0.025, 0.975))))
}

# ── MC loop ───────────────────────────────────────────────────────────────────

run_scenario <- function(sc_label, overlap_strength) {
  cat(sprintf("=== %s ===\n", sc_label)); .flush()
  t0   <- Sys.time()
  rows <- vector("list", N_MC)

  for (i in seq_len(N_MC)) {
    dat <- generate_data(N_OBS, overlap_strength)

    iptw_if   <- tryCatch(.iptw_fit(dat),         error=function(e) NULL)
    iptw_boot <- tryCatch(.boot_se(dat,.iptw_fit,seed=i), error=function(e) NULL)
    tmle_if   <- tryCatch(.tmle_fit(dat),         error=function(e) NULL)
    tmle_boot <- tryCatch(.boot_se(dat,.tmle_fit,seed=i), error=function(e) NULL)
    tcf_if    <- tryCatch(.tmle_cf_fit(dat),      error=function(e) NULL)
    tcf_boot  <- tryCatch(.boot_se(dat,.tmle_cf_fit,seed=i), error=function(e) NULL)
    mt_if     <- tryCatch(.match_tmle_fit(dat),   error=function(e) NULL)
    mt_boot   <- tryCatch(.boot_se(dat,.match_tmle_fit,seed=i), error=function(e) NULL)

    make_row <- function(method, if_fit, boot_fit) {
      if (is.null(if_fit)) return(NULL)
      cov_if   <- as.integer(!is.na(if_fit$ci_lower) &&
                               if_fit$ci_lower <= truth_rd &&
                               truth_rd <= if_fit$ci_upper)
      cov_boot <- if (!is.null(boot_fit) && all(is.finite(boot_fit$ci)))
        as.integer(boot_fit$ci[1] <= truth_rd && truth_rd <= boot_fit$ci[2])
      else NA_integer_
      data.frame(rep=i, method=method,
                 estimate=if_fit$estimate, se_if=if_fit$se,
                 se_boot=if(!is.null(boot_fit)) boot_fit$se else NA_real_,
                 covers_if=cov_if, covers_boot=cov_boot,
                 stringsAsFactors=FALSE)
    }

    rows[[i]] <- rbind(make_row("IPTW",       iptw_if, iptw_boot),
                       make_row("TMLE",       tmle_if, tmle_boot),
                       make_row("TMLE_CF",    tcf_if,  tcf_boot),
                       make_row("Match_TMLE", mt_if,   mt_boot))

    if (i %% 10 == 0) {
      el  <- as.numeric(difftime(Sys.time(), t0, units="mins"))
      rem <- if (i>1) el/i*(N_MC-i) else NA
      cat(sprintf("  rep %d/%d  (%.1f min", i, N_MC, el))
      if (!is.na(rem)) cat(sprintf(", ~%.0f min remaining", rem))
      cat(")\n"); .flush()
    }
  }
  do.call(rbind, Filter(Negate(is.null), rows))
}

set.seed(SEED + 99L)
results_C <- run_scenario("Scenario C: Very Good Overlap", 0.25)
results_A <- run_scenario("Scenario A: Good Overlap",      0.5)
results_B <- run_scenario("Scenario B: Marginal Overlap",  1.5)

# ── Summary ───────────────────────────────────────────────────────────────────

build_summary <- function(df, truth) {
  do.call(rbind, lapply(unique(df$method), function(m) {
    d  <- df[df$method==m & !is.na(df$estimate), ]
    if (nrow(d)==0) return(NULL)
    emp <- sd(d$estimate)
    mif  <- mean(d$se_if,   na.rm=TRUE)
    mbt  <- mean(d$se_boot, na.rm=TRUE)
    data.frame(method=m, n_reps=nrow(d),
               bias=round(mean(d$estimate)-truth, 5),
               emp_sd=round(emp,4),
               mean_se_if=round(mif,4),  se_sd_ratio_if=round(mif/emp,3),
               coverage_if=round(mean(d$covers_if,na.rm=TRUE),3),
               mean_se_boot=round(mbt,4), se_sd_ratio_boot=round(mbt/emp,3),
               coverage_boot=round(mean(d$covers_boot,na.rm=TRUE),3),
               stringsAsFactors=FALSE)
  }))
}

smC <- build_summary(results_C, truth_rd)
smA <- build_summary(results_A, truth_rd)
smB <- build_summary(results_B, truth_rd)

cat("\n=== SUMMARY: Scenario C (Very Good Overlap) ===\n"); print(smC)
cat("\n=== SUMMARY: Scenario A (Good Overlap) ===\n"); print(smA)
cat("\n=== SUMMARY: Scenario B (Marginal Overlap) ===\n"); print(smB)
.flush()

combined <- rbind(cbind(scenario="C: Very Good Overlap", smC),
                  cbind(scenario="A: Good Overlap",      smA),
                  cbind(scenario="B: Marginal Overlap",  smB))

saveRDS(list(results_C=results_C, results_A=results_A, results_B=results_B,
             summary_C=smC, summary_A=smA, summary_B=smB, combined=combined,
             truth=truth_rd,
             config=list(N_MC=N_MC,B_BOOT=B_BOOT,N_OBS=N_OBS,
                         TRUNC=TRUNC,SEED=SEED)),
        file.path(out_dir, "bootstrap_variance.rds"))
write.csv(combined, file.path(out_dir, "bootstrap_variance.csv"), row.names=FALSE)

cat(sprintf("\nCOMPLETED: %s\n", format(Sys.time())))
cat("BOOTSTRAP_VARIANCE_DONE\n")
