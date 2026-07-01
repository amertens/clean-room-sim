#!/usr/bin/env Rscript
# Diagnostic 5: FORWARD direction via an engineered near-positivity threat.
# Threat amplifies the covariate->treatment slopes so a subgroup approaches
# deterministic treatment (estimated PS -> 0/1). Light truncation inherits the
# exploding IPW tail; heavy truncation caps it. Mild outcome misspecification
# (nonlinear true surface + interaction-aware plasmode Q0, candidates keep a
# main-effects SL.glm Q) gives heavy truncation a small baseline bias cost.
# Candidate fits use the package's own .plasmode_fit_one_candidate for fidelity.
suppressWarnings(suppressMessages(library(pkgload)))
.this_dir <- dirname(normalizePath(sub("^--file=", "",
              commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkgload::load_all(file.path(repo_root, "cleanTMLE"), quiet = TRUE)
fit_one <- cleanTMLE:::.plasmode_fit_one_candidate

generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL, nl = 0) {
  if (!is.null(seed)) set.seed(seed)
  age <- rnorm(n, 55, 10); sex <- rbinom(n, 1, 0.55)
  biomarker <- rnorm(n, 0, 1); comorbidity <- sample(0:2, n, TRUE, c(.5,.3,.2))
  ckd <- rbinom(n, 1, 0.12)
  lp_trt <- -0.5 + overlap_strength * (0.03*(age-55) + 0.8*sex + 0.6*biomarker +
            0.5*ckd + 0.3*comorbidity)
  treatment <- rbinom(n, 1, plogis(lp_trt))
  inter <- nl * (1.2*biomarker*ckd + 0.9*sex*comorbidity + 0.03*(age-55)*biomarker)
  lp_out <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*biomarker + 0.6*ckd +
            0.25*comorbidity + inter + effect_size/0.15*treatment
  event_24 <- rbinom(n, 1, plogis(lp_out))
  nc_outcome <- rbinom(n, 1, plogis(-1 + 0.01*(age-55) + 0.1*sex + 0.15*biomarker))
  data.frame(age = round(age,1), sex = sex, biomarker = round(biomarker,3),
             comorbidity = comorbidity, ckd = ckd, treatment = treatment,
             event_24 = event_24, nc_outcome = nc_outcome)
}

cov <- c("age","sex","biomarker","comorbidity","ckd")
SEED <- 20260530L
Q0LIB <- "SL.glm.interaction"

# Custom near-positivity plasmode-DQ loop. Returns metrics rows per candidate
# per severity. severity s scales the centered PS linear predictor (s=1 baseline).
positivity_stress <- function(ref, cands, severities, reps, es = 0.05,
                              q0_lib = Q0LIB) {
  n <- nrow(ref); A <- ref$treatment
  # Q0 (covariate-only) capturing nonlinearity:
  Q0 <- SuperLearner::SuperLearner(Y = ref$event_24,
          X = ref[, cov, drop = FALSE], family = binomial(),
          SL.library = q0_lib, env = asNamespace("SuperLearner"))
  p_base <- pmin(pmax(as.numeric(Q0$SL.predict), 0.001), 0.999)
  ps_mod <- glm(reformulate(cov, "treatment"), ref, family = binomial())
  lp0 <- as.numeric(predict(ps_mod, type = "link"))
  lp_bar <- mean(lp0); lpc <- lp0 - lp_bar
  out <- list()
  for (s in severities) {
    ps_s <- plogis(lp_bar + s * lpc)
    rr <- vector("list", reps)
    for (r in seq_len(reps)) {
      set.seed(SEED + r + round(s * 1000))
      A_rep <- rbinom(n, 1L, ps_s)
      p1 <- pmin(pmax(p_base + es, 0.001), 0.999); p0 <- p_base
      p_obs <- ifelse(A_rep == 1, p1, p0)
      Y_sim <- rbinom(n, 1L, p_obs)
      truth <- mean(p1) - mean(p0)
      res <- lapply(cands, function(cd) tryCatch(
        fit_one(Y_sim = Y_sim, A = A_rep, W_data = ref, treatment = "treatment",
                covariates = cov, cand = cd, n = n),
        error = function(e) list(est = NA, se = NA, ci_lower = NA, ci_upper = NA)))
      names(res) <- vapply(cands, function(x) x$candidate_id, "")
      rr[[r]] <- c(res, list(.truth = truth))
    }
    truth_v <- vapply(rr, function(x) x$.truth, numeric(1))
    for (cid in names(res)) {
      est <- vapply(rr, function(x) x[[cid]]$est, numeric(1))
      se  <- vapply(rr, function(x) x[[cid]]$se,  numeric(1))
      lo  <- vapply(rr, function(x) x[[cid]]$ci_lower, numeric(1))
      hi  <- vapply(rr, function(x) x[[cid]]$ci_upper, numeric(1))
      ok <- !is.na(est)
      out[[length(out)+1]] <- data.frame(
        scenario = if (s == 1) "none" else "near_positivity",
        level = sprintf("slope_x%.1f", s), candidate = cid,
        bias = mean(est[ok]-truth_v[ok]),
        rmse = sqrt(mean((est[ok]-truth_v[ok])^2)),
        coverage = mean(lo[ok] <= truth_v[ok] & truth_v[ok] <= hi[ok]),
        maxwt = NA, n_converged = sum(ok), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

for (ov in c(1.2, 1.4)) for (nl in c(0.6, 1.0)) {
  ref <- generate_data(2000L, ov, -0.05, seed = SEED, nl = nl)
  ps <- predict(glm(reformulate(cov,"treatment"), ref, family=binomial()), type="response")
  cands <- list(
    tmle_candidate("aggressive", g_library="SL.glm", truncation=0.001),
    tmle_candidate("middle",     g_library="SL.glm", truncation=0.025),
    tmle_candidate("robust",     g_library="SL.glm", truncation=0.10))
  m <- positivity_stress(ref, cands, severities = c(1.0, 2.0, 3.0, 4.0), reps = 60L)
  cat(sprintf("\n=== ov=%.1f nl=%.1f  PS[%.4f,%.4f] frac<.10=%.1f%% ===\n",
              ov, nl, min(ps), max(ps), 100*mean(ps<.10)))
  w <- reshape(m[, c("level","candidate","rmse")], idvar="level",
               timevar="candidate", direction="wide")
  names(w) <- sub("rmse.","",names(w))
  print(w[, c("level","aggressive","middle","robust")], row.names = FALSE)
  base <- m[m$scenario=="none",]; deg <- m[m$scenario!="none",]
  bwin <- base$candidate[which.min(base$rmse)]
  wc <- tapply(deg$rmse, deg$candidate, max)
  mmwin <- names(wc)[which.min(wc)]
  cat(sprintf("  baseline winner=%s   worstcase: aggr=%.5f mid=%.5f rob=%.5f  minmax=%s  %s\n",
      bwin, wc["aggressive"], wc["middle"], wc["robust"], mmwin,
      if (bwin=="aggressive" && mmwin=="robust") "**FORWARD!**"
      else if (bwin!=mmwin) "diverge(other)" else "agree"))
}
cat("\nDONE\n")
