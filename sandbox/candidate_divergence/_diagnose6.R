#!/usr/bin/env Rscript
# Diagnostic 6: concentrate the Q-misspecification in the low-PS tail so heavy
# truncation pays a real BASELINE bias cost (aggressive wins baseline), while
# the near-positivity threat still blows up aggressive (robust wins worst-case).
# Lever: biomarker is a strong PS driver AND carries a strong age:biomarker
# outcome interaction that the candidates' main-effects SL.glm Q cannot
# represent; the interaction-aware plasmode Q0 puts it into the synthetic Y.
suppressWarnings(suppressMessages(library(pkgload)))
.this_dir <- dirname(normalizePath(sub("^--file=", "",
              commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])))
repo_root <- normalizePath(file.path(.this_dir, "..", ".."))
pkgload::load_all(file.path(repo_root, "cleanTMLE"), quiet = TRUE)
fit_one <- cleanTMLE:::.plasmode_fit_one_candidate

generate_data <- function(n, overlap_strength = 1.0, effect_size = -0.05,
                          seed = NULL, ps_bio = 1.2, out_ab = 0.15) {
  if (!is.null(seed)) set.seed(seed)
  age <- rnorm(n, 55, 10); sex <- rbinom(n, 1, 0.55)
  biomarker <- rnorm(n, 0, 1); comorbidity <- sample(0:2, n, TRUE, c(.5,.3,.2))
  ckd <- rbinom(n, 1, 0.12)
  lp_trt <- -0.5 + overlap_strength * (0.03*(age-55) + 0.8*sex + ps_bio*biomarker +
            0.5*ckd + 0.3*comorbidity)
  treatment <- rbinom(n, 1, plogis(lp_trt))
  # Outcome: strong age:biomarker interaction (captured by SL.glm.interaction,
  # missed by main-effects glm) concentrated where biomarker is extreme = where
  # PS is extreme. Biomarker also has a real main outcome effect (confounder).
  lp_out <- -2.0 + 0.015*(age-55) + 0.3*sex + 0.5*biomarker + 0.6*ckd +
            0.25*comorbidity + out_ab*(age-55)*biomarker +
            effect_size/0.15*treatment
  event_24 <- rbinom(n, 1, plogis(lp_out))
  nc_outcome <- rbinom(n, 1, plogis(-1 + 0.01*(age-55) + 0.1*sex + 0.15*biomarker))
  data.frame(age = round(age,1), sex = sex, biomarker = round(biomarker,3),
             comorbidity = comorbidity, ckd = ckd, treatment = treatment,
             event_24 = event_24, nc_outcome = nc_outcome)
}

cov <- c("age","sex","biomarker","comorbidity","ckd")
SEED <- 20260530L
Q0LIB <- "SL.glm.interaction"

positivity_stress <- function(ref, cands, severities, reps, es = 0.05) {
  n <- nrow(ref)
  Q0 <- SuperLearner::SuperLearner(Y = ref$event_24, X = ref[, cov, drop = FALSE],
          family = binomial(), SL.library = Q0LIB, env = asNamespace("SuperLearner"))
  p_base <- pmin(pmax(as.numeric(Q0$SL.predict), 0.001), 0.999)
  ps_mod <- glm(reformulate(cov, "treatment"), ref, family = binomial())
  lp0 <- as.numeric(predict(ps_mod, type = "link")); lp_bar <- mean(lp0); lpc <- lp0 - lp_bar
  out <- list()
  for (s in severities) {
    ps_s <- plogis(lp_bar + s * lpc)
    rr <- vector("list", reps)
    for (r in seq_len(reps)) {
      set.seed(SEED + r + round(s * 1000))
      A_rep <- rbinom(n, 1L, ps_s)
      p1 <- pmin(pmax(p_base + es, 0.001), 0.999); p0 <- p_base
      Y_sim <- rbinom(n, 1L, ifelse(A_rep == 1, p1, p0))
      truth <- mean(p1) - mean(p0)
      res <- lapply(cands, function(cd) tryCatch(
        fit_one(Y_sim = Y_sim, A = A_rep, W_data = ref, treatment = "treatment",
                covariates = cov, cand = cd, n = n),
        error = function(e) list(est = NA, se = NA, ci_lower = NA, ci_upper = NA)))
      names(res) <- vapply(cands, function(x) x$candidate_id, "")
      rr[[r]] <- c(res, list(.truth = truth))
    }
    tv <- vapply(rr, function(x) x$.truth, numeric(1))
    for (cid in names(res)) {
      est <- vapply(rr, function(x) x[[cid]]$est, numeric(1)); ok <- !is.na(est)
      out[[length(out)+1]] <- data.frame(
        scenario = if (s == 1) "none" else "near_positivity",
        level = sprintf("x%.1f", s), candidate = cid,
        bias = mean(est[ok]-tv[ok]), rmse = sqrt(mean((est[ok]-tv[ok])^2)),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, out); out
}

grid <- expand.grid(ov = c(0.9, 1.2), ps_bio = c(1.2, 1.8), out_ab = c(0.12, 0.20))
for (i in seq_len(nrow(grid))) {
  ov <- grid$ov[i]; pb <- grid$ps_bio[i]; oa <- grid$out_ab[i]
  ref <- generate_data(2000L, ov, -0.05, seed = SEED, ps_bio = pb, out_ab = oa)
  ps <- predict(glm(reformulate(cov,"treatment"), ref, family=binomial()), type="response")
  cands <- list(
    tmle_candidate("aggressive", g_library="SL.glm", truncation=0.001),
    tmle_candidate("middle",     g_library="SL.glm", truncation=0.025),
    tmle_candidate("robust",     g_library="SL.glm", truncation=0.10))
  m <- positivity_stress(ref, cands, severities = c(1.0, 2.5, 4.0), reps = 70L)
  base <- m[m$scenario=="none",]; deg <- m[m$scenario!="none",]
  br <- setNames(base$rmse, base$candidate); bb <- setNames(base$bias, base$candidate)
  wc <- tapply(deg$rmse, deg$candidate, max)
  bwin <- base$candidate[which.min(base$rmse)]; mmwin <- names(wc)[which.min(wc)]
  fwd <- bwin=="aggressive" && mmwin=="robust"
  cat(sprintf("ov=%.1f psbio=%.1f outab=%.2f frac<.10=%.1f%% | base rmse a=%.5f m=%.5f r=%.5f (bias a=%.4f r=%.4f) | wc a=%.4f r=%.4f | bwin=%s mm=%s %s\n",
      ov, pb, oa, 100*mean(ps<.10), br["aggressive"], br["middle"], br["robust"],
      bb["aggressive"], bb["robust"], wc["aggressive"], wc["robust"], bwin, mmwin,
      if (fwd) "**FORWARD**" else if (bwin!=mmwin) "diverge" else "agree"))
}
cat("DONE\n")
