# Shared helpers for tuning (subset of divergence_study.R).
generate_data <- function(n, overlap_strength = 0.5, effect_size = -0.05,
                          seed = NULL, U_prevalence = 0, U_trt_OR = 1, U_out_OR = 1) {
  if (!is.null(seed)) set.seed(seed)
  age <- rnorm(n, 55, 10); sex <- rbinom(n, 1, 0.55)
  biomarker <- rnorm(n, 0, 1); comorbidity <- sample(0:2, n, TRUE, c(.5,.3,.2))
  ckd <- rbinom(n, 1, 0.12)
  U <- if (U_prevalence > 0) rbinom(n, 1, U_prevalence) else rep(0L, n)
  lp_trt <- -0.5 + overlap_strength * (0.03*(age-55) + 0.8*sex + 0.6*biomarker +
            0.5*ckd + 0.3*comorbidity) + log(U_trt_OR)*U
  treatment <- rbinom(n, 1, plogis(lp_trt))
  lp_out <- -2.5 + 0.015*(age-55) + 0.3*sex + 0.2*biomarker + 0.6*ckd +
            0.25*comorbidity + effect_size/0.15*treatment + log(U_out_OR)*U
  event_24 <- rbinom(n, 1, plogis(lp_out))
  nc_outcome <- rbinom(n, 1, plogis(-1 + 0.01*(age-55) + 0.1*sex + 0.15*biomarker))
  data.frame(age = round(age,1), sex = sex, biomarker = round(biomarker,3),
             comorbidity = comorbidity, ckd = ckd, treatment = treatment,
             event_24 = event_24, nc_outcome = nc_outcome)
}

near_positivity_stress <- function(lock, candidates, slopes, reps,
                                   effect_sizes = c(0.05)) {
  fit_one <- cleanTMLE:::.plasmode_fit_one_candidate
  data <- lock$data; n <- nrow(data); cov <- lock$covariates
  ids <- vapply(candidates, function(x) x$candidate_id, "")
  Q0 <- stats::glm(stats::reformulate(cov, response = "event_24"), data = data,
                   family = stats::binomial())
  p_base <- pmin(pmax(as.numeric(stats::predict(Q0, type = "response")), 0.001), 0.999)
  ps_mod <- stats::glm(stats::reformulate(cov, response = "treatment"), data = data,
                       family = stats::binomial())
  lp0 <- as.numeric(stats::predict(ps_mod, type = "link")); lp_bar <- mean(lp0); lpc <- lp0 - lp_bar
  rows <- list(); wdiag <- list()
  for (si in seq_along(slopes)) {
    s <- slopes[si]; ps_s <- plogis(lp_bar + s * lpc)
    for (es in effect_sizes) {
      rr <- vector("list", reps); maxwt <- numeric(reps)
      for (r in seq_len(reps)) {
        set.seed(lock$seed + r + 500000L + si * 1000L)
        A_rep <- stats::rbinom(n, 1L, ps_s)
        p1 <- pmin(pmax(p_base + es, 0.001), 0.999); p0 <- p_base
        Y_sim <- stats::rbinom(n, 1L, ifelse(A_rep == 1, p1, p0)); truth <- mean(p1) - mean(p0)
        ps_hat <- as.numeric(stats::predict(suppressWarnings(stats::glm(
          stats::reformulate(cov, response = ".A."),
          data = cbind(data, .A. = A_rep), family = stats::binomial())), type = "response"))
        maxwt[r] <- max(1/pmax(ps_hat,0.001), 1/pmax(1-ps_hat,0.001))
        res <- lapply(candidates, function(cd) tryCatch(
          fit_one(Y_sim = Y_sim, A = A_rep, W_data = data, treatment = "treatment",
                  covariates = cov, cand = cd, n = n),
          error = function(e) list(est=NA_real_, se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_)))
        names(res) <- ids; rr[[r]] <- c(res, list(.truth = truth))
      }
      tv <- vapply(rr, function(x) x$.truth, numeric(1))
      for (cid in ids) {
        est <- vapply(rr, function(x) x[[cid]]$est, numeric(1))
        se  <- vapply(rr, function(x) x[[cid]]$se,  numeric(1))
        lo  <- vapply(rr, function(x) x[[cid]]$ci_lower, numeric(1))
        hi  <- vapply(rr, function(x) x[[cid]]$ci_upper, numeric(1))
        ok <- !is.na(est); e <- est[ok]; t <- tv[ok]
        emp_sd <- stats::sd(e); mean_se <- mean(se[ok])
        rows[[length(rows)+1]] <- data.frame(
          scenario="near_positivity", level=sprintf("slope_x%.1f", s),
          effect_size=es, candidate=cid, bias=round(mean(e-t),5),
          rmse=round(sqrt(mean((e-t)^2)),5),
          coverage=round(mean(lo[ok] <= t & t <= hi[ok]),3),
          emp_sd=round(emp_sd,5), mean_se=round(mean_se,5),
          se_cal=round(if (emp_sd>0) mean_se/emp_sd else NA,3),
          n_converged=sum(ok), stringsAsFactors=FALSE)
      }
      wdiag[[length(wdiag)+1]] <- data.frame(level=sprintf("slope_x%.1f", s),
        effect_size=es, mean_max_weight=round(mean(maxwt),1), stringsAsFactors=FALSE)
    }
  }
  list(metrics = do.call(rbind, rows), weights = do.call(rbind, wdiag))
}
