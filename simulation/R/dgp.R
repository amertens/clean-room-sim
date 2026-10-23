# Data-generating processes for the estimability simulation (PLAN.md section 4).
#
# A scenario is one row of scenario_table(): an overlap level (gamma of the
# treatment model), a true outcome model (inside, near or outside the declared
# family of outcome surfaces) and an effect (non-null or null). Everything a
# scenario needs to generate data (gamma, b0, the standardizing constants of the
# inside model) is stored on its row, so a dataset depends on the row and its
# name only.

POP_N <- 1e6          # size of the population draw for calibration and truths
POP_SEED <- 20261007L # seed of the population draw
BETA_SEED <- 7L       # seed of the fixed direction X beta of the inside model
BASE_RISK <- 0.15     # control-arm risk that b0 is calibrated to
TRIM_BAND <- c(0.05, 0.95)
ESTIMANDS <- c("ATE", "trimmed_ATE", "ATT", "ATO")

# Every seeded step uses R's default generators, whatever the caller has set.
with_seed <- function(seed, code) {
  withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
                   .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

# Seed from a name: a hash of its parts, never a loop position. Kept below
# 2^28 so that seed + repetition and seed + 1000 * fold stay integers.
seed_from <- function(...) {
  h <- digest::digest(paste(..., sep = "|"), algo = "xxhash32", serialize = FALSE)
  strtoi(substr(h, 1, 7), base = 16L)
}

gen_covariates <- function(n, seed) {
  with_seed(seed, data.frame(
    X1 = stats::runif(n, -2, 2),
    X2 = stats::rbinom(n, 1L, 0.5),
    X3 = stats::runif(n, -1, 1),
    X4 = stats::runif(n, -1, 1)))
}

true_pi <- function(X, gamma) stats::plogis(-0.2 + gamma * X$X1 + 0.3 * X$X2 + 0.3 * X$X3)

# Population moments of X1 to X4 (exact), for the standardized direction X beta.
X_MEAN <- c(X1 = 0, X2 = 0.5, X3 = 0, X4 = 0)
X_SD <- c(X1 = 4 / sqrt(12), X2 = 0.5, X3 = 2 / sqrt(12), X4 = 2 / sqrt(12))
BETA <- with_seed(BETA_SEED, stats::rnorm(4))

xbeta <- function(X) {
  Z <- sweep(sweep(as.matrix(X[c("X1", "X2", "X3", "X4")]), 2, X_MEAN), 2, X_SD, "/")
  drop(Z %*% BETA)
}

# Linear predictor of the inside model's prognostic index and effect modifier.
# `pars` holds the population means and SDs used by z() (lp_mean, lp_sd,
# xb_mean, xb_sd) and c_inside, the scale that makes the prognostic index have
# SD 1, as the package's linear surface does with prognostic = 1.
inside_parts <- function(X, pi, pars) {
  zlp <- (stats::qlogis(pi) - pars$lp_mean) / pars$lp_sd
  zxb <- (xbeta(X) - pars$xb_mean) / pars$xb_sd
  list(prog = pars$c_inside * (0.5 * zlp + 0.5 * zxb), zlp = zlp)
}

# True outcome risk mu(a, X). `a` is 0, 1 or a vector; null = TRUE drops every
# term in A.
true_mu <- function(model, a, X, pi, b0, null = FALSE, pars = NULL) {
  eff <- if (null) 0 else 1
  eta <- switch(model,
    inside = {
      p <- inside_parts(X, pi, pars)
      b0 + p$prog + eff * a * (0.5 + 0.5 * p$zlp)
    },
    near = b0 + 0.8 * X$X1 + 0.6 * X$X1^2 + 0.4 * X$X3 + eff * a * (0.5 + 0.6 * X$X1),
    outside = b0 + 0.8 * X$X1 + 0.4 * X$X3 +
      eff * a * (0.5 + 1.0 * X$X4 + 0.8 * (X$X1 > 1.5)),
    stop("Unknown model: ", model))
  stats::plogis(eta)
}

pop_covariates <- function() gen_covariates(POP_N, POP_SEED)

# Standardizing constants of the inside model at one gamma, on the population draw.
inside_pars <- function(Xp, gamma) {
  lp <- stats::qlogis(true_pi(Xp, gamma))
  xb <- xbeta(Xp)
  pars <- list(lp_mean = mean(lp), lp_sd = stats::sd(lp), xb_mean = mean(xb),
               xb_sd = stats::sd(xb), c_inside = 1)
  zsum <- 0.5 * (lp - pars$lp_mean) / pars$lp_sd + 0.5 * (xb - pars$xb_mean) / pars$xb_sd
  pars$c_inside <- 1 / stats::sd(zsum)
  pars
}

calib_b0 <- function(model, Xp, pi, pars) {
  stats::uniroot(function(b) mean(true_mu(model, 0, Xp, pi, b, pars = pars)) - BASE_RISK,
                 c(-15, 10), tol = 1e-10)$root
}

# For each gamma: share of the population with pi outside the band, and b0 for
# the model (b0 depends on gamma only for the inside model, whose prognostic
# index contains logit pi).
calibrate <- function(model, gamma_grid, Xp = pop_covariates()) {
  do.call(rbind, lapply(gamma_grid, function(gm) {
    pi <- true_pi(Xp, gm)
    pars <- if (model == "inside") inside_pars(Xp, gm) else NULL
    data.frame(model = model, gamma = gm,
               share_outside = mean(pi < TRIM_BAND[1] | pi > TRIM_BAND[2]),
               b0 = calib_b0(model, Xp, pi, pars),
               c_inside = if (is.null(pars)) NA_real_ else pars$c_inside)
  }))
}

# Overlap levels (PLAN.md, revisions after pilot steps 1 and 2): gamma 1.25,
# 1.5, 1.75, 2.0 and 2.5. Each level is named by its population share of
# patients with true pi outside the band, in whole percent ("out04" is about
# 4% outside), and labelled with the share itself.
OVERLAP_GAMMAS <- c(1.25, 1.5, 1.75, 2.0, 2.5)

overlap_levels <- function(gammas = OVERLAP_GAMMAS, Xp = pop_covariates()) {
  share <- vapply(gammas, function(gm) {
    pi <- true_pi(Xp, gm)
    mean(pi < TRIM_BAND[1] | pi > TRIM_BAND[2])
  }, numeric(1))
  data.frame(id = sprintf("out%02d", round(100 * share)), gamma = gammas,
             share_outside = share,
             label = sprintf("%s%% outside", formatC(100 * share, format = "fg", digits = 2)),
             stringsAsFactors = FALSE)
}

# The scenario table: overlap level x model x effect, with everything needed
# to draw. `levels` is overlap_levels().
scenario_table <- function(levels, Xp = pop_covariates()) {
  g <- expand.grid(effect = c("effect", "null"), model = c("inside", "near", "outside"),
                   overlap = levels$id, stringsAsFactors = FALSE)
  g <- g[c("overlap", "model", "effect")]
  g$id <- paste(g$overlap, g$model, g$effect, sep = "_")
  g$gamma <- levels$gamma[match(g$overlap, levels$id)]
  rows <- lapply(seq_len(nrow(g)), function(i) {
    r <- g[i, ]
    pi <- true_pi(Xp, r$gamma)
    pars <- if (r$model == "inside") inside_pars(Xp, r$gamma) else
      list(lp_mean = NA_real_, lp_sd = NA_real_, xb_mean = NA_real_, xb_sd = NA_real_,
           c_inside = NA_real_)
    cbind(r, b0 = calib_b0(r$model, Xp, pi, if (r$model == "inside") pars else NULL),
          as.data.frame(pars))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

scen_pars <- function(scenario) {
  if (scenario$model != "inside") return(NULL)
  as.list(scenario[c("lp_mean", "lp_sd", "xb_mean", "xb_sd", "c_inside")])
}

scen_mu <- function(scenario, a, X, pi) {
  true_mu(scenario$model, a, X, pi, scenario$b0, null = scenario$effect == "null",
          pars = scen_pars(scenario))
}

# True estimands on the risk-difference scale, from the population draw.
population_truth <- function(scenario, Xp = pop_covariates()) {
  pi <- true_pi(Xp, scenario$gamma)
  tau <- scen_mu(scenario, 1, Xp, pi) - scen_mu(scenario, 0, Xp, pi)
  k <- pi >= TRIM_BAND[1] & pi <= TRIM_BAND[2]
  h <- pi * (1 - pi)
  c(ATE = mean(tau), trimmed_ATE = mean(tau[k]), ATT = sum(pi * tau) / sum(pi),
    ATO = sum(h * tau) / sum(h))
}

# Area under the ROC curve of score s for binary A (rank formula).
auc <- function(s, A) {
  r <- rank(s)
  n1 <- as.numeric(sum(A == 1))  # numeric: n1 * n0 overflows an integer at n = 1e6
  n0 <- as.numeric(sum(A == 0))
  (sum(r[A == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Table S2 diagnostics for one scenario on the population draw.
dgp_diagnostics <- function(scenario, Xp = pop_covariates()) {
  pi <- true_pi(Xp, scenario$gamma)
  A <- with_seed(POP_SEED + 1L, stats::rbinom(length(pi), 1L, pi))
  mu <- scen_mu(scenario, A, Xp, pi)
  prev <- mean(mu)
  tv <- population_truth(scenario, Xp)
  data.frame(id = scenario$id, overlap = scenario$overlap, model = scenario$model,
             effect = scenario$effect, gamma = scenario$gamma, b0 = scenario$b0,
             t(tv), share_outside = mean(pi < TRIM_BAND[1] | pi > TRIM_BAND[2]),
             c_stat_true_pi = auc(pi, A), share_treated = mean(A),
             prevalence = prev, control_risk = mean(scen_mu(scenario, 0, Xp, pi)),
             # Var(mu(A, X)) / Var(Y): the share of outcome variance explained.
             var_explained = stats::var(mu) / (prev * (1 - prev)),
             row.names = NULL)
}

# Datasets. Covariates and treatment depend on the overlap level only, seeded
# from ("data", overlap, n, idx): the map, which never sees the outcome, then
# depends on (overlap, n, idx) only and is shared by every outcome model and
# effect at that level. Each scenario draws its outcome on those covariates
# and treatment with its own seed ("outcome", scenario id, n, idx). Seeds come
# from names, never from loop positions.
design_seed <- function(overlap, n, idx) seed_from("data", overlap, n, idx)

# `level` is a row with `id` (the overlap id) and `gamma`; a scenario row works
# too, through its `overlap` column.
level_id <- function(level) if (!is.null(level$overlap)) level$overlap else level$id

draw_design <- function(level, n, idx) {
  s <- design_seed(level_id(level), n, idx)
  X <- gen_covariates(n, s)
  pi <- true_pi(X, level$gamma)
  A <- with_seed(s + 1L, stats::rbinom(n, 1L, pi))
  data.frame(id = seq_len(n), a = A, X, pi = pi)
}

draw_data <- function(scenario, n, idx) {
  d <- draw_design(scenario, n, idx)
  mu <- scen_mu(scenario, d$a, d, d$pi)
  d$y <- with_seed(seed_from("outcome", scenario$id, n, idx), stats::rbinom(n, 1L, mu))
  d
}
