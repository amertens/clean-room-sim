# ============================================================
# Negative-Control Reconciliation
# ============================================================
# Resolves the contradiction between:
#   reports/archive/manuscript_draft.qmd  — RD=0.033 / -0.019 SIGNIFICANT
#   rescueCo/results/negative_control_results.csv  — both NCs PASS
#
# Strategy: run all 5 NCs in 4 ways on the same fresh load of the data
#   (a) Unadjusted regression on FULL ambulance cohort (n=2,197)
#   (b) Unadjusted regression on PRIMARY cohort (n=1,703, transfers excluded)
#   (c) IPTW-adjusted regression on primary cohort
#   (d) TMLE on primary cohort (the cleanTMLE pipeline's method)
# Output: rescueCo/results/nc_reconciliation.csv
# ============================================================

suppressPackageStartupMessages({
  library(haven); library(data.table); library(janitor)
})

source(file.path(if (dir.exists("rescueCo")) "rescueCo/R" else "R",
                  "bootstrap.R"))
source("rescueCo/R/utils.R")
cfg <- load_cr_config()
cr_log("=== NC reconciliation: 5 NCs × 4 methods ===")

# Load registry directly (avoid relying on stage1_cohort.rds which excludes
# transfers and may have different NC encoding)
reg <- as.data.table(read_dta(cfg$data$trauma_registry_dta))
fup <- as.data.table(read_dta(cfg$data$sixmo_followup_dta))
names(reg) <- make_clean_names(names(reg))
names(fup) <- make_clean_names(names(fup))
as_num <- function(x) suppressWarnings(as.numeric(x))

# Treatment indicator (matches the modernised analysis/02 logic)
reg[, A := fifelse(as_num(ems_ambulance) == 1, 1L,
          fifelse(as_num(ems_ambulance) == 0, 0L, NA_integer_))]
reg[is.na(A) & as_num(hosp_transport) == 1, A := 0L]
reg[, is_transfer := as.character(haven::as_factor(inter_facility_transfer)) == "Yes*"]
reg[is.na(is_transfer), is_transfer := FALSE]

# 5 negative-control candidates from cfg
ncs <- c("chronic_hypertension", "chronic_diabetes_insulin",
         "chronic_hiv_art", "household_urban", "fuel_wood")

# Build NC indicators on registry
reg[, chronic_hypertension     := as.integer(as_num(chronic_illness_2))]
reg[, chronic_diabetes_insulin := as.integer(as_num(chronic_illness_0))]
reg[, chronic_hiv_art          := as.integer(as_num(chronic_illness_3))]
reg[, household_urban := as.integer(
  as.character(haven::as_factor(household_area)) == "Urban")]
reg[, fuel_wood := as.integer(as_num(cooking_fuel) == 0)]

# Cohorts
amb_full     <- reg[!is.na(A)]                       # all ambulance (n ≈ 2,197)
amb_primary  <- reg[!is.na(A) & is_transfer == FALSE] # primary cohort (n ≈ 1,703)

cr_log(paste("Full ambulance:", nrow(amb_full),
              " | Primary (no transfers):", nrow(amb_primary)))

# Per-cohort prevalence of each NC
prev_tab <- rbind(
  data.frame(cohort = "Full ambulance", n = nrow(amb_full),
             t(sapply(ncs, function(v) round(mean(amb_full[[v]], na.rm = TRUE), 4)))),
  data.frame(cohort = "Primary", n = nrow(amb_primary),
             t(sapply(ncs, function(v) round(mean(amb_primary[[v]], na.rm = TRUE), 4))))
)
print(prev_tab)
write.csv(prev_tab, "rescueCo/results/nc_prevalence_by_cohort.csv",
          row.names = FALSE)

# ============================================================
# Method (a): Unadjusted on full ambulance
# ============================================================
fit_unadj <- function(d, nc) {
  v <- d[[nc]]; A <- d$A
  ok <- !is.na(v) & !is.na(A)
  if (sum(ok) < 30 || length(unique(v[ok])) < 2) {
    return(data.frame(estimate = NA_real_, se = NA_real_,
                       p_value = NA_real_, n = sum(ok)))
  }
  m <- lm(v[ok] ~ A[ok])
  s <- summary(m)$coefficients["A[ok]", ]
  data.frame(estimate = s[1], se = s[2], p_value = s[4], n = sum(ok))
}

# ============================================================
# Method (c): IPTW-adjusted (PS reweighting)
# ============================================================
# Use a simple GLM PS on the primary cohort with the same covariate set
# as the locked analysis (read from analysis_outputs).
fit_iptw <- function(d, nc, ps_covs) {
  ok <- !is.na(d[[nc]]) & !is.na(d$A)
  d_ok <- d[ok]
  cov_avail <- intersect(ps_covs, names(d_ok))
  if (length(cov_avail) == 0) return(data.frame(estimate=NA, se=NA, p_value=NA, n=nrow(d_ok)))
  fmla <- as.formula(paste("A ~", paste(cov_avail, collapse = " + ")))
  ps <- tryCatch(
    glm(fmla, data = d_ok, family = binomial())$fitted.values,
    error = function(e) NULL)
  if (is.null(ps)) return(data.frame(estimate=NA, se=NA, p_value=NA, n=nrow(d_ok)))
  ps <- pmax(pmin(ps, 0.99), 0.01)
  w <- ifelse(d_ok$A == 1, 1/ps, 1/(1-ps))
  fit <- lm(d_ok[[nc]] ~ d_ok$A, weights = w)
  s <- summary(fit)$coefficients["d_ok$A", ]
  data.frame(estimate = s[1], se = s[2], p_value = s[4], n = nrow(d_ok))
}

# ============================================================
# Method (d): TMLE (matches cleanTMLE NC IPTW behaviour)
# ============================================================
fit_tmle <- function(d, nc, ps_covs) {
  cov_avail <- intersect(ps_covs, names(d))
  ok <- complete.cases(d[, c("A", nc, cov_avail), with = FALSE])
  d_ok <- d[ok]
  if (nrow(d_ok) < 50 || length(unique(d_ok[[nc]])) < 2)
    return(data.frame(estimate=NA, se=NA, p_value=NA, n=nrow(d_ok)))
  Y <- d_ok[[nc]]
  A <- d_ok$A
  W <- as.matrix(d_ok[, ..cov_avail])
  res <- tryCatch(tmle::tmle(Y = Y, A = A, W = W,
                              Q.SL.library = c("SL.glm", "SL.mean"),
                              g.SL.library = c("SL.glm", "SL.mean"),
                              family = if (length(unique(Y)) == 2) "binomial" else "gaussian"),
                  error = function(e) NULL)
  if (is.null(res)) return(data.frame(estimate=NA, se=NA, p_value=NA, n=nrow(d_ok)))
  ate <- res$estimates$ATE
  data.frame(estimate = ate$psi, se = sqrt(ate$var.psi),
             p_value = ate$pvalue, n = nrow(d_ok))
}

# Covariates available on registry directly (same as W_matrix construction
# in script 01 minus the dropped NZV columns)
ps_covs <- c("age", "sex", "marital_status", "education", "inj_distance",
             "occupation", "intent")
ps_covs <- intersect(ps_covs, names(reg))
# Convert haven_labelled to numeric for these
for (cv in ps_covs) {
  amb_full[[cv]]    <- as_num(amb_full[[cv]])
  amb_primary[[cv]] <- as_num(amb_primary[[cv]])
}

# Run all 5 NCs × 4 methods
out <- list()
for (nc in ncs) {
  cr_log(paste("Running NC:", nc))
  out[[length(out)+1]] <- cbind(method = "(a) Unadj, full ambulance",
    nc = nc, fit_unadj(amb_full, nc))
  out[[length(out)+1]] <- cbind(method = "(b) Unadj, primary cohort",
    nc = nc, fit_unadj(amb_primary, nc))
  out[[length(out)+1]] <- cbind(method = "(c) IPTW-adjusted, primary",
    nc = nc, fit_iptw(amb_primary, nc, ps_covs))
  out[[length(out)+1]] <- cbind(method = "(d) TMLE, primary",
    nc = nc, fit_tmle(amb_primary, nc, ps_covs))
}
recon <- do.call(rbind, out)
recon <- recon[order(recon$nc, recon$method), ]
recon$flagged <- !is.na(recon$p_value) & recon$p_value < 0.05
recon[, c("estimate","se","p_value")] <- lapply(recon[, c("estimate","se","p_value")], round, 5)

cr_log("Reconciliation table (5 NCs × 4 methods):")
print(recon, row.names = FALSE)
write.csv(recon, "rescueCo/results/nc_reconciliation.csv", row.names = FALSE)
cr_log("Saved nc_reconciliation.csv")
