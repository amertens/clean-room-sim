#!/usr/bin/env Rscript
# ===========================================================================
# 99_smoke_test.R
# Quick validation: runs a small N simulation through all stages and
# asserts that key outputs are non-NA with reasonable uncertainty.
#
# Usage: Rscript scripts/99_smoke_test.R
# ===========================================================================

message("=== Smoke Test ===")

# Set working directory to project root
if (file.exists("scripts/99_smoke_test.R")) {
  # Already in project root
} else if (file.exists("../scripts/99_smoke_test.R")) {
  setwd("..")
}

# Source all R modules
r_files <- list.files("R", pattern = "\\.R$", recursive = TRUE,
                      full.names = TRUE)
for (f in r_files) source(f)

# Use a temp directory for outputs
tmp_dir <- tempdir()
s1_dir <- file.path(tmp_dir, "stage1")
s2_dir <- file.path(tmp_dir, "stage2")
s3_dir <- file.path(tmp_dir, "stage3")
for (d in c(s1_dir, s2_dir, s3_dir)) dir.create(d, recursive = TRUE,
                                                   showWarnings = FALSE)

cfg <- load_config("config/default.yml")
# Override for speed
cfg$stage1$min_events <- 10

# Small N simulation
set_deterministic_seed(99999)
dat <- generate_hcv_data(
  N = 2000, h0 = 8e-4, HR_early = 1.30, HR_late = 0.75,
  tau = 45, max_follow = 365,
  np_hazard = FALSE, dep_censor = FALSE, complexity = FALSE,
  switch_on = FALSE, seed = 99999
)

n_events <- sum(dat$event)
message("Events: ", n_events, " / ", nrow(dat),
        " (", round(mean(dat$event) * 100, 1), "%)")
stopifnot(n_events >= 10)

# Stage 1
s1 <- stage1_build_cohort(dat, cfg = cfg, output_dir = s1_dir)
stopifnot(s1$checkpoint == "PASS")
message("Stage 1: PASS")

# Stage 2
s2 <- stage2_design_checks(s1$cohort, cfg = cfg, output_dir = s2_dir)
stopifnot(s2$checkpoint == "PASS")
message("Stage 2: PASS")

# Stage 3 (small bootstrap)
s3 <- stage3_estimation(s1$cohort, stage2_result = s2, cfg = cfg,
                        output_dir = s3_dir, B = 50)

# Validate TMLE at 180 days
tmle_180 <- s3$results$t180$tmle
stopifnot(!is.na(tmle_180$RD))
stopifnot(!is.na(tmle_180$se_RD))
stopifnot(tmle_180$se_RD > 0)
message("TMLE 180d RD: ", round(tmle_180$RD, 5),
        " (SE: ", round(tmle_180$se_RD, 5), ")")

# Validate G-comp SE is not zero
gcomp_180 <- s3$results$t180$gcomp
stopifnot(!is.na(gcomp_180$RD))
stopifnot(!is.na(gcomp_180$se_RD))
stopifnot(gcomp_180$se_RD > 0)
message("G-comp 180d RD: ", round(gcomp_180$RD, 5),
        " (SE: ", round(gcomp_180$se_RD, 5), ")")

# Validate CIs are non-degenerate
stopifnot(diff(tmle_180$ci_RD) > 0)
stopifnot(diff(gcomp_180$ci_RD) > 0)

message("\n=== Smoke Test: OK ===")
