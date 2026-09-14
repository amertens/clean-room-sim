# build_stage_manifest / summarize_stage_path read a cleanTMLE audit, so these
# tests need cleanTMLE to construct one.

make_audit <- function(seed = 1) {
  dat   <- cleanTMLE::sim_func1(n = 120, seed = seed)
  lock  <- cleanTMLE::create_analysis_lock(dat, "treatment", "event_24",
                                           c("age", "sex", "biomarker"),
                                           seed = seed)
  audit <- cleanTMLE::create_audit_log(lock)
  audit <- cleanTMLE::record_stage(audit, "Stage 1a", "Lock created",
                                   decision = NA)
  audit <- cleanTMLE::record_stage(audit, "Check Point 1", "CP1 evaluated",
                                   decision = "GO")
  audit
}

test_that("build_stage_manifest prints the recorded stage path", {
  skip_if_not_installed("cleanTMLE")
  audit <- make_audit(1)
  out <- capture.output(res <- build_stage_manifest(audit))
  expect_true(any(grepl("Clean-Room Stage Path", out)))
  expect_true(any(grepl("Check Point 1", out)))
  expect_type(res, "character")
})

test_that("build_stage_manifest handles an empty audit", {
  skip_if_not_installed("cleanTMLE")
  dat   <- cleanTMLE::sim_func1(n = 60, seed = 2)
  lock  <- cleanTMLE::create_analysis_lock(dat, "treatment", "event_24",
                                           c("age", "sex"), seed = 2)
  audit <- cleanTMLE::create_audit_log(lock)
  out <- capture.output(build_stage_manifest(audit))
  expect_true(any(grepl("No stage entries", out)))
})

test_that("summarize_stage_path renders a compact narrative", {
  skip_if_not_installed("cleanTMLE")
  audit <- make_audit(3)
  out <- capture.output(res <- summarize_stage_path(audit))
  expect_true(any(grepl("Stage 1a", out)))
  expect_true(any(grepl("Check Point 1: GO", out)))
  expect_true(grepl("->", res))
})

test_that("summarize_stage_path rejects a non-audit", {
  expect_error(summarize_stage_path(list()), "cleantmle_audit")
})
