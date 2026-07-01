# attach_estimand / declare_sensitivity_plan operate on a cleanTMLE lock, so
# these tests need cleanTMLE to construct one.

test_that("attach_estimand stashes estimand metadata without touching the hash", {
  skip_if_not_installed("cleanTMLE")
  dat  <- cleanTMLE::sim_func1(n = 120, seed = 1)
  lock <- cleanTMLE::create_analysis_lock(dat, "treatment", "event_24",
                                          c("age", "sex", "biomarker"), seed = 1)
  h0 <- lock$lock_hash
  lock2 <- attach_estimand(lock,
    description = "Effect on 24-month event risk",
    treatment_strategies = c("Treatment", "Control"),
    contrast = "risk_difference")
  expect_equal(lock2$estimand$contrast, "risk_difference")
  expect_equal(lock2$estimand$treatment_strategies, c("Treatment", "Control"))
  expect_identical(lock2$lock_hash, h0)  # metadata does not change the hash
})

test_that("attach_estimand rejects a non-lock", {
  expect_error(attach_estimand(list()), "cleanroom_lock")
})

test_that("declare_sensitivity_plan appends named plans", {
  skip_if_not_installed("cleanTMLE")
  dat  <- cleanTMLE::sim_func1(n = 120, seed = 2)
  lock <- cleanTMLE::create_analysis_lock(dat, "treatment", "event_24",
                                          c("age", "sex", "biomarker"), seed = 2)
  lock <- declare_sensitivity_plan(lock, label = "trunc",
    description = "Vary PS truncation",
    settings = list(truncation = c(0.01, 0.05, 0.10)))
  expect_true("trunc" %in% names(lock$sensitivity_plans))
  expect_equal(lock$sensitivity_plans$trunc$settings$truncation,
               c(0.01, 0.05, 0.10))
})

test_that("declare_sensitivity_plan requires a single-string label", {
  skip_if_not_installed("cleanTMLE")
  dat  <- cleanTMLE::sim_func1(n = 80, seed = 3)
  lock <- cleanTMLE::create_analysis_lock(dat, "treatment", "event_24",
                                          c("age", "sex"), seed = 3)
  expect_error(declare_sensitivity_plan(lock, label = c("a", "b")),
               "single character string")
})
