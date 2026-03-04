test_that("DGP is reproducible with same seed", {
  d1 <- generate_hcv_data(N = 1000, seed = 42, complexity = FALSE,
                          np_hazard = FALSE, dep_censor = FALSE,
                          switch_on = FALSE)
  d2 <- generate_hcv_data(N = 1000, seed = 42, complexity = FALSE,
                          np_hazard = FALSE, dep_censor = FALSE,
                          switch_on = FALSE)
  expect_equal(nrow(d1), nrow(d2))
  expect_equal(mean(d1$age), mean(d2$age))
  expect_equal(mean(d1$event), mean(d2$event))
  expect_equal(mean(d1$follow_time), mean(d2$follow_time))
  expect_equal(sum(d1$treatment), sum(d2$treatment))
})

test_that("DGP returns expected columns", {
  d <- generate_hcv_data(N = 500, seed = 1, complexity = FALSE,
                         np_hazard = FALSE, dep_censor = FALSE,
                         switch_on = FALSE)
  required_cols <- c("id", "age", "sex_male", "treatment", "event",
                     "follow_time", "switch", "ckd", "cirrhosis")
  expect_true(all(required_cols %in% names(d)))
})

test_that("DGP produces valid data ranges", {
  d <- generate_hcv_data(N = 2000, seed = 7, complexity = FALSE,
                         np_hazard = FALSE, dep_censor = FALSE,
                         switch_on = FALSE)
  expect_true(all(d$age >= 18))
  expect_true(all(d$treatment %in% c(0, 1)))
  expect_true(all(d$event %in% c(0, 1)))
  expect_true(all(d$follow_time > 0))
})

test_that("DGP counterfactual overrides work", {
  d1 <- generate_hcv_data(N = 500, seed = 10, treat_override = "all_treated",
                          complexity = FALSE, np_hazard = FALSE,
                          dep_censor = FALSE, switch_on = FALSE)
  expect_true(all(d1$treatment == 1))

  d0 <- generate_hcv_data(N = 500, seed = 10, treat_override = "all_control",
                          complexity = FALSE, np_hazard = FALSE,
                          dep_censor = FALSE, switch_on = FALSE)
  expect_true(all(d0$treatment == 0))
})

test_that("compute_true_risk returns finite values", {
  tr <- compute_true_risk(N_truth = 5000, t_eval = 180, seed = 42,
                          h0 = 3e-4, complexity = FALSE, np_hazard = FALSE,
                          dep_censor = FALSE, switch_on = FALSE)
  expect_true(is.finite(tr$risk_1))
  expect_true(is.finite(tr$risk_0))
  expect_true(is.finite(tr$RD))
  expect_true(tr$risk_1 >= 0 && tr$risk_1 <= 1)
  expect_true(tr$risk_0 >= 0 && tr$risk_0 <= 1)
})
