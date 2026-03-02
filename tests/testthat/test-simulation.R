test_that("Simulation driver produces expected output columns", {
  skip_if_not_installed("SuperLearner")

  # Run a minimal simulation: 1 scenario, 2 replicates, 1 time point
  d <- generate_hcv_data(N = 500, seed = 99, h0 = 3e-4,
                         complexity = FALSE, np_hazard = FALSE,
                         dep_censor = FALSE, switch_on = FALSE)

  # Manually test one replicate of each estimator
  tmle_res <- tmle_survival_risk(
    data = d, t_eval = 90,
    sl_lib_Q = c("SL.glm", "SL.mean"),
    sl_lib_g = c("SL.glm", "SL.mean"),
    sl_lib_cens = c("SL.glm", "SL.mean")
  )

  iptw_res <- iptw_survival(data = d, t_eval = 90,
                            sl_lib = c("SL.glm", "SL.mean"))

  gcomp_res <- gcomp_risk(data = d, t_eval = 90,
                          sl_lib = c("SL.glm", "SL.mean"))

  cox_res <- cox_ph_estimator(d)

  # Check all return finite values

  expect_true(is.finite(tmle_res$RD))
  expect_true(is.finite(iptw_res$RD))
  expect_true(is.finite(gcomp_res$RD))
  expect_true(is.finite(cox_res$HR))
})

test_that("Config loads and has all scenario names", {
  cfg <- load_config()
  expected_scenarios <- c("simple", "nonlinear", "dep_censor", "switching",
                          "np_hazard")
  expect_true(all(expected_scenarios %in% names(cfg$scenarios)))
})

test_that("Comparator estimators return expected structure", {
  d <- generate_hcv_data(N = 500, seed = 88, h0 = 3e-4,
                         complexity = FALSE, np_hazard = FALSE,
                         dep_censor = FALSE, switch_on = FALSE)

  # Cox
  cox <- cox_ph_estimator(d)
  expect_true("HR" %in% names(cox))
  expect_true("ph_test" %in% names(cox))
  expect_true(is.finite(cox$HR))

  # Cox risk at t
  cox_r <- cox_risk_at_t(cox, d, t_eval = 90)
  expect_true(is.finite(cox_r$RD))

  # IPTW
  iptw <- iptw_survival(d, t_eval = 90, sl_lib = c("SL.glm", "SL.mean"))
  expect_true("RD" %in% names(iptw))
  expect_true(is.finite(iptw$RD))

  # G-comp
  gc <- gcomp_risk(d, t_eval = 90, sl_lib = c("SL.glm", "SL.mean"))
  expect_true("RD" %in% names(gc))
  expect_true(is.finite(gc$RD))
})
