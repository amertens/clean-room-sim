test_that("tmle_point_binary returns finite estimates on small data", {
  skip_if_not_installed("tmle")
  skip_if_not_installed("SuperLearner")

  set.seed(123)
  n <- 500
  W <- data.frame(W1 = rnorm(n), W2 = rbinom(n, 1, 0.5))
  A <- rbinom(n, 1, plogis(0.5 * W$W1 + 0.3 * W$W2))
  Y <- rbinom(n, 1, plogis(-1 + 0.5 * A + 0.3 * W$W1))

  res <- tmle_point_binary(Y, A, W,
                           sl_lib_Q = c("SL.glm", "SL.mean"),
                           sl_lib_g = c("SL.glm", "SL.mean"))

  expect_true(is.finite(res$psi))
  expect_true(is.finite(res$se))
  expect_true(all(is.finite(res$ci)))
  expect_true(res$se > 0)
  expect_true(res$ci[1] < res$ci[2])
})

test_that("tmle_survival_risk returns finite estimates", {
  skip_if_not_installed("SuperLearner")

  d <- generate_hcv_data(N = 1000, seed = 77, h0 = 3e-4,
                         complexity = FALSE, np_hazard = FALSE,
                         dep_censor = FALSE, switch_on = FALSE)

  res <- tmle_survival_risk(
    data = d, t_eval = 180,
    sl_lib_Q = c("SL.glm", "SL.mean"),
    sl_lib_g = c("SL.glm", "SL.mean"),
    sl_lib_cens = c("SL.glm", "SL.mean")
  )

  expect_true(is.finite(res$RD))
  expect_true(is.finite(res$se_RD))
  expect_true(res$se_RD > 0)
  expect_true(all(is.finite(res$ci_RD)))
  expect_true(res$risk_1 >= 0 && res$risk_1 <= 1)
  expect_true(res$risk_0 >= 0 && res$risk_0 <= 1)
})

test_that("filter_sl_libraries always returns at least two libraries", {
  result <- filter_sl_libraries(c("SL.nonexistent_package_xyz"))
  expect_true(length(result) >= 1)
  expect_true("SL.glm" %in% result || "SL.mean" %in% result)
})
