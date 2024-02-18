test_that("estim works for single-stage designs", {
  design <- setupOneStageBasket(k = 3, p0 = 0.2)

  # No bias when p1 is equal for all baskets and prior parameters
  # correspond to p1
  res1 <- estim(design = design, p1 = 0.5, n = 15,
    weight_fun = weights_fujikawa)
  res_loop1 <- estim(design = design, p1 = c(0.2, 0.5, 0.5), n = 15,
    weight_fun = weights_fujikawa)

  expect_equal(res1$Mean, rep(0.5, 3))
  expect_true(all(res_loop1$MSE > res1$MSE))

  # Calculate posterior means for a design with equal p1
  res2 <- estim(design = design, p1 = c(0.2, 0.2, 0.2), n = 15,
    weight_fun = weights_mml)
  res_loop2 <- estim_loop(design = design, p1 = c(0.2, 0.2, 0.2), n = 15,
    weight_fun = weights_mml, weight_params = list())

  expect_equal(res2, res_loop2)

  # Calculate posterior means for a design with unequal p1
  res3 <- estim(design = design, p1 = c(0.2, 0.4, 0.5), n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0))
  res_loop3 <- estim_loop(design = design, p1 = c(0.2, 0.4, 0.5), n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0))

  expect_equal(res3, res_loop3)
})

test_that("estim works for two-stage designs", {
  design <- setupTwoStageBasket(k = 3, p0 = 0.2)

  # Calculate posterior means for a design with equal p1
  res1 <- estim(design = design, p1 = c(0.2, 0.2, 0.2), n = 16, n1 = 8,
    lambda = 0.95, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_cpp, weight_params = list(a = 2, b = 2))

  res_loop1 <- estim_twostage_loop(design = design, p1 = c(0.2, 0.2, 0.2),
    n = 16, n1 = 8, lambda = 0.95, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_cpp, weight_params = list(a = 2, b = 2))

  expect_equal(res1, res_loop1)

  # Calculate posterior means for a design with unequal p1
  res2 <- estim(design = design, p1 = c(0.4, 0.5, 0.2), n = 16, n1 = 8,
    lambda = 0.95, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.2, prob_effstop = 0.8),
    weight_fun = weights_cpp, weight_params = list(a = 1, b = 1))

  res_loop2 <- estim_twostage_loop(design = design, p1 = c(0.4, 0.5, 0.2),
    n = 16, n1 = 8, lambda = 0.95, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.2, prob_effstop = 0.8),
    weight_fun = weights_cpp, weight_params = list(a = 1, b = 1))

  expect_equal(res2, res_loop2)
})
