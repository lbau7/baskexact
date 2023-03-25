test_that("ecd works", {
  design <- setupOneStageBasket(k = 3, theta0 = 0.2)
  # Calculate expected number of correction decisions for a design with equal
  # theta1
  ecd1 <- ecd(design, theta1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)
  ecd_loop1 <- ecd_loop(design, theta1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)

  expect_equal(ecd1, ecd_loop1)

  # Calculate expected number of correction decisions for a design with unequal
  # theta1
  ecd2 <- ecd(design, theta1 = c(0.2, 0.4, 0.6), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)
  ecd_loop2 <- ecd_loop(design, theta1 = c(0.2, 0.4, 0.6), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)

  expect_equal(ecd2, ecd_loop2)

  # Calculate expected number of correct decisions with global weights with
  # equal theta1
  ecd3 <- ecd(design, theta1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))
  ecd_loop3 <- ecd_loop(design, theta1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))

  expect_equal(ecd3, ecd_loop3)

  # Calculate expected number of correct decisions with global weights with
  # unequal theta1
  ecd4 <- ecd(design, theta1 = c(0.2, 0.3, 0.4), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))
  ecd_loop4 <- ecd_loop(design, theta1 = c(0.2, 0.3, 0.4), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))

  expect_equal(ecd4, ecd_loop4)
})

