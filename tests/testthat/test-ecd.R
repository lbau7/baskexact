test_that("ecd works", {
  design <- setupOneStageBasket(k = 3, theta0 = 0.2)
  # Calculate expected number of correction decisions for a design with equal
  # theta1
  ecd1 <- ecd(design, theta1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)
  ecd_loop1 <- ecd_loop(design, theta1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)

  # Calculate expected number of correction decisions for a design with unequal
  # theta1
  ecd2 <- ecd(design, theta1 = c(0.2, 0.4, 0.6), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)
  ecd_loop2 <- ecd_loop(design, theta1 = c(0.2, 0.4, 0.6), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)

  expect_equal(ecd1, ecd_loop1)
  expect_equal(ecd2, ecd_loop2)
})

