test_that("estim works", {
  design <- setupOneStageBasket(k = 3, theta0 = 0.2)
  res1 <- estim(design = design, theta1 = 0.5, n = 15,
    weight_fun = weights_fujikawa)
  res2 <- estim(design = design, theta1 = c(0.2, 0.5, 0.5), n = 15,
    weight_fun = weights_fujikawa)

  # No bias when theta1 is equal for all baskets and prior parameters
  # correspond to theta1
  expect_equal(res1$Mean, rep(0.5, 3))
  expect_true(all(res2$MSE > res1$MSE))
})
