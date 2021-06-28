test_that("toer works", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.5, 0.5, 0.5))
  toer_group <- suppressMessages(toer(design = design, n = 24, epsilon = 2,
    tau = 0.5, logbase = exp(1), lambda = 0.99, results = "group"))
  toer_fwer <- suppressMessages(toer(design = design, n = 24, epsilon = 2,
    tau = 0.5, logbase = exp(1), lambda = 0.99, results = "fwer"))
  rej_expect <- c(0.03239555, 0.03239555, 0.03239555)

  expect_equal(toer_group$rejection_probabilities, rej_expect,
    tolerance = 10e-8)
  expect_equal(toer_fwer, 0.06315308, tolerance = 10e-8)
  expect_equal(toer_group$fwer, toer_fwer)
})

test_that("toer stops when it's supposed to ", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5))
  expect_error(toer(design = design, n = c(10, 15, 20), lambda = 0.99,
    epsilon = 2, tau = 0, logbase = 2, results = "fwer"))
  expect_error(toer(design = design, n = 20, lambda = 1.1,
    epsilon = 2, tau = 0, logbase = 2, results = "fwer"))
  expect_error(toer(design = design, n = 20, lambda = 0.99,
    epsilon = -2, tau = 0, logbase = 2, results = "fwer"))
  expect_error(toer(design = design, n = 20, lambda = 0.99,
    epsilon = 2, tau = 1, logbase = 2, results = "fwer"))
  expect_error(toer(design = design, n = 20, lambda = 0.99,
    epsilon = 2, tau = 0, logbase = -2, results = "fwer"))
})
