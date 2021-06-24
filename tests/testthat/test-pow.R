test_that("pow calculates gwp equal to ewp for one active basket", {
  design <- setupBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.5, 0.2, 0.2))
  pow_ewp <- pow(design = design, n = 15, lambda = 0.99, epsilon = 2,
    tau = 0, results = "ewp")
  pow_group <- pow(design = design, n = 15, lambda = 0.99, epsilon = 2,
    tau = 0, results = "group")$ewp

  expect_equal(pow_ewp, pow_group)
})

test_that("rejection probabilities with pow are equal to toer", {
  design <- setupBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.3,
    theta1 = c(0.6, 0.6, 0.3))
  toer_group <- toer(design = design, n = 20, lambda = 0.95, epsilon = 1,
    tau = 0.2, results = "group")
  pow_group <- pow(design = design, n = 20, lambda = 0.95, epsilon = 1,
    tau = 0.2, results = "group")

  expect_equal(toer_group$rejection_probabilities,
    pow_group$rejection_probabilities)
  expect_false(toer_group$fwer == pow_group$ewp)
})

test_that("pow stops when it's supposed to ", {
  design <- setupBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5))
  expect_error(pow(design = design, n = c(10, 15, 20), lambda = 0.99,
    epsilon = 2, tau = 0, logbase = 2, results = "fwer"))
  expect_error(pow(design = design, n = 20, lambda = 1.1,
    epsilon = 2, tau = 0, logbase = 2, results = "fwer"))
  expect_error(pow(design = design, n = 20, lambda = 0.99,
    epsilon = -2, tau = 0, logbase = 2, results = "fwer"))
  expect_error(pow(design = design, n = 20, lambda = 0.99,
    epsilon = 2, tau = 1, logbase = 2, results = "fwer"))
  expect_error(pow(design = design, n = 20, lambda = 0.99,
    epsilon = 2, tau = 0, logbase = -2, results = "fwer"))
})
