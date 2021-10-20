test_that("pow works without pruning", {
  # Compare groupwise power with experimentwise power with one active basket
  # and compare with reject_single_loop
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  pow_ewp <- pow(design = design1, theta1 = c(0.5, 0.2, 0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE), results = "ewp")
  pow_group <- pow(design = design1, theta1 = c(0.5, 0.2, 0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE), results = "group")
  pow_loop <- reject_single_loop(design = design1, theta = c(0.5, 0.2, 0.2),
    n = 15, lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE), prob = "pwr")

  expect_equal(pow_ewp, pow_group$ewp)
  expect_equal(pow_loop$ewp, pow_group$ewp)
  expect_equal(pow_loop$rejection_probabilities,
    pow_group$rejection_probabilities)

  # Compare rejection probabilities of pow and toer
  design2 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.3)
  toer_probs <- toer(design = design2, theta1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = FALSE), results = "group")
  pow_probs <- pow(design = design2, theta1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = FALSE), results = "group")

  expect_equal(toer_probs$rejection_probabilities,
    pow_probs$rejection_probabilities)
  expect_false(toer_probs$fwer == pow_probs$ewp)
})

test_that("pow works with pruning", {
  # Compare groupwise power with experimentwise power with one active basket
  # and compare with reject_single_loop
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  pow_ewp <- pow(design = design1, theta1 = c(0.5, 0.2,0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE), results = "ewp")
  pow_group <- pow(design = design1, theta1 = c(0.5, 0.2, 0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE), results = "group")
  pow_loop <- reject_single_loop(design = design1, theta1 = c(0.5, 0.2, 0.2),
    n = 15, lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE), prob = "pwr")

  expect_equal(pow_ewp, pow_group$ewp)
  expect_equal(pow_loop$ewp, pow_group$ewp)
  expect_equal(pow_loop$rejection_probabilities,
    pow_group$rejection_probabilities)

  # Compare rejection probabilities of pow and toer
  design2 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.3)
  toer_probs <- toer(design = design2, theta1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = TRUE), results = "group")
  pow_probs <- pow(design = design2, theta1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = TRUE), results = "group")

  expect_equal(toer_probs$rejection_probabilities,
    pow_probs$rejection_probabilities)
  expect_false(toer_probs$fwer == pow_probs$ewp)
})
