test_that("basket_test works", {
  ## Without Pruning
  # Reproduced from Fujikawa et al., 2020, Supplement R Code
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  res1 <- basket_test(design = design, n = 24, r = c(5, 3, 8), lambda = 0.99,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
    logbase = exp(1), prune = FALSE))

  # Test if weights are correct
  weights_exp <- c(0.7750081, 0.6661536, 0.3430768)
  weights <- res1$weights[rbind(c(1, 2), c(3, 1), c(3, 2))]
  expect_equal(weights, weights_exp, tolerance = 10e-7)

  # Test if posterior distributions are correct
  shape_exp <- c(15.09542, 48.37479, 11.73774, 43.33247,
    14.36923, 37.87076)
  shape <- as.vector(res1$post_dist_borrow)
  expect_equal(shape, shape_exp, tolerance = 10e-7)

  # Test if posterior probabilities are correct
  prob_exp <- c(0.7524369, 0.5703997, 0.8941389)
  prob <- as.vector(res1$post_prob_borrow)
  expect_equal(prob, prob_exp, tolerance = 10e-7)

  ## With Pruning
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  res2 <- basket_test(design = design, n = 24, r = c(4, 4, 5), lambda = 0.99,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
    logbase = exp(1), prune = TRUE))

  # Results are equal when all baskets are pruned
  expect_equal(res2$post_dist_noborrow, res2$post_dist_borrow)
  expect_equal(res2$post_prob_noborrow, res2$post_prob_borrow)

  # With Global Weight
  res3 <- basket_test(design = design, n = 24, r = c(5, 3, 8), lambda = 0.99,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
      logbase = exp(1), prune = FALSE), globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 2, w = 0.5))

  # Weights are strictly smaller when global weight is used
  expect_true(all(res3$weights[rbind(c(1, 2), c(3, 1), c(3, 2))] <
    res1$weights[rbind(c(1, 2), c(3, 1), c(3, 2))]))
})

test_that("errors in basket_test work", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)

  expect_error(basket_test(design = design, n = 20, r = c(-1, 10, 10),
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE)))
  expect_error(basket_test(design = design, n = 20, r = c(1, 25, 10),
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE)))
  expect_error(basket_test(design = design, n = 20, r = c(2, 3),
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
      epsilon = 2, tau = 0, logbase = 2, prune = FALSE)))
})
