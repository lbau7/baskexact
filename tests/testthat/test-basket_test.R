test_that("test works", {
  ## Without Pruning
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  res <- basket_test(design = design, n = 24, r = c(5, 6, 9), lambda = 0.99,
    epsilon = 2, tau = 0, logbase = exp(1), prune = FALSE)

  # Test if weights are correct
  weights_exp <- c(0.947424, 0.5192557, 0.6848343)
  weights <- res$weights[rbind(c(1, 2), c(3, 1), c(3, 2))]
  expect_equal(weights, weights_exp, tolerance = 10e-7)

  # Test if posterior distributions are correct
  shape_exp <- c(17.82453, 46.30915, 19.53289, 48.90583,
    17.90937, 39.39697)
  shape <- as.vector(res$post_dist_borrow)
  expect_equal(shape, shape_exp, tolerance = 10e-7)

  # Test ist posterior probabilities are correct
  prob_exp <- c(0.9262286, 0.9497929, 0.9753582)
  prob <- as.vector(res$post_prob_borrow)
  expect_equal(prob, prob_exp, tolerance = 10e-7)

  ## With Pruning
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  res <- basket_test(design = design, n = 24, r = c(4, 4, 5), lambda = 0.99,
    epsilon = 2, tau = 0, logbase = exp(1), prune = TRUE)

  # Results are Equal when all Baskets are Pruned
  expect_equal(res$post_dist_noborrow, res$post_dist_borrow)
  expect_equal(res$post_prob_noborrow, res$post_prob_borrow)
})

test_that("errors in basket_test work", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)

  expect_error(basket_test(design = design, n = 20, r = c(-1, 10, 10),
    lambda = 0.99, epsilon = 2, tau = 0, logbase = 2, prune = FALSE))
  expect_error(basket_test(design = design, n = 20, r = c(1, 25, 10),
    lambda = 0.99, epsilon = 2, tau = 0, logbase = 2, prune = FALSE))
})
