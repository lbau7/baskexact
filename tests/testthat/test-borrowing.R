test_that("get_weight works", {
  # Reproduced from Fujikawa et al. 2019, Supplement R code
  design <- setupBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.5, 0.5, 0.2))
  weight_mat <- get_weights(design = design, n = 15, epsilon = 2, tau = 0,
    logbase = exp(1))
  r <- c(5, 1, 3)
  elmnts <- all_combs <- t(utils::combn(r, 2)) + 1
  weights <- as.vector(weight_mat[elmnts])
  weights_exp <- c(0.3206983, 0.7493639, 0.6509846)
  expect_equal(weights, weights_exp, tolerance = 1e-7)
})
