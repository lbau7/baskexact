test_that("get_weight works", {
  # Reproduced from Fujikawa et al., 2020, Supplement R code
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.5, 0.5, 0.2))
  weight_mat <- get_weights(design = design, n = 15, epsilon = 2, tau = 0,
    logbase = exp(1))
  r <- c(5, 1, 3)
  elmnts <- all_combs <- t(utils::combn(r, 2)) + 1
  weights <- as.vector(weight_mat[elmnts])
  weights_exp <- c(0.3206983, 0.7493639, 0.6509846)

  expect_equal(weights, weights_exp, tolerance = 10e-7)
})

test_that("beta_borrow works", {
  # Reproduced from Fujikawa et al., 2020, Supplement R Code
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.2, 0.2, 0.2))
  weight_mat <- get_weights(design = design, n = 24, epsilon = 2, tau = 0.5,
    logbase = exp(1))
  r <- c(7, 2, 5)
  shape_post <- matrix(c(design@shape1 + r, design@shape2 + 24 - r),
    byrow = TRUE, ncol = design@k)
  shape_borrow <- beta_borrow(k = 3, r = r, weight_mat = weight_mat,
    shape = shape_post)
  shape_expect <- matrix(c(12.9215409, 34.4051363, 6.33262523, 34.1087508,
    14.2283671, 47.5396861), nrow = 2)

  expect_equal(shape_borrow, shape_expect, tolerance = 10e-7)
})
