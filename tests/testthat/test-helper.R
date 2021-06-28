test_that("get_crit works", {
  # Reproduced from Fujikawa et al. 2020
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.5, 0.5, 0.2))
  crit <- get_crit(design = design, n = 24, lambda = 0.99)
  expect_equal(crit, 10)
})

test_that("get_crit returns NA if sample size is too small", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.7,
    theta1 = c(0.9, 0.9, 0.7))
  crit <- get_crit(design = design, n = 11, lambda = 0.99)
  expect_equal(crit, NA_integer_)
})

test_that("get_crit_pool works", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.5, 0.5, 0.5))
  n <- 20
  crit <- get_crit_pool(design = design, n = n, lambda = 0.99)
  nocrit <- crit - 1
  weight_mat <- get_weights(design = design, n = n, epsilon = 2, tau = 0,
    logbase = 2)

  shape_crit <- matrix(c(design@shape1 + rep(crit, design@k),
    design@shape2 + n - rep(crit, design@k)), byrow = TRUE, ncol = design@k)
  shape_nocrit <- matrix(c(design@shape1 + rep(nocrit, design@k),
    design@shape2 + n - rep(nocrit, design@k)), byrow = TRUE, ncol = design@k)

  shape_crit <- beta_borrow(k = design@k, r = rep(crit, design@k),
    weight_mat = weight_mat, shape = shape_crit)
  shape_nocrit <- beta_borrow(k = design@k, r = rep(crit, design@k),
    weight_mat = weight_mat, shape = shape_nocrit)

  post_prob_crit <- post_beta(shape = shape_crit, theta0 = design@theta0)
  post_prob_nocrit <- post_beta(shape = shape_nocrit, theta0 = design@theta0)

  expect_true(all(post_prob_crit > 0.99))
  expect_true(all(post_prob_nocrit <= 0.99))
})

test_that("get_crit_pool returns NA if sample size is too small", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.8,
    theta1 = c(0.9, 0.9, 0.9))
  crit <- get_crit_pool(design = design, n = 10, lambda = 0.99)
  expect_true(is.na(crit))
})

test_that("get_targ works", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.5, 0.2, 0.2))
  targ_toer <- get_targ(design = design, prob = "toer")
  targ_pwr <- get_targ(design = design, prob = "pwr")

  expect_equal(targ_toer, c(FALSE, TRUE, TRUE))
  expect_equal(targ_pwr, c(TRUE, FALSE, FALSE))
})

test_that("prune_weights works", {
  design <- setupOneStageBasket(k = 6, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = rep(0.5, 6))
  weight_mat <- get_weights(design = design, n = 15, epsilon = 2, tau = 0,
    logbase = 2)
  weight_mat <- prune_weights(weight_mat, cut = 8)

  r <- c(5, 6, 7, 8, 9, 10)
  shape_post <- matrix(c(design@shape1 + r, design@shape2 + 15 - r),
    byrow = TRUE, ncol = design@k)
  shape_borrow <- beta_borrow(k = design@k, r = r, weight_mat = weight_mat,
    shape = shape_post)

  expect_equal(shape_post[, 1:3], shape_borrow[, 1:3])
  expect_false(any(shape_post[, 4:6] == shape_borrow[, 4:6]))
})
