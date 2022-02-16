test_that("get_crit works", {
  # Reproduced from Fujikawa et al. 2020
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  crit <- get_crit(design = design, n = 24, lambda = 0.99)

  expect_equal(crit, 10)
})

test_that("get_crit returns NA if sample size is too small", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.7)
  crit <- get_crit(design = design, n = 11, lambda = 0.99)

  expect_equal(crit, NA_integer_)
})

test_that("get_crit_pool works", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  n <- 20
  crit <- get_crit_pool(design = design, n = n, lambda = 0.99)
  nocrit <- crit - 1
  weight_mat <- weights_fujikawa(design = design, n = n, epsilon = 2, tau = 0,
    logbase = 2, prune = FALSE)

  shape_crit <- matrix(c(design@shape1 + rep(crit, design@k),
    design@shape2 + n - rep(crit, design@k)), byrow = TRUE, ncol = design@k)
  shape_nocrit <- matrix(c(design@shape1 + rep(nocrit, design@k),
    design@shape2 + n - rep(nocrit, design@k)), byrow = TRUE, ncol = design@k)

  shape_crit <- beta_borrow(weight_mat = weight_mat, design = design, n = n,
    r = rep(crit, design@k))
  shape_nocrit <- beta_borrow(weight_mat = weight_mat, design = design, n = n,
    r = rep(nocrit, design@k))

  post_prob_crit <- post_beta(shape = shape_crit, theta0 = design@theta0)
  post_prob_nocrit <- post_beta(shape = shape_nocrit, theta0 = design@theta0)

  expect_true(all(post_prob_crit > 0.99))
  expect_true(all(post_prob_nocrit <= 0.99))
})

test_that("get_crit and get_crit_pool are identical with one basket", {
  design <- setupOneStageBasket(k = 1, shape1 = 1, shape2 = 1, theta0 = 0.2)
  crit_pool <- get_crit_pool(design = design, n = 20, lambda = 0.99)
  crit <- get_crit(design = design, n = 20, lambda = 0.99)

  expect_equal(crit_pool, crit)
})

test_that("get_crit_pool returns NA if sample size is too small", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.8)
  crit <- get_crit_pool(design = design, n = 10, lambda = 0.99)

  expect_true(is.na(crit))
})

test_that("get_targ works", {
  targ_toer <- get_targ(theta0 = 0.2, theta1 = c(0.5, 0.2, 0.2), prob = "toer")
  targ_pwr <- get_targ(theta0 = 0.2, theta1 = c(0.5, 0.2, 0.2), prob = "pwr")

  expect_equal(targ_toer, c(FALSE, TRUE, TRUE))
  expect_equal(targ_pwr, c(TRUE, FALSE, FALSE))
})

test_that("prune_weights works", {
  design <- setupOneStageBasket(k = 6, shape1 = 1, shape2 = 1, theta0 = 0.2)
  weight_mat <- weights_fujikawa(design = design, n = 15, epsilon = 2, tau = 0,
    logbase = 2, prune = FALSE)
  weight_mat <- prune_weights(weight_mat, cut = 8)

  r <- c(5, 6, 7, 8, 9, 10)
  shape_post <- matrix(c(design@shape1 + r, design@shape2 + 15 - r),
    byrow = TRUE, ncol = design@k)
  shape_borrow <- beta_borrow(weight_mat = weight_mat, design = design, n = 15,
    r = r)

  expect_equal(shape_post[, 1:3], shape_borrow[, 1:3])
  expect_false(any(shape_post[, 4:6] == shape_borrow[, 4:6]))
})

test_that("vectorization of get_prob works", {
  prob1 <- get_prob(n = 5, r = 2, theta = 0.2)
  prob2 <- get_prob(n = 10, r = 3, theta = 0.4)
  prob3 <- get_prob(n = 15, r = 4, theta = 0.5)
  prob_prod <- prob1 * prob2 * prob3
  prob_all <- get_prob(n = c(5, 10, 15), r = c(2, 3, 4),
    theta = c(0.2, 0.4, 0.5))

  expect_equal(prob_prod, prob_all)
})

test_that("beta_borrow works", {
  # Reproduced from Fujikawa et al., 2020, Supplement R Code
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
  weight_mat <- weights_fujikawa(design = design, n = 24, epsilon = 2,
    tau = 0.5, logbase = exp(1), prune = FALSE)
  r <- c(7, 2, 5)
  shape_borrow <- beta_borrow(weight_mat = weight_mat, design = design, n = 24,
    r = r)
  shape_expect <- matrix(c(12.9215409, 34.4051363, 6.33262523, 34.1087508,
    14.2283671, 47.5396861), nrow = 2)

  expect_equal(shape_borrow, shape_expect, tolerance = 10e-7)
})

