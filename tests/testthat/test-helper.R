test_that("get_crit works", {
  # Reproduced from Fujikawa et al. 2020
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  crit <- get_crit(design = design, n = 24, lambda = 0.99)

  expect_equal(crit, 10)
})

test_that("get_crit returns NA if sample size is too small", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.7)
  crit <- get_crit(design = design, n = 11, lambda = 0.99)

  expect_equal(crit, NA_integer_)
})

test_that("get_crit_pool works", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
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

  post_prob_crit <- post_beta(shape = shape_crit, p0 = design@p0)
  post_prob_nocrit <- post_beta(shape = shape_nocrit, p0 = design@p0)

  expect_true(all(post_prob_crit > 0.99))
  expect_true(all(post_prob_nocrit <= 0.99))
})

test_that("get_crit and get_crit_pool are identical with one basket", {
  design <- setupOneStageBasket(k = 1, shape1 = 1, shape2 = 1, p0 = 0.2)
  crit_pool <- get_crit_pool(design = design, n = 20, lambda = 0.99)
  crit <- get_crit(design = design, n = 20, lambda = 0.99)

  expect_equal(crit_pool, crit)
})

test_that("get_crit_pool returns NA if sample size is too small", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.8)
  crit <- get_crit_pool(design = design, n = 10, lambda = 0.99)

  expect_true(is.na(crit))
})

test_that("get_targ works", {
  targ_toer <- get_targ(p0 = 0.2, p1 = c(0.5, 0.2, 0.2), prob = "toer")
  targ_pwr <- get_targ(p0 = 0.2, p1 = c(0.5, 0.2, 0.2), prob = "pwr")

  expect_equal(targ_toer, c(FALSE, TRUE, TRUE))
  expect_equal(targ_pwr, c(TRUE, FALSE, FALSE))
})

test_that("prune_weights works", {
  design <- setupOneStageBasket(k = 6, shape1 = 1, shape2 = 1, p0 = 0.2)
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
  prob1 <- get_prob(n = 5, r = 2, p = 0.2)
  prob2 <- get_prob(n = 10, r = 3, p = 0.4)
  prob3 <- get_prob(n = 15, r = 4, p = 0.5)
  prob_prod <- prob1 * prob2 * prob3
  prob_all <- get_prob(n = c(5, 10, 15), r = c(2, 3, 4),
    p = c(0.2, 0.4, 0.5))

  expect_equal(prob_prod, prob_all)
})

test_that("mean_beta works", {
  shape <- matrix(rep(1, 6), ncol = 3)
  res <- mean_beta(shape)

  expect_equal(res, rep(0.5, 3))
})

test_that("post_pred works", {
  # Reproduced from Fujikawa et al., 2020, Supplement R Code
  design <- setupTwoStageBasket(k = 3, p0 = 0.2)
  crit <- get_crit(design = design, n = 24, lambda = 0.99)
  shape <- matrix(c(6, 11, 2, 15, 4, 13), nrow = 2)

  weights <- weights_fujikawa(design = design, n = 25, n1 = 15,
    epsilon = 2, tau = 0.5, logbase = exp(1))
  shape_post <- beta_borrow(weight_mat = weights, design = design, n = 15,
    r = c(5, 1, 3))
  res <- post_pred(n = 24, n1 = 15, r1 = c(5, 1, 3), shape = shape_post,
    crit = crit)
  prob_expect <- c(0.1309378, 4.769149e-06, 0.002926508)

  expect_equal(res, prob_expect, tolerance = 1e-6)
})
