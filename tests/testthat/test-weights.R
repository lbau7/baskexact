test_that("weight_fujikawa works", {
  # Single-stage design
  # Reproduced from Fujikawa et al., 2020, Supplement R code
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  weight_fuj1 <- weights_fujikawa(design = design1, n = 15, epsilon = 2,
    tau = 0, logbase = exp(1), prune = FALSE)
  r <- c(5, 1, 3)
  elmnts <- all_combs <- t(utils::combn(r, 2)) + 1
  weights <- as.vector(weight_fuj1[elmnts])
  weights_exp <- c(0.3206983, 0.7493639, 0.6509846)

  expect_equal(weights, weights_exp, tolerance = 10e-7)
  expect_s3_class(weight_fuj1, "fujikawa")
  expect_true(isSymmetric(unclass(weight_fuj1)))

  # Two-stage design
  design2 <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  weight_fuj2 <- weights_fujikawa(design = design2, n = 15, n1 = 7, epsilon = 2,
    tau = 0, logbase = exp(1), prune = FALSE)

  expect_s3_class(weight_fuj2, "fujikawa")
  expect_true(isSymmetric(unclass(weight_fuj2)))

  # Compare single-stage and two-stage weight matrices
  expect_equal(unclass(weight_fuj1), weight_fuj2[-(1:8), -(1:8)])
})

test_that("weight_jsd works", {
  # Single-stage design
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  weight_jsd1 <- weights_jsd(design = design1, n = 15, epsilon = 2, tau = 0,
    logbase = 2, prune = FALSE)
  weight_fujikawa1 <- weights_fujikawa(design = design1, n = 15, epsilon = 2,
    tau = 0, logbase = 2, prune = FALSE)

  # Weight matrix for weight_fujikawa and weight_jsd is identical,
  # only the class differs
  expect_equal(unclass(weight_jsd1), unclass(weight_fujikawa1))
  expect_s3_class(weight_jsd1, "pp")

  # Two-stage design
  design2 <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  weight_jsd3 <- weights_jsd(design = design2, n = 15, n1 = 7, epsilon = 2,
    tau = 0, logbase = 2, prune = FALSE)
  weight_fujikawa3 <- weights_fujikawa(design = design2, n = 15, n1 = 7,
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE)

  expect_equal(unclass(weight_jsd3), unclass(weight_fujikawa3))
  expect_s3_class(weight_jsd3, "pp")

  # Compare single-stage and two-stage weight matrices
  expect_equal(unclass(weight_jsd1), weight_jsd3[-(1:8), -(1:8)])
})

test_that("weight_cpp works", {
  # Single-stage design
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  weight_cpp1 <- weights_cpp(design = design1, n = 20, a = 1, b = 1)

  x11 <- c(rep(0, 7), rep(1, 13))
  x21 <- c(rep(0, 3), rep(1, 17))
  sks1 <- as.numeric(ks.test(x11, x21)$statistic)
  s1 <- 20^(1 / 4) * sks1
  w1 <- 1 / (1 + exp(1 + 1 * log(s1)))

  expect_equal(w1, weight_cpp1[14, 18])
  expect_s3_class(weight_cpp1, "pp")
  expect_true(isSymmetric(unclass(weight_cpp1)))

  # Two-stage design
  design2 <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  weight_cpp2 <- weights_cpp(design = design2, n = 20, n1 = 10, a = 1, b = 1)

  x12 <- c(rep(0, 6), rep(1, 14))
  x22 <- c(rep(0, 8), rep(1, 2))
  sks2 <- as.numeric(ks.test(x12, x22)$statistic)
  s2 <- 20^(1 / 4) * sks2
  w2 <- 1 / (1 + exp(1 + 1 * log(s2)))

  x13 <- c(rep(0, 5), rep(1, 5))
  x23 <- c(rep(0, 8), rep(1, 2))
  sks3 <- as.numeric(ks.test(x13, x23)$statistic)
  s3 <- 10^(1 / 4) * sks3
  w3 <- 1 / (1 + exp(1 + 1 * log(s3)))

  expect_equal(w2, weight_cpp2[26, 3])
  expect_equal(w3, weight_cpp2[6, 3])
  expect_s3_class(weight_cpp2, "pp")
  expect_true(isSymmetric(unclass(weight_cpp2)))

  # Compare single-stage and two-stage weight matrices
  weight_cpp3 <- weights_cpp(design = design1, n = 10, a = 1, b = 1)

  expect_equal(unclass(weight_cpp1), weight_cpp2[-(1:11), -(1:11)])
  expect_equal(unclass(weight_cpp3), weight_cpp2[1:11, 1:11])
})

test_that("weight_mml works", {
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  weights_mml1 <- weights_mml(design = design1, n = 20)

  design2 <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  weights_mml2 <- weights_mml(design = design2, n = 20, n1 = 10)

  weights_mml3 <- weights_mml(design = design1, n = 10)

  expect_equal(unclass(weights_mml1), weights_mml2[-(1:11), -(1:11)],
    tolerance = 1e-6)
  expect_equal(unclass(weights_mml3), weights_mml2[1:11, 1:11],
    tolerance = 1e-6)
})

test_that("weight_separate works", {
  # Single-stage design
  design <- setupOneStageBasket(k = 3, p0 = 0.2)

  toer1 <- toer(
    design = design,
    n = 20,
    lambda = 0.99,
    weight_fun = weights_separate,
    results = "group"
  )

  toer2 <- 0
  for (i in 0:20) {
    shape <- data.frame(shape = c(1 + i, 1 + 20 - i))
    rej <- post_beta(shape = shape, p0 = 0.2) >= 0.99
    if (rej) toer2 <- toer2 + get_prob(n = 20, r = i, p = 0.2)
  }

  expect_equal(toer1$rejection_probabilities[1], toer2)

  pow1 <- pow(
    design = design,
    p1 = c(0.5, 0.5, 0.5),
    n = 20,
    lambda = 0.99,
    weight_fun = weights_separate,
    results = "group",
  )

  pow2 <- 0
  for (i in 0:20) {
    shape <- data.frame(shape = c(1 + i, 1 + 20 - i))
    rej <- post_beta(shape = shape, p0 = 0.2) >= 0.99
    if (rej) pow2 <- pow2 + get_prob(n = 20, r = i, p = 0.5)
  }

  expect_equal(pow1$rejection_probabilities[1], pow2)

  ecd <- ecd(
    design = design,
    p1 = c(0.5, 0.5, 0.5),
    n = 20,
    lambda = 0.99,
    weight_fun = weights_separate
  )

  expect_equal(ecd, 3 * pow2)

  estim1 <- estim(
    design = design,
    p1 = c(0.4, 0.4, 0.4),
    n = 20,
    weight_fun = weights_separate
  )

  estim2 <- 0
  mse2 <- 0
  for (i in 0:20) {
    shape <- data.frame(shape = c(1 + i, 1 + 20 - i))
    prob <- get_prob(n = 20, r = i, p = 0.4)
    estim2 <- estim2 + mean_beta(shape) * prob
    mse2 <- mse2 + (mean_beta(shape) - 0.4)^2 * prob
  }

  expect_equal(estim1$Mean[1], as.numeric(estim2))
  expect_equal(estim1$MSE[1], as.numeric(mse2))

  # Two-stage design
  design2 <- setupTwoStageBasket(k = 3, p0 = 0.2)
  toer_2stage1 <- toer(
    design = design2,
    n = 14,
    n1 = 7,
    lambda = 0.99,
    interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_separate,
    results = "group"
  )

  toer_2stage2 <- 0
  for (i in 0:7) {
    shape <- data.frame(shape = c(1 + i, 1 + 7 - i))
    pbeta_int <- post_beta(shape = shape, p0 = 0.2)
    rej_interim <- pbeta_int > 0.9
    stop_interim <- pbeta_int < 0.1
    stop_interim <-
    if (rej_interim) {
      toer_2stage2 <- toer_2stage2 + get_prob(n = 7, r = i, p = 0.2)
    } else if (!rej_interim & !stop_interim) {
      for (j in 0:7) {
        shape <- data.frame(shape = c(1 + i + j, 1 + 14 - i - j))
        rej <- post_beta(shape = shape, p0 = 0.2) >= 0.99
        if (rej) {
          toer_2stage2 <- toer_2stage2 +
          get_prob(n = 7, r = i, p = 0.2) * get_prob(n = 7, r = j, p = 0.2)
        }
      }
    }
  }

  expect_equal(toer_2stage1$rejection_probabilities[1], toer_2stage2)
})

test_that("weight_separate works", {
  design1 <- setupOneStageBasket(k = 3, p0 = 0.2)
  weights1 <- weights_pool(design = design1, n = 10)

  expect_true(all(weights1 == 1))

  design2 <- setupTwoStageBasket(k = 3, p0 = 0.2)
  weights2 <- weights_pool(design = design2, n = 10, n1 = 5)

  expect_true(all(weights2 == 1))
})

