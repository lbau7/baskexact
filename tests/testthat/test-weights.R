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

  # Single-stage design with pruning
  weight_jsd2 <- weights_jsd(design = design1, n = 15, lambda = 0.95,
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE)
  weight_fujikawa2 <- weights_fujikawa(design = design1, n = 15, lambda = 0.95,
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE)

  expect_equal(unclass(weight_jsd2), unclass(weight_fujikawa2))

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


test_that("weight_spline works", {
  # Single-stage design

  n_test <- 20
  clamplim_test <- c(0,1)
  diffknots_test <- c(1,0.5,0)
  weightknots_test <- c(0,0,0.8)

  design1 <- setupOneStageBasket(k = 4, shape1 = 1, shape2 = 1, p0 = 0.2)
  weight_mspline1 <- weights_spline(design = design1,
                                    n = n_test,
                                    diffknots  = diffknots_test,
                                    weightknots = weightknots_test,
                                    splinemethod = "monoH.FC",
                                    clamplim = clamplim_test)
  r <- c(15, 5, 8, 16)

  elmnts <- all_combs <- t(utils::combn(r, 2)) + 1
  weights <- as.vector(weight_mspline1[elmnts])
  weights_expected <- c(0.0000, 0.1224, 0.7128, 0.5096, 0.0000, 0.0576)

  expect_equal(weights, weights_expected)


  #Right S3 class
  expect_s3_class(weight_mspline1, "pp")

  #Symmetric matrix in output
  expect_true(isSymmetric(unclass(weight_mspline1)))

  #Main diagonale (same
  upper_weightlim <- max(weightknots_test)
  expect_true(all.equal(rep(upper_weightlim,n_test+1) , diag(unclass(weight_mspline1))))
})
