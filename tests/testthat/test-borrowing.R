test_that("beta_borrow works", {
  design <- setupOneStageBasket(k = 3, p0 = 0.2)

  ## Without Global Weight
  # Reproduced from Fujikawa et al., 2020, Supplement R Code
  weights_fuj1 <- weights_fujikawa(design = design, n = 24, epsilon = 2,
    tau = 0.5, logbase = exp(1), prune = FALSE)
  postshape_comp <- beta_borrow(weight_mat = weights_fuj1, design = design,
    n = 24, r = c(7, 2, 5))
  shape_expect <- matrix(c(12.9215409, 34.4051363, 6.33262523, 34.1087508,
    14.2283671, 47.5396861), nrow = 2)

  expect_equal(postshape_comp, shape_expect, tolerance = 10e-7)

  # When information is fully shared, posterior parameters only differ
  # by the amount of prior information that is shared
  weights_fuj <- weights_fujikawa(design = design, n = 10, lambda = 0.95,
    epsilon = 2, tau = 0.99, logbase = 2)
  weights_jsd <- weights_jsd(design = design, n = 10, lambda = 0.95,
    epsilon = 2, tau = 0.99, logbase = 2)

  postshape_fuj1 <- beta_borrow(weights_fuj, design = design, n = 10,
    r = c(3, 3, 3))
  postshape_jsd1 <- beta_borrow(weights_jsd, design = design, n = 10,
    r = c(3, 3, 3))

  expect_equal(postshape_fuj1, postshape_jsd1 + 2)

  # When no information is shared, posterior paramaters are identical
  postshape_fuj2 <- beta_borrow(weights_fuj, design = design, n = 10,
    r = c(3, 4, 5))
  postshape_jsd2 <- beta_borrow(weights_jsd, design = design, n = 10,
    r = c(3, 4, 5))

  expect_equal(postshape_fuj2, postshape_jsd2)

  ## With Global Weight
  postshape_fuj3 <- beta_borrow(weights_fuj, globalweight_fun =
      globalweights_diff, globalweight_params = list(eps_global = 2, w = 0.5),
    design = design, n = 10, r = c(3, 3, 3))
  postshape_jsd3 <- beta_borrow(weights_jsd, globalweight_fun =
      globalweights_diff, globalweight_params = list(eps_global = 2, w = 0.5),
    design = design, n = 10, r = c(3, 3, 3))

  expect_equal(postshape_fuj3, postshape_jsd3 + 1)

  # Only half of the information is shared if w = 0.5
  expect_equal(postshape_jsd1[1, ], postshape_jsd3[1, ] + 3)
  expect_equal(postshape_jsd1[2, ], postshape_jsd3[2, ] + 7)
  expect_equal(postshape_fuj1[1, ], postshape_fuj3[1, ] + 4)
  expect_equal(postshape_fuj1[2, ], postshape_fuj3[2, ] + 8)

  ## Additional tests with more baskets
  design4 <- setupOneStageBasket(k = 4, p0 = 0.2)
  r <- c(3, 3, 6, 8)
  n <- 15
  cpp_fun <- function(x, n, a, b) {
    1 / (1 + exp(a + b * log(n^(1 / 4) * abs(x[1]/n - x[2]/n))))
  }
  cpp_w1 <- c(cpp_fun(r[c(1, 2)], 15, 1.5, 2), cpp_fun(r[c(1, 3)], 15, 1.5, 2),
    cpp_fun(r[c(1, 4)], 15, 1.5, 2))
  s11 <- 1 + r[1] + r[2] * cpp_w1[1] + r[3] * cpp_w1[2] + r[4] * cpp_w1[3]
  s21 <- 1 + (n - r[1]) + (n - r[2]) * cpp_w1[1] + (n - r[3]) * cpp_w1[2] +
    (n - r[4]) * cpp_w1[3]

  cpp_w2 <- c(cpp_fun(r[c(2, 1)], 15, 1.5, 2), cpp_fun(r[c(2, 3)], 15, 1.5, 2),
    cpp_fun(r[c(2, 4)], 15, 1.5, 2))
  s12 <- 1 + r[2] + r[1] * cpp_w2[1] + r[3] * cpp_w2[2] + r[4] * cpp_w2[3]
  s22 <- 1 + (n - r[2]) + (n - r[1]) * cpp_w2[1] + (n - r[3]) * cpp_w2[2] +
    (n - r[4]) * cpp_w2[3]

  cpp_w3 <- c(cpp_fun(r[c(3, 1)], 15, 1.5, 2), cpp_fun(r[c(3, 2)], 15, 1.5, 2),
    cpp_fun(r[c(3, 4)], 15, 1.5, 2))
  s13 <- 1 + r[3] + r[1] * cpp_w3[1] + r[2] * cpp_w3[2] + r[4] * cpp_w3[3]
  s23 <- 1 + (n - r[3]) + (n - r[1]) * cpp_w3[1] + (n - r[2]) * cpp_w3[2] +
    (n - r[4]) * cpp_w3[3]

  cpp_w4 <- c(cpp_fun(r[c(4, 1)], 15, 1.5, 2), cpp_fun(r[c(4, 2)], 15, 1.5, 2),
    cpp_fun(r[c(4, 3)], 15, 1.5, 2))
  s14 <- 1 + r[4] + r[1] * cpp_w4[1] + r[2] * cpp_w4[2] + r[3] * cpp_w4[3]
  s24 <- 1 + (n - r[4]) + (n - r[1]) * cpp_w4[1] + (n - r[2]) * cpp_w4[2] +
    (n - r[3]) * cpp_w4[3]

  wmat_cpp <- weights_cpp(design = design4, n = 15, a = 1.5, b = 2)
  postshape_val1 <- beta_borrow(weight_mat = wmat_cpp, design = design4,
    n = 15, r = c(3, 3, 6, 8))
  postshape_val2 <- unname(val_borrow_cpp(design = design4, n = 15,
    r = c(3, 3, 6, 8), a = 1.5, b = 2))

  expect_equal(postshape_val1, postshape_val2)
  expect_equal(postshape_val1[, 1], c(s11, s21))
  expect_equal(postshape_val1[, 2], c(s12, s22))
  expect_equal(postshape_val1[, 3], c(s13, s23))
  expect_equal(postshape_val1[, 4], c(s14, s24))

  wmat <- basket_test(design = design4, n = 15, r = c(3, 3, 6, 8),
    lambda = 0.95, weight_fun = weights_cpp,
    weight_params = list(a = 1.5, b = 2))$weights

  expect_equal(unname(wmat[1, -1]), cpp_w1)
  expect_equal(unname(wmat[2, -2]), cpp_w2)
  expect_equal(unname(wmat[3, -3]), cpp_w3)
  expect_equal(unname(wmat[4, -4]), cpp_w4)

  postshape_val3 <- beta_borrow(weight_mat = wmat_cpp,
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 2), design = design4,
    n = 15, r = c(4, 9, 3, 1))
  postshape_val4 <- unname(val_borrow_cpp(design = design4, n = 15,
    r = c(4, 9, 3, 1), a = 1.5, b = 2, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 2)))

  expect_equal(postshape_val3, postshape_val4)

  wmat_fuj <- weights_fujikawa(design = design4, n = 15, epsilon = 1.5,
    tau = 0.2, logbase = 2)
  postshape_val5 <- beta_borrow(weight_mat = wmat_fuj, design = design4,
    n = 15, r = c(3, 3, 6, 8))
  postshape_val6 <- unname(val_borrow_fujikawa(design = design4, n = 15,
    r = c(3, 3, 6, 8), epsilon = 1.5, tau = 0.2, logbase = 2))

  expect_equal(postshape_val5, postshape_val6)

  design5 <- setupOneStageBasket(k = 5, p0 = 0.15)
  postshape_val7 <- beta_borrow(weight_mat = wmat_cpp, design = design5,
    n = 15, r = c(3, 5, 7, 9, 12))
  postshape_val8 <- unname(val_borrow_cpp(design = design5, n = 15,
    r = c(3, 5, 7, 9, 12), a = 1.5, b = 2))

  expect_equal(postshape_val7, postshape_val8)

  # Test cpp function
  wmat_valid <- weights_cpp(design = design5, n = 17, a = 2, b = 1.5)
  res <- val_borrow_mat(design5, 20, r = c(7, 9, 12, 7, 3), wmat_valid)
  expect_equal(res[[1]], res[[2]])
  expect_true(isSymmetric(res[[1]]))
})

test_that("beta_borrow_int works", {
  design <- setupTwoStageBasket(k = 3, p0 = 0.2)

  # Reproduced from Fujikawa et al., 2020 Supplement R code
  weights_fuj1 <- weights_fujikawa(design = design, n = 24, n1 = 15,
    lambda = 0.99, epsilon = 2, tau = 0.5, logbase = exp(1))
  postshape_comp <- beta_borrow_int(weights_fuj1, design = design, n = 24,
    n1 = 15, r = c(7, 1, 3), res_int = c(0, -1, -1))
  shape_expect <- c(11.3577, 28.91253, 4.603939, 23.4628,
    12.01737, 37.87443)

  expect_equal(as.vector(postshape_comp), shape_expect, tolerance = 1e-7)

  # When information is fully shared, posterior parameters only differ
  # by the amount of prior information that is shared
  weights_fuj <- weights_fujikawa(design = design, n = 15, n1 = 7,
    lambda = 0.95, epsilon = 2, tau = 0.99, logbase = 2)
  weights_jsd <- weights_jsd(design = design, n = 15, n1 = 7, lambda = 0.95,
    epsilon = 2, tau = 0.99, logbase = 2)

  postshape_fuj1 <- beta_borrow_int(weights_fuj, design = design, n = 15,
    n1 = 7, r = c(3, 3, 3), res_int = c(-1, -1, -1))
  postshape_jsd1 <- beta_borrow_int(weights_jsd, design = design, n = 15,
    n1 = 7, r = c(3, 3, 3), res_int = c(-1, -1, -1))

  expect_equal(postshape_fuj1, postshape_jsd1 + 2)

  # When no information is shared, posterior paramaters are identical
  postshape_fuj2 <- beta_borrow_int(weights_fuj, design = design, n = 15,
    n1 = 7, r = c(3, 4, 5), res_int = c(-1, -1, -1))
  postshape_jsd2 <- beta_borrow_int(weights_jsd, design = design, n = 15,
    n1 = 7, r = c(3, 4, 5), res_int = c(-1, -1, -1))

  expect_equal(postshape_fuj2, postshape_jsd2)
})
