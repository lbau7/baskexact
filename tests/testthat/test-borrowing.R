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
