test_that("toer works without pruning", {
  # Compare Fujikawa et al., 2020
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)

  # Proposed design (i) in Fujikawa et al.
  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_single_loop
  toer_group1 <- toer(design = design, n = 24, lambda = 0.99,
    weight_fun = weights_fujikawa, tuning_params = list(epsilon = 2, tau = 0,
    logbase = exp(1), prune = FALSE), results = "group")
  toer_fwer1 <- toer(design = design, n = 24, lambda = 0.99,
    weight_fun = weights_fujikawa, tuning_params = list(epsilon = 2, tau = 0,
    logbase = exp(1), prune = FALSE), results = "fwer")
  toer_loop1 <- reject_single_loop(design = design, theta1 = rep(0.2, 3),
    n = 24, lambda = 0.99, weight_fun = weights_fujikawa,
    tuning_params = list(epsilon = 2, tau = 0, logbase = exp(1), prune = FALSE),
    prob = "toer")

  # In Fujikawa et al., based on simulation:
  # Basketwise 0.019, 0.020, 0.022
  # Experimentwise: 0.035
  rej_expect1 <- c(0.02158174, 0.02158174, 0.02158174)
  fwer_expect1 <- 0.03600149

  expect_equal(toer_group1$rejection_probabilities, rej_expect1,
    tolerance = 10e-7)
  expect_equal(toer_fwer1, fwer_expect1, tolerance = 10e-7)
  expect_equal(toer_fwer1, toer_group1$fwer)
  expect_equal(toer_fwer1, toer_loop1$fwer)
  expect_equal(toer_group1$rejection_probabilities,
    toer_loop1$rejection_probabilities)

  # Proposed design (ii) in Fujikawa et al.
  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_single_loop
  toer_group2 <- toer(design = design, n = 24, lambda = 0.99,
    weight_fun = weights_fujikawa, tuning_params = list(epsilon = 2, tau = 0.5,
      logbase = exp(1), prune = FALSE), results = "group")
  toer_fwer2 <- toer(design = design, n = 24, lambda = 0.99,
    weight_fun = weights_fujikawa, tuning_params = list(epsilon = 2, tau = 0.5,
      logbase = exp(1), prune = FALSE), results = "fwer")
  toer_loop2 <- reject_single_loop(design = design, theta1 = rep(0.2, 3),
    n = 24, lambda = 0.99, weight_fun = weights_fujikawa,
    tuning_params = list(epsilon = 2, tau = 0.5, logbase = exp(1),
    prune = FALSE), prob = "toer")

  # In Fujikawa et al., based on simulation:
  # Basketwise: 0.029, 0.032, 0.034
  # Experimentwise: 0.063
  rej_expect2 <- c(0.03239555, 0.03239555, 0.03239555)
  fwer_expect2 <- 0.06315308

  expect_equal(toer_group2$rejection_probabilities, rej_expect2,
    tolerance = 10e-7)
  expect_equal(toer_fwer2, fwer_expect2, tolerance = 10e-7)
  expect_equal(toer_fwer2, toer_group2$fwer)
  expect_equal(toer_fwer2, toer_loop2$fwer)
  expect_equal(toer_group2$rejection_probabilities,
    toer_loop2$rejection_probabilities)

  # Compare the results of "fwer" and "group" when null hypothesis is not
  # global null
  toer_group3 <- toer(design = design, theta1 = c(0.2, 0.4, 0.5), n = 24,
    lambda = 0.99, weight_fun = weights_fujikawa,
    tuning_params = list(epsilon = 1, tau = 0, logbase = 2, prune = FALSE),
    results = "group")
  toer_fwer3 <- toer(design = design, theta1 = c(0.2, 0.4, 0.5), n = 24,
    lambda = 0.99, weight_fun = weights_fujikawa,
    tuning_params = list(epsilon = 1, tau = 0, logbase = 2, prune = FALSE),
    results = "fwer")
  toer_loop3 <- reject_single_loop(design = design, theta1 = c(0.2, 0.4, 0.5),
    n = 24, lambda = 0.99, weight_fun = weights_fujikawa,
    tuning_params = list(epsilon = 1, tau = 0, logbase = 2,
    prune = FALSE), prob = "toer")

  expect_equal(toer_fwer3, toer_group3$fwer)
  expect_equal(toer_fwer3, toer_loop3$fwer)
  expect_equal(toer_group3$rejection_probabilities,
    toer_loop3$rejection_probabilities)
})

test_that("toer works with pruning", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)

  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_single_loop
  toer_group1 <- toer(design = design, n = 15, lambda = 0.95,
    weight_fun = weights_fujikawa, tuning_params = list(epsilon = 1, tau = 0.2,
      logbase = 2, prune = TRUE), results = "group")
  toer_fwer1 <- toer(design = design, n = 15, lambda = 0.95,
    weight_fun = weights_fujikawa, tuning_params = list(epsilon = 1, tau = 0.2,
      logbase = 2, prune = TRUE), results = "fwer")
  toer_loop1 <- reject_single_loop(design = design, theta1 = rep(0.2, 3),
    n = 15, lambda = 0.95, weight_fun = weights_fujikawa,
    tuning_params = list(epsilon = 1, tau = 0.2, logbase = 2, prune = TRUE),
    prob = "toer")

  expect_equal(toer_fwer1, toer_group1$fwer)
  expect_equal(toer_fwer1, toer_loop1$fwer)
  expect_equal(toer_group1$rejection_probabilities,
    toer_loop1$rejection_probabilities)
})
