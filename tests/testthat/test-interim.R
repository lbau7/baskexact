test_that("analyses with interim_posterior work", {
  design <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  ess1 <- ess(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_posterior,
    interim_params = list(
      prob_futstop = 0.3,
      prob_effstop = 0.95
    ),
    weight_fun = weights_fujikawa,
    weight_params = list(
      epsilon = 2,
      tau = 0,
      logbase = exp(1)
    ),
    results = "group")

  ess2 <- ess(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_posterior,
    interim_params = list(
      prob_futstop = 0.3,
      prob_effstop = 0.9
    ),
    weight_fun = weights_fujikawa,
    weight_params = list(
      epsilon = 2,
      tau = 0,
      logbase = exp(1)
    ),
    results = "group")

  ess3 <- ess(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_posterior,
    interim_params = list(
      prob_futstop = 0.5,
      prob_effstop = 0.95
    ),
    weight_fun = weights_fujikawa,
    weight_params = list(
      epsilon = 2,
      tau = 0,
      logbase = exp(1)
    ),
    results = "group")

  # Expected sample sizes drecrease when boundary for efficacy stop
  # is lowered and boundary for futility stop is raised.
  expect_true(all(ess1 > ess2))
  expect_true(all(ess1 > ess3))
})

test_that("analyses with interim_postpred work", {
  design <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  ess1 <- ess(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(
      prob_futstop = 0.1,
      prob_effstop = 0.9
    ),
    weight_fun = weights_fujikawa,
    weight_params = list(
      epsilon = 2,
      tau = 0,
      logbase = exp(1)
    ),
    results = "group")

  ess2 <- ess(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(
      prob_futstop = 0.1,
      prob_effstop = 0.8
    ),
    weight_fun = weights_fujikawa,
    weight_params = list(
      epsilon = 2,
      tau = 0,
      logbase = exp(1)
    ),
    results = "group")

  ess3 <- ess(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(
      prob_futstop = 0.3,
      prob_effstop = 0.9
    ),
    weight_fun = weights_fujikawa,
    weight_params = list(
      epsilon = 2,
      tau = 0,
      logbase = exp(1)
    ),
    results = "group")

  # Expected sample sizes drecrease when boundary for efficacy stop
  # is lowered and boundary for futility stop is raised.
  expect_true(all(ess1 > ess2))
  expect_true(all(ess1 > ess3))
})
