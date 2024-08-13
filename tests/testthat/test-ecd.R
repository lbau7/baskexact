test_that("ecd works for a single-stage design", {
  design <- setupOneStageBasket(k = 3, p0 = 0.2)
  # Calculate expected number of correction decisions for a design with equal
  # p1
  ecd1 <- ecd(design, p1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)
  ecd_loop1 <- ecd_loop(design, p1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)

  expect_equal(ecd1, ecd_loop1)

  # Calculate expected number of correction decisions for a design with unequal
  # p1
  ecd2 <- ecd(design, p1 = c(0.2, 0.4, 0.6), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)
  ecd_loop2 <- ecd_loop(design, p1 = c(0.2, 0.4, 0.6), n = 20, lambda = 0.9,
    weight_fun = weights_fujikawa)

  expect_equal(ecd2, ecd_loop2)

  # Calculate expected number of correct decisions with global weights with
  # equal p1
  ecd3 <- ecd(design, p1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))
  ecd_loop3 <- ecd_loop(design, p1 = c(0.2, 0.2, 0.2), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))

  expect_equal(ecd3, ecd_loop3)

  # Calculate expected number of correct decisions with global weights with
  # unequal p1
  ecd4 <- ecd(design, p1 = c(0.2, 0.3, 0.4), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))
  ecd_loop4 <- ecd_loop(design, p1 = c(0.2, 0.3, 0.4), n = 20, lambda = 0.9,
    weight_fun = weights_cpp, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1))

  expect_equal(ecd4, ecd_loop4)
})

test_that("ecd works for a two-stage design", {
  design <- setupTwoStageBasket(k = 3, p0 = 0.2)
  # Calculate expected number of correction decisions for a design with equal
  # p1
  ecd1 <- ecd(design, p1 = c(0.2, 0.2, 0.2), n = 15, n1 = 7, lambda = 0.95,
    weight_fun = weights_mml, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.15, prob_effstop = 0.85))
  ecd_loop1 <- ecd_twostage_loop(design, p1 = c(0.2, 0.2, 0.2), n = 15, n1 = 7,
    lambda = 0.95, weight_fun = weights_mml, weight_params = list(),
    interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.15, prob_effstop = 0.85))

  expect_equal(ecd1, ecd_loop1)

  # Calculate expected number of correction decisions for a design with unequal
  # p1
  ecd2 <- ecd(design, p1 = c(0.2, 0.4, 0.5), n = 15, n1 = 7, lambda = 0.9,
    weight_fun = weights_mml, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9))
  ecd_loop2 <- ecd_twostage_loop(design, p1 = c(0.2, 0.4, 0.5), n = 15, n1 = 7,
    lambda = 0.9, weight_fun = weights_mml, weight_params = list(),
    interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9))

  expect_equal(ecd2, ecd_loop2)

  # Calculate expected number of correct decisions with global weights with
  # equal p1
  ecd3 <- ecd(design, p1 = c(0.2, 0.2, 0.2), n = 15, n1 = 7, lambda = 0.9,
    weight_fun = weights_cpp, weight_params = list(a = 1, b = 2),
    interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 2))

  ecd_loop3 <- ecd_twostage_loop(design, p1 = c(0.2, 0.2, 0.2), n = 15, n1 = 7,
    lambda = 0.9, weight_fun = weights_cpp, weight_params = list(a = 1, b = 2),
    interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 2))

  expect_equal(ecd3, ecd_loop3)

  # Calculate expected number of correct decisions with global weights with
  # unequal p1
  ecd4 <- ecd(design, p1 = c(0.2, 0.4, 0.5), n = 15, n1 = 7, lambda = 0.9,
    weight_fun = weights_cpp, weight_params = list(a = 1, b = 2),
    interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.7))

  ecd_loop4 <- ecd_twostage_loop(design, p1 = c(0.2, 0.4, 0.5), n = 15, n1 = 7,
    lambda = 0.9, weight_fun = weights_cpp, weight_params = list(a = 1, b = 2),
    interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.7))

  expect_equal(ecd4, ecd_loop4)
})

test_that("ecd works for a single-stage design with small sample size and
          high lambda", {
    design <- setupOneStageBasket(k = 4, shape1 = 1, shape2 = 1, p0 = 0.15)
    ecd1 <- ecd(design, p1 = c(0.15, 0.15, 0.15, 0.15), n = 10,
                lambda = 0.9999999999999998889777,
                weight_fun = weights_fujikawa)
    expect_equal(ecd1, 4)
})
