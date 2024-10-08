test_that("toer works for a single-stage design without pruning", {
  # Compare Fujikawa et al., 2020
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  # Proposed design (i) in Fujikawa et al.
  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_single_loop
  toer_group1 <- toer(design = design, n = 24, lambda = 0.99,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
    logbase = exp(1), prune = FALSE), results = "group")
  toer_fwer1 <- toer(design = design, n = 24, lambda = 0.99,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
    logbase = exp(1), prune = FALSE), results = "fwer")
  toer_loop1 <- reject_single_loop(design = design, p1 = rep(0.2, 3),
    n = 24, lambda = 0.99, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0, logbase = exp(1), prune = FALSE),
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
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0.5,
      logbase = exp(1), prune = FALSE), results = "group")
  toer_fwer2 <- toer(design = design, n = 24, lambda = 0.99,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0.5,
      logbase = exp(1), prune = FALSE), results = "fwer")
  toer_loop2 <- reject_single_loop(design = design, p1 = rep(0.2, 3),
    n = 24, lambda = 0.99, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0.5, logbase = exp(1),
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
  toer_group3 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 24,
    lambda = 0.99, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0, logbase = 2, prune = FALSE),
    results = "group")
  toer_fwer3 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 24,
    lambda = 0.99, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0, logbase = 2, prune = FALSE),
    results = "fwer")
  toer_loop3 <- reject_single_loop(design = design, p1 = c(0.2, 0.4, 0.5),
    n = 24, lambda = 0.99, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0, logbase = 2,
    prune = FALSE), prob = "toer")

  expect_equal(toer_fwer3, toer_group3$fwer)
  expect_equal(toer_fwer3, toer_loop3$fwer)
  expect_equal(toer_group3$rejection_probabilities,
    toer_loop3$rejection_probabilities)

  # Compare the results of "fwer" and "group" when a global weight is used
  toer_group4 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 20,
    lambda = 0.95, weight_fun = weights_cpp,
    weight_params = list(a = 1, b = 1), globalweight_fun =
      globalweights_diff, globalweight_params = list(eps_global = 1),
    results = "group")
  toer_fwer4 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 20,
    lambda = 0.95, weight_fun = weights_cpp,
    weight_params = list(a = 1, b = 1), globalweight_fun =
      globalweights_diff, globalweight_params = list(eps_global = 1),
    results = "fwer")
  toer_loop4 <- reject_single_loop(design = design, p1 = c(0.2, 0.4, 0.5),
    n = 20, lambda = 0.95, weight_fun = weights_cpp,
    weight_params = list(a = 1, b = 1), globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1), prob = "toer")

  expect_equal(toer_fwer4, toer_group4$fwer)
  expect_equal(toer_fwer4, toer_loop4$fwer)
  expect_equal(toer_group4$rejection_probabilities,
    toer_loop4$rejection_probabilities)

  # Compare then results when a global weight and pruning is used
  toer_group5 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 15,
    lambda = 0.98, weight_fun = weights_jsd,
    weight_params = list(tau = 0, prune = TRUE),
    globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3),
    results = "group")
  mat_jsd <- weights_jsd(design = design, n = 15, lambda = 0.98,
    tau = 0, prune = TRUE, globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3))
  toer_prob5 <- reject_prob_group(design, p1 = c(0.2, 0.4, 0.5), n = 15,
    lambda = 0.98, weight_mat = mat_jsd, globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3), prob = "toer")
  toer_fwer5 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 15,
    lambda = 0.98, weight_fun = weights_jsd,
    weight_params = list(tau = 0, prune = TRUE),
    globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3),
    results = "fwer")
  toer_loop5 <- reject_single_loop(design = design, p1 = c(0.2, 0.4, 0.5),
    n = 15, lambda = 0.98, weight_fun = weights_jsd,
    weight_params = list(tau = 0, prune = TRUE),
    globalweight_fun = globalweights_fix, globalweight_params = list(w = 0.3),
    prob = "toer")

  expect_equal(toer_fwer5, toer_group5$fwer)
  expect_equal(toer_fwer5, toer_loop5$fwer)
  expect_equal(toer_group5$rejection_probabilities,
    toer_loop5$rejection_probabilities)
  expect_equal(toer_group5$rejection_probabilities,
    toer_prob5$rejection_probabilities)
})

test_that("toer works for a single-stage design with pruning", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_single_loop
  toer_group1 <- toer(design = design, n = 15, lambda = 0.95,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0.2,
      logbase = 2, prune = TRUE), results = "group")
  toer_fwer1 <- toer(design = design, n = 15, lambda = 0.95,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0.2,
      logbase = 2, prune = TRUE), results = "fwer")
  toer_loop1 <- reject_single_loop(design = design, p1 = rep(0.2, 3),
    n = 15, lambda = 0.95, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0.2, logbase = 2, prune = TRUE),
    prob = "toer")

  expect_equal(toer_fwer1, toer_group1$fwer)
  expect_equal(toer_fwer1, toer_loop1$fwer)
  expect_equal(toer_group1$rejection_probabilities,
    toer_loop1$rejection_probabilities)
})
test_that("toer works for a single-stage design with extremal tau and lambda", {
  design <- baskexact::setupOneStageBasket(k = 4, shape1 = 1, shape2 = 1,
                                           p0 = 0.15)
  toer1 <- baskexact::toer(design, p1 = c(0.15, 0.15, 0.15, 0.15), n = 10,
                  lambda = 0.9999999999,
                  weight_fun = baskexact::weights_fujikawa,
                  weight_params = list(epsilon = 2, tau = 1, logbase = 2.72),
                  results = "group")
  expect_equal(toer1$rejection_probabilities,
               c(0, 0, 0, 0))
  expect_equal(toer1$fwer,
               0)
})
test_that("toer works for a two-stage design", {
  # Compare Fujikawa et al., 2020
  design <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  # Proposed design (i) in Fujikawa et al.
  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_twostage_loop
  toer_group1 <- toer(design = design, n = 24, n1 = 15, lambda = 0.99,
    interim_fun = interim_postpred, interim_params = list(prob_futstop = 0.1,
      prob_effstop = 0.9), weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0, logbase = exp(1)),
    results = "group")
  toer_fwer1 <- toer(design = design, n = 24, n1 = 15, lambda = 0.99,
    interim_fun = interim_postpred, interim_params = list(prob_futstop = 0.1,
      prob_effstop = 0.9), weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0, logbase = exp(1)),
    results = "fwer")
  toer_loop1 <- reject_twostage_loop(design = design, p1 = c(0.2, 0.2, 0.2),
    n = 24, n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
      logbase = exp(1)), prob = "toer")

  # In Fujikawa et al., based on simulation:
  # Basketwise 0.013, 0.016, 0.016
  # Experimentwise: 0.035
  rej_expect1 <- c(0.01703198, 0.01703198, 0.01703198)
  fwer_expect1 <- 0.03722851

  expect_equal(toer_group1$rejection_probabilities, rej_expect1,
    tolerance = 10e-7)
  expect_equal(toer_fwer1, fwer_expect1, tolerance = 10e-7)
  expect_equal(toer_fwer1, toer_group1$fwer)
  expect_equal(toer_fwer1, toer_loop1$fwer)
  expect_equal(toer_group1$rejection_probabilities,
    toer_loop1$rejection_probabilities)

  # Proposed design (ii) in Fujikawa et al.
  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_twostage_loop
  toer_group2 <- toer(design = design, n = 24, n1 = 15, lambda = 0.99,
    interim_fun = interim_postpred, interim_params = list(prob_futstop = 0.1,
      prob_effstop = 0.9), weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0.5, logbase = exp(1)),
    results = "group")
  toer_fwer2 <- toer(design = design, n = 24, n1 = 15, lambda = 0.99,
    interim_fun = interim_postpred, interim_params = list(prob_futstop = 0.1,
      prob_effstop = 0.9), weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0.5, logbase = exp(1)),
    results = "fwer")
  toer_loop2 <- reject_twostage_loop(design = design, p1 = c(0.2, 0.2, 0.2),
    n = 24, n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0.5,
      logbase = exp(1)), prob = "toer")

  # In Fujikawa et al., based on simulation:
  # Basketwise 0.017, 0.021, 0.021
  # Experimentwise: 0.047
  rej_expect2 <- c(0.02175429, 0.02175429, 0.02175429)
  fwer_expect2 <- 0.04955128

  expect_equal(toer_group2$rejection_probabilities, rej_expect2,
    tolerance = 10e-7)
  expect_equal(toer_fwer2, fwer_expect2, tolerance = 10e-7)
  expect_equal(toer_fwer2, toer_group2$fwer)
  expect_equal(toer_fwer2, toer_loop2$fwer)
  expect_equal(toer_group2$rejection_probabilities,
    toer_loop2$rejection_probabilities)

  # Compare the results of "fwer" and "group" when null hypothesis is not
  # global null
  # Proposed design (i) in Fujikawa et al.
  # Compare the results of reject_prob_ew, reject_prob_group and
  # reject_twostage_loop
  toer_group3 <- toer(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
      logbase = exp(1)), results = "group")
  toer_fwer3 <- toer(design = design, p1 = c(0.5, 0.2, 0.2), n = 24,
    n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
      logbase = exp(1)), results = "fwer")
  toer_loop3 <- reject_twostage_loop(design = design, p1 = c(0.5, 0.2, 0.2),
    n = 24, n1 = 15, lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
      logbase = exp(1)), prob = "toer")

  # In Fujikawa et al., based on simulation:
  # Basketwise 0.806, 0.058, 0.068
  # Experimentwise: 0.808 (different definition)
  rej_expect3 <- c(0.79791970, 0.06210063, 0.06210063)
  fwer_expect3 <- 0.1079397

  expect_equal(toer_group3$rejection_probabilities, rej_expect3,
    tolerance = 10e-7)
  expect_equal(toer_fwer3, fwer_expect3, tolerance = 10e-7)
  expect_equal(toer_fwer3, toer_group3$fwer)
  expect_equal(toer_fwer3, toer_loop3$fwer)
  expect_equal(toer_group3$rejection_probabilities,
    toer_loop3$rejection_probabilities)

  # Compare the results of "fwer" and "group" when a global weight is used
  toer_group4 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 15, n1 = 7,
    lambda = 0.95, interim_fun = interim_postpred, weight_fun = weights_cpp,
    weight_params = list(a = 1, b = 1), globalweight_fun =
      globalweights_diff, globalweight_params = list(eps_global = 1),
    results = "group")
  toer_fwer4 <- toer(design = design, p1 = c(0.2, 0.4, 0.5), n = 15, n1 = 7,
    lambda = 0.95, interim_fun = interim_postpred, weight_fun = weights_cpp,
    weight_params = list(a = 1, b = 1), globalweight_fun =
      globalweights_diff, globalweight_params = list(eps_global = 1),
    results = "fwer")
  toer_loop4 <- reject_twostage_loop(design = design, p1 = c(0.2, 0.4, 0.5),
    n = 15, n1 = 7, lambda = 0.95, interim_fun = interim_postpred,
    weight_fun = weights_cpp, weight_params = list(a = 1, b = 1),
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1), prob = "toer")

  expect_equal(toer_fwer4, toer_group4$fwer)
  expect_equal(toer_fwer4, toer_loop4$fwer)
  expect_equal(toer_group4$rejection_probabilities,
    toer_loop4$rejection_probabilities)
})
