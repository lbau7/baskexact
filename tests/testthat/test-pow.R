test_that("pow works for a single-stage design without pruning", {
  # Compare groupwise power with experimentwise power with one active basket
  # and compare with reject_single_loop
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  pow_ewp <- pow(design = design1, p1 = c(0.5, 0.2, 0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE), results = "ewp")
  pow_group <- pow(design = design1, p1 = c(0.5, 0.2, 0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE), results = "group")
  pow_loop <- reject_single_loop(design = design1, p1 = c(0.5, 0.2, 0.2),
    n = 15, lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE), prob = "pwr")

  expect_equal(pow_ewp, pow_group$ewp)
  expect_equal(pow_loop$ewp, pow_group$ewp)
  expect_equal(pow_loop$rejection_probabilities,
    pow_group$rejection_probabilities)

  # Compare rejection probabilities of pow and toer
  design2 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.3)
  toer_probs <- toer(design = design2, p1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = FALSE), results = "group")
  pow_probs <- pow(design = design2, p1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = FALSE), results = "group")

  expect_equal(toer_probs$rejection_probabilities,
    pow_probs$rejection_probabilities)
  expect_false(toer_probs$fwer == pow_probs$ewp)
})

test_that("pow works for a single-stage design with pruning", {
  # Compare groupwise power with experimentwise power with one active basket
  # and compare with reject_single_loop
  design1 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
  pow_ewp <- pow(design = design1, p1 = c(0.5, 0.2,0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE), results = "ewp")
  pow_group <- pow(design = design1, p1 = c(0.5, 0.2, 0.2), n = 15,
    lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE), results = "group")
  pow_loop <- reject_single_loop(design = design1, p1 = c(0.5, 0.2, 0.2),
    n = 15, lambda = 0.99, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE), prob = "pwr")

  expect_equal(pow_ewp, pow_group$ewp)
  expect_equal(pow_loop$ewp, pow_group$ewp)
  expect_equal(pow_loop$rejection_probabilities,
    pow_group$rejection_probabilities)

  # Compare rejection probabilities of pow and toer
  design2 <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.3)
  toer_probs <- toer(design = design2, p1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = TRUE), results = "group")
  pow_probs <- pow(design = design2, p1 = c(0.6, 0.6, 0.3), n = 20,
    lambda = 0.95, weight_fun = weights_fujikawa, weight_params = list(
    epsilon = 1, tau = 0.2, logbase = 2, prune = TRUE), results = "group")

  expect_equal(toer_probs$rejection_probabilities,
    pow_probs$rejection_probabilities)
  expect_false(toer_probs$fwer == pow_probs$ewp)

  # Compare then results when a global weight and pruning is used
  pow_group2 <- pow(design = design1, p1 = c(0.2, 0.4, 0.5), n = 15,
    lambda = 0.98, weight_fun = weights_jsd,
    weight_params = list(tau = 0, prune = TRUE),
    globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3), results = "group")
  pow_ewp2 <- pow(design = design1, p1 = c(0.2, 0.4, 0.5), n = 15,
    lambda = 0.98, weight_fun = weights_jsd,
    weight_params = list(tau = 0, prune = TRUE),
    globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3), results = "ewp")
  mat_jsd <- weights_jsd(design = design1, n = 15, lambda = 0.98,
    tau = 0, prune = TRUE, globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3))
  pow_prob2 <- reject_prob_group(design, p1 = c(0.2, 0.4, 0.5), n = 15,
    lambda = 0.98, weight_mat = mat_jsd, globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.3), prob = "pow")
  pow_loop2 <- reject_single_loop(design = design1, p1 = c(0.2, 0.4, 0.5),
    n = 15, lambda = 0.98, weight_fun = weights_jsd,
    weight_params = list(tau = 0, prune = TRUE),
    globalweight_fun = globalweights_fix, globalweight_params = list(w = 0.3),
    prob = "pow")
  pow_false2 <- pow(design = design1, p1 = c(0.2, 0.4, 0.5), n = 15,
    lambda = 0.98, weight_fun = weights_jsd,
    weight_params = list(tau = 0, prune = TRUE),
    globalweight_fun = globalweights_fix,
    globalweight_params = list(w = 0.9), results = "ewp")

  expect_equal(pow_group2$rejection_probabilities,
    pow_loop2$rejection_probabilities)
  expect_equal(pow_group2$ewp, pow_ewp2)
  expect_equal(pow_group2$rejection_probabilities,
    pow_loop2$rejection_probabilities)
  expect_equal(pow_group2$rejection_probabilities,
    pow_prob2$rejection_probabilities)
  expect_equal(pow_prob2$ewp, pow_ewp2)
  expect_true(pow_ewp2 != pow_false2)
})

test_that("pow works for a two-stage design", {
  design <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  # Compare the results of toer, pow and loop
  rej_toer <- toer(design = design, p1 = c(0.2, 0.5, 0.5), n = 18,
    n1 = 9, lambda = 0.99, interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_cpp, weight_params = list(a = 1.5, b = 1.5),
    results = "group")
  rej_pow <- pow(design = design, p1 = c(0.2, 0.5, 0.5), n = 18,
    n1 = 9, lambda = 0.99, interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_cpp, weight_params = list(a = 1.5, b = 1.5),
    results = "group")
  rej_loop <- reject_twostage_loop(design = design, p1 = c(0.2, 0.5, 0.5),
    n = 18, n1 = 9, lambda = 0.99, interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_cpp, weight_params = list(a = 1.5, b = 1.5),
    prob = "pow")
  ewp_pow <- pow(design = design, p1 = c(0.2, 0.5, 0.5), n = 18,
    n1 = 9, lambda = 0.99, interim_fun = interim_posterior,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_cpp, weight_params = list(a = 1.5, b = 1.5),
    results = "ewp")

  expect_equal(rej_toer$rejection_probabilities,
    rej_pow$rejection_probabilities)
  expect_equal(rej_loop$rejection_probabilities,
    rej_pow$rejection_probabilities)
  expect_equal(rej_loop$ewp, rej_pow$ewp)
  expect_equal(ewp_pow, rej_pow$ewp)

  # With a global weight
  rej_toer2 <- toer(design = design, p1 = c(0.2, 0.2, 0.5), n = 15,
    n1 = 7, lambda = 0.97, interim_fun = interim_postpred,
    weight_fun = weights_fujikawa,
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1), results = "group")
  rej_pow2 <- pow(design = design, p1 = c(0.2, 0.2, 0.5), n = 15,
    n1 = 7, lambda = 0.97, interim_fun = interim_postpred,
    weight_fun = weights_fujikawa,
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1), results = "group")
  rej_loop2 <- reject_twostage_loop(design = design, p1 = c(0.2, 0.2, 0.5),
    n = 15, n1 = 7, lambda = 0.97, interim_fun = interim_postpred,
    weight_fun = weights_fujikawa, globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1), prob = "pow")
  ewp_pow2 <- pow(design = design, p1 = c(0.2, 0.2, 0.5), n = 15,
    n1 = 7, lambda = 0.97, interim_fun = interim_postpred,
    weight_fun = weights_fujikawa,
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 1), results = "ewp")
  ewp_false2 <- pow(design = design, p1 = c(0.2, 0.2, 0.5), n = 15,
    n1 = 7, lambda = 0.97, interim_fun = interim_postpred,
    weight_fun = weights_fujikawa,
    globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = 2), results = "ewp")

  expect_equal(rej_toer2$rejection_probabilities,
    rej_pow2$rejection_probabilities)
  expect_equal(rej_loop2$rejection_probabilities,
    rej_pow2$rejection_probabilities)
  expect_equal(rej_loop2$ewp, rej_pow2$ewp)
  expect_equal(ewp_pow2, rej_pow2$ewp)
  expect_true(ewp_pow2 != ewp_false2)
})
