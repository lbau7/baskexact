test_that("ess works", {
  # Compare Fujikawa et al., 2020
  design <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  # Proposed design (i) in Fujikawa et al. - p = (0.2, 0.2, 0.5)
  ess1 <- ess(design = design, p1 = c(0.2, 0.2, 0.5), n = 24, n1 = 15,
    lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0, logbase = exp(1)))

  # In Fujikawa et al., based on simulation: 16.44, 16.46, 18.58
  ess_expect1 <- c(16.48675, 16.48675, 18.63614)
  expect_equal(ess1, ess_expect1, tolerance = 10e-7)

  # Proposed design (ii) in Fujikawa et al. - p = (0.2, 0.5, 0.5)
  ess2 <- ess(design = design, p1 = c(0.2, 0.5, 0.5), n = 24, n1 = 15,
    lambda = 0.99, interim_fun = interim_postpred,
    interim_params = list(prob_futstop = 0.1, prob_effstop = 0.9),
    weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2, tau = 0.5, logbase = exp(1)))

  # In Fujikawa et al., based on simulation: 16.49, 18.65, 18.67
  ess_expect2 <- c(16.47148, 18.67549, 18.67549)
  expect_equal(ess2, ess_expect2, tolerance = 10e-7)
})
