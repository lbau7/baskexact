test_that("opt_design works", {
  design <- setupOneStageBasket(k = 3, p0 = 0.2)
  optres1 <- opt_design(design = design, n = 15, alpha = 0.05,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = c(1, 2),
    tau = c(0)), globalweight_fun = globalweights_diff, globalweight_params =
      list(eps_global = c(1, 2)), scenarios = get_scenarios(design, 0.5),
    prec_digits = 4)

  lambdres1 <- adjust_lambda(design = design, alpha = 0.05, n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon =
        optres1[1, ]$epsilon, tau = optres1[1, ]$tau),
    globalweight_fun = globalweights_diff, globalweight_params =
      list(eps_global = optres1[1, ]$eps_global), prec_digits = 4)
  ecdres11 <- ecd(design = design, n = 15, weight_fun = weights_fujikawa,
    lambda = lambdres1$lambda, weight_params = list(epsilon =
        optres1[1, ]$epsilon, tau = optres1[1, ]$tau), globalweight_fun =
      globalweights_diff, globalweight_params = list(eps_global =
          optres1[1, ]$eps_global))
  ecdres21 <- ecd(design = design, n = 15, p1 = c(0.2, 0.2, 0.5),
    lambda = lambdres1$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = optres1[1, ]$epsilon,
      tau = optres1[1, ]$tau), globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = optres1[1, ]$eps_global))
  ecdres31 <- ecd(design = design, n = 15, p1 = c(0.2, 0.5, 0.5),
    lambda = lambdres1$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = optres1[1, ]$epsilon,
      tau = optres1[1, ]$tau), globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = optres1[1, ]$eps_global))
  ecdres41 <- ecd(design = design, n = 15, p1 = c(0.5, 0.5, 0.5),
    lambda = lambdres1$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = optres1[1, ]$epsilon,
      tau = optres1[1, ]$tau), globalweight_fun = globalweights_diff,
    globalweight_params = list(eps_global = optres1[1, ]$eps_global))

  expect_equal(optres1[1, 4], lambdres1$lambda)
  expect_equal(optres1[1, 5], ecdres11)
  expect_equal(optres1[1, 6], ecdres21)
  expect_equal(optres1[1, 7], ecdres31)
  expect_equal(optres1[1, 8], ecdres41)
  expect_equal(optres1[1, 9], mean(c(ecdres11, ecdres21, ecdres31, ecdres41)))

  # Border case with no tuning parameter values
  optres2 <- opt_design(design = design, n = 15, alpha = 0.05,
    weight_fun = weights_fujikawa, globalweight_fun = globalweights_diff,
    scenarios = get_scenarios(design, 0.5), prec_digits = 4)

  lambdres2 <- adjust_lambda(design = design, alpha = 0.05, n = 15,
    weight_fun = weights_fujikawa, globalweight_fun = globalweights_diff,
    prec_digits = 4)
  ecdres12 <- ecd(design = design, n = 15, weight_fun = weights_fujikawa,
    lambda = lambdres2$lambda, globalweight_fun = globalweights_diff)
  ecdres22 <- ecd(design = design, n = 15, p1 = c(0.2, 0.2, 0.5),
    lambda = lambdres2$lambda, weight_fun = weights_fujikawa,
    globalweight_fun = globalweights_diff)
  ecdres32 <- ecd(design = design, n = 15, p1 = c(0.2, 0.5, 0.5),
    lambda = lambdres2$lambda, weight_fun = weights_fujikawa,
    globalweight_fun = globalweights_diff)
  ecdres42 <- ecd(design = design, n = 15, p1 = c(0.5, 0.5, 0.5),
    lambda = lambdres2$lambda, weight_fun = weights_fujikawa,
    globalweight_fun = globalweights_diff)

  expect_equal(unname(optres2[1]), lambdres2$lambda)
  expect_equal(unname(optres2[2]), ecdres12)
  expect_equal(unname(optres2[3]), ecdres22)
  expect_equal(unname(optres2[4]), ecdres32)
  expect_equal(unname(optres2[5]), ecdres42)
  expect_equal(unname(optres2[6]), mean(c(ecdres12, ecdres22, ecdres32, ecdres42)))
})
