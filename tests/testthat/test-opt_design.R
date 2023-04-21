test_that("opt_design works", {
  design <- setupOneStageBasket(k = 3, theta0 = 0.2)
  optres <- opt_design(design = design, n = 15, alpha = 0.05,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = c(1, 2),
    tau = c(0, 0.5)), scenarios = get_scenarios(design, 0.5), prec_digits = 4)

  lambdres <- adjust_lambda(design = design, alpha = 0.05, n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon =
        optres[1, ]$epsilon, tau = optres[1, ]$tau), prec_digits = 4)
  ecdres1 <- ecd(design = design, n = 15, weight_fun = weights_fujikawa,
    lambda = lambdres$lambda, weight_params = list(epsilon =
        optres[1, ]$epsilon, tau = optres[1, ]$tau))
  ecdres2 <- ecd(design = design, n = 15, theta1 = c(0.2, 0.2, 0.5),
    lambda = lambdres$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = optres[1, ]$epsilon, tau = optres[1, ]$tau))
  ecdres3 <- ecd(design = design, n = 15, theta1 = c(0.2, 0.5, 0.5),
    lambda = lambdres$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = optres[1, ]$epsilon, tau = optres[1, ]$tau))
  ecdres4 <- ecd(design = design, n = 15, theta1 = c(0.5, 0.5, 0.5),
    lambda = lambdres$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = optres[1, ]$epsilon, tau = optres[1, ]$tau))

  expect_equal(optres[1, 3], lambdres$lambda)
  expect_equal(optres[1, 4], ecdres1)
  expect_equal(optres[1, 5], ecdres2)
  expect_equal(optres[1, 6], ecdres3)
  expect_equal(optres[1, 7], ecdres4)
  expect_equal(optres[1, 8], mean(c(ecdres1, ecdres2, ecdres3, ecdres4)))
})
