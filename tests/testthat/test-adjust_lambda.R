test_that("adjust_lambda works for a single-stage design", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  # Without Pruning
  adj_res1 <- adjust_lambda(design = design, alpha = 0.025, n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
    logbase = 2, prune = FALSE), prec_digits = 4)
  lambda_high1 <- adj_res1$lambda - 0.0001
  toer_adj1 <- toer(design = design, n = 15, lambda = adj_res1$lambda,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
    logbase = 2, prune = FALSE), results = "fwer")
  toer_high1 <- toer(design = design, n = 15, lambda = lambda_high1,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
    logbase = 2, prune = FALSE), results = "fwer")

  expect_equal(adj_res1$toer, toer_adj1)
  expect_lte(adj_res1$toer, 0.025)
  expect_gt(toer_high1, 0.025)

  # With Pruning
  adj_res2 <- adjust_lambda(design = design, alpha = 0.025, n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
    logbase = 2, prune = TRUE), prec_digits = 4)
  lambda_high2 <- adj_res2$lambda - 0.0001
  toer_adj2 <- toer(design = design, n = 15, lambda = adj_res2$lambda,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
    logbase = 2, prune = TRUE), results = "fwer")
  toer_high2 <- toer(design = design, n = 15, lambda = lambda_high2,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
    logbase = 2, prune = TRUE), results = "fwer")

  expect_equal(adj_res2$toer, toer_adj2)
  expect_lte(adj_res2$toer, 0.025)
  expect_gt(toer_high2, 0.025)

  # Cases with an additional step after uniroot
  adj_res3 <- adjust_lambda(design = design, alpha = 0.025, n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 4,
    tau = 0, logbase = 2, prune = FALSE), prec_digits = 4)
  lambda_high3 <- adj_res3$lambda - 0.0001
  toer_adj3 <- toer(design = design, n = 15, lambda = adj_res3$lambda,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 4, tau = 0,
    logbase = 2, prune = FALSE), results = "fwer")
  toer_high3 <- toer(design = design, n = 15, lambda = lambda_high3,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 4, tau = 0,
    logbase = 2, prune = FALSE), results = "fwer")

  adj_res4 <- adjust_lambda(design = design, alpha = 0.025, n = 15,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0.2,
    logbase = 2, prune = FALSE), prec_digits = 4)
  lambda_high4 <- adj_res4$lambda - 0.0001
  toer_adj4 <- toer(design = design, n = 15, lambda = adj_res4$lambda,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0.2,
    logbase = 2, prune = FALSE), results = "fwer")
  toer_high4 <- toer(design = design, n = 15, lambda = lambda_high4,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0.2,
    logbase = 2, prune = FALSE), results = "fwer")

  expect_equal(adj_res3$toer, toer_adj3)
  expect_lte(adj_res3$toer, 0.025)
  expect_gt(toer_high3, 0.025)
  expect_equal(adj_res3$toer, toer_adj3)
  expect_lte(adj_res3$toer, 0.025)
  expect_gt(toer_high3, 0.025)

  expect_equal(adj_res4$toer, toer_adj4)
  expect_lte(adj_res4$toer, 0.025)
  expect_gt(toer_high4, 0.025)
  expect_equal(adj_res4$toer, toer_adj4)
  expect_lte(adj_res4$toer, 0.025)
  expect_gt(toer_high4, 0.025)
})

test_that("errors in adjust_lambda work for a single-stage design", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  expect_error(adjust_lambda(design = design, alpha = 1.1, n = 20,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
    logbase = 2, prune = TRUE), prec_digits = 3))
  expect_error(adjust_lambda(design = design, alpha = 0.025, n = c(20, 10, 10),
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
    logbase = 2, prune = TRUE), prec_digits = 3))
})

test_that("adjust_lambda works for a two-stage design", {
  design <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  adj_res1 <- adjust_lambda(design = design, alpha = 0.025, n = 20, n1 = 10,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
      logbase = 2), interim_fun = interim_postpred,
    interim_params = list(futstop = 0.1, effstop = 0.1), prec_digits = 4)
  lambda_high1 <- adj_res1$lambda - 0.0001
  toer_adj1 <- toer(design = design, n = 20, n1 = 10,
    lambda = adj_res1$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1))
  toer_high1 <- toer(design = design, n = 20, n1 = 10,
    lambda = lambda_high1, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1))

  expect_equal(adj_res1$toer, toer_adj1)
  expect_lte(adj_res1$toer, 0.025)
  expect_gt(toer_high1, 0.025)

  adj_res2 <- adjust_lambda(design = design, alpha = 0.025, n = 20, n1 = 10,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1.5, tau = 0,
      logbase = 2), interim_fun = interim_postpred,
    interim_params = list(futstop = 0.1, effstop = 0.1), prec_digits = 4)
  lambda_high2 <- adj_res2$lambda - 0.0001
  toer_adj2 <- toer(design = design, n = 20, n1 = 10,
    lambda = adj_res2$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1.5, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1))
  toer_high2 <- toer(design = design, n = 20, n1 = 10,
    lambda = lambda_high2, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1.5, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1))

  expect_equal(adj_res2$toer, toer_adj2)
  expect_lte(adj_res2$toer, 0.025)
  expect_gt(toer_high2, 0.025)

  # Cases with an additional step after uniroot
  adj_res3 <- adjust_lambda(design = design, alpha = 0.025, n = 20, n1 = 10,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2.5, tau = 0,
      logbase = 2), interim_fun = interim_postpred,
    interim_params = list(futstop = 0.1, effstop = 0.1), prec_digits = 4)
  lambda_high3 <- adj_res3$lambda - 0.0001
  toer_adj3 <- toer(design = design, n = 20, n1 = 10,
    lambda = adj_res3$lambda, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2.5, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1))
  toer_high3 <- toer(design = design, n = 20, n1 = 10,
    lambda = lambda_high3, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 2.5, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1))

  expect_equal(adj_res3$toer, toer_adj3)
  expect_lte(adj_res3$toer, 0.025)
  expect_gt(toer_high3, 0.025)

  adj_res4 <- adjust_lambda(design = design, alpha = 0.025, n = 23, n1 = 10,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
      logbase = 2), interim_fun = interim_postpred,
    interim_params = list(futstop = 0.1, effstop = 0.1), prec_digits = 4)
  lambda_high4 <- adj_res4$lambda - 0.0001
  toer_adj4 <- toer(design = design, n = 23, n1 = 10, lambda = adj_res4$lambda,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0,
      logbase = 2), interim_fun = interim_postpred, interim_params =
      list(futstop = 0.1, effstop = 0.1))
  toer_high4 <- toer(design = design, n = 23, n1 = 10, lambda = lambda_high4,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 2.5, tau = 0,
      logbase = 2), interim_fun = interim_postpred, interim_params =
      list(futstop = 0.1, effstop = 0.1))

  expect_equal(adj_res4$toer, toer_adj4)
  expect_lte(adj_res4$toer, 0.025)
  expect_gt(toer_high4, 0.025)
})

test_that("errors in adjust_lambda work for a two-stage design", {
  design <- setupTwoStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)

  expect_error(adjust_lambda(design = design, alpha = 1.1, n = 20, n1 = 10,
    weight_fun = weights_fujikawa, weight_params = list(epsilon = 1, tau = 0,
      logbase = 2), interim_fun = interim_postpred,
    interim_params = list(futstop = 0.1, effstop = 0.1), prec_digits = 4))
  expect_error(adjust_lambda(design = design, alpha = 0.025, n = c(20, 10, 10),
    n1 = 10, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1), prec_digits = 4))
  expect_error(adjust_lambda(design = design, alpha = 0.025, n = 20,
    n1 = c(10, 5, 5), weight_fun = weights_fujikawa,
    weight_params = list(epsilon = 1, tau = 0, logbase = 2),
    interim_fun = interim_postpred, interim_params = list(futstop = 0.1,
      effstop = 0.1), prec_digits = 4))
})

