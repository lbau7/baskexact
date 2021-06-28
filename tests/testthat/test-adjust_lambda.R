test_that("adjust_lambda works", {
  design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = c(0.2, 0.2, 0.2))

  # Without Pruning
  adj_res1 <- adjust_lambda(design = design, alpha = 0.025, n = 15,
    epsilon = 1, tau = 0, logbase = 2, prune = FALSE, prec_digits = 4)
  lambda_high1 <- adj_res1$lambda - 0.0001
  toer_adj1 = toer(design = design, n = 15, lambda = adj_res1$lambda,
    epsilon = 1, tau = 0, logbase = 2, prune = FALSE, results = "fwer")
  toer_high1 <- toer(design = design, n = 15, lambda = lambda_high1,
    epsilon = 1, tau = 0, logbase = 2, prune = FALSE, results = "fwer")

  # With Pruning
  adj_res2 <- adjust_lambda(design = design, alpha = 0.025, n = 15,
    epsilon = 1, tau = 0, logbase = 2, prune = TRUE, prec_digits = 4)
  lambda_high2 <- adj_res2$lambda - 0.0001
  toer_adj2 = toer(design = design, n = 15, lambda = adj_res2$lambda,
    epsilon = 1, tau = 0, logbase = 2, prune = TRUE, results = "fwer")
  toer_high2 <- toer(design = design, n = 15, lambda = lambda_high2,
    epsilon = 1, tau = 0, logbase = 2, prune = TRUE, results = "fwer")

  expect_equal(adj_res1$toer, toer_adj1)
  expect_lte(adj_res1$toer, 0.025)
  expect_gt(toer_high1, 0.025)
  expect_equal(adj_res2$toer, toer_adj2)
  expect_lte(adj_res2$toer, 0.025)
  expect_gt(toer_high2, 0.025)
})
