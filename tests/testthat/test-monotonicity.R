test_that("check_mon_within works", {
  # Without Pruning
  design <- setupBasket(k = 4, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = rep(0.5, 4))
  r1 <- check_mon_within(design = design, n = 24, lambda = 0.99, epsilon = 0.5,
    tau = 0, logbase = 2, prune = FALSE, details = TRUE)$Events

  r2 <- r3 <- r4 <- r5 <- r1
  r2[1] <- r2[1] - 1
  r3[2] <- r3[2] + 1
  r4[3] <- r4[3] - 1
  r5[4] <- r5[4] + 1

  weights1 <- get_weights(design = design, n = 24, epsilon = 0.5, tau = 0,
    logbase = 2)
  res1 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r1,
    weight_mat = weights1)
  res2 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r2,
    weight_mat = weights1)
  res3 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r3,
    weight_mat = weights1)
  res4 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r4,
    weight_mat = weights1)
  res5 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r5,
    weight_mat = weights1)

  monres1 <- check_mon_within(design = design, n = 24, lambda = 0.99,
    epsilon = 0.5, tau = 0, logbase = 2, prune = FALSE, details = FALSE)

  # With Pruning
  r6 <- check_mon_within(design = design, n = 24, lambda = 0.99, epsilon = 7,
    tau = 0, logbase = 2, prune = TRUE, details = TRUE)$Events[1, ]
  r7 <- r8 <- r9 <- r6
  r7[3] <- r7[3] + 1
  r8[4] <- r8[4] - 1
  r9[4] <- r9[4] + 1

  weights2 <- get_weights(design = design, n = 24, epsilon = 7, tau = 0,
    logbase = 2)
  crit_pool <- get_crit_pool(design = design, n = 24, lambda = 0.99)
  weights2 <- prune_weights(weights2, cut = crit_pool)

  res6 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r6,
    weight_mat = weights2)
  res7 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r7,
    weight_mat = weights2)
  res8 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r8,
    weight_mat = weights2)
  res9 <- bskt_final(design = design, n = 24, lambda = 0.99, r = r9,
    weight_mat = weights2)

  monres2 <- check_mon_within(design = design, n = 24, lambda = 0.99,
    epsilon = 7, tau = 0, logbase = 2, prune = TRUE, details = FALSE)

  expect_true(any(res1 != cummax(res1)))
  expect_false(any(res2 != cummax(res2)))
  expect_false(any(res3 != cummax(res3)))
  expect_false(any(res4 != cummax(res4)))
  expect_false(any(res5 != cummax(res5)))
  expect_false(monres1)
  expect_true(any(res6 != cummax(res6)))
  expect_false(any(res7 != cummax(res7)))
  expect_false(any(res8 != cummax(res8)))
  expect_false(any(res9 != cummax(res9)))
  expect_false(monres2)
})

test_that("check_mon_between works", {
  # Without Pruning
  design <- setupBasket(k = 4, shape1 = 1, shape2 = 1, theta0 = 0.2,
    theta1 = rep(0.5, 4))
  ev1 <- check_mon_between(design = design, n = 15, lambda = 0.99,
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE, details = TRUE)[[1]]

  weights1 <- get_weights(design = design, n = 15, epsilon = 2, tau = 0,
    logbase = 2)
  res1 <- bskt_final(design = design, n = 15, lambda = 0.99,
    r = ev1$Events[1, ], weight_mat = weights1)
  res2 <- bskt_final(design = design, n = 15, lambda = 0.99,
    r = ev1$Events[2, ], weight_mat = weights1)

  res_nodet1 <- check_mon_between(design = design, n = 15, lambda = 0.99,
    epsilon = 2, tau = 0, logbase = 2, prune = FALSE, details = FALSE)

  # With Pruning
  res_nodet2 <- check_mon_between(design = design, n = 15, lambda = 0.99,
    epsilon = 2, tau = 0, logbase = 2, prune = TRUE, details = FALSE)

  crit_pool <- get_crit_pool(design = design, n = 15, lambda = 0.99)
  weights2 <- prune_weights(weights1, cut = crit_pool)
  res3 <- bskt_final(design = design, n = 15, lambda = 0.99,
    r = ev1$Events[1, ], weight_mat = weights2)
  res4 <- bskt_final(design = design, n = 15, lambda = 0.99,
    r = ev1$Events[2, ], weight_mat = weights2)

  expect_equal(res1, ev1$Results[1, ])
  expect_true(any(res1 == 1))
  expect_equal(res2, ev1$Results[2, ])
  expect_true(all(res2 == 0))
  expect_false(res_nodet1)
  expect_equal(res3, res4)
  expect_true(res_nodet2)
})
