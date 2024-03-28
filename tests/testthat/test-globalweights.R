test_that("globalweights_diff works", {
  res1 <- globalweights_diff(n = 20, r = c(0, 10, 20), eps_global = 1)
  res2 <- globalweights_diff(n = 40, r = c(0, 10, 20, 30, 40), eps_global = 1)

  # Weight is 0 when heterogeneity is maximal
  expect_equal(res1, 0)
  expect_equal(res2, 0)

  res3 <- globalweights_diff(n = 20, r = c(5, 5, 5), eps_global = 1)
  res4 <- globalweights_diff(n = 20, r = c(5, 5, 5, 5), eps_global = 1)

  # Weight is 1 when there is no heterogeneity
  expect_equal(res3, 1)
  expect_equal(res4, 1)

  res5 <- globalweights_diff(n = 20, r = c(1, 3, 7, 9), eps_global = 1)
  res6 <- globalweights_diff(n = 20, r = c(9, 7, 3, 1), eps_global = 1)

  # Weight is permutation invariant
  expect_equal(res5, res6)

  res7 <- globalweights_diff(n = 20, r = c(1, 2, 3, 4), eps_global = 1)
  res8 <- globalweights_diff(n = 20, r = c(1, 2, 3, 5), eps_global = 1)

  # Weight is smaller when heterogeneity is larger
  expect_gt(res7, res8)

  res9 <- globalweights_diff(n = 20, r = c(1, 2, 3, 4), eps_global = 1)
  res10 <- globalweights_diff(n = 20, r = c(1, 2, 3, 4), eps_global = 2)

  # Weight is smaller when epsilon is larger
  expect_gt(res9, res10)
})

test_that("globalweights_fix works", {
  res <- globalweights_fix(n = 20, r = c(1, 3, 5), w = 0.3)
  expect_equal(res, 0.3)
})

test_that("globalweights_maxdiff works", {
  res1 <- globalweights_maxdiff(n = 20, r = c(5, 10, 15), eps_global = 1)
  res2 <- globalweights_maxdiff(n = 20, r = c(5, 10, 15), eps_global = 1,
    w = 0.5)

  expect_equal(res1, 0.5)
  expect_equal(res2, 0.25)
})
