test_that("globalweights_diff works", {
  res1 <- globalweights_diff(n = 20, r = c(0, 10, 20), eps_global = 1)
  res2 <- globalweights_diff(n = 20, r = c(5, 5, 5), eps_global = 1)
  res3 <- globalweights_diff(n = 20, r = c(5, 5, 5), eps_global = 1, w = 0.5)

  expect_equal(res1, 0)
  expect_equal(res2, 1)
  expect_equal(res3, 0.5)
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
