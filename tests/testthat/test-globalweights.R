test_that("global_weights_diff works", {
  res1 <- globalweights_diff(n = 20, r = c(0, 10, 20), epsilon = 1)
  res2 <- globalweights_diff(n = 20, r = c(5, 5, 5), epsilon = 1)

  expect_equal(res1, 0)
  expect_equal(res2, 1)
})
