test_that("setupBasket works", {
  design <- setupOneStageBasket(k = 3, p0 = 0.2)

  expect_s4_class(design, "OneStageBasket")
})

test_that("validity method works", {
  expect_error(setupOneStageBasket(k = c(1, 2), p0 = 0.2))
  expect_error(setupOneStageBasket(k = 2.3, p0 = 0.2))
  expect_error(setupOneStageBasket(k = -3, p0 = 0.2))
  expect_error(setupOneStageBasket(k = 3, shape1 = c(1, 2, 3), p0 = 0.2))
  expect_error(setupOneStageBasket(k = 3, shape2 = c(1, 2, 3), p0 = 0.2))
  expect_error(setupOneStageBasket(k = 3, shape1 = 0, p0 = 0.2))
  expect_error(setupOneStageBasket(k = 3, shape2 = 0, p0 = 0.2))
  expect_error(setupOneStageBasket(k = 3, p0 = c(0.2, 0.3, 0.3)))
  expect_error(setupBasket(k = 3, p0 = -0.1))
  expect_error(setupBasket(k = 3, p0 = 1.1))
})
