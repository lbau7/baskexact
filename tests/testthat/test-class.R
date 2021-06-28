test_that("setupBasket works", {
  design <- setupOneStageBasket(k = 3, theta0 = 0.2, theta1 = c(0.2, 0.5, 0.5))

  expect_s4_class(design, "OneStageBasket")
})

test_that("validity method works", {
  expect_error(setupOneStageBasket(k = c(1, 2), theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupOneStageBasket(k = 2.3, theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupOneStageBasket(k = -3, theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupOneStageBasket(k = 3, shape1 = c(1, 2, 3), theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupOneStageBasket(k = 3, shape2 = c(1, 2, 3), theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupOneStageBasket(k = 3, shape1 = 0, theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupOneStageBasket(k = 3, shape2 = 0, theta0 = 0.2,
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupOneStageBasket(k = 3, theta0 = c(0.2, 0.3, 0.3),
    theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupBasket(k = 3, theta0 = 0.2, theta1 = c(0.2, 0.5)))
  expect_error(setupBasket(k = 3, theta0 = -0.1,theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupBasket(k = 3, theta0 = 1.1,theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupBasket(k = 3, theta0 = 0.2,theta1 = c(0.2, 0.5, 1.1)))
  expect_error(setupBasket(k = 3, theta0 = 0.2,theta1 = c(0.2, 0.5, -0.1)))
  expect_error(setupBasket(k = 3, theta0 = 0.3,theta1 = c(0.2, 0.5, 0.5)))
  expect_error(setupBasket(k = 3, theta0 = 0.2,theta1 = c(0.2, 0.2, 0.2)))
})
