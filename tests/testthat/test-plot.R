test_that("plot_weights_work", {
  plot1 <- plot_weights(
    design = setupOneStageBasket(k = 3, p0 = 0.2),
    n = 20, r1 = 10,
    weight_fun = weights_cpp,
    weight_params = list(a = 1, b = 1)
  )

  expect_s3_class(plot1, "ggplot")
  expect_true(all(plot1$data$r == 0:20))

  plot2 <- plot_weights(
    design = setupOneStageBasket(k = 3, p0 = 0.2),
    n = 20, r1 = 10,
    weight_fun = weights_cpp,
    weight_params = list(a = c(1, 2, 3), b = 1)
  )

  expect_s3_class(plot2, "ggplot")
  expect_true(all(unique(plot2$param) %in% 1:3))

  plot3 <- plot_weights(
    design = setupOneStageBasket(k = 3, p0 = 0.2),
    n = 20, r1 = 10,
    weight_fun = weights_cpp,
    weight_params = list(a = c(1, 2, 3), b = c(1, 2, 3))
  )

  expect_s3_class(plot3, "ggplot")
  expect_true(all(unique(plot3$param1) %in% 1:3))
  expect_true(all(unique(plot3$param2) %in% 1:3))

  expect_error(plot_weights(n = 20, r1 = 10, weight_fun = weights_fujikawa,
    weight_params = list(epsilon = c(1, 2, 3), tau = c(0, 0.5),
      logbase = c(2, 3))))
})
