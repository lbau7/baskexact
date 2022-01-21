setGeneric("ecd",
  function(design, ...) standardGeneric("ecd")
)

setMethod("ecd", "OneStageBasket",
  function(design, theta1 = NULL, n, lambda, weight_fun, weight_params = list(),
           ...) {
    check_params(n = n, lambda = lambda)

    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda))

    ecd_calc(design = design, theta1 = theta1, n = n, lambda = lambda,
      weight_mat = weight_mat)
  })
