#' @include class.R
NULL

#' Expected Sample Size
setGeneric("ess",
  function(design, ...) standardGeneric("ess")
)

#' @describeIn ess Expected sample size for two-stage basket design.
setMethod("ess", "TwoStageBasket",
  function(design, theta1 = NULL, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, weight_params = list(), ...)  {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "pwr")

    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, n1, lambda = lambda))

    events_int <- arrangements::permutations(x = 0:n1, k = design@k,
      replace = TRUE)
    int_fun <- function(x) do.call(interim_fun, args = c(interim_params,
      design = design, n = n, n1 = n1, r1 = list(x), lambda = lambda,
      weight_mat = list(weight_mat)))
    res_int <- t(apply(events_int, 1, int_fun))
    sampsize <- ifelse(res_int == 0, n, n1)

    probs <- apply(events_int, 1, function(x)
      get_prob(n = n1, r = x, theta = theta1))
    colSums(apply(sampsize, 2, function(x) x * probs))
  })
