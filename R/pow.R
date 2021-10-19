#' @include generics.R
NULL

#' @describeIn pow Power for a single-stage basket design.
setMethod("pow", "OneStageBasket",
  function(design, theta1, n, lambda, weight_fun, tuning_params,
           results = c("ewp", "group"), ...) {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "pwr")
    check_params(n = n, lambda = lambda)

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(tuning_params, design = design,
      n = n, lambda = lambda))

    if (results == "ewp") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "pwr")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat, prob = "pwr")
    }
  })

#' @describeIn pow Power for a two-stage basket design.
setMethod("pow", "TwoStageBasket",
  function(design, theta1, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, tuning_params = list(),
           results = c("ewp", "group"), ...)  {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "pwr")
    check_params(n = n, lambda = lambda)

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(tuning_params, design = design,
      n = n, n1, lambda = lambda))

    if (results == "ewp") {
      reject_prob_ew2(design = design, theta1 = theta1, n = n, n1 = n1,
        lambda = lambda, interim_fun = interim_fun,
        interim_params = interim_params, weight_mat = weight_mat, prob = "pwr")
    } else {
      reject_prob_group2(design = design, theta1 = theta1, n = n, n1 = n1,
        lambda = lambda, interim_fun = interim_fun,
        interim_params = interim_params, weight_mat, prob = "pwr")
    }
  })
