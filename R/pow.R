#' @include generics.R
NULL

#' @describeIn pow Power for a single-stage basket design.
setMethod("pow", "OneStageBasket",
  function(design, theta1, n, lambda, weight_fun, tuning_params,
           results = c("ewp", "group"), ...) {
    results <- match.arg(results)
    weight_mat <- do.call(weights, args = c(tuning_params, design = design,
      n = n, lambda = lambda))

    if (results == "ewp") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "pwr")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat, prob = "pwr")
    }
  })
