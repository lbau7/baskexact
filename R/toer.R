#' @include generics.R
NULL

#' @describeIn toer Type 1 error rate for a single-stage basket design.
setMethod("toer", "OneStageBasket",
  function(design, theta1 = NULL, n, lambda, weight_fun, tuning_params,
           results = c("fwer", "group"), ...) {
    results <- match.arg(results)
    weight_mat <- do.call(weights, args = c(tuning_params, design = design,
      n = n, lambda = lambda))

    if (results == "fwer") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "toer")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat, prob = "toer")
    }
  })
