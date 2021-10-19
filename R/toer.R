#' @include generics.R
NULL

#' @describeIn toer Type 1 error rate for a single-stage basket design.
#' @template theta1_toer
#' @template n
#' @template weights
#' @param results Whether only the family wise error rate (option \code{fwer})
#'   or also the rejection probabilities per group (option \code{group}) should
#'   be returned.
#' @template dotdotdot
setMethod("toer", "OneStageBasket",
  function(design, theta1 = NULL, n, lambda, weight_fun, tuning_params,
           results = c("fwer", "group"), ...) {
    check_params(n = n, lambda = lambda)
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "toer")

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(tuning_params, design = design,
      n = n, lambda = lambda))

    if (results == "fwer") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "toer")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat, prob = "toer")
    }
  })

#' @describeIn toer Type 1 error rate for two-stage basket design.
setMethod("toer", "TwoStageBasket",
  function(design, theta1 = NULL, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, tuning_params = list(),
           results = c("fwer", "group"), ...)  {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "toer")

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(tuning_params, design = design,
      n = n, n1, lambda = lambda))

    if (results == "fwer") {
      reject_prob_ew2(design = design, theta1 = theta1, n = n, n1 = n1,
        lambda = lambda, interim_fun = interim_fun,
        interim_params = interim_params, weight_mat = weight_mat, prob = "toer")
    } else {
      reject_prob_group2(design = design, theta1 = theta1, n = n, n1 = n1,
        lambda = lambda, interim_fun = interim_fun,
        interim_params = interim_params, weight_mat, prob = "toer")
    }
  })
