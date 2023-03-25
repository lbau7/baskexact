#' @include class.R
NULL

#' Type 1 Error Rate
#'
#' Computes the exact family wise type 1 error rate of a basket trial .
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{toer} computes the exact family wise type 1 error rate and the
#' exact rejection probabilities per group. The family wise type 1 error rate
#' is the probability to reject at least one null hypothesis for a basket with
#' theta1 = theta0. If all theta1 > theta0 then the family wise type 1 error
#' rate under the global null hypothesis is computed. The rejection
#' probabilities correspond to the type 1 error rate for baskets with theta1 =
#' theta 0 and to the power for baskets with theta1 > theta 0.
#'
#' @return If \code{results = "fwer"} then the family wise type 1 error rate is
#' returned as a numeric value. If \code{results = "group"} then a list with
#' the rejection probabilities per group and the family wise type 1 error rate
#' is returned. If all theta1 > theta0 then the family wise type 1 error rate
#' is calculated under the global null hypothesis. For baskets with theta1 =
#' theta0 the rejection probabilities corresponds to the type 1 error rate, for
#' baskets with theta1 > theta0 the rejection probabilities corresponds to the
#' power.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_fujikawa)
setGeneric("toer",
  function(design, ...) standardGeneric("toer")
)

#' @describeIn toer Type 1 error rate for a single-stage basket design.
#'
#' @template design
#' @template theta1_toer
#' @template n
#' @template lambda
#' @template weights
#' @template results_toer
#' @template dotdotdot
setMethod("toer", "OneStageBasket",
  function(design, theta1 = NULL, n, lambda, weight_fun, weight_params = list(),
           globalweight_fun = NULL, globalweight_params,
           results = c("fwer", "group"), ...) {
    check_params(n = n, lambda = lambda)
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "toer")

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda))

    if (results == "fwer") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, globalweight_fun,
        globalweight_params = globalweight_params, prob = "toer")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat,
        globalweight_fun = globalweight_fun,
        globalweight_params = globalweight_params, prob = "toer")
    }
  })

#' @describeIn toer Type 1 error rate for two-stage basket design.
#'
#' @template design
#' @template theta1_toer
#' @template n
#' @template n1
#' @template lambda
#' @template interim
#' @template weights
#' @template results_toer
#' @template dotdotdot
setMethod("toer", "TwoStageBasket",
  function(design, theta1 = NULL, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, weight_params = list(),
           results = c("fwer", "group"), ...)  {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "toer")

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
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
