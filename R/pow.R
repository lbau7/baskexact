#' @include class.R
NULL

#' Power
#'
#' Computes the exact power for a basket trial.
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{pow} computes the exact experimentwise power and the
#' exact rejection probabilities per group. The experimentwise power
#' is the probability to reject at least one null hypothesis for a basket with
#' theta1 > theta0. The rejection probabilities correspond to the type 1 error
#' rate for baskets with theta1 = theta 0 and to the power for baskets with
#' theta1 > theta 0.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
#'
#' This method is implemented for the class \code{\link{OneStageBasket}}.
#'
#' @return If \code{results = "ewp"} then the experimentwise power is
#' returned as a numeric value. If \code{results = "group"} then a list with
#' the rejection probabilities per group and the experimentwise power
#' is returned. For baskets with theta1 = theta0 the rejection probabilities
#' corresponds to the type 1 error rate, for baskets with theta1 > theta0 the
#' rejection probabilities corresponds to the power.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' pow(design, theta1 = c(0.2, 0.5, 0.5), n = 15, lambda = 0.99,
#'   weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0))
setGeneric("pow",
  function(design, ...) standardGeneric("pow")
)

#' @describeIn pow Power for a single-stage basket design.
#'
#' @template design
#' @template theta1_pow
#' @template n
#' @template lambda
#' @template weights
#' @template results_pow
setMethod("pow", "OneStageBasket",
  function(design, theta1, n, lambda, weight_fun, weight_params = list(),
           globalweight_fun = NULL, globalweight_params = list(),
           results = c("ewp", "group"), ...) {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "pwr")
    check_params(n = n, lambda = lambda)

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda))

    if (results == "ewp") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, globalweight_fun = globalweight_fun,
        globalweight_params = globalweight_params, prob = "pwr")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat,
        globalweight_fun = globalweight_fun,
        globalweight_params = globalweight_params, prob = "pwr")
    }
  })

#' @describeIn pow Power for a two-stage basket design.
#'
#' @template design
#' @template theta1_pow
#' @template n
#' @template n1
#' @template lambda
#' @template interim
#' @template weights
#' @template results_pow
setMethod("pow", "TwoStageBasket",
  function(design, theta1, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, weight_params = list(),
           results = c("ewp", "group"), ...)  {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "pwr")
    check_params(n = n, lambda = lambda)

    results <- match.arg(results)
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
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
