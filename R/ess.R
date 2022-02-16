#' @include class.R
NULL

#' Expected Sample Size
#'
#' Computes the expected sample size of a two-stage basket trial.
#'
#' @template design
#' @template dotdotdot
#'
#' @export
#'
#' @examples
#' design <- setupTwoStageBasket(k = 3, theta0 = 0.2)
#' ess(design, n = 20, n1 = 10, lambda = 0.99, weight_fun = weights_fujikawa,
#'   interim_fun = interim_postpred)
setGeneric("ess",
  function(design, ...) standardGeneric("ess")
)

#' @describeIn ess Expected sample size for two-stage basket design.
#'
#' @template design
#' @template theta1_toer
#' @template n
#' @template n1
#' @template lambda
#' @template interim
#' @template weights
#' @template dotdotdot
setMethod("ess", "TwoStageBasket",
  function(design, theta1 = NULL, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, weight_params = list(), ...)  {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "ess")

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
