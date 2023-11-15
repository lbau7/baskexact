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
#' design <- setupTwoStageBasket(k = 3, p0 = 0.2)
#' ess(design, n = 20, n1 = 10, lambda = 0.99, weight_fun = weights_fujikawa,
#'   interim_fun = interim_postpred)
setGeneric("ess",
  function(design, ...) standardGeneric("ess")
)

#' @describeIn ess Expected sample size for two-stage basket design.
#'
#' @template design
#' @template p1_toer
#' @template n
#' @template n1
#' @template lambda
#' @template interim
#' @template weights
#' @template globalweights
#' @template dotdotdot
setMethod("ess", "TwoStageBasket",
  function(design, p1 = NULL, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, weight_params = list(),
           globalweight_fun = NULL, globalweight_params = list(), ...)  {
    p1 <- check_p1(design = design, p1 = p1, type = "ess")

    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, n1 = n1, lambda = lambda, globalweight_fun = globalweight_fun,
      globalweight_params = globalweight_params))

    events_int <- arrangements::permutations(x = 0:n1, k = design@k,
      replace = TRUE)
    int_fun <- function(x) do.call(interim_fun, args = c(interim_params,
      design = design, n = n, n1 = n1, r1 = list(x), lambda = lambda,
      weight_mat = list(weight_mat), globalweight_fun = globalweight_fun,
      globalweight_params = list(globalweight_params)))
    res_int <- t(apply(events_int, 1, int_fun))
    sampsize <- ifelse(res_int == 0, n, n1)

    probs <- apply(events_int, 1, function(x)
      get_prob(n = n1, r = x, p = p1))
    colSums(apply(sampsize, 2, function(x) x * probs))
  })
