#' @include class.R
NULL

#' Expected number of correct decisions
#'
#' Computes the expected number of correct decisions of a basket trial.
#'
#' @template design
#' @template dotdotdot
#'
#' @details Computes the expected number of correction decisions, i.e. the
#' expected number of actually active baskets that are declared active and
#' actually inactive baskets that are declared inactive.
#'
#' @return A numeric value.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
#' ecd(design = design, p1 = c(0.5, 0.2, 0.2), n = 20, lambda = 0.99,
#' weight_fun = weights_fujikawa)
setGeneric("ecd",
  function(design, ...) standardGeneric("ecd")
)

#' @describeIn ecd Expected number of correction decisions for a single-stage
#'   basket design.
#'
#' @template design
#' @template p1_pow
#' @template n
#' @template lambda
#' @template weights
#' @template globalweights
#' @template dotdotdot
setMethod("ecd", "OneStageBasket",
  function(design, p1 = NULL, n, lambda, weight_fun, weight_params = list(),
           globalweight_fun = NULL, globalweight_params = list(), ...) {
    check_params(n = n, lambda = lambda)
    if (is.null(p1)) p1 <- rep(design@p0, design@k)

    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda, globalweight_fun = globalweight_fun,
      globalweight_params = list(globalweight_params)))

    ecd_calc(design = design, p1 = p1, n = n, lambda = lambda,
      weight_mat = weight_mat, globalweight_fun = globalweight_fun,
      globalweight_params = globalweight_params)
  })

#' @describeIn ecd Expected number of correction decisions for a two-stage
#'   basket design.
#'
#' @template design
#' @template p1_pow
#' @template n
#' @template n1
#' @template lambda
#' @template interim
#' @template weights
#' @template globalweights
#' @template dotdotdot
setMethod("ecd", "TwoStageBasket",
  function(design, p1 = NULL, n, n1, lambda, interim_fun,
           interim_params = list(), weight_fun, weight_params = list(),
           globalweight_fun = NULL, globalweight_params = list(), ...) {
    check_params(n = n, lambda = lambda)
    if (is.null(p1)) p1 <- rep(design@p0, design@k)

    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, n1 = n1, lambda = lambda, globalweight_fun = globalweight_fun,
      globalweight_params = list(globalweight_params)))

    ecd_calc2(design = design, p1 = p1, n = n, n1 = n1, lambda = lambda,
      weight_mat = weight_mat, interim_fun = interim_fun,
      interim_params = interim_params, globalweight_fun = globalweight_fun,
      globalweight_params = globalweight_params)
  })
