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
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' ecd(design = design, theta1 = c(0.5, 0.2, 0.2), n = 20, lambda = 0.99,
#' weight_fun = weights_fujikawa)
setGeneric("ecd",
  function(design, ...) standardGeneric("ecd")
)

#' @describeIn ecd Expected number of correction decisions for a single-stage
#'   basket design.
#'
#' @template design
#' @template theta1_pow
#' @template n
#' @template lambda
#' @template weights
#' @template dotdotdot
setMethod("ecd", "OneStageBasket",
  function(design, theta1 = NULL, n, lambda, weight_fun, weight_params = list(),
           ...) {
    check_params(n = n, lambda = lambda)

    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda))

    ecd_calc(design = design, theta1 = theta1, n = n, lambda = lambda,
      weight_mat = weight_mat)
  })
