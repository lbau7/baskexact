#' Class OneStageBasket
#'
#' OneStageBasket is an S4 class. An object of this class contains the most
#' important design features of a single-stage basket trial.
#'
#' @slot k The number of baskets.
#' @slot shape1 First common shape parameter of the beta prior.
#' @slot shape2 Second common shape parameter of the beta prior.
#' @slot theta0 A common probability under the null hypothesis.
#'
#' @details
#' This class implements a single-stage basket trial based on the design
#' proposed by Fujikawa et al. In this design, at first separate posterior
#' distributions are calculated for each basket based on a beta-binomial model.
#' Information is then borrowed between baskets by calculating weights that
#' reflect the similarity between the basket and computing posterior
#' distributions for each basket where the parameters of the beta posterior are
#' calculated as a weighted sum of the individual posterior distributions.
#' The weight between two baskets i and j is found as (1 - JSD(i, j))^epsilon
#' where JSD(i, j) is the Jensen-Shannon divergence between basket i and j.
#' A small value of epsilon results in stronger borrowing also across baskets
#' with heterogenous results. If epsilon is large then information is only
#' borrowed between baskets with similar results. If a weight is smaller than
#' tau it is set to 0, which results in no borrowing. If for a basket
#' the posterior probability that \eqn{\theta} > theta0 is greater than
#' lambda, then the null hypothesis is rejected.
#'
#' Currently only common prior distributions and a common null
#' hypothesis are supported.
#'
#' @references Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
#' A Bayesian basket trial design that borrows information across strata based
#' on the similarity between the posterior distributions of the response
#' probability. Biometrical Journal, 62(2), 330-338.
#'
#' @aliases OneStageBasket
setClass("OneStageBasket",
  slots = c(
    k = "numeric",
    shape1 = "numeric",
    shape2 = "numeric",
    theta0 = "numeric"
  ))

#' Setup OneStageBasket
#'
#' Creates an object of class \code{\link{OneStageBasket}}.
#'
#' @param k The number of baskets.
#' @param shape1 First common shape parameter of the beta prior.
#' @param shape2 Second common shape parameter of the beta prior.
#' @param theta0 A common probability under the null hypothesis.
#'
#' @details A \code{\link{OneStageBasket}} object contains the most important
#' design features of a basket trial. Currently only common prior distributions
#' and a common null hypothesis are supported.
#'
#' @return An S4 object of class \code{\link{OneStageBasket}}.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
setupOneStageBasket <- function(k, shape1 = 1, shape2 = 1, theta0) {
  methods::new("OneStageBasket",
    k = k,
    shape1 = shape1,
    shape2 = shape2,
    theta0 = theta0
  )}

setValidity("OneStageBasket", function(object) {
  if (length(object@k) != 1) {
    "k must have length 1"
  } else if ((object@k %% 1 != 0) | (object@k <= 0)) {
    "k must be a positive integer"
  } else if (length(object@shape1) != 1) {
    "shape1 must have length 1"
  } else if (length(object@shape2) != 1) {
    "shape2 must have length 1"
  } else if (object@shape1 <= 0) {
    "shape1 must be positive"
  } else if (object@shape2 <= 0) {
    "shape2 must be positive"
  } else if (length(object@theta0) != 1) {
    "theta0 must have length 1"
  } else if ((object@theta0 <= 0) | (object@theta0 >= 1)) {
    "theta0 must be between 0 and 1"
  } else {
    TRUE
  }})
