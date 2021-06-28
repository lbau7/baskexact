#' Class OneStageBasket
#'
#' OneStageBasket is an S4 class. An object of this class conatains the most
#' important design features of a single-stage basket trial.
#'
#' @slot k The number of baskets.
#' @slot shape1 First common shape parameter of the beta prior.
#' @slot shape2 Second common shape parameter of the beta prior.
#' @slot theta0 A common probability under the null hypothesis.
#' @slot theta1 A probability under the alternative hypothesis for each basket.
#'
#' @details Currently only common prior distributions and a common null
#' hypothesis are supported.
#'
#' @aliases OneStageBasket
setClass("OneStageBasket",
  slots = c(
    k = "numeric",
    shape1 = "numeric",
    shape2 = "numeric",
    theta0 = "numeric",
    theta1 = "numeric"
  ))

#' Setup OneStageBasket
#'
#' Creates an object of class \code{OneStageBasket}.
#'
#' @param k The number of baskets.
#' @param shape1 First common shape parameter of the beta prior.
#' @param shape2 Second common shape parameter of the beta prior.
#' @param theta0 A common probability under the null hypothesis.
#' @param theta1 A probability under the alternative hypothesis for each basket.
#'
#' @details A \code{OneStageBasket} object contains the most important design
#' features of a basket trial. Currently only common prior distributions and a
#' common null hypothesis are supported.
#'
#' @return An S4 object of class \code{OneStageBasket}.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2, theta1 = c(0.2, 0.5, 0.5))
setupOneStageBasket <- function(k, shape1 = 1, shape2 = 1, theta0, theta1) {
  methods::new("OneStageBasket",
    k = k,
    shape1 = shape1,
    shape2 = shape2,
    theta0 = theta0,
    theta1 = theta1
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
  } else if (length(object@theta1) != object@k) {
    "theta1 must have length k"
  } else if ((object@theta0 <= 0) | (object@theta0 >= 1)) {
    "theta0 must be between 0 and 1"
  } else if (any(object@theta1 <= 0) | any(object@theta1 >= 1)) {
    "theta1 must be between 0 and 1"
  } else if (any(object@theta0 > object@theta1)) {
    "theta1 must be greater than or equal to theta0"
  } else {
    TRUE
  }})
