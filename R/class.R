# Parent-class Basket for internal use only.
setClass("Basket",
  slots = c(
    k = "numeric",
    shape1 = "numeric",
    shape2 = "numeric",
    p0 = "numeric"
  ))

#' Class OneStageBasket
#'
#' OneStageBasket is an S4 class. An object of this class contains the most
#' important design features of a single-stage basket trial.
#'
#' @slot k The number of baskets.
#' @slot shape1 First common shape parameter of the beta prior.
#' @slot shape2 Second common shape parameter of the beta prior.
#' @slot p0 A common probability under the null hypothesis.
#'
#' @details
#' This class implements a single-stage basket trial based on the power prior
#' design or the design proposed by Fujikawa et al. In these designs,
#' information is borrowed between baskets by calculating weights that
#' reflect the similarity between the baskets (and optionally the overall
#' heterogeneity). Posterior distributions for each basket are beta
#' distributions where the parameters are found by adding weighted sums
#' of the observed responses and non-responses in each basket to the
#' prior parameters (or in case of Fujikawa's design by calculating
#' weighted sums of the individual posterior distributions).
#'
#' Currently only common prior distributions and a common null
#' hypothesis are supported.
#'
#' @references
#' Baumann, L., Sauer, L., & Kieser, M. (2024). A basket trial design based on
#' power priors. arXiv:2309.06988.
#'
#' Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
#' A Bayesian basket trial design that borrows information across strata based
#' on the similarity between the posterior distributions of the response
#' probability. Biometrical Journal, 62(2), 330-338.
#'
#' @aliases OneStageBasket
setClass("OneStageBasket", contains = "Basket")

#' Class TwoStageBasket
#'
#' TwoStageBasket is an S4 class. An object of this class contains the most
#' important design features of a two-stage basket trial.
#'
#' @slot k The number of baskets.
#' @slot shape1 First common shape parameter of the beta prior.
#' @slot shape2 Second common shape parameter of the beta prior.
#' @slot p0 A common probability under the null hypothesis.
#'
#' @details
#' This class implements a two-stage basket trial based on the power prior
#' design or the design proposed by Fujikawa et al. In these designs,
#' information is borrowed between baskets by calculating weights that
#' reflect the similarity between the baskets (and optionally the overall
#' heterogeneity). Posterior distributions for each basket are beta
#' distributions where the parameters are found by adding weighted sums
#' of the observed responses and non-responses in each basket to the
#' prior parameters (or in case of Fujikawa's design by calculating
#' weighted sums of the individual posterior distributions).
#'
#' @references
#' Baumann, L., Sauer, L., & Kieser, M. (2024). A basket trial design based on
#' power priors. arXiv:2309.06988.
#'
#' Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
#' A Bayesian basket trial design that borrows information across strata based
#' on the similarity between the posterior distributions of the response
#' probability. Biometrical Journal, 62(2), 330-338.
#'
#' @aliases TwoStageBasket
setClass("TwoStageBasket", contains = "Basket")

#' Setup OneStageBasket
#'
#' Creates an object of class \code{\link{OneStageBasket}}.
#'
#' @param k The number of baskets.
#' @param shape1 First common shape parameter of the beta prior.
#' @param shape2 Second common shape parameter of the beta prior.
#' @param p0 A common probability under the null hypothesis.
#'
#' @details A \code{\link{OneStageBasket}} object contains the most important
#' design features of a basket trial. Currently only common prior distributions
#' and a common null hypothesis are supported.
#'
#' @return An S4 object of class \code{\link{OneStageBasket}}.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
setupOneStageBasket <- function(k, shape1 = 1, shape2 = 1, p0) {
  methods::new("OneStageBasket",
    k = k,
    shape1 = shape1,
    shape2 = shape2,
    p0 = p0
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
  } else if (length(object@p0) != 1) {
    "p0 must have length 1"
  } else if ((object@p0 <= 0) | (object@p0 >= 1)) {
    "p0 must be between 0 and 1"
  } else {
    TRUE
  }})

#' Setup TwoStageBasket
#'
#' Creates an object of class \code{\link{TwoStageBasket}}.
#'
#' @param k The number of baskets.
#' @param shape1 First common shape parameter of the beta prior.
#' @param shape2 Second common shape parameter of the beta prior.
#' @param p0 A common probability under the null hypothesis.
#'
#' @details A \code{\link{TwoStageBasket}} object contains the most important
#' design features of a basket trial. Currently only common prior distributions
#' and a common null hypothesis are supported.
#'
#' @return An S4 object of class \code{\link{TwoStageBasket}}.
#' @export
#'
#' @examples
#' design <- setupTwoStageBasket(k = 3, p0 = 0.2)
setupTwoStageBasket <- function(k, shape1 = 1, shape2 = 1, p0) {
  methods::new("TwoStageBasket",
    k = k,
    shape1 = shape1,
    shape2 = shape2,
    p0 = p0
  )}

setValidity("TwoStageBasket", function(object) {
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
  } else if (length(object@p0) != 1) {
    "p0 must have length 1"
  } else if ((object@p0 <= 0) | (object@p0 >= 1)) {
    "p0 must be between 0 and 1"
  } else {
    TRUE
  }})

