#' @include class.R
NULL

#' Adjust Lambda
#'
#' Finds the value for \code{lambda} such that the family wise error
#' rate is protected at level \code{alpha}.
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{adjust_alpha} finds the greatest value with
#' \code{prec_digits} for \code{lambda} which controls the family wise error
#' rate at level \code{alpha} (one-sided). A combination of the uniroot
#' function followed by a grid search is used to finde the correct value
#' for \code{lambda}.
#'
#' @return The greatest value with \code{prec_digits} decimal places for
#' \code{lambda} which controls the family wise error rate at level
#' \code{alpha} (one-sided) and the exact family wise error rate for this
#' value of \code{lambda}.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
#' adjust_lambda(design = design, alpha = 0.025, n = 15,
#'   weight_fun = weights_fujikawa, prec_digits = 4)
setGeneric("adjust_lambda",
  function(design, ...) standardGeneric("adjust_lambda")
)

#' @describeIn adjust_lambda Adjust lambda for a single-stage design.
#'
#' @template design
#' @template alpha
#' @template p1_toer
#' @template n
#' @template weights
#' @template globalweights
#' @template prec_digits
#' @template dotdotdot
setMethod("adjust_lambda", "OneStageBasket",
  function(design, alpha = 0.025, p1 = NULL, n, weight_fun,
           weight_params = list(), globalweight_fun = NULL,
           globalweight_params = list(), prec_digits, ...) {
    p1 <- check_p1(design = design, p1 = p1, type = "toer")
    if (alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1")
    if (length(n) != 1) stop("n must have length 1")

    upper_lim <- 1 - 10^(-prec_digits)
    root_fun <- function(x) toer(design = design, p1 = p1, n = n,
      lambda = x, weight_fun = weight_fun, weight_params = weight_params,
      globalweight_fun = globalweight_fun,
      globalweight_params = globalweight_params, results = "fwer") - alpha

    # Use uniroot to find lambda close to the smallest lambda that protects
    # the significance level at alpha
    uni_root <- stats::uniroot(root_fun, interval = c(0.5, upper_lim),
      tol = 10^(-prec_digits))

    if (uni_root$f.root > 0 ) {
      # If rejection prob is greater than alpha after uniroot, round lambda up
      root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
      val <- root_fun(root)
      if (val > 0) {
        # If rejection prob is still above alpha, increase lambda
        while (val > 0) {
          root <- root + 10^(-prec_digits)
          val <- root_fun(root)
        }
      }
    } else {
      # If rejection prob is less than alpha after uniroot, round lambda down
      root <- floor(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
      val <- root_fun(root)
      if (val > 0) {
        # If the rejection prob is greater now than alpha with the rounded-down
        # lambda, then round lambda up
        root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
        val <- root_fun(root)
      } else {
        # If the rejection prob is still below alpha, decrease lambda
        repeat {
          root_old <- root
          val_old <- val
          root <- root - 10^(-prec_digits)
          val <- root_fun(root)
          if (val > 0) {
            root <- root_old
            val <- val_old
            break
          }
        }
      }
    }
    return(list(
      lambda = root,
      toer = val + alpha
    ))
  })

#' @describeIn adjust_lambda Adjust lambda for a two-stage design.
#'
#' @template design
#' @template alpha
#' @template p1_toer
#' @template n
#' @template n1
#' @template interim
#' @template weights
#' @template prec_digits
#' @template dotdotdot
setMethod("adjust_lambda", "TwoStageBasket",
  function(design, alpha = 0.025, p1 = NULL, n, n1, interim_fun,
           interim_params = list(), weight_fun, weight_params = list(),
           prec_digits, ...) {
    p1 <- check_p1(design = design, p1 = p1, type = "toer")
    if (alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1")
    if (length(n) != 1) stop("n must have length 1")
    if (length(n1) != 1) stop("n1 must have length 1")

    upper_lim <- 1 - 10^(-prec_digits)
    root_fun <- function(x) toer(design = design, p1 = p1, n = n,
      n1 = n1, lambda = x, interim_fun = interim_fun,
      interim_params = interim_params, weight_fun = weight_fun,
      weight_params = weight_params, results = "fwer") - alpha

    # Use uniroot to find lambda close to the smallest lambda that protects
    # the significance level at alpha
    uni_root <- stats::uniroot(root_fun, interval = c(0.5, upper_lim),
      tol = 10^(-prec_digits))

    if (uni_root$f.root > 0 ) {
      # If rejection prob is greater than alpha after uniroot, round lambda up
      root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
      val <- root_fun(root)
      if (val > 0) {
        # If rejection prob is still above alpha, increase lambda
        while (val > 0) {
          root <- root + 10^(-prec_digits)
          val <- root_fun(root)
        }
      }
    } else {
      # If rejection prob is less than alpha after uniroot, round lambda down
      root <- floor(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
      val <- root_fun(root)
      if (val > 0) {
        # If the rejection prob is greater now than alpha with the rounded-down
        # lambda, then round lambda up
        root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
        val <- root_fun(root)
      } else {
        # If the rejection prob is still below alpha, decrease lambda
        repeat {
          root_old <- root
          val_old <- val
          root <- root - 10^(-prec_digits)
          val <- root_fun(root)
          if (val > 0) {
            root <- root_old
            val <- val_old
            break
          }
        }
      }
    }
    return(list(
      lambda = root,
      toer = val + alpha
    ))
  })


