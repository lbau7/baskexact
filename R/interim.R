#' @include class.R
NULL

#' Interim analysis based on the posterior predictive probability
#'
#' Conducts an interim analysis based on the posterior predictive probability.
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{interim_postpred} conducts an interim analysis with possible
#' stop for efficacy and futility based on the posterior predictive probability.
#' If the posterior predictive probability is less than \code{prob_fustop} the
#' basket is  stopped for futility, if the posterior predictive probability is
#' greater than \code{prob_effstop} the basket is stopped for efficacy. If
#' \code{prob_fustop = 0} or \code{prob_effstop = 1} then no futility-stop and
#' no efficacy stop is possible, respectively.
#'
#' The function is generally not called by the user but passed to another
#' function such as \code{\link{toer}} and \code{\link{pow}} to specify which
#' interim analysis is conducted.
#'
#' @return A vector with a length equal to the number of baskets with
#' elements -1, 0 or 1 where -1 means stop for futility, 0 means continuation
#' and 1 means stop for efficacy.
#' @export
#'
#' @examples
#' design <- setupTwoStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 20, n1 = 10, lambda = 0.99, interim_fun = interim_postpred,
#'   weight_fun = weights_fujikawa)
setGeneric("interim_postpred",
  function(design, ...) standardGeneric("interim_postpred")
)

#' @describeIn interim_postpred Interim analysis based on the posterior
#'   predictive probabilty for two-stage basket designs.
#'
#' @template design
#' @template n
#' @template n1
#' @template r1
#' @template lambda
#' @param weight_mat The matrix with all weights. Automatically calculated
#'   in the functions to which \code{interim_postpred} is passed.
#' @param prob_futstop Probability cut-off for stopping for futility.
#' @param prob_effstop Probability cut-off for stopping for efficacy.
setMethod("interim_postpred", "TwoStageBasket",
  function(design, n, n1, r1, lambda, weight_mat, prob_futstop = 0.1,
           prob_effstop = 0.9, ...) {
    crit <- get_crit(design = design, n = n, lambda = lambda)
    shape_borrow <- beta_borrow(weight_mat = weight_mat, design = design,
      n = n1, r = r1)
    pred_prob <- post_pred(n = n, n1 = n1, r1 = r1, shape = shape_borrow,
      crit = crit)
    ifelse(pred_prob < prob_futstop, -1, ifelse(pred_prob > prob_effstop, 1, 0))
  })

#' Interim analysis based on the posterior probability
#'
#' Conducts an interim analysis based on the posterior probability.
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{interim_posterior} conducts an interim analysis with possible
#' stop for efficacy and futility based on the posterior probability. If the
#' posterior probability is less than \code{prob_fustop} the basket is stopped
#' for futility, if the posterior probability is greater than
#' \code{prob_effstop} the basket is stopped for efficacy. If
#' \code{prob_fustop = 0} or \code{prob_effstop = 1} then no futility-stop and
#' no efficacy stop is possible, respectively.
#'
#' The function is generally not called by the user but passed to another
#' function such as \code{\link{toer}} and \code{\link{pow}} to specify which
#' interim analysis is conducted.
#'
#' @return A vector with a length equal to the number of baskets with
#' elements -1, 0 or 1 where -1 means stop for futility, 0 means continuation
#' and 1 means stop for efficacy.
#' @export
#'
#' @examples
#' design <- setupTwoStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 20, n1 = 10, lambda = 0.99, weight_fun = weights_fujikawa,
#'   interim_fun = interim_posterior, interim_params = list(prob_futstop = 0.05,
#'     prob_effstop = 0.95))
setGeneric("interim_posterior",
  function(design, ...) standardGeneric("interim_posterior")
)

#' @describeIn interim_posterior Interim analysis based on the posterior
#'   probabilty for two-stage basket designs.
#'
#' @template design
#' @template n1
#' @template r1
#' @param weight_mat The matrix with all weights. Automatically calculated
#'   in the functions to which \code{interim_postpred} is passed.
#' @param prob_futstop Probability cut-off for stopping for futility.
#' @param prob_effstop Probability cut-off for stopping for efficacy.
setMethod("interim_posterior", "TwoStageBasket",
  function(design, n1, r1, weight_mat, prob_futstop, prob_effstop, ...) {
    shape_borrow <- beta_borrow(weight_mat = weight_mat, design = design,
      n = n1, r = r1)
    post_prob <- post_beta(shape = shape_borrow, theta0 = design@theta0)
    ifelse(post_prob < prob_futstop, -1, ifelse(post_prob > prob_effstop, 1, 0))
  })
