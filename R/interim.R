#' @include class.R
NULL

#' Interim analysis based on the posterior predictive probability
setGeneric("interim_postpred",
  function(design, ...) standardGeneric("interim_postpred")
)

#' @describeIn interim_postpred Interim analysis based on the posterior
#'   predictive probabilty for two-stage basket designs.
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
setGeneric("interim_posterior",
  function(design, ...) standardGeneric("interim_posterior")
)

#' @describeIn interim_posterior Interim analysis based on the posterior
#'   probabilty for two-stage basket designs.
setMethod("interim_posterior", "TwoStageBasket",
  function(design, n1, r1, weight_mat, prob_futstop, prob_effstop, ...) {
    shape_borrow <- beta_borrow(weight_mat = weight_mat, design = design,
      n = n1, r = r)
    post_prob <- post_beta(shape = shape_borrow, theta0 = design@theta0)
    ifelse(post_prob < prob_futstop, -1, ifelse(post_prob > prob_effstop, 1, 0))
  })
