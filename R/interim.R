#' @include class.R
NULL

#' Interim Analysis Based on the Posterior Predictive Probability
setGeneric("interim_postpred",
  function(design, ...) standardGeneric("interim_postpred")
)

#' @describeIn interim_postpred Interim analysis based on the posterior
#'   predictive probabilty for two-stage basket designs.
setMethod("interim_postpred", "TwoStageBasket",
  function(design, n, n1, r1, lambda, weight_mat, prob_futstop = 0.1,
           prob_effstop = 0.9, ...) {
    crit <- get_crit(design = design, n = n, lambda = lambda)
    shape_post <- matrix(c(design@shape1 + r1, design@shape2 + n1 - r1),
      byrow = TRUE, ncol = design@k)
    shape_borrow <- beta_borrow(k = design@k, r = r1, weight_mat = weight_mat,
      shape = shape_post)
    pred_prob <- post_pred(n = n, n1 = n1, r1 = r1, shape = shape_borrow,
      crit = crit)
    ifelse(pred_prob < prob_futstop, -1, ifelse(pred_prob > prob_effstop, 1, 0))
  })
