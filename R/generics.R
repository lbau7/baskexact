#' @include class.R
NULL

#' Type 1 Error Rate
#'
#' Computes the exact family wise type 1 error rate of a basket trial .
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{toer} computes the exact family wise type 1 error rate and the
#' exact rejection probabilities per group. The family wise type 1 error rate
#' is the probability to reject at least one null hypothesis for a basket with
#' theta1 = theta0. If all theta1 > theta0 then the family wise type 1 error
#' rate under the global null hypothesis is computed. The rejection
#' probabilities correspond to the type 1 error rate for baskets with theta1 =
#' theta 0 and to the power for baskets with theta1 > theta 0.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
#'
#' This method is implemented for the class \code{\link{OneStageBasket}}.
#'
#' @return If \code{results = "fwer"} then the family wise type 1 error rate is
#' returned as a numeric value. If \code{results = "group"} then a list with
#' the rejection probabilities per group and the family wise type 1 error rate
#' is returned. If all theta1 > theta0 then the family wise type 1 error rate
#' is calculated under the global null hypothesis. For baskets with theta1 =
#' theta0 the rejection probabilities corresponds to the type 1 error rate, for
#' baskets with theta1 > theta0 the rejection probabilities corresponds to the
#' power.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_fujikawa,
#'   weight_params = list(epsilon = 2, tau = 0))
setGeneric("toer",
  function(design, ...) standardGeneric("toer")
)

#' Power
#'
#' Computes the exact power for a basket trial.
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{pow} computes the exact experimentwise power and the
#' exact rejection probabilities per group. The experimentwise power
#' is the probability to reject at least one null hypothesis for a basket with
#' theta1 > theta0. The rejection probabilities correspond to the type 1 error
#' rate for baskets with theta1 = theta 0 and to the power for baskets with
#' theta1 > theta 0.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
#'
#' This method is implemented for the class \code{\link{OneStageBasket}}.
#'
#' @return If \code{results = "ewp"} then the experimentwise power is
#' returned as a numeric value. If \code{results = "group"} then a list with
#' the rejection probabilities per group and the experimentwise power
#' is returned. For baskets with theta1 = theta0 the rejection probabilities
#' corresponds to the type 1 error rate, for baskets with theta1 > theta0 the
#' rejection probabilities corresponds to the power.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' pow(design, theta1 = c(0.2, 0.5, 0.5), n = 15, lambda = 0.99,
#'   weight_fun = weights_fujikawa, weight_params = list(epsilon = 2, tau = 0))
setGeneric("pow",
  function(design, ...) standardGeneric("pow")
)

#' Check Within-Trial Monotonicity
#'
#' Checks whether the within-trial monotonicity condition holds.
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{check_mon_within} checks whether the within-trial
#' monotonicity condition holds. For a single-stage design with equal
#' prior distributions and equal sample sizes for each basket this condition
#' states that there are no cases where the null hypothesis of a basket is
#' rejected when there is at least one other basket with more observed
#' responses for which the null hypothesis cannot be rejected.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
#'
#' This method is implemented for the class \code{\link{OneStageBasket}}.
#'
#' @return If \code{details = FALSE} then only a logical value is returned.
#' If \code{details = TRUE} then if there are any cases where the
#' within-trial monotonicity condition is violated, a list of these cases and
#' their results are returned.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 4, shape1 = 1, shape2 = 1, theta0 = 0.2)
#' check_mon_within(design = design, n = 24, lambda = 0.99, epsilon = 0.5,
#'   tau = 0, prune = FALSE, details = TRUE)
setGeneric("check_mon_within",
  function(design, ...) standardGeneric("check_mon_within")
)

#' Check Between-Trial Monotonicity
#'
#' Checks whether the between-trial monotonicity condition holds.
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{check_mon_between} checks whether the between-trial
#' monotonicity condition holds. For a single-stage design with equal prior
#' distributions and equal sample sizes for each basket this condition states
#' that there are no cases where at least one null hypothesis is rejected when
#' when there is a case with an equal or higher number of responses in each
#' basket for which no null hypothesis is rejected.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
#'
#' This method is implemented for the class \code{\link{OneStageBasket}}.
#'
#' @return If \code{details = FALSE} then only a logical value is returned.
#' If \code{details = TRUE} then if there are any cases where the
#' between-trial monotonicity condition is violated, a list of theses cases
#' and their results are returned.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 4, shape1 = 1, shape2 = 1, theta0 = 0.2)
#' check_mon_between(design = design, n = 24, lambda = 0.99, epsilon = 3,
#'   tau = 0, prune = FALSE, details = TRUE)
setGeneric("check_mon_between",
  function(design, ...) standardGeneric("check_mon_between")
)

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
#' This method is implemented for the class \code{\link{OneStageBasket}}.
#'
#' @return The greatest value with \code{prec_digits} decimal places for
#' \code{lambda} which controls the family wise error rate at level
#' \code{alpha} (one-sided) and the exact family wise error rate for this
#' value of \code{lambda}.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
#' adjust_lambda(design = design, alpha = 0.025, n = 15, epsilon = 1, tau = 0,
#'   logbase = 2, prune = FALSE, prec_digits = 4)
setGeneric("adjust_lambda",
  function(design, ...) standardGeneric("adjust_lambda")
)

#' Test for the Results of a Basket Trial
#'
#' \code{basket_test} evaluates the results of a basket trial and calculates
#' the posterior distributions with and without borrowing.
#'
#' @template design
#' @template n
#' @param r The vector of observed responses.
#' @template lambda
#' @template tuning
#' @template prune
#' @template dotdotdot
#'
#' @return A list, including matrices of the weights that are used for
#' borrowing, posterior distribution parameters for all baskets without and
#' with borrowing, as well as the posterior probabilities for all baskets
#' without and with borrowing.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
#' basket_test(design = design, n = 24, r = c(5, 9, 10), lambda = 0.99,
#'   epsilon = 1, tau = 0, logbase = 2, prune = FALSE)
setGeneric("basket_test",
  function(design, ...)
    standardGeneric("basket_test")
)

#' Weights Based on Fujikawa et al.'s Design
setGeneric("weights_fujikawa",
  function(design, ...) standardGeneric("weights_fujikawa")
)

#' Interim Analysis Based on the Posterior Predictive Probability
setGeneric("interim_postpred",
  function(design, ...) standardGeneric("interim_postpred")
)
