#' Type 1 Error Rate
#'
#' Computes the exact family wise type 1 error rate of a basket trial .
#'
#' @template design
#' @template n
#' @template lambda
#' @template tuning
#' @template prune
#' @param results Whether only the family wise error rate (option \code{fwer})
#'   or also the rejection probabilities per group (option \code{group}) should
#'   be returned.
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
#' Calculations are based on the design of Fujikawa et al.:
#' At first separate posterior distributions are calculated for each basket
#' based on a beta-binomial model. Information is then borrowed between baskets
#' by calculating weights that reflect the similarity between the basket and
#' computing posterior distributions for each basket where the parameters
#' of the beta posterior are calculated as a weighted sum of the individual
#' posterior distributions. The weight between two baskets i and j is found as
#' (1 - JSD(i, j))^epsilon where JSD(i, j) is the Jensen-Shannon divergence
#' between basket i and j. A small value of epsilon results in stronger
#' borrowing also across baskets with heterogenous results. If epsilon is
#' large then information is only borrowed between baskets with similar results.
#' If a weight is smaller than tau it is set to 0, which results in no
#' borrowing. If for a basket the posterior probability that \eqn{\theta} >
#' theta0 is greater than lambda, then the null hypothesis is rejected.
#' @references Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
#' A Bayesian basket trial design that borrows information across strata based
#' on the similarity between the posterior distributions of the response
#' probability. Biometrical Journal, 62(2), 330-338.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
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
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2, theta1 = c(0.2, 0.5, 0.5))
#' toer(design, n = 15, lambda = 0.99, epsilon = 2, tau = 0)
setGeneric("toer",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune = FALSE,
           results = c("fwer", "group"), ...) standardGeneric("toer")
)

#' Power
#'
#' Computes the exact power for a basket trial.
#'
#' @template design
#' @template n
#' @template lambda
#' @template tuning
#' @template prune
#' @param results Whether only the experimentwise power (option \code{ewp})
#'   or also the rejection probabilities per group (option \code{group}) should
#'   be returned.
#' @template dotdotdot
#'
#' @details \code{pow} computes the exact experimentwise power and the
#' exact rejection probabilities per group. The experimentwise power
#' is the probability to reject at least one null hypothesis for a basket with
#' theta1 > theta0. The rejection probabilities correspond to the type 1 error
#' rate for baskets with theta1 = theta 0 and to the power for baskets with
#' theta1 > theta 0.
#'
#' Calculations are based on the design of Fujikawa et al.:
#' At first separate posterior distributions are calculated for each basket
#' based on a beta-binomial model. Information is then borrowed between baskets
#' by calculating weights that reflect the similarity between the basket and
#' computing posterior distributions for each basket where the parameters
#' of the beta posterior are calculated as a weighted sum of the individual
#' posterior distributions. The weight between two baskets i and j is found as
#' (1 - JSD(i, j))^epsilon where JSD(i, j) is the Jensen-Shannon divergence
#' between basket i and j. A small value of epsilon results in stronger
#' borrowing also across baskets with heterogenous results. If epsilon is
#' large then information is only borrowed between baskets with similar results.
#' If a weight is smaller than tau it is set to 0, which results in no
#' borrowing. If for a basket the posterior probability that \eqn{\theta} >
#' theta0 is greater than lambda, then the null hypothesis is rejected.
#' @references Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
#' A Bayesian basket trial design that borrows information across strata based
#' on the similarity between the posterior distributions of the response
#' probability. Biometrical Journal, 62(2), 330-338.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
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
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2, theta1 = c(0.2, 0.5, 0.5))
#' pow(design, n = 15, lambda = 0.99, epsilon = 2, tau = 0)
setGeneric("pow",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune = FALSE,
    results = c("ewp", "group"), ...) standardGeneric("pow")
)

#' Check Within-Trial Monotonicity
#'
#' Checks whether the within-trial monotonicity condition holds.
#'
#' @template design
#' @template n
#' @template lambda
#' @template tuning
#' @template prune
#' @template details
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
#' @return If \code{details = FALSE} then only a logical value is returned.
#' If \code{details = TRUE} then if there are any cases where the
#' within-trial monotonicity condition is violated, a list of these cases and
#' their results are returned.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 4, shape1 = 1, shape2 = 1, theta0 = 0.2,
#'   theta1 = c(0.2, 0.5, 0.5, 0.5))
#' check_mon_within(design = design, n = 24, lambda = 0.99, epsilon = 0.5,
#'   tau = 0, prune = FALSE, details = TRUE)
setGeneric("check_mon_within",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune, details, ...)
    standardGeneric("check_mon_within")
)

#' Check Between-Trial Monotonicity
#'
#' Checks whether the between-trial monotonicity condition holds.
#'
#' @template design
#' @template n
#' @template lambda
#' @template tuning
#' @template prune
#' @template details
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
#' @return If \code{details = FALSE} then only a logical value is returned.
#' If \code{details = TRUE} then if there are any cases where the
#' between-trial monotonicity condition is violated, a list of theses cases
#' and their results are returned.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 4, shape1 = 1, shape2 = 1, theta0 = 0.2,
#'   theta1 = c(0.2, 0.5, 0.5, 0.5))
#' check_mon_between(design = design, n = 24, lambda = 0.99, epsilon = 3,
#'   tau = 0, prune = FALSE, details = TRUE)
setGeneric("check_mon_between",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune, details, ...)
    standardGeneric("check_mon_between")
)

#' Adjust Lambda
#'
#' Finds the value for \code{lambda} such that the family wise error
#' rate is protected at level \code{alpha}.
#'
#' @template design
#' @param alpha The one-sided signifance level.
#' @template n
#' @template tuning
#' @template prune
#' @param prec_digits Number of decimal places that are considered when
#'   adjusting lambda
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
#' design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2,
#'   theta1 = c(0.2, 0.2, 0.2))
#' adjust_lambda(design = design, alpha = 0.025, n = 15, epsilon = 1, tau = 0,
#'   logbase = 2, prune = FALSE, prec_digits = 4)
setGeneric("adjust_lambda",
  function(design, alpha = 0.025, n, epsilon, tau, logbase, prune,
           prec_digits, ...)
    standardGeneric("adjust_lambda")
)
