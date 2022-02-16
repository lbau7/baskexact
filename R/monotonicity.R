#' @include class.R
NULL

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
#' check_mon_within(design = design, n = 24, lambda = 0.99,
#'   weight_fun = weights_fujikawa, weight_params = list(epsilon = 0.5,
#'    tau = 0), details = TRUE)
setGeneric("check_mon_within",
  function(design, ...) standardGeneric("check_mon_within")
)

#' @describeIn check_mon_within Within-trial monotonicity condition for a
#'   single-stage design.
#'
#' @template design
#' @template n
#' @template lambda
#' @template weights
#' @template details
#' @template dotdotdot
setMethod("check_mon_within", "OneStageBasket",
  function(design, n, lambda, weight_fun, weight_params = list(), details,
           ...) {
    # Not working with different priors and different n!
    check_params(n = n, lambda = lambda)

    crit <- get_crit(design = design, n = n, lambda = lambda)
    crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda))

    # Create matrix with all possible outcomes (without permutations)
    events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
    # Discard events where all or no baskets are significanct
    sel_events <- apply(events, 1,
      function(x) all(x >= crit) | all(x < crit_pool))
    events <- events[!sel_events, ]

    # If pruning is used and no list of violating outcomes is desired,
    # then outcomes with a different number of responses in baskets that
    # are pruned can be ignored -
    ### funktioniert aktuell nicht mehr -- verallgemeinern?
    # if (prune & !details) {
    #   events[which(events < crit_pool)] <- 0
    #   events <- events[!duplicated(events), ]
    # }

    func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
      r = x, weight_mat = weight_mat)

    # Conduct test for all remaining outcomes
    res <- t(apply(events, 1, func))
    res_sig <- rowSums(res)
    # Investigate outcomes where some (but not all) baskets are significant
    res_inv <- res[which(res_sig > 0 & res_sig < design@k), ]
    no_mon <- apply(res_inv, 1, function(x) any(x != cummax(x)))
    check <- sum(no_mon) == 0

    if (check) {
      check
    } else {
      if (details) {
        events_inv <- events[which(res_sig > 0 & res_sig < design@k), ]
        events_nomon <- events_inv[no_mon, ]
        res_nomon <- res_inv[no_mon, ]
        list(
          Events = events_nomon,
          Results = res_nomon
        )
      } else {
        check
      }
    }
  })

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
#' check_mon_between(design = design, n = 24, lambda = 0.99,
#'   weight_fun = weights_fujikawa, weight_params = list(epsilon = 3,
#'     tau = 0), details = TRUE)
setGeneric("check_mon_between",
  function(design, ...) standardGeneric("check_mon_between")
)

#' @describeIn check_mon_between Between-trial monotonicity condition for a
#'   single-stage design.
#'
#' @template design
#' @template n
#' @template lambda
#' @template weights
#' @template details
#' @template dotdotdot
setMethod("check_mon_between", "OneStageBasket",
  function(design, n, lambda, weight_fun, weight_params = list(), details,
           ...) {
    check_params(n = n, lambda = lambda)

    crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda))

    # Create matrix with all possible outcomes (without permutations)
    events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
    func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
      r = x, weight_mat = weight_mat)

    # Conduct test for all outcomes
    res <- t(apply(events, 1, func))
    # Select outcomes with at least one rejected null hypothesis
    res_sig <- apply(res, 1, function(x) any(x == 1))
    res_test <- res_sig

    if (details) detlist <- list()
    checkout <- TRUE
    # For every outcome with at least one rejected null hypothesis:
    # check whether there are any outcomes with at least as many responses
    # in each baskets but no rejected null hypotheses
    for (i in 1:nrow(events)) {
      if (res_test[i]) {
        # Check condition
        events_sel <- apply(events, 1, function(x) all(x >= events[i, ]))
        res_sel <- res_sig[events_sel]
        check <- sum(res_sel) == length(res_sel)
        if (check) {
          # If the condition is fulfilled for this outcome, then it is also
          # fulfilled for all outcomes with at least as many responses in
          # each basket
          res_test[events_sel] <- FALSE
        } else {
          if (details) {
            if (checkout) checkout <- FALSE
            detlist <- c(detlist, list(list(
              Events = rbind(events[i, ], events[events_sel & !res_sig, ]),
              Results = rbind(res[i, ], res[events_sel & !res_sig, ])
            )))
          } else {
            checkout <- FALSE
            break
          }
        }
      }
    }

    if (checkout) {
      TRUE
    } else {
      if (details) {
        detlist
      } else {
        FALSE
      }
    }
  })
