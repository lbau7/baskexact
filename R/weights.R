#' @include class.R
NULL

#' Weights Based on Fujikawa et al.'s Design
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{weights_fujikawa} calculates the weights used for sharing
#' information between baskets based on the proposal by Fujikawa et al. (2020).
#' The weight between two baskets i and j is found as (1 - JSD(i, j))^epsilon
#' where JSD(i, j) is the Jensen-Shannon divergence between the individual
#' posterior distributions of the response probabilities of basket i and j.
#' Note that Fujikawa's weights also share the prior information between the
#' baskets.
#'
#' A small value of epsilon results in stronger borrowing also across baskets
#' with heterogenous results. If epsilon is large then information is only
#' borrowed between baskets with similar results. If a weight is smaller than
#' tau it is set to 0, which results in no borrowing.
#'
#' If \code{prune = TRUE} then the baskets with an observed number of baskets
#' smaller than the pooled critical value are not borrowed from. The
#' pooled critical value is the smallest integer c for which all null
#' hypotheses can be rejected if the number of responses is exactly c for
#' all baskets.
#'
#' The function is generally not called by the user but passed to another
#' function such as \code{\link{toer}} and \code{\link{pow}} to specificy
#' how the weights are calculated.
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_fujikawa)
setGeneric("weights_fujikawa",
  function(design, ...) standardGeneric("weights_fujikawa")
)

#' @describeIn weights_fujikawa Fujikawa-weights for a single-stage basket
#'   design.
#'
#' @template design
#' @template n
#' @template lambda
#' @template tuning_fujikawa
#' @template prune
setMethod("weights_fujikawa", "OneStageBasket",
  function(design, n, lambda, epsilon = 1.25, tau = 0.5, logbase = 2,
           prune = FALSE, ...) {
    shape1_post <- design@shape1 + c(0:n)
    shape2_post <- design@shape2 + c(n:0)
    n_sum <- n + 1

    p <- function(x) stats::dbeta(x, shape1_post[i], shape2_post[i])
    q <- function(x) stats::dbeta(x, shape1_post[j], shape2_post[j])
    m <- function(x) 0.5 * (p(x) + q(x))
    f <- function(x) p(x) * log(p(x) / m(x), base = logbase)
    g <- function(x) q(x) * log(q(x) / m(x), base = logbase)
    h <- function(x) 0.5 * f(x) + 0.5 * g(x)
    jsd_mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
          next
        } else {
          p <- function(x) stats::dbeta(x, shape1_post[i], shape2_post[i])
          q <- function(x) stats::dbeta(x, shape1_post[j], shape2_post[j])
          m <- function(x) 0.5 * (p(x) + q(x))
          f <- function(x) p(x) * log(p(x) / m(x), base = logbase)
          g <- function(x) q(x) * log(q(x) / m(x), base = logbase)
          kl_f <- stats::integrate(f, 0, 1)$value
          kl_g <- stats::integrate(g, 0, 1)$value
          jsd_mat[i, j] <- 0.5 * kl_f + 0.5 * kl_g
        }
      }
    }
    jsd_mat <- jsd_mat + t(jsd_mat)
    weight_mat <- (1 - jsd_mat)^epsilon
    weight_mat[weight_mat <= tau] <- 0

    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    class(weight_mat) <- "fujikawa"
    weight_mat
  })

#' @describeIn weights_fujikawa Fujikawa-weights for a two-stage basket design.
#'
#' @template design
#' @template n
#' @template n1
#' @template lambda
#' @template tuning_fujikawa
#' @template prune
setMethod("weights_fujikawa", "TwoStageBasket",
  function(design, n, n1, lambda, epsilon = 1.25, tau = 0, logbase = 2,
           prune = FALSE, ...) {
    # prune is currently not used - pruning-method not easily extendable
    # to a two-stage design
    shape1_post <- design@shape1 + c(0:n1, 0:n)
    shape2_post <- design@shape2 + c(n1:0, n:0)
    n_sum <- n + n1 + 2

    p <- function(x) stats::dbeta(x, shape1_post[i], shape2_post[i])
    q <- function(x) stats::dbeta(x, shape1_post[j], shape2_post[j])
    m <- function(x) 0.5 * (p(x) + q(x))
    f <- function(x) p(x) * log(p(x) / m(x), base = logbase)
    g <- function(x) q(x) * log(q(x) / m(x), base = logbase)
    h <- function(x) 0.5 * f(x) + 0.5 * g(x)
    jsd_mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
          next
        } else {
          p <- function(x) stats::dbeta(x, shape1_post[i], shape2_post[i])
          q <- function(x) stats::dbeta(x, shape1_post[j], shape2_post[j])
          m <- function(x) 0.5 * (p(x) + q(x))
          f <- function(x) p(x) * log(p(x) / m(x), base = logbase)
          g <- function(x) q(x) * log(q(x) / m(x), base = logbase)
          kl_f <- stats::integrate(f, 0, 1)$value
          kl_g <- stats::integrate(g, 0, 1)$value
          jsd_mat[i, j] <- 0.5 * kl_f + 0.5 * kl_g
        }
      }
    }
    jsd_mat <- jsd_mat + t(jsd_mat)
    weight_mat <- (1 - jsd_mat)^epsilon
    weight_mat[weight_mat <= tau] <- 0

    class(weight_mat) <- "fujikawa"
    weight_mat
  })

#' Weights Based on Maximising the Marginal Likelihood
#'
#' @template design
#' @template dotdotdot
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_mml)
setGeneric("weights_mml",
  function(design, ...) standardGeneric("weights_mml")
)

#' @describeIn weights_mml Maximum marginal likelihood weights for a
#'   single-stage basket design.
#'
#' @template design
#' @template n
#' @template lambda
#' @template dotdotdot
setMethod("weights_mml", "OneStageBasket",
  function(design, n, lambda, ...) {
    x1 <- x2 <- c(0:n)
    n_sum <- n + 1

    weight_mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
        } else {
          l_x1 <- function(x) stats::dbinom(x1[i], n, prob = x)
          l_x2 <- function(x) stats::dbinom(x2[j], n, prob = x)
          prior <- function(x) stats::dbeta(x = x, shape1 = design@shape1,
            shape2 = design@shape2)
          l_marg <- function(delta) {
            a <- stats::integrate(function(x) l_x1(x) * l_x2(x)^delta *
              prior(x), lower = 0, upper = 1)$value
            b <- stats::integrate(function(x) l_x1(x)^delta * prior(x),
              lower = 0, upper = 1)$value
            - a / b
          }
          weight_mat[i, j] <- stats::optim(0.5, l_marg, method = "Brent",
            lower = 0, upper = 1)$par
        }
      }
    }
    # Achtung: Borrowing funktioniert anders als mit Fujikawa
    weight_mat <- weight_mat + t(weight_mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"
    weight_mat
  })

#' @describeIn weights_mml Maximum marginal likelihood weights for a
#'   two-stage basket design.
#'
#' @template design
#' @template n
#' @template n1
#' @template lambda
#' @template dotdotdot
setMethod("weights_mml", "TwoStageBasket",
  function(design, n, n1, lambda, ...) {
    x1 <- x2 <- c(0:n1, 0:n)
    n_sum <- n + n1 + 2

    weight_mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
          next
        } else {
          l_x1 <- function(x) stats::dbinom(x1[i], n, prob = x)
          l_x2 <- function(x) stats::dbinom(x2[j], n, prob = x)
          prior <- function(x) stats::dbeta(x = x, shape1 = design@shape1,
            shape2 = design@shape2)
          l_marg <- function(delta) {
            a <- stats::integrate(function(x) l_x1(x) * l_x2(x)^delta *
              prior(x), lower = 0, upper = 1)$value
            b <- stats::integrate(function(x) l_x1(x)^delta * prior(x),
              lower = 0, upper = 1)$value
            - a / b
          }
          weight_mat[i, j] <- stats::optim(0.5, l_marg, method = "Brent",
            lower = 0, upper = 1)$par
        }
      }
    }
    # Achtung: Borrowing funktioniert anders als mit Fujikawa
    weight_mat <- weight_mat + t(weight_mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"
    weight_mat
  })
