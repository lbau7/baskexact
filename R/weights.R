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
#' @template tuning_fujikawa
#' @template prune
setMethod("weights_fujikawa", "TwoStageBasket",
  function(design, n, n1, epsilon = 1.25, tau = 0, logbase = 2,
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


#' Weights Based on the Equivalence Probability Weight
#'
#' @template design
#' @template dotdotdot
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
setGeneric("weights_eqprob",
  function(design, ...) standardGeneric("weights_eqprob")
)

#' @describeIn weights_eqprob Equivalence probability weights for a
#'   single-stage basket design.
#'
#' @template design
#' @template n
#' @param bound Equivalence bound.
#' @template dotdotdot
setMethod("weights_eqprob", "OneStageBasket",
  function(design, n, bound, ...) {
    n_sum <- n + 1
    weight_mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    p1 <- p2 <- 0:n / n

    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
        } else {
          denom <- sqrt(((p1[i] * (1 - p1[i])) / n) +
              ((p2[j] * (1 - p2[j])) / n))
          a <- pnorm(bound - (p1[i] - p2[j]) / denom)
          b <- pnorm(-bound - (p1[i] - p2[j]) / denom)
          weight_mat[i, j] <- a - b
        }
      }
    }
    weight_mat <- weight_mat + t(weight_mat)

    # Calculate diagonal weight using arbitrary element
    # Caution: If 0 responses in both arms the result is NaN
    denom <- sqrt(((p1[2] * (1 - p1[2])) / n) +
        ((p2[2] * (1 - p2[2])) / n))
    a <- pnorm(bound - (p1[2] - p2[2]) / denom)
    b <- pnorm(-bound - (p1[2] - p2[2]) / denom)
    diag_weight <- a - b

    diag(weight_mat) <- diag_weight
    class(weight_mat) <- "pp"
    weight_mat
  })

#' Weights Based on the Probability Weight
#'
#' @template design
#' @template dotdotdot
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_prob)
setGeneric("weights_prob",
  function(design, ...) standardGeneric("weights_prob")
)

#' @describeIn weights_eqprob Probability weights for a single-stage basket
#'   design.
#'
#' @template design
#' @template n
#' @template dotdotdot
setMethod("weights_prob", "OneStageBasket",
  function(design, n, ...) {
    shape1_post <- design@shape1 + c(0:n)
    shape2_post <- design@shape2 + c(n:0)
    n_sum <- n + 1
    weight_mat <- matrix(0, nrow = n_sum, ncol = n_sum)

    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
        } else {
          f <- function(x) {
            stats::dbeta(
              x = x,
              shape1 = shape1_post[i],
              shape2 = shape2_post[i]
            ) *
              stats::pbeta(
                q = x,
                shape1 = shape1_post[j],
                shape2 = shape2_post[j]
              )
          }
          temp <- stats::integrate(f, lower = 0, upper = 1)$value
          weight_mat[i, j] <- 2 * min(temp, 1 - temp)
        }
      }
    }

    weight_mat <- weight_mat + t(weight_mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"
    weight_mat
  })

#' Weights Based on the Calibrated Power Prior
#'
#' @template design
#' @template dotdotdot
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, theta0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_cpp)
setGeneric("weights_cpp",
  function(design, ...) standardGeneric("weights_cpp")
)

#' @describeIn weights_cpp Calibrated power prior weights for a single-stage
#'   basket design.
#'
#' @template design
#' @template n
#' @param a first tuning parameter
#' @param b second tuning parameter
#' @template dotdotdot
setMethod("weights_cpp", "OneStageBasket",
  function(design, n, a = 1, b = 1, ...) {
    n_sum <- n + 1
    weight_mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    r1 <- r2 <- 0:n

    g <- function(s, a, b) {
      1 / (1 + exp(a + b * log(s)))
    }

    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
          next
        } else {
          vec1 <- rep(0:1, c(n - r1[i], r1[i]))
          vec2 <- rep(0:1, c(n - r2[j], r2[j]))
          ks <- suppressWarnings(stats::ks.test(vec1, vec2)$statistic)
          s <- n^(1/4) * ks
          weight_mat[i, j] <- g(s = s, a = a, b = b)
        }
      }
    }
    weight_mat <- weight_mat + t(weight_mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"
    weight_mat
  })

#' @describeIn weights_cpp Calibrated power prior weights for a two-stage
#'   basket design.
#'
#' @template design
#' @template n
#' @template n1
#' @param a first tuning parameter
#' @param b second tuning parameter
#' @template dotdotdot
setMethod("weights_cpp", "TwoStageBasket",
  function(design, n, n1, a = 1, b = 1, ...) {
    r <-  c(0:n1, 0:n)
    nvec <- c(rep(n1, n1 + 1), rep(n, n + 1))
    n_sum <- n + n1 + 2
    weight_mat <- matrix(0, nrow = n_sum, ncol = n_sum)

    g <- function(s, a, b) {
      1 / (1 + exp(a + b * log(s)))
    }

    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
          next
        } else {
          vec1 <- rep(0:1, c(nvec[i] - r[i], r[i]))
          vec2 <- rep(0:1, c(nvec[j] - r[j], r[j]))
          ks <- suppressWarnings(stats::ks.test(vec1, vec2)$statistic)
          s <- n^(1/4) * ks
          weight_mat[i, j] <- g(s = s, a = a, b = b)
        }
      }
    }
    weight_mat <- weight_mat + t(weight_mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"
    weight_mat
  })
