#' @include generics.R
NULL

#' @describeIn weights_fujikawa Fujikawa-weights for a single-stage basket
#'   design.
setMethod("weights_fujikawa", "OneStageBasket",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune = FALSE, ...) {
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

    weight_mat
  })

#' @describeIn weights_fujikawa Fujikawa-weights for a two-stage basket design.
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

    weight_mat
  })
