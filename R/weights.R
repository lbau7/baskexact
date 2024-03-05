#' @include class.R
NULL

#' Weights Based on Fujikawa et al.'s Design
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{weights_fujikawa} calculates the weights used for sharing
#' information between baskets based on the proposal by Fujikawa et al. (2020).
#' The weight for two baskets i and j is found as
#' \eqn{(1 - JSD(i, j))^\varepsilon} where \eqn{JSD(i, j)} is the Jensen-Shannon
#' divergence between the individual posterior distributions of the response
#' probabilities of basket i and j. Note that Fujikawa's weights also share the
#' prior information between the baskets.
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
#' @references Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
#' A Bayesian basket trial design that borrows information across strata based
#' on the similarity between the posterior distributions of the response
#' probability. Biometrical Journal, 62(2), 330-338.
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
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
#' @template tuning_jsd
#' @template prune
#' @template globalweights
setMethod("weights_fujikawa", "OneStageBasket",
  function(design, n, lambda, epsilon = 1.25, tau = 0.5, logbase = 2,
           prune = FALSE, globalweight_fun = NULL,
           globalweight_params = list(), ...) {
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

    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda,
        weight_mat = weight_mat, globalweight_fun = globalweight_fun,
        globalweight_params = globalweight_params)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    weight_mat
  })

#' @describeIn weights_fujikawa Fujikawa-weights for a two-stage basket design.
#'
#' @template design
#' @template n
#' @template n1
#' @template tuning_jsd
setMethod("weights_fujikawa", "TwoStageBasket",
  function(design, n, n1, epsilon = 1.25, tau = 0, logbase = 2, ...) {
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

#' Weights Based on the Jensen-Shannon Divergence
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{weights_jsd} calculates the weights used for sharing
#' information between baskets based on the Jensen-Shannon divergence (JSD).
#' The weight for two baskets i and j is found as
#' \eqn{(1 - JSD(i, j))^\varepsilon} where \eqn{JSD(i, j)} is the Jensen-Shannon
#' divergence between the individual posterior distributions of the response
#' probabilities of basket i and j. This is identical to how the weights are
#' calculated in \code{\link{weights_fujikawa}}, however when Fujikawa's weights
#' are used the prior information is also shared.
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
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_jsd)
setGeneric("weights_jsd",
  function(design, ...) standardGeneric("weights_jsd")
)

#' @describeIn weights_jsd Jensen-Shannon Divergence weights for a
#'   single-stage basket design.
#'
#' @template design
#' @template n
#' @template lambda
#' @template tuning_jsd
#' @template globalweights
#' @template prune
setMethod("weights_jsd", "OneStageBasket",
  function(design, n, lambda, epsilon = 1.25, tau = 0.5, logbase = 2,
           prune = FALSE, globalweight_fun = NULL,
           globalweight_params = list(), ...) {
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
          kl_f <- stats::integrate(f, 0, 1)$value
          kl_g <- stats::integrate(g, 0, 1)$value
          jsd_mat[i, j] <- 0.5 * kl_f + 0.5 * kl_g
        }
      }
    }
    jsd_mat <- jsd_mat + t(jsd_mat)
    weight_mat <- (1 - jsd_mat)^epsilon
    weight_mat[weight_mat <= tau] <- 0
    class(weight_mat) <- "pp"

    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda,
        weight_mat = weight_mat, globalweight_fun = globalweight_fun,
        globalweight_params = globalweight_params)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    weight_mat
  })

#' @describeIn weights_jsd Jensen-Shannon Divergence weights for a two-stage
#'   basket design.
#'
#' @template design
#' @template n
#' @template n1
#' @template tuning_jsd
#' @template prune
setMethod("weights_jsd", "TwoStageBasket",
  function(design, n, n1, epsilon = 1.25, tau = 0, logbase = 2, ...) {
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
          kl_f <- stats::integrate(f, 0, 1)$value
          kl_g <- stats::integrate(g, 0, 1)$value
          jsd_mat[i, j] <- 0.5 * kl_f + 0.5 * kl_g
        }
      }
    }
    jsd_mat <- jsd_mat + t(jsd_mat)
    weight_mat <- (1 - jsd_mat)^epsilon
    weight_mat[weight_mat <= tau] <- 0

    class(weight_mat) <- "pp"
    weight_mat
  })

#' Weights Based on the Calibrated Power Prior
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{weights_cpp} calculates the weights based on an approach
#' by Pan & Yuan (2017). The weight for two baskets i and j is found by at
#' first calculating \eqn{S_{KS;i,j}} as the Kolmogorov-Smirnov statistic,
#' which is equal to the difference in response rates for binary variables.
#' \eqn{S_{KS;i,j}} is then transformed to \eqn{S_{i,j} = n^{1/4}S_{KS;i,j}}.
#' Then the weight is found as \eqn{1 / (1 + exp(a + b * log(S_{i,j})))}, where
#' a and b are tuning parameters.
#'
#' The function is generally not called by the user but passed to another
#' function such as \code{\link{toer}} and \code{\link{pow}} to specificy
#' how the weights are calculated.
#'
#' @references Pan, H., Yuan, Y., & Xia, J. (2017). A calibrated power prior
#' approach to borrow information from historical data with application to
#' biosimilar clinical trials. Journal of the Royal Statistical Society Series
#' C: Applied Statistics, 66(5), 979-996.
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_cpp)
setGeneric("weights_cpp",
  function(design, ...) standardGeneric("weights_cpp")
)

#' @describeIn weights_cpp Calibrated power prior weights for a single-stage
#'   basket design.
#'
#' @template design
#' @template n
#' @template tuning_cpp
#' @template prune
#' @template lambda
#' @template globalweights
#' @template dotdotdot
setMethod("weights_cpp", "OneStageBasket",
  function(design, n, a = 1, b = 1, prune = FALSE, lambda,
    globalweight_fun = NULL, globalweight_params = list(), ...) {
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
          s <- n^(1 / 4) * ks
          weight_mat[i, j] <- g(s = s, a = a, b = b)
        }
      }
    }
    weight_mat <- weight_mat + t(weight_mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"

    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda,
        weight_mat = weight_mat, globalweight_fun = globalweight_fun,
        globalweight_params = globalweight_params)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    weight_mat
  })

#' @describeIn weights_cpp Calibrated power prior weights for a two-stage
#'   basket design.
#'
#' @template design
#' @template n
#' @template n1
#' @template tuning_cpp
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
          s <- max(nvec[i], nvec[j])^(1 / 4) * ks
          weight_mat[i, j] <- g(s = s, a = a, b = b)
        }
      }
    }
    weight_mat <- weight_mat + t(weight_mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"
    weight_mat
  })

#' Weights Based on the Marginal Maximum Likelihood
#'
#' @template design
#' @template dotdotdot
#'
#' @details \code{weights_mml} calculates the weights based on the marginal
#' maximum likelihood approach by Gravestock & Held (2017). In this approach,
#' the weight is found as the maximum of the marginal likelihood of the
#' weight-parameter given the dataset that information should be borrowed
#' from. However, since this can lead to non-symmetric weights (meaning that
#' the amount of information that data set 1 borrows from data set 2 is
#' generally not identical to the information data set 2 borrows from data set
#' 1), a symmetrised version is used here: For the sharing-weight of
#' Basket 1 and Basket 2 the MML is calculted two times - once conditional
#' on the data of Basket 1 and once conditional on the data of Basket 2.
#' The mean of these two weights is then used, resulting in symmetrical
#' sharing.
#'
#' @references Gravestock, I., & Held, L. (2017). Adaptive power priors with
#' empirical Bayes for clinical trials. Pharmaceutical statistics, 16(5),
#' 349-360.
#'
#' @return A matrix including the weights of all possible pairwise outcomes.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_mml)
setGeneric("weights_mml",
  function(design, ...) standardGeneric("weights_mml")
)

#' @describeIn weights_mml Maximum marginal likelihood weights for a
#'   single-stage basket design
#'
#' @template design
#' @template n
#' @template dotdotdot
setMethod("weights_mml", "OneStageBasket",
  function(design, n, ...) {
    n_sum <- n + 1
    mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    r <- 0:n
    for (i in 1:n_sum) {
      for (j in 1:n_sum) {
        f <- function(delta) -extraDistr::dbbinom(
          x = r[i],
          size = n,
          alpha = design@shape1 + delta * r[j],
          beta = design@shape2 + delta * (n - r[j])
        )

        l <- stats::optim(0.5, fn = f, lower = 0, upper = 1,
          method = "Brent")$par
        mat[i, j] <- ifelse(l <= 6.474096e-09, 0, l)
      }
    }
    mat <- (mat + t(mat)) / 2
    diag(mat) <- 1
    class(mat) <- "pp"
    mat
  })

#' @describeIn weights_mml Maximum marginal likelihood weights for a
#'   two-stage basket design
#'
#' @template design
#' @template n
#' @template n1
#' @template dotdotdot
setMethod("weights_mml", "TwoStageBasket",
  function(design, n, n1, ...) {
    r <-  c(0:n1, 0:n)
    nvec <- c(rep(n1, n1 + 1), rep(n, n + 1))
    n_sum <- n + n1 + 2
    mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    for (i in 1:n_sum) {
      for (j in i:n_sum) {
        if (i == j) {
          next
        } else {
          f1 <- function(delta) -extraDistr::dbbinom(
            x = r[i],
            size = nvec[i],
            alpha = design@shape1 + delta * r[j],
            beta = design@shape2 + delta * (nvec[j] - r[j])
          )

          f2 <- function(delta) -extraDistr::dbbinom(
            x = r[j],
            size = nvec[j],
            alpha = design@shape1 + delta * r[i],
            beta = design@shape2 + delta * (nvec[i] - r[i])
          )

          l1 <- stats::optim(0.5, fn = f1, lower = 0, upper = 1,
            method = "Brent")$par
          l2 <- stats::optim(0.5, fn = f2, lower = 0, upper = 1,
            method = "Brent")$par
          l <- mean(c(l1, l2))
          mat[i, j] <- ifelse(l <= 6.474096e-09, 0, l)
        }
      }
    }
    weight_mat <- mat + t(mat)
    diag(weight_mat) <- 1
    class(weight_mat) <- "pp"
    weight_mat
  })

#' Separate Analysis in Each Basket
#'
#' @template design
#' @template dotdotdot
#'
#' @details When \code{weights_separate} is used as a weight function, a
#' separate analysis performed in each basket.
#'
#' @return A weight matrix where all weights are 0.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_separate)
setGeneric("weights_separate",
  function(design, ...) standardGeneric("weights_separate")
)

#' @describeIn weights_separate Separate analysis for a single-stage basket
#'   design
#'
#' @template design
#' @template n
#' @template dotdotdot
setMethod("weights_separate", "OneStageBasket",
  function(design, n, ...) {
    n_sum <- n + 1
    mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    class(mat) <- "pp"
    mat
  })

#' @describeIn weights_separate Separate analysis for a two-stage basket
#'   design
#'
#' @template design
#' @template n
#' @template n1
#' @template dotdotdot
setMethod("weights_separate", "TwoStageBasket",
  function(design, n, n1, ...) {
    n_sum <- n + n1 + 2
    mat <- matrix(0, nrow = n_sum, ncol = n_sum)
    class(mat) <- "pp"
    mat
  })

#' Pooled Analysis
#'
#' @template design
#' @template dotdotdot
#'
#' @details When \code{weights_pool} is used as a weight function, all data
#' are pooled.
#'
#' @return A weight matrix where all weights are 1.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
#' toer(design, n = 15, lambda = 0.99, weight_fun = weights_pool)
setGeneric("weights_pool",
  function(design, ...) standardGeneric("weights_pool")
)

#' @describeIn weights_pool Pooled analysis for a single-stage basket design
#'
#' @template design
#' @template n
#' @template n1
#' @template dotdotdot
setMethod("weights_pool", "OneStageBasket",
  function(design, n, ...) {
    n_sum <- n + 1
    mat <- matrix(1, nrow = n_sum, ncol = n_sum)
    class(mat) <- "pp"
    mat
  })

#' @describeIn weights_pool Pooled analysis for a single-stage basket design
#'
#' @template design
#' @template n
#' @template dotdotdot
setMethod("weights_pool", "OneStageBasket",
  function(design, n, ...) {
    n_sum <- n + 1
    mat <- matrix(1, nrow = n_sum, ncol = n_sum)
    class(mat) <- "pp"
    mat
  })

#' @describeIn weights_pool Pooled analysis for a two-stage basket design
#'
#' @template design
#' @template n
#' @template n1
#' @template dotdotdot
setMethod("weights_pool", "TwoStageBasket",
  function(design, n, n1, ...) {
    n_sum <- n + n1 + 2
    mat <- matrix(1, nrow = n_sum, ncol = n_sum)
    class(mat) <- "pp"
    mat
  })

