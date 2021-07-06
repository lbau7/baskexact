# Calculates the weight-matrix for all possible outcomes in a single-stage
# design
get_weights <- function(design, n, epsilon, tau, logbase) {
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
  weight_mat
}

# Computes the posterior distribution with borrowing
beta_borrow <- function(k, r, weight_mat, shape) {
  all_combs <- arrangements::combinations(r, 2) + 1
  weights_vec <- weight_mat[all_combs]
  weight_beta(k = k, weights = weights_vec, shape = shape)
}
