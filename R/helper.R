# Get the critical value for individual baskets
get_crit <- function(design, n, lambda) {
  shape1post <- design@shape1 + 1:n
  shape2post <- design@shape1 + n - 1:n
  betafun <- function(x, y) 1 - stats::pbeta(design@theta0, x, y)
  prob <- mapply(betafun, shape1post, shape2post)
  which(prob >= lambda)[1]
}

# Get the critical value with borrowing when the number of responses
# is equal in all baskets
get_crit_pool <- function(design, n, lambda) {
  shape1pool <- (design@shape1 + 1:n) * design@k
  shape2pool <- (design@shape1 + n - 1:n) * design@k
  betafun <- function(x, y) 1 - stats::pbeta(design@theta0, x, y)
  prob <- mapply(betafun, shape1pool, shape2pool)
  which(prob >= lambda)[1]
}

get_targ <- function(design, prob) {
  if (prob == "toer") {
    design@theta0 == design@theta1
  } else {
    design@theta0 != design@theta1
  }
}

get_prob <- function(n, r, theta) {
  prod(stats::dbinom(x = r, size = n, prob = theta))
}

post_beta <- function(shape, theta0) {
  stats::pbeta(theta0, shape1 = shape[1, ], shape2 = shape[2, ],
    lower.tail = FALSE)
}

# Prevents that information is borrowed from baskets with
# less than cut events
prune_weights <- function(weight_mat, cut) {
  weight_mat[0:cut, ] <- 0
  weight_mat[, 0:cut] <- 0
  weight_mat
}
