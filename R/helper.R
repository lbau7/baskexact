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

# Returns a vector that determines which baskets are of interest
# to compute the type 1 error rate or the power
get_targ <- function(theta0, theta1, prob) {
  if (prob == "toer") {
    theta0 == theta1
  } else {
    theta0 != theta1
  }
}

# Calculate probability for an event to occur
get_prob <- function(n, r, theta) {
  prod(stats::dbinom(x = r, size = n, prob = theta))
}

# Calculate the posterior probability
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

# Computes the posterior predictive probability
post_pred <- function(n, n1, r1, shape, crit) {
  extraDistr::pbbinom(
    q = crit - r1 - 1,
    size = n - n1,
    alpha = shape[1, ],
    beta = shape[2, ],
    lower.tail = FALSE
  )
}

# Calculates the number of permutations of results
get_permutations <- function(x) {
  tab <- tabulate(x + 1)
  tab <- tab[tab != 0]
  ifelse(length(unique(x)) == 1, 1,
    arrangements::npermutations(x = unique(x), freq = tab))
}
