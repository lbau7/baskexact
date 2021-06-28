# These functions cannot be accessed by the user. They are only used to
# validate the results of the main functions, by using other (mostly
# slower but easier to follow) algorithms to calculate the same values.

# Loop-based calculation of the rejection probabilities of a single-stage
# basket design with 3 baskets
reject_single_loop <- function(design, n, lambda, epsilon, tau,
                               logbase = exp(1), prune = FALSE,
                               prob = c("toer", "pwr")) {
  if (all(design@theta1 != design@theta0)) {
    design@theta1 <- rep(design@theta0, design@k)
  }

  targ <- get_targ(design = design, prob = prob)
  rej_ew <- 0
  rej_group <- c(0, 0, 0)
  weights <- get_weights(design = design, n = n, epsilon = epsilon, tau = tau,
    logbase = logbase)

  if (prune) {
    crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
    weights <- prune_weights(weight_mat = weights, cut = crit_pool)
  }

  for (i1 in 0:n) {
    for (i2 in 0:n) {
      for (i3 in 0:n) {
        events <- c(i1, i2, i3)
        res <- bskt_final(design = design, n = n, lambda = lambda, r = events,
          weight_mat = weights)

        if (any(res[targ] == 1)) {
          prob_temp <- get_prob(n = n, r = events, theta = design@theta1)
          rej_ew <- rej_ew + prob_temp
          rej_group[which(res[targ] == 1)] <- rej_group[which(res[targ] == 1)] +
            prob_temp
        }
      }
    }
  }

  if (prob == "toer") {
    list(
      rejection_probabilities = rej_group,
      fwer = rej_ew
    )
  } else {
    list(
      rejection_probabilities = rej_group,
      ewp = rej_ew
    )
  }
}


