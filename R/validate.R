# These functions cannot be accessed by the user. They are only used to
# validate the results of the main functions, by using other (mostly
# slower but easier to follow) algorithms to calculate the same values.

# Loop-based calculation of the rejection probabilities of a single-stage
# basket design with 3 baskets
reject_single_loop <- function(design, n, lambda, epsilon, tau,
                               logbase = exp(1), prune = FALSE,
                               prob = c("toer", "pwr")) {
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

        if (any(res == 1)) {
          prob_temp <- get_prob(n = n, r = events, theta = design@theta1)
          rej_group[which(res == 1)] <- rej_group[which(res == 1)] +
            prob_temp

          if (any(res[targ] == 1)) {
            rej_ew <- rej_ew + prob_temp
          }
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

# Loop-based version of check_mon_within
mon_within_loop <- function(design, n, lambda, epsilon, tau, logbase = 2,
                            prune, ...) {
  weights <- get_weights(design = design, n = n, epsilon = epsilon, tau = tau,
    logbase = logbase)

  if (prune) {
    crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
    weights <- prune_weights(weight_mat = weights, cut = crit_pool)
  }

  events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
  func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
    r = x, weight_mat = weights)

  viol <- c()
  for (i in 1:nrow(events)) {
    res_loop <- func(events[i, ])
    if (any(res_loop != cummax(res_loop))) viol <- rbind(viol, events[i, ])
  }

  if (length(viol) == 0) {
    TRUE
  } else {
    viol
  }
}



