# These functions cannot be accessed by the user. They are only used to
# validate the results of the main functions, by using other (mostly
# slower but easier to follow) algorithms to calculate the same values.

# Loop-based calculation of the rejection probabilities of a single-stage
# basket design with 3 baskets
reject_single_loop <- function(design, p1, n, lambda, weight_fun,
                               weight_params, globalweight_fun = NULL,
                               globalweight_params = list(),
                               prob = c("toer", "pwr")) {
  targ <- get_targ(p0 = design@p0, p1 = p1, prob = prob)
  rej_ew <- 0
  rej_group <- c(0, 0, 0)
  weights <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda))

  for (i1 in 0:n) {
    for (i2 in 0:n) {
      for (i3 in 0:n) {
        events <- c(i1, i2, i3)
        res <- bskt_final(design = design, n = n, lambda = lambda, r = events,
          weight_mat = weights, globalweight_fun = globalweight_fun,
          globalweight_params = globalweight_params)

        if (any(res == 1)) {
          prob_temp <- get_prob(n = n, r = events, p = p1)
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

# Loop-based calculation of the rejection probabilities of a two-stage
# basket design with 3 baskets
reject_twostage_loop <- function(design, p1, n, n1, lambda, interim_fun,
                                 interim_params = list(), weight_fun,
                                 weight_params = list(),
                                 prob = c("toer", "pwr")) {
  targ <- get_targ(p0 = design@p0, p1 = p1, prob = prob)
  rej_ew <- 0
  rej_group <- c(0, 0, 0)
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, n1 = n1, lambda = lambda))
  for (i1 in 0:n1) {
    for (i2 in 0:n1) {
      for (i3 in 0:n1) {
        events1 <- c(i1, i2, i3)
        res_int <- do.call(interim_fun, args = c(interim_params,
          design = design, n = n, n1 = n1, r1 = list(events1), lambda = lambda,
          weight_mat = list(weight_mat)))
        if (sum(res_int) == -design@k) {
          next
        } else if (all(res_int %in% c(-1, 1))) {
          prob_temp <- get_prob(n = n1, r = events1, p = p1)
          rej_group[which(res_int == 1)] <- rej_group[which(res_int == 1)] +
            prob_temp
          if (any(res_int[targ] == 1)) {
            rej_ew <- rej_ew + prob_temp
          }
        } else {
          if (any(res_int == 1)) {
            prob_temp <- get_prob(n = n1, r = events1, p = p1)
            rej_group[which(res_int == 1)] <- rej_group[which(res_int == 1)] +
              prob_temp
            if (any(res_int[targ] == 1)) {
              rej_ew <- rej_ew + prob_temp
            }
          }
          # Muss signifikante Interimsergebnise berücksichtigen
          if (sum(res_int == 0) == 1) {
            for (i4 in 0:(n - n1)) {
              events2 <- numeric(3)
              events2[which(res_int == 0)] <- i4
              events_sum <- events1 + events2
              res_fin <- bskt_final_int(design = design, n = n, n1 = n1,
                r = events_sum, res_int = res_int, lambda = lambda,
                weight_mat = weight_mat)
              if (any(res_fin == 1)) {
                prob_temp <- get_prob(n = n1, r = events1, p = p1) *
                  get_prob(n = (n - n1), r = i4,
                    p = p1[which(res_int == 0)])
                rej_group[which(res_fin == 1)] <-
                  rej_group[which(res_fin == 1)] + prob_temp
                if (any(res_fin[targ] == 1) & all(res_int[targ] != 1)) {
                  rej_ew <- rej_ew + prob_temp
                }
              }
            }
          } else if (sum(res_int == 0) == 2) {
            for (i4 in 0:(n - n1)) {
              for (i5 in 0:(n - n1)) {
                events2 <- numeric(3)
                events2[which(res_int == 0)] <- c(i4, i5)
                events_sum <- events1 + events2
                res_fin <- bskt_final_int(design = design, n = n, n1 = n1,
                  r = events_sum, res_int = res_int, lambda = lambda,
                  weight_mat = weight_mat)
                if (any(res_fin == 1)) {
                  prob_temp <- get_prob(n = n1, r = events1, p = p1) *
                    get_prob(n = (n - n1), r = c(i4, i5),
                      p = p1[which(res_int == 0)])
                  rej_group[which(res_fin == 1)] <-
                    rej_group[which(res_fin == 1)] + prob_temp
                  if (any(res_fin[targ] == 1) & all(res_int[targ] != 1)) {
                    rej_ew <- rej_ew + prob_temp
                  }
                }
              }
            }
          } else {
            for (i4 in 0:(n - n1)) {
              for (i5 in 0:(n - n1)) {
                for (i6 in 0:(n - n1)) {
                  events2 <- c(i4, i5, i6)
                  events_sum <- events1 + events2
                  res_fin <- bskt_final_int(design = design, n = n, n1 = n1,
                    r = events_sum, res_int = res_int, lambda = lambda,
                    weight_mat = weight_mat)
                  if (any(res_fin == 1)) {
                    prob_temp <- get_prob(n = n1, r = events1, p = p1) *
                      get_prob(n = (n - n1), r = c(i4, i5, i6),
                        p = p1)
                    rej_group[which(res_fin == 1)] <-
                      rej_group[which(res_fin == 1)] + prob_temp
                    if (any(res_fin[targ] == 1) & all(res_int[targ] != 1)) {
                      rej_ew <- rej_ew + prob_temp
                    }
                  }
                }
              }
            }
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
mon_within_loop <- function(design, n, lambda, weight_fun, weight_params,
                            globalweight_fun = NULL,
                            globalweight_params = list()) {
  weights <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda))

  events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
  func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
    r = x, weight_mat = weights, globalweight_fun = globalweight_fun,
    globalweight_params = globalweight_params)

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

# Loop-based version of check_mon_between without shortcuts
mon_between_loop <- function(design, n, lambda, weight_fun, weight_params,
                             globalweight_fun = NULL,
                             globalweight_params = list()) {
  weights <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda))

  events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
  func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
    r = x, weight_mat = weights, globalweight_fun = globalweight_fun,
    globalweight_params = globalweight_params)

  res <- numeric(nrow(events))
  for (i in 1:nrow(events)) {
    res_loop <- func(events[i, ])
    res[i] <- any(res_loop == 1)
  }

  viol <- c()
  for (i in 1:nrow(events)) {
    if (res[i]) {
      events_sel <- apply(events, 1, function(x) all(x >= events[i, ]))
      res_sel <- res[events_sel]
      check <- sum(res_sel) == length(res_sel)
      if (!check) viol <- rbind(viol, events[i, ])
    }
  }

  if (length(viol) == 0) {
    TRUE
  } else {
    viol
  }
}

# Loop-based version of ecd
ecd_loop <- function(design, p1, n, lambda, weight_fun,
                     weight_params = list(), globalweight_fun = NULL,
                     globalweight_params = list()) {
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda))
  events <- arrangements::permutations(0:n, k = design@k, replace = TRUE)
  targ <- get_targ(p0 = design@p0, p1 = p1, prob = "pwr")

  cd <- prob <- numeric(nrow(events))
  for (i in 1:nrow(events)) {
    res_loop <- bskt_final(design = design, n = n, lambda = lambda,
      r = events[i, ], weight_mat = weight_mat,
      globalweight_fun = globalweight_fun,
      globalweight_params = globalweight_params)
    cd[i] <- sum(res_loop == targ)
    prob[i] <- get_prob(n = n, r = events[i, ], p = p1)
  }

  sum(cd * prob)
}






