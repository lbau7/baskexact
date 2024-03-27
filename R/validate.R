# These functions cannot be accessed by the user. They are only used to
# validate the results of the main functions, by using other (mostly
# slower but easier to follow) algorithms to calculate the same values.

# Validate borrowing
val_borrow_cpp <- function(design, n, r, a, b,
                           globalweight_fun = NULL,
                           globalweight_params = list()) {
  cpp_fun <- function(x, n, a, b) {
    1 / (1 + exp(a + b * log(n^(1 / 4) * abs(x[1]/n - x[2]/n))))
  }

  shape1 <- shape2 <- numeric(length(r))
  if (!is.null(globalweight_fun)) {
    gw <- do.call(globalweight_fun, args = c(globalweight_params, n = n,
      r = list(r)))
  }

  for (i in 1:length(r)) {
    w <- sapply((1:length(r))[-i],
      function(x) cpp_fun(x = r[c(i, x)], n = n, a = a, b = b))

    if (!is.null(globalweight_fun)) {
      w <- w * gw
    }

    shape1[i] <- design@shape1 + r[i] + sum(r[-i] * w)
    shape2[i] <- design@shape2 + (n - r[i]) + sum((n - r[-i]) * w)
  }
  rbind(shape1, shape2)
}

val_borrow_fujikawa <- function(design, n, r, epsilon, tau, logbase) {
  kl_fun <- function(x, y) {
    f <- function(z) x(z) * log(x(z) / y(z), base = logbase)
    integrate(f, lower = 0, upper = 1)$value
  }

  jsd_fun <- function(sp1, sp2, n, epsilon, tau, logbase) {
    j1 <- function(x) dbeta(x, shape1 = sp1[1], shape2 = sp2[1])
    j2 <- function(x) dbeta(x, shape1 = sp1[2], shape2 = sp2[2])
    m <- function(x) (1 / 2) * (j1(x) + j2(x))
    jsd <- (1 / 2) * kl_fun(j1, m) + (1 / 2) * kl_fun(j2, m)
    w <- (1 - jsd)^epsilon
    ifelse(w <= tau, 0, w)
  }

  shape1 <- shape2 <- numeric(length(r))
  shape_prior1 <- design@shape1 + r
  shape_prior2 <- design@shape2 + (n - r)
  for (i in 1:length(r)) {
    w <- sapply((1:length(r))[-i],
      function(x) jsd_fun(sp1 = shape_prior1[c(i, x)],
        sp2 = shape_prior2[c(i, x)], n = n, epsilon = epsilon, tau = tau,
        logbase = logbase))
    shape1[i] <- shape_prior1[i] + sum(shape_prior1[-i] * w)
    shape2[i] <- shape_prior2[i] + sum(shape_prior2[-i] * w)
  }
  rbind(shape1, shape2)
}

# Loop-based calculation of the rejection probabilities of a single-stage
# basket design with 3 baskets
reject_single_loop <- function(design, p1, n, lambda, weight_fun,
                               weight_params, globalweight_fun = NULL,
                               globalweight_params = list(),
                               prob = c("toer", "pwr")) {
  targ <- get_targ(p0 = design@p0, p1 = p1, prob = prob)
  rej_ew <- 0
  rej_group <- c(0, 0, 0)
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda, globalweight_fun = globalweight_fun,
    globalweight_params = list(globalweight_params)))

  for (i1 in 0:n) {
    for (i2 in 0:n) {
      for (i3 in 0:n) {
        events <- c(i1, i2, i3)
        res <- bskt_final(design = design, n = n, lambda = lambda, r = events,
          weight_mat = weight_mat, globalweight_fun = globalweight_fun,
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

# Loop-based calculation of the rejection probabilities of a single-stage
# basket design with 4 baskets with additional use of validation functions
reject_single_loop4 <- function(design, p1, n, lambda, weightval_fun,
                                weightval_params, globalweight_fun = NULL,
                                globalweight_params = list(),
                                prob = c("toer", "pwr")) {

  targ <- get_targ(p0 = design@p0, p1 = p1, prob = prob)
  targ_ecd <- design@p0 != p1

  rej_ew <- 0
  rej_group <- c(0, 0, 0, 0)
  ecd <- 0

  for (i1 in 0:n) {
    for (i2 in 0:n) {
      for (i3 in 0:n) {
        for (i4 in 0:n) {
          events <- c(i1, i2, i3, i4)
          shape_post <- do.call(weightval_fun, args = c(weightval_params,
            globalweight_fun = globalweight_fun,
            globalweight_params = list(globalweight_params),
            design = design, n = n, r = list(events)))
          postprob <- stats::pbeta(q = design@p0, shape1 = shape_post[1, ],
            shape2 = shape_post[2, ], lower.tail = FALSE)
          res <- ifelse(postprob >= lambda, 1, 0)
          prob_temp <- get_prob(n = n, r = events, p = p1)
          ecd <- ecd + sum(res == targ_ecd) * prob_temp

          if (any(res == 1)) {
            rej_group[which(res == 1)] <- rej_group[which(res == 1)] +
              prob_temp
            if (any(res[targ] == 1)) {
              rej_ew <- rej_ew + prob_temp
          }
        }
      }
    }
  }
  }
  if (prob == "toer") {
    list(
      rejection_probabilities = rej_group,
      fwer = rej_ew,
      ecd = ecd
    )
  } else {
    list(
      rejection_probabilities = rej_group,
      ewp = rej_ew,
      ecd = ecd
    )
  }
}

# Loop-based calculation of the rejection probabilities of a two-stage
# basket design with 3 baskets
reject_twostage_loop <- function(design, p1, n, n1, lambda, interim_fun,
                                 interim_params = list(), weight_fun,
                                 weight_params = list(),
                                 globalweight_fun = NULL,
                                 globalweight_params = list(),
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
          weight_mat = list(weight_mat), globalweight_fun = globalweight_fun,
          globalweight_params = list(globalweight_params)))
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
          # Consider significant interim results
          if (sum(res_int == 0) == 1) {
            for (i4 in 0:(n - n1)) {
              events2 <- numeric(3)
              events2[which(res_int == 0)] <- i4
              events_sum <- events1 + events2
              res_fin <- bskt_final_int(design = design, n = n, n1 = n1,
                r = events_sum, res_int = res_int, lambda = lambda,
                weight_mat = weight_mat, globalweight_fun = globalweight_fun,
                globalweight_params = globalweight_params)
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
                  weight_mat = weight_mat, globalweight_fun = globalweight_fun,
                  globalweight_params = globalweight_params)
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
                    weight_mat = weight_mat,
                    globalweight_fun = globalweight_fun,
                    globalweight_params = globalweight_params)
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

# Loop-based calculation of the mean posterior means of a single-stage
# basket design with 3 baskets
estim_loop <- function(design, p1, n, lambda = NULL, weight_fun,
                       weight_params, globalweight_fun = NULL,
                       globalweight_params = list()) {
  post_means_w <- c(0, 0, 0)
  mse_w <- c(0, 0, 0)
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda, globalweight_fun = globalweight_fun,
    globalweight_params = list(globalweight_params)))

  for (i1 in 0:n) {
    for (i2 in 0:n) {
      for (i3 in 0:n) {
        events_loop <- c(i1, i2, i3)
        prob <- get_prob(n = n, r = events_loop, p = p1)
        bshape <- beta_borrow(weight_mat = weight_mat,
          globalweight_fun = globalweight_fun,
          globalweight_params = globalweight_params, design = design,
          n = n, r = events_loop)
        bmean <- mean_beta(bshape)
        mse <- (bmean - p1)^2
        post_means_w <- post_means_w + bmean * prob
        mse_w <- mse_w + mse * prob
      }
    }
  }
  list(
    Mean = post_means_w,
    MSE = mse_w
  )
}

# Loop-based calculation of the mean posterior means of a two-stage
# basket design with 3 baskets
estim_twostage_loop <- function(design, p1, n, n1, lambda, interim_fun,
                                interim_params = list(), weight_fun,
                                weight_params = list(), globalweight_fun = NULL,
                                globalweight_params = list()) {
  post_means_w <- c(0, 0, 0)
  mse_w <- c(0, 0, 0)
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, n1 = n1, lambda = lambda))
  for (i1 in 0:n1) {
    for (i2 in 0:n1) {
      for (i3 in 0:n1) {
        events1 <- c(i1, i2, i3)
        res_int <- do.call(interim_fun, args = c(interim_params,
          design = design, n = n, n1 = n1, r1 = list(events1), lambda = lambda,
          weight_mat = list(weight_mat), globalweight_fun = globalweight_fun,
          globalweight_params = list(globalweight_params)))
        if (all(res_int %in% c(-1, 1))) {
          prob <- get_prob(n = n1, r = events1, p = p1)
          bshape <- beta_borrow(weight_mat = weight_mat,
            globalweight_fun = globalweight_fun,
            globalweight_params = globalweight_params, design = design,
            n = n1, r = events1)
          bmean <- mean_beta(bshape)
          mse <- (bmean - p1)^2
          post_means_w <- post_means_w + bmean * prob
          mse_w <- mse_w + mse * prob
        } else {
          if (sum(res_int == 0) == 1) {
            for (i4 in 0:(n - n1)) {
              events2 <- numeric(3)
              events2[which(res_int == 0)] <- i4
              events_sum <- events1 + events2

              prob1 <- get_prob(n = n1, r = events1, p = p1)
              prob2 <- get_prob(n = n - n1, r = i4,
                p = p1[which(res_int == 0)])
              bshape <- beta_borrow_int(weight_mat = weight_mat,
                globalweight_fun = globalweight_fun,
                globalweight_params = globalweight_params, design = design,
                n = n, n1 = n1, r = events_sum, res_int = res_int)
              bmean <- mean_beta(bshape)
              mse <- (bmean - p1)^2
              post_means_w <- post_means_w + bmean * prob1 * prob2
              mse_w <- mse_w + mse * prob1 * prob2
            }
          } else if (sum(res_int == 0) == 2) {
            for (i4 in 0:(n - n1)) {
              for (i5 in 0:(n - n1)) {
                events2 <- numeric(3)
                events2[which(res_int == 0)] <- c(i4, i5)
                events_sum <- events1 + events2

                prob1 <- get_prob(n = n1, r = events1, p = p1)
                prob2 <- get_prob(n = n - n1, r = c(i4, i5),
                  p = p1[which(res_int == 0)])
                bshape <- beta_borrow_int(weight_mat = weight_mat,
                  globalweight_fun = globalweight_fun,
                  globalweight_params = globalweight_params, design = design,
                  n = n, n1 = n1, r = events_sum, res_int = res_int)
                bmean <- mean_beta(bshape)
                mse <- (bmean - p1)^2
                post_means_w <- post_means_w + bmean * prob1 * prob2
                mse_w <- mse_w + mse * prob1 * prob2
              }
            }
          } else {
            for (i4 in 0:(n - n1)) {
              for (i5 in 0:(n - n1)) {
                for (i6 in 0:(n - n1)) {
                  events2 <- c(i4, i5, i6)
                  events_sum <- events1 + events2

                  prob1 <- get_prob(n = n1, r = events1, p = p1)
                  prob2 <- get_prob(n = n - n1, r = events2, p = p1)
                  bshape <- beta_borrow_int(weight_mat = weight_mat,
                    globalweight_fun = globalweight_fun,
                    globalweight_params = globalweight_params, design = design,
                    n = n, n1 = n1, r = events_sum, res_int = res_int)
                  bmean <- mean_beta(bshape)
                  mse <- (bmean - p1)^2
                  post_means_w <- post_means_w + bmean * prob1 * prob2
                  mse_w <- mse_w + mse * prob1 * prob2
                }
              }
            }
          }
        }
      }
    }
  }
  list(
    Mean = post_means_w,
    MSE = mse_w
  )
}

# Loop-based version of check_mon_within
mon_within_loop <- function(design, n, lambda, weight_fun, weight_params,
                            globalweight_fun = NULL,
                            globalweight_params = list()) {
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda, globalweight_fun = globalweight_fun,
    globalweight_params = list(globalweight_params)))

  events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
  func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
    r = x, weight_mat = weight_mat, globalweight_fun = globalweight_fun,
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
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda, globalweight_fun = globalweight_fun,
    globalweight_params = list(globalweight_params)))

  events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
  func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
    r = x, weight_mat = weight_mat, globalweight_fun = globalweight_fun,
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

# Loop-based version of ecd for single-stage design
ecd_loop <- function(design, p1, n, lambda, weight_fun,
                     weight_params = list(), globalweight_fun = NULL,
                     globalweight_params = list()) {
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, lambda = lambda, globalweight_fun = globalweight_fun,
    globalweight_params = list(globalweight_params)))
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

# Loop-based version of ecd for two-stage design
ecd_twostage_loop <- function(design, p1, n, n1, lambda, interim_fun,
                              interim_params, weight_fun, weight_params,
                              globalweight_fun = NULL,
                              globalweight_params = list()) {
  weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
    n = n, n1 = n1, lambda = lambda, globalweight_fun = globalweight_fun,
    globalweight_params = list(globalweight_params)))

  events <- arrangements::permutations(0:n, k = design@k, replace = TRUE)
  targ <- get_targ(p0 = design@p0, p1 = p1, prob = "pwr")
  cd <- 0

  for (i1 in 0:n1) {
    for (i2 in 0:n1) {
      for (i3 in 0:n1) {
        events1 <- c(i1, i2, i3)
        res_int <- do.call(interim_fun, args = c(interim_params,
          design = design, n = n, n1 = n1, r1 = list(events1), lambda = lambda,
          weight_mat = list(weight_mat), globalweight_fun = globalweight_fun,
          globalweight_params = list(globalweight_params)))
        if (sum(res_int) == -design@k) {
          cd <- cd + sum(c(0, 0, 0) == targ) *
            get_prob(n = n1, r = events1, p = p1)
        } else if (all(res_int %in% c(-1, 1))) {
          res_int <- ifelse(res_int == -1, 0, res_int)
          cd <- cd + sum(res_int == targ) * get_prob(n = n1, r = events1,
            p = p1)
        } else {
          if (sum(res_int == 0) == 1) {
            for (i4 in 0:(n - n1)) {
              events2 <- numeric(3)
              events2[which(res_int == 0)] <- i4
              events_sum <- events1 + events2
              res_fin <- bskt_final_int(design = design, n = n, n1 = n1,
                r = events_sum, res_int = res_int, lambda = lambda,
                weight_mat = weight_mat, globalweight_fun = globalweight_fun,
                globalweight_params = globalweight_params)
              res_fin <- ifelse(res_int == 1, 1, res_fin)
              res_fin <- ifelse(res_fin == -1, 0, res_fin)
              cd <- cd + sum(res_fin == targ) * get_prob(n = n1, r = events1,
                p = p1) * get_prob(n = (n - n1), r = i4,
                  p = p1[which(res_int == 0)])
            }
          } else if (sum(res_int == 0) == 2) {
            for (i4 in 0:(n - n1)) {
              for (i5 in 0:(n - n1)) {
                events2 <- numeric(3)
                events2[which(res_int == 0)] <- c(i4, i5)
                events_sum <- events1 + events2
                res_fin <- bskt_final_int(design = design, n = n, n1 = n1,
                  r = events_sum, res_int = res_int, lambda = lambda,
                  weight_mat = weight_mat, globalweight_fun = globalweight_fun,
                  globalweight_params = globalweight_params)
                res_fin <- ifelse(res_int == 1, 1, res_fin)
                res_fin <- ifelse(res_fin == -1, 0, res_fin)
                cd <- cd + sum(res_fin == targ) * get_prob(n = n1, r = events1,
                  p = p1) * get_prob(n = (n - n1), r = c(i4, i5),
                    p = p1[which(res_int == 0)])
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
                    weight_mat = weight_mat,
                    globalweight_fun = globalweight_fun,
                    globalweight_params = globalweight_params)
                  res_fin <- ifelse(res_int == 1, 1, res_fin)
                  res_fin <- ifelse(res_fin == -1, 0, res_fin)
                  cd <- cd + sum(res_fin == targ) * get_prob(n = n1, r = events1,
                    p = p1) * get_prob(n = (n - n1), r = c(i4, i5, i6), p = p1)
                }
              }
            }
          }
        }
      }
    }
  }
  cd
}
