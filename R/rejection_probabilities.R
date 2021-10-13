# Calculates the experimentwise rejection probability for a single-stage design
reject_prob_ew <- function(design, theta1, n, lambda, weight_mat,
                           prob = c("toer", "pwr")) {
  # Computational shortcuts don't work with unequal priors or n!
  targ <- get_targ(theta0 = design@theta0, theta1 = theta1, prob = prob)
  # Create matrix with all possible outcomes (without permutations)
  events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)

  # Remove outcomes when no significant results are possible
  crit <- get_crit(design = design, n = n, lambda = lambda)
  crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
  sel_nosig <- apply(events, 1, function(x) all(x < crit_pool))
  events <- events[!sel_nosig, ]

  # Remove outcomes where all results are significant and calculate the
  # probabilities later
  sel_sig <- apply(events, 1, function(x) all(x >= crit))
  events_sig <- events[sel_sig, ]
  events <- events[!sel_sig, ]

  # Conduct test for the remaining outcomes
  fun <- function(x) bskt_final(design = design, n = n, lambda = lambda, r = x,
    weight_mat = weight_mat)
  res <- t(apply(events, 1, fun))

  # Select outcomes with at least one rejected null hypothesis
  # and the corresponding results (including events_sig saved before)
  eff_vec <- apply(res, 1, function(x) any(x == 1))
  events_eff <- rbind(events[eff_vec, ], events_sig)
  res_eff <- rbind(
    res[eff_vec, ],
    matrix(1, nrow = nrow(events_sig), ncol = design@k)
  )

  # If all theta1 are equal each permutation has the same probability
  if ((sum(targ) == design@k) & (length(unique(theta1)) == 1)) {
    probs_eff <- apply(events_eff, 1,
      function(x) get_prob(n = n, r = x, theta = theta1))
    # Helper function that calculates the number of permutations
    perm_fun <- function(x) {
      tab <- tabulate(x + 1)
      tab <- tab[tab != 0]
      ifelse(length(unique(x)) == 1, 1,
        arrangements::npermutations(x = unique(x), freq = tab))
    }
    eff_perm <- apply(events_eff, 1, perm_fun)
    rej_prob <- sum(probs_eff * eff_perm)
  } else {
    # If not all theta1 are equal calculate probability for each permutation
    rej_prob <- numeric(nrow(res_eff))
    for (i in 1:nrow(res_eff)) {
      # If number of responses is equal in each basket, each permutation
      # has the same probability even when not all theta1 are equal
      if (length(unique(events_eff[i, ])) == 1) {
        rej_prob[i] <- get_prob(n = n, r = events_eff[i, ],
          theta = theta1)
      } else {
      events_loop <- arrangements::permutations(events_eff[i, ])
      res_loop <- arrangements::permutations(
        res_eff[i, ])[!duplicated(events_loop), ]
      events_loop <- events_loop[!duplicated(events_loop), ]
      eff_loop <- apply(res_loop, 1, function(x) any(x[targ] == 1))
      events_loop <- events_loop[eff_loop, , drop = FALSE]
      rej_prob[i] <- sum(apply(events_loop, 1, function(x) get_prob(n = n,
        r = x, theta = theta1)))
      }
    }
    rej_prob <- sum(rej_prob)
  }
  rej_prob
}

# Calculates the groupwise rejection probabilities for a two-stage design
reject_prob_group <- function(design, theta1, n, lambda, weight_mat,
  prob = c("toer", "pwr")) {
  targ <- get_targ(theta0 = design@theta0, theta1 = theta1, prob = prob)
  # Create matrix with all possible outcomes
  events <- arrangements::permutations(0:n, k = design@k, replace = TRUE)

  # Conduct test for all possible outcomes
  fun <- function(x) bskt_final(design = design, n = n, lambda = lambda, r = x,
    weight_mat = weight_mat)
  res <- t(apply(events, 1, fun))

  eff_vec <- apply(res, 1, function(x) any(x == 1))
  eff_vec_targ <- apply(res[eff_vec, ], 1, function(x) any(x[targ] == 1))
  events_eff <- events[eff_vec, ]
  # Calculate probability of ouctomes where any null hypothesis was rejected
  probs_eff <- apply(events_eff, 1,
    function(x) get_prob(n = n, r = x, theta = theta1))
  res_eff <- res[eff_vec,]
  rej <- colSums(apply(res_eff == 1, 2, function(x) x * probs_eff))
  # Use only the probabilities of outcomes with a rejected null hypothesis
  # where a targeted basket was significant to calculate experimentwise
  # rejection probability
  rej_ew <- sum(probs_eff[eff_vec_targ])

  if (prob == "toer") {
    list(
      rejection_probabilities = rej,
      fwer = rej_ew
    )
  } else {
    list(
      rejection_probabilities = rej,
      ewp = rej_ew
    )
  }
}

# Calculates the experimentwise rejection probability for a two-stage design
reject_prob_ew2 <- function(design, theta1, n, n1, lambda, weight_mat,
                            prob = c("toer", "pwr")) {
  targ <- get_targ(theta0 = theta0, theta1 = theta1, prob = prob)
  if (sum(targ) == k) {
    events_int <- arrangements::combinations(0:n1, k = k, replace = TRUE)
  } else {
    events_int <- arrangements::permutations(x = 0:n1, k = k, replace = TRUE)
  }

  int_fun <- function(x) do.call(interim_fun, args = c(interim_params, n = n,
    n1 = n1, r1 = x, weight_mat = weight_mat))
  res_int <- t(apply(events_int, 1, int_fun))

  # Calculate probability for events where at least one group
  # was stopped for efficacy
  eff_vec <- apply(res_int, 1, function(x) any(x[targ] == 1))
  events_eff <- events_int[eff_vec, ]
  probs_eff <- apply(events_eff, 1,
    function(x) get_prob(n = n1, r = x, theta = theta1))

  if (sum(targ) == k) {
    # Multiply probability by number of permutations for each row
    perm_fun <- function(x) ifelse(length(unique(x)) == 1, 1,
      arrangements::npermutations(x = unique(x), freq = table(x)))
    eff_perm <- apply(events_eff, 1, perm_fun)
    rej_prob <- sum(probs_eff * eff_perm)
  } else {
    rej_prob <- sum(probs_eff)
  }

  if (interim) {
    # Continue only when there is no relevant significant result and
    # there is at least one basket left in the second stage
    cont_vec <- !eff_vec
    cont_vec[which(cont_vec)] <- apply(res_int[cont_vec, ], 1,
      function(x) any(x == 0))
    events_cont <- events_int[cont_vec, , drop = FALSE]
    res_cont <- res_int[cont_vec, , drop = FALSE]

    all_events_stg2 <- lapply(1:k, function(x)
      arrangements::permutations(x = 0:(n - n1), k = x, replace = TRUE))

    for (i in 1:nrow(events_cont)) {
      no_cont <- sum(res_cont[i, ] == 0)
      # Event-Matrix für die Gruppen die in Stage 2 gehen
      events_stg2 <- all_events_stg2[[no_cont]]
      if (no_cont < k) {
        events_loop <- matrix(0, nrow = nrow(events_stg2), ncol = k)
        events_loop[, res_cont[i, ] == 0] <- events_stg2
      } else {
        events_loop <- events_stg2
      }

      events_fin <- t(t(events_loop) + events_cont[i, ])
      fin_func <- function(x) bskt_final(k = k, n = n, n1 = n1, r = x,
        shape1 = shape1, shape2 = shape2, weight_mat = weight_mat,
        res_int = res_cont[i, ], theta_0 = theta_0, lambda = lambda)

      fin_res <- t(apply(events_fin, 1, fin_func))
      sig_vec <- apply(fin_res, 1, function(x) any(x[targ] == 1))

      prob_cont <- get_prob(n = n1, r = events_cont[i, ], theta = theta_1)
      prob_sig <- apply(events_stg2[sig_vec, , drop = FALSE], 1, function(x)
        get_prob(n = n - n1, r = x, theta = theta_1[res_cont[i, ] == 0]))
      rej_prob_temp <- prob_cont * sum(prob_sig)

      if ((rej_prob_temp > 0) & (length(unique(events_cont[i, ])) > 1) &
          (sum(targ) == k)) {
        perm_temp <- arrangements::npermutations(
          x = unique(events_cont[i, ]), freq = table(events_cont[i, ]))
        rej_prob_temp <- rej_prob_temp * perm_temp
      }
      rej_prob <- rej_prob + rej_prob_temp
    }
  }
  rej_prob
}

reject_prob_group <- function(k, n, n1, shape1 = 1, shape2 = 1, theta_0,
  theta_1 = NULL, lambda, interim,
  weight_mat = NULL, crit,
  prob = c("toer", "pwr")) {
  if (is.null(n1)) n1 <- n
  targ <- get_targ(theta_0 = theta_0, theta_1 = theta_1, prob = prob)
  events_int <- arrangements::permutations(x = 0:n1, k = k, replace = TRUE)

  if (interim) {
    int_fun <- function(x) bskt_interim(k = k, n = n, n1 = n1, r1 = x,
      shape1 = shape1, shape2 = shape2, crit = crit, weight_mat = weight_mat)
  } else {
    int_fun <- function(x) bskt_final_noint(k = k, n = n, r = x,
      shape1 = shape1, shape2 = shape2, weight_mat = weight_mat,
      theta_0 = theta_0, lambda = lambda)
  }

  res_int <- t(apply(events_int, 1, int_fun))
  eff_vec <- apply(res_int, 1, function(x) any(x == 1))
  eff_vec_targ <- apply(res_int[eff_vec, ], 1, function(x) any(x[targ] == 1))
  events_eff <- events_int[eff_vec, ]
  probs_eff <- apply(events_eff, 1,
    function(x) get_prob(n = n1, r = x, theta = theta_1))
  res_eff <- res_int[eff_vec,]
  rej <- colSums(apply(res_eff == 1, 2, function(x) x * probs_eff))
  rej_ew <- sum(probs_eff[eff_vec_targ])

  if (interim) {
    cont_vec <- apply(res_int, 1, function(x) any(x == 0))
    events_cont <- events_int[cont_vec, , drop = FALSE]
    res_cont <- res_int[cont_vec, , drop = FALSE]

    all_events_stg2 <- lapply(1:k, function(x)
      arrangements::permutations(x = 0:(n - n1), k = x, replace = TRUE))

    for (i in 1:nrow(events_cont)) {
      no_cont <- sum(res_cont[i, ] == 0)
      events_stg2 <- all_events_stg2[[no_cont]]
      if (no_cont < k) {
        events_loop <- matrix(0, nrow = nrow(events_stg2), ncol = k)
        events_loop[, res_cont[i,] == 0] <- events_stg2
      } else {
        events_loop <- events_stg2
      }

      events_fin <- t(t(events_loop) + events_cont[i, ])
      fin_func <- function(x) bskt_final(k = k, n = n, n1 = n1, r = x,
        shape1 = shape1, shape2 = shape2, weight_mat = weight_mat,
        res_int = res_cont[i, ], theta_0 = theta_0, lambda = lambda)
      res_fin <- t(apply(events_fin, 1, fin_func))

      prob_cont <- get_prob(n = n1, r = events_cont[i, ], theta = theta_1)
      prob_stg2 <- apply(events_stg2, 1, function(x)
        get_prob(n = n - n1, r = x, theta = theta_1[res_cont[i, ] == 0]))
      rej <- rej + colSums(apply(res_fin == 1, 2, function(x) x * prob_stg2)) *
        prob_cont

      # Ignore results which were already counted for rej_ew after the
      # interim analysis
      if (all(res_cont[i, ][targ] != 1)) {
        rej_ew <- rej_ew + sum(apply(res_fin, 1, function(x) any(x[targ] == 1)) *
            prob_stg2) * prob_cont
      }
    }
  }

  if (prob == "toer") {
    list(
      rejection_probabilities = rej,
      fwer = rej_ew
    )
  } else {
    list(
      rejection_probabilities = rej,
      ewp = rej_ew
    )
  }
}

reject_prob_sim <- function(k, n, n1, shape1 = 1, shape2 = 1, theta_0, theta_1,
  lambda, interim, weight_mat, crit,
  prob = c("toer", "pwr"), results = c("ew", "group"),
  data = data) {
  targ <- get_targ(theta_0 = theta_0,
    theta_1 = theta_1, prob = prob)
  if (results == "group") {
    rej <- matrix(0, nrow = nrow(data), ncol = k)
  }
  rej_ew <- numeric(nrow(data))

  for (i in 1:nrow(data)) {
    if (interim) {
      events_int <- data[i, 1:k]
      res_int <- bskt_interim(k = k, n = n, n1 = n1, r1 = events_int,
        shape1 = shape1, shape2 = shape2, crit = crit,
        weight_mat = weight_mat)
      if (all(res_int == -1)) {
        next
      }
      if (any(res_int == 1)) {
        if (results == "group") {
          rej[i, ] <- res_int == 1
        }
        rej_ew[i] <- any(res_int[targ] == 1)
        if ((results == "ew") & (rej_ew[i] == 1)) {
          next
        }
      }
      select_vec <- ((k + 1):(2 * k))[which(res_int == 0)]
      events_stg2 <- rep(0, k)
      events_stg2[which(res_int == 0)] <- data[i, select_vec]
      events_fin <- events_int + events_stg2
      res_fin <- bskt_final(k = k, n = n, n1 = n1, r = events_fin,
        shape1 = shape1, shape2 = shape2, weight_mat = weight_mat,
        res_int = res_int, theta_0 = theta_0, lambda = lambda)
      if (any(res_fin == 1)) {
        if (results == "group") {
          rej[i, ] <- pmax(rej[i, ], res_fin == 1)
        }
        rej_ew[i] <- max(rej_ew[i], any(res_fin[targ] == 1))
      }
    } else {
      events <- data[i, 1:k]
      res <- bskt_final_noint(k = k, n = n, r = events, shape1 = shape1,
        shape2 = shape2, weight_mat = weight_mat, theta_0 = theta_0,
        lambda = lambda)
      if (results == "group") {
        rej[i, ] <- res == 1
      }
      rej_ew[i] <- any(res[targ] == 1)
    }
  }

  if (prob == "toer") {
    if (results == "group") {
      list(
        rejection_probabilities = colMeans(rej),
        fwer = mean(rej_ew)
      )
    } else {
      mean(rej_ew)
    }
  } else {
    if (results == "group") {
      list(
        rejection_probabilities = colMeans(rej),
        ewp = mean(rej_ew)
      )
    } else {
      mean(rej_ew)
    }
  }
}

