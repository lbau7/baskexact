# Calculates expected number of correct decisions for a single-stage design
ecd_calc <- function(design, theta1, n, lambda, weight_mat, globalweight_fun,
                     globalweight_params) {
  targ <- get_targ(theta0 = design@theta0, theta1 = theta1, prob = "pwr")

  # Create matrix with all possible outcomes (without permutations)
  events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)

  # Remove outcomes when no significant results are possible
  crit <- get_crit(design = design, n = n, lambda = lambda)
  crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
  sel_nosig <- apply(events, 1, function(x) all(x < crit_pool))
  events_nosig <- events[sel_nosig, ]
  events_sel <- events[!sel_nosig, ]

  # Remove outcomes where all results are significant and calculate the
  # probabilities later
  sel_sig <- apply(events_sel, 1, function(x) all(x >= crit))
  events_sig <- events_sel[sel_sig, ]
  events_sel <- events_sel[!sel_sig, ]

  # Conduct test for the remaining outcomes
  fun <- function(x) bskt_final(design = design, n = n, lambda = lambda, r = x,
    weight_mat = weight_mat, globalweight_fun = globalweight_fun,
    globalweight_params = globalweight_params)
  res_sel <- t(apply(events_sel, 1, fun))

  # Add results for outcomes with all or no significant baskets
  res_nosig <- matrix(rep(0, times = design@k * sum(sel_nosig)),
    ncol = design@k)
  res_allsig <- matrix(rep(1, times = design@k * sum(sel_sig)),
    ncol = design@k)
  res <- rbind(res_nosig, res_sel, res_allsig)

  # Reorder events to allign with res
  events <- rbind(events_nosig, events_sel, events_sig)

  # If all theta1 are equal each permutation has the same probability
  if (length(unique(theta1)) == 1) {
    # Compute number of correct decisions for all outcomes
    cd <- apply(res_sel, 1, function(x) sum(x == targ))

    # Add number of correct decisions for outcomes with no or all significant
    cd_nosig <- sum(rep(0, times = design@k) == targ)
    cd_allsig <- sum(rep(1, times = design@k) == targ)
    cd <- c(rep(cd_nosig, times = sum(sel_nosig)), cd, rep(cd_allsig,
      times = sum(sel_sig)))

    probs <- apply(events, 1,
      function(x) get_prob(n = n, r = x, theta = theta1))
    nperm <- apply(events, 1, get_permutations)
    exp_cd <- sum(cd * probs * nperm)
  } else {
    # If not all theta1 are equal calculate probability for each permutation
    exp_cd <- numeric(nrow(events))
    for (i in 1:nrow(events)) {
      # If number of responses is equal in each basket, there is only
      # one permutation (when n is equal)
      if (length(unique(events[i, ])) == 1) {
        probs_loop <- get_prob(n = n, r = events[i, ], theta = theta1)
        cd_loop <- sum(res[i, ] == targ)
        exp_cd[i] <- cd_loop * probs_loop
      } else {
        events_loop <- arrangements::permutations(events[i, ])
        res_loop <- arrangements::permutations(
          res[i, ])[!duplicated(events_loop), ]
        events_loop <- events_loop[!duplicated(events_loop), ]
        probs_loop <- apply(events_loop, 1, function(x) get_prob(n = n,
          r = x, theta = theta1))
        cd_loop <- apply(res_loop, 1, function(x) sum(x == targ))
        exp_cd[i] <- sum(cd_loop * probs_loop)
      }
    }
    exp_cd <- sum(exp_cd)
  }
  exp_cd
}

# Calculates the experimentwise rejection probability for a single-stage design
reject_prob_ew <- function(design, theta1, n, lambda, weight_mat,
                           globalweight_fun, globalweight_params,
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
    weight_mat = weight_mat, globalweight_fun = globalweight_fun,
    globalweight_params = globalweight_params)
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
    eff_perm <- apply(events_eff, 1, get_permutations)
    rej_prob <- sum(probs_eff * eff_perm)
  } else {
    # If not all theta1 are equal calculate probability for each permutation
    rej_prob <- numeric(nrow(res_eff))
    for (i in 1:nrow(res_eff)) {
      # If number of responses is equal in each basket, there is only
      # one permutation (when n is equal)
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

# Calculates the groupwise rejection probabilities for a single-stage design
reject_prob_group <- function(design, theta1, n, lambda, weight_mat,
                              globalweight_fun = globalweight_fun,
                              globalweight_params = globalweight_params,
                              prob = c("toer", "pwr")) {
  targ <- get_targ(theta0 = design@theta0, theta1 = theta1, prob = prob)
  # Create matrix with all possible outcomes
  events <- arrangements::permutations(0:n, k = design@k, replace = TRUE)

  # Conduct test for all possible outcomes
  fun <- function(x) bskt_final(design = design, n = n, lambda = lambda, r = x,
    weight_mat = weight_mat, globalweight_fun = globalweight_fun,
    globalweight_params = globalweight_params)
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
reject_prob_ew2 <- function(design, theta1, n, n1, lambda, interim_fun,
                            interim_params, weight_mat,
                            prob = c("toer", "pwr")) {
  targ <- get_targ(theta0 = design@theta0, theta1 = theta1, prob = prob)

  # Überprüfen, ob diese Unterscheidung notwendig ist
  if (sum(targ) == design@k) {
    events_int <- arrangements::combinations(0:n1, k = design@k, replace = TRUE)
  } else {
    events_int <- arrangements::permutations(x = 0:n1, k = design@k,
      replace = TRUE)
  }

  int_fun <- function(x) do.call(interim_fun, args = c(interim_params,
    design = design, n = n, n1 = n1, r1 = list(x), lambda = lambda,
    weight_mat = list(weight_mat)))
  res_int <- t(apply(events_int, 1, int_fun))

  # Calculate probability for events where at least one group
  # was stopped for efficacy
  eff_vec <- apply(res_int, 1, function(x) any(x[targ] == 1))
  events_eff <- events_int[eff_vec, ]
  probs_eff <- apply(events_eff, 1,
    function(x) get_prob(n = n1, r = x, theta = theta1))

  if (sum(targ) == design@k) {
    # Multiply probability by number of permutations for each row
    eff_perm <- apply(events_eff, 1, get_permutations)
    rej_prob <- sum(probs_eff * eff_perm)
  } else {
    rej_prob <- sum(probs_eff)
  }

  # Continue only when there is no relevant significant result and
  # there is at least one basket left in the second stage
  cont_vec <- !eff_vec
  cont_vec[which(cont_vec)] <- apply(res_int[cont_vec, ], 1,
    function(x) any(x == 0))
  events_cont <- events_int[cont_vec, , drop = FALSE]
  res_cont <- res_int[cont_vec, , drop = FALSE]

  all_events_stg2 <- lapply(1:design@k, function(x)
    arrangements::permutations(x = 0:(n - n1), k = x, replace = TRUE))

  for (i in 1:nrow(events_cont)) {
    no_cont <- sum(res_cont[i, ] == 0)
    # Event-Matrix für die Gruppen die in Stage 2 gehen
    events_stg2 <- all_events_stg2[[no_cont]]
    if (no_cont < design@k) {
      events_loop <- matrix(0, nrow = nrow(events_stg2), ncol = design@k)
      events_loop[, res_cont[i, ] == 0] <- events_stg2
    } else {
      events_loop <- events_stg2
    }

    events_fin <- t(t(events_loop) + events_cont[i, ])
    fin_func <- function(x) bskt_final_int(design = design, n = n, n1 = n1,
      r = x, res_int = res_cont[i, ], lambda = lambda, weight_mat = weight_mat)

    fin_res <- t(apply(events_fin, 1, fin_func))
    sig_vec <- apply(fin_res, 1, function(x) any(x[targ] == 1))

    prob_cont <- get_prob(n = n1, r = events_cont[i, ], theta = theta1)
    prob_sig <- apply(events_stg2[sig_vec, , drop = FALSE], 1, function(x)
      get_prob(n = n - n1, r = x, theta = theta1[res_cont[i, ] == 0]))
    rej_prob_temp <- prob_cont * sum(prob_sig)

    if ((rej_prob_temp > 0) & (length(unique(events_cont[i, ])) > 1) &
        (sum(targ) == design@k)) {
      perm_temp <- arrangements::npermutations(
        x = unique(events_cont[i, ]), freq = table(events_cont[i, ]))
      rej_prob_temp <- rej_prob_temp * perm_temp
    }
    rej_prob <- rej_prob + rej_prob_temp
  }
  rej_prob
}

# Calculates the groupwise rejection probability for a two-stage design
reject_prob_group2 <- function(design, theta1, n, n1, lambda, interim_fun,
                               interim_params, weight_mat,
                               prob = c("toer", "pwr")) {
  targ <- get_targ(theta0 = design@theta0, theta1 = theta1, prob = prob)
  events_int <- arrangements::permutations(x = 0:n1, k = design@k,
    replace = TRUE)

  int_fun <- function(x) do.call(interim_fun, args = c(interim_params,
    design = design, n = n, n1 = n1, r1 = list(x), lambda = lambda,
    weight_mat = list(weight_mat)))

  res_int <- t(apply(events_int, 1, int_fun))
  eff_vec <- apply(res_int, 1, function(x) any(x == 1))
  eff_vec_targ <- apply(res_int[eff_vec, ], 1, function(x) any(x[targ] == 1))
  events_eff <- events_int[eff_vec, ]
  probs_eff <- apply(events_eff, 1,
    function(x) get_prob(n = n1, r = x, theta = theta1))
  res_eff <- res_int[eff_vec,]
  rej <- colSums(apply(res_eff == 1, 2, function(x) x * probs_eff))
  rej_ew <- sum(probs_eff[eff_vec_targ])

  cont_vec <- apply(res_int, 1, function(x) any(x == 0))
  events_cont <- events_int[cont_vec, , drop = FALSE]
  res_cont <- res_int[cont_vec, , drop = FALSE]

  all_events_stg2 <- lapply(1:design@k, function(x)
    arrangements::permutations(x = 0:(n - n1), k = x, replace = TRUE))

  for (i in 1:nrow(events_cont)) {
    no_cont <- sum(res_cont[i, ] == 0)
    events_stg2 <- all_events_stg2[[no_cont]]
    if (no_cont < design@k) {
      events_loop <- matrix(0, nrow = nrow(events_stg2), ncol = design@k)
      events_loop[, res_cont[i,] == 0] <- events_stg2
    } else {
      events_loop <- events_stg2
    }

    events_fin <- t(t(events_loop) + events_cont[i, ])
    fin_func <- function(x) bskt_final_int(design = design, n = n, n1 = n1,
      r = x, res_int = res_cont[i, ], lambda = lambda, weight_mat = weight_mat)
    res_fin <- t(apply(events_fin, 1, fin_func))

    prob_cont <- get_prob(n = n1, r = events_cont[i, ], theta = theta1)
    prob_stg2 <- apply(events_stg2, 1, function(x)
      get_prob(n = n - n1, r = x, theta = theta1[res_cont[i, ] == 0]))
    rej <- rej + colSums(apply(res_fin == 1, 2, function(x) x * prob_stg2)) *
      prob_cont

    # Ignore results which were already counted for rej_ew after the
    # interim analysis
    if (all(res_cont[i, ][targ] != 1)) {
      rej_ew <- rej_ew + sum(apply(res_fin, 1, function(x) any(x[targ] == 1)) *
          prob_stg2) * prob_cont
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
