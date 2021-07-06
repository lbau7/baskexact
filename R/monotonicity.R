#' @include generics.R
NULL

#' @describeIn check_mon_within Within-trial monotonicity condition for a
#' single-stage design.
setMethod("check_mon_within", "OneStageBasket",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune, details, ...) {
    # Not working with different priors and different n!
    if (length(n) != 1) stop("n must have length 1")
    if (lambda <= 0 | lambda >= 1) stop("lambda must be between 0 and 1")
    if (epsilon < 0) stop("epsilon must non-negative positive")
    if (tau < 0 | tau >= 1) stop("tau must be in [0, 1)")
    if (logbase <= 0) stop("logbase must be positive")

    crit <- get_crit(design = design, n = n, lambda = lambda)
    crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
    weight_mat <- get_weights(design = design, n = n, epsilon = epsilon,
      tau = tau, logbase = logbase)
    if (prune) {
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
    # Discard events where all or no baskets are significanct
    sel_events <- apply(events, 1,
      function(x) all(x >= crit) | all(x < crit_pool))
    events <- events[!sel_events, ]

    # If pruning is used and no list of violating outcomes is desired,
    # then outcomes with different number of responses in baskets that
    # are pruned can be ignored
    if (prune & !details) {
      events[which(events < crit_pool)] <- 0
      events <- events[!duplicated(events), ]
    }

    func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
      r = x, weight_mat = weight_mat)

    # Conduct test for all remaining outcomes
    res <- t(apply(events, 1, func))
    res_sig <- rowSums(res)
    # Investigate outcomes where some (but not all) baskets are significant
    res_inv <- res[which(res_sig > 0 & res_sig < design@k), ]
    no_mon <- apply(res_inv, 1, function(x) any(x != cummax(x)))
    check <- sum(no_mon) == 0

    if (check) {
      check
    } else {
      if (details) {
        events_inv <- events[which(res_sig > 0 & res_sig < design@k), ]
        events_nomon <- events_inv[no_mon, ]
        res_nomon <- res_inv[no_mon, ]
        list(
          Events = events_nomon,
          Results = res_nomon
        )
      } else {
        check
      }
    }
  })

#' @describeIn check_mon_between Between-trial monotonicity condition for a
#' single-stage design.
setMethod("check_mon_between", "OneStageBasket",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune, details, ...) {
    if (length(n) != 1) stop("n must have length 1")
    if (lambda <= 0 | lambda >= 1) stop("lambda must be between 0 and 1")
    if (epsilon < 0) stop("epsilon must non-negative positive")
    if (tau < 0 | tau >= 1) stop("tau must be in [0, 1)")
    if (logbase <= 0) stop("logbase must be positive")

    crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
    weight_mat <- get_weights(design = design, n = n, epsilon = epsilon,
      tau = tau, logbase = logbase)
    if (prune) {
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    events <- arrangements::combinations(0:n, k = design@k, replace = TRUE)
    func <- function(x) bskt_final(design = design, n = n, lambda = lambda,
      r = x, weight_mat = weight_mat)

    res <- t(apply(events, 1, func))
    res_sig <- apply(res, 1, function(x) any(x == 1))
    res_test <- res_sig

    if (details) detlist <- list()
    checkout <- TRUE
    for (i in 1:length(res_sig)) {
      if (res_test[i]) {
        events_sel <- apply(events, 1, function(x) all(x >= events[i, ]))
        res_sel <- res_sig[events_sel]
        check <- sum(res_sel) == length(res_sel)
        if (check) {
          res_test[events_sel] <- FALSE
        } else {
          if (details) {
            if (checkout) checkout <- FALSE
            detlist <- c(detlist, list(list(
              Events = rbind(events[i, ], events[events_sel & !res_sig, ]),
              Results = rbind(res[i, ], res[events_sel & !res_sig, ])
            )))
          } else {
            checkout <- FALSE
            break
          }
        }
      }
    }

    if (checkout) {
      TRUE
    } else {
      if (details) {
        detlist
      } else {
        FALSE
      }
    }
  })
