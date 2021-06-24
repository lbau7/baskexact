# These functions cannot be accessed by the user. They are only used to
# validate the results of the main functions, by using other (mostly
# slower but easier to follow) algorithms to calculate the same values.

reject_prob_ew_slow <- function(design, n, lambda, weight_mat,
                                prob = c("toer", "pwr")) {
  targ <- get_targ(design = design, prob = prob)
  events <- arrangements::permutations(0:n, k = design@k, replace = TRUE)

  fun <- function(x) bskt_final(design = design, n = n, lambda = lambda, r = x,
    weight_mat = weight_mat)
  res <- t(apply(events, 1, fun))

  eff_vec <- apply(res, 1, function(x) any(x[targ] == 1))
  events_eff <- events[eff_vec, ]
  probs_eff <- apply(events_eff, 1,
    function(x) get_prob(n = n, r = x, theta = design@theta1))
  sum(probs_eff)
}
