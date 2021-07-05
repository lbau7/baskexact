#' @describeIn pow Power for a single-stage design.
setMethod("pow", "OneStageBasket",
  function(design, n, lambda, epsilon, tau, logbase = 2, prune = FALSE,
           results = c("ewp", "group"), ...) {
    if (length(n) != 1) stop("n must have length 1")
    if (lambda <= 0 | lambda >= 1) stop("lambda must be between 0 and 1")
    if (epsilon < 0) stop("epsilon must be positive")
    if (tau < 0 | tau >= 1) stop("tau must be in [0, 1)")
    if (logbase <= 0) stop("logbase must be positive")

    results <- match.arg(results)
    weight_mat <- get_weights(design = design, n = n, epsilon = epsilon,
      tau = tau, logbase = logbase)

    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    if (results == "ewp") {
      reject_prob_ew(design = design, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "pwr")
    } else {
      reject_prob_group(design = design, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "pwr")
    }
  })
