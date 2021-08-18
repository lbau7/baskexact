#' @include generics.R
NULL

#' @describeIn toer Type 1 error rate for a single-stage basket design.
setMethod("toer", "OneStageBasket",
  function(design, theta1 = NULL, n, lambda, epsilon, tau, logbase = 2,
           prune = FALSE, results = c("fwer", "group"), ...) {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "toer")
    check_tuning(epsilon = epsilon, tau = tau, logbase = logbase)
    check_params(n = n, lambda = lambda)

    results <- match.arg(results)
    weight_mat <- get_weights(design = design, n = n, epsilon = epsilon,
      tau = tau, logbase = logbase)

    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    if (results == "fwer") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "toer")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat, prob = "toer")
    }
  })
