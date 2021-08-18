#' @include generics.R
NULL

#' @describeIn pow Power for a single-stage basket design.
setMethod("pow", "OneStageBasket",
  function(design, theta1, n, lambda, epsilon, tau, logbase = 2, prune = FALSE,
           results = c("ewp", "group"), ...) {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "pwr")
    check_tuning(epsilon = epsilon, tau = tau, logbase = logbase)
    check_params(n = n, lambda = lambda)

    results <- match.arg(results)
    weight_mat <- get_weights(design = design, n = n, epsilon = epsilon,
      tau = tau, logbase = logbase)

    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    if (results == "ewp") {
      reject_prob_ew(design = design, theta1 = theta1, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "pwr")
    } else {
      reject_prob_group(design = design, theta1 = theta1, n = n,
        lambda = lambda, weight_mat = weight_mat, prob = "pwr")
    }
  })
