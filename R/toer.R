#' @describeIn toer Type 1 error rate for a single stage design.
setMethod("toer", "Basket",
  function(design, n, lambda, epsilon, tau, logbase = 2,
           results = c("fwer", "group"), ...) {
    if (all(design@theta1 != design@theta0)) {
      message("all theta1 > theta0, computing type 1 error rates under global null")
      design@theta1 <- rep(design@theta0, design@k)
    }
    if (length(n) != 1) stop("n must have length 1")
    if (lambda <= 0 | lambda >= 1) stop("lambda must be between 0 and 1")
    if (epsilon < 0) stop("epsilon must be positive")
    if (tau < 0 | tau >= 1) stop("tau must be in [0, 1)")
    if (logbase <= 0) stop("logbase must be positive")

    results <- match.arg(results)
    weight_mat <- get_weights(design = design, n = n, epsilon = epsilon,
      tau = tau, logbase = logbase)

    if (results == "fwer") {
      reject_prob_ew(design = design, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "toer")
    } else {
      reject_prob_group(design = design, n = n, lambda = lambda,
        weight_mat = weight_mat, prob = "toer")
    }
  })
