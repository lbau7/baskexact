#' @describeIn adjust_lambda Adjust lambda for a single-stage design.
adjust_lambda <- function(design, alpha = 0.025, n, epsilon, tau, logbase,
                          prune, prec_digits, ...) {
  upper_lim <- 1 - 10^(-prec_digits)
  root_fun <- function(x) toer(design = design, n = n, lambda = x,
    epsilon = epsilon, tau = tau, logbase = logbase, prune = prune,
    results = "fwer") - alpha

  # Use uniroot to find lambda close to the smallest lambda that protects
  # the significance level at alpha
  uni_root <- stats::uniroot(root_fun, interval = c(0.5, upper_lim),
    tol = 10^(-prec_digits))

  if (uni_root$f.root > 0 ) {
    # If rejection prob is greater than alpha after uniroot, round lambda up
    root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
    val <- root_fun(root)
    if (val > 0) {
      # If rejection prob is still above alpha, increase lambda
      while (val > 0) {
        root <- root + 10^(-prec_digits)
        val <- root_fun(root)
      }
    }
  } else {
    # If rejection prob is less than alpha after uniroot, round lambda down
    root <- floor(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
    val <- root_fun(root)
    if (val > 0) {
      # If the rejection prob is greater now than alpha with the rounded-down
      # lambda, then round lambda up
      root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
      val <- root_fun(root)
    } else {
      # If the rejection prob is still below alpha, decrease lambda
      repeat {
        root_old <- root
        val_old <- val
        root <- root - 10^(-prec_digits)
        val <- root_fun(root)
        if (val > 0) {
          root <- root_old
          val <- val_old
          break
        }
      }
    }
  }
  return(list(
    lambda = root,
    toer = val + alpha
  ))
}
