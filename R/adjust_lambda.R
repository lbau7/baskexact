#' @include generics.R
NULL

#' @describeIn adjust_lambda Adjust lambda for a single-stage design.
setMethod("adjust_lambda", "OneStageBasket",
  function(design, alpha = 0.025, theta1 = NULL, n, epsilon, tau, logbase,
           prune = FALSE, prec_digits, ...) {
    theta1 <- check_theta1(design = design, theta1 = theta1, type = "toer")
    check_tuning(epsilon = epsilon, tau = tau, logbase = logbase)
    if (alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1")
    if (length(n) != 1) stop("n must have length 1")

    upper_lim <- 1 - 10^(-prec_digits)
    root_fun <- function(x) toer(design = design, theta1 = theta1, n = n,
      lambda = x, epsilon = epsilon, tau = tau, logbase = logbase,
      prune = prune, results = "fwer") - alpha

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
  })
