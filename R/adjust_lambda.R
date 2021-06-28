#' @describeIn adjust_lambda Adjust lambda for a single-stage design.
adjust_lambda <- function(design, alpha = 0.025, n, epsilon, tau, logbase,
                          prune, prec_digits, ...) {
  upper_lim <- 1
  crit_max <- get_crit(design = design, n = n, lambda = upper_lim)
  if (is.na(crit_max)) upper_lim <- 1 - 10^(-prec_digits)

  root_fun <- function(x) toer(design = design, n = n, lambda = x,
    epsilon = epsilon, tau = tau, logbase = logbase, prune = prune,
    results = "fwer") - alpha
  uni_root <- stats::uniroot(root_fun, interval = c(0.5, upper_lim),
    tol = 10^(-prec_digits))

  if (uni_root$f.root > 0 ) {
    root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
    val <- root_fun(root)
    if (val > 0) {
      while (val > 0) {
        root <- root + 10^(-prec_digits)
        val <- root_fun(root)
      }
    }
  } else {
    root <- floor(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
    val <- root_fun(root)
    if (val > 0) {
      root <- ceiling(uni_root$root * 10^(prec_digits)) / 10^(prec_digits)
      val <- root_fun(root)
    } else {
      repeat {
        root_temp <- root - 10^(-prec_digits)
        val_temp <- root_fun(root)
        if (val_temp < 0) break
        root <- root_temp
        val <- val_temp
      }
    }
  }
  return(list(
    lambda = root,
    toer = val + alpha
  ))
}
