#' @include class.R
NULL

check_tuning <- function(epsilon, tau, logbase) {
  if (epsilon < 0) stop("epsilon must be non-negative")
  if (tau < 0 | tau >= 1) stop("tau must be in [0, 1)")
  if (logbase <= 0) stop("logbase must be positive")
}

check_theta1 <- function(design, theta1, type) {
  if (is.null(theta1)) theta1 <- rep(design@theta0, design@k)
  if (any(theta1 < design@theta0)) {
    stop("all theta1 must be greater than or equal to theta0")
  }

  if ((length(theta1) != design@k) & length(theta1) != 1) {
    stop("theta1 must have length k or length 1")
  }
  if (length(theta1) == 1) theta1 <- rep(theta1, design@k)
    if ((type == "toer") & all(theta1 != design@theta0)) {
    stop("no true null hyoptheses, cannot compute type 1 error rate")
  }
  if ((type == "pwr") & all(theta1 == design@theta0)) {
    stop("no true alternative hyoptheses, cannot compute power")
  }

  theta1
}

check_params <- function(n, lambda) {
  if (length(n) != 1) stop("n must have length 1")
  if (lambda <= 0 | lambda >= 1) stop("lambda must be between 0 and 1")
}

