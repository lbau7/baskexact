#' @include class.R
NULL

check_tuning <- function(epsilon, tau, logbase) {
  if (epsilon < 0) stop("epsilon must be non-negative")
  if (tau < 0 | tau >= 1) stop("tau must be in [0, 1)")
  if (logbase <= 0) stop("logbase must be positive")
}

check_p1 <- function(design, p1, type) {
  if (is.null(p1)) p1 <- rep(design@p0, design@k)
  if (any(p1 < design@p0)) {
    stop("all p1 must be greater than or equal to p0")
  }

  if ((length(p1) != design@k) & length(p1) != 1) {
    stop("p1 must have length k or length 1")
  }
  if (length(p1) == 1) p1 <- rep(p1, design@k)
  if ((type == "toer") & all(p1 != design@p0)) {
    stop("no true null hypotheses, cannot compute type 1 error rate")
  }
  if ((type == "pwr") & all(p1 == design@p0)) {
    stop("no true alternative hyoptheses, cannot compute power")
  }

  p1
}

check_params <- function(n, lambda) {
  if (length(n) != 1) stop("n must have length 1")
  if (lambda <= 0 | lambda >= 1) stop("lambda must be between 0 and 1")
}

