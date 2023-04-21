#' Global Weights Based on Response Rate Differences
#'
#' @template n
#' @template r
#' @template tuning_globalweights
#'
#' @details \code{globalweights_diff} calculates a weight based on the
#' heterogeneity of the response rates of all baskets that is multiplied
#' to the pairwise weights calculated with the function that is passed to
#' \code{weight_fun}. The weight is 1 when the number of responses is identical
#' in all baskets and 0 if the response rates are an equidistant sequence
#' from 0 to 1. If the maximum weight should be smaller than 1, \code{w}
#' can  be set to a smaller value.
#'
#' @return A numeric value.
#' @export
#'
#' @examples
#' globalweights_diff(n = 20, r = c(1, 3, 5), eps_global = 2)
globalweights_diff <- function(n, r, eps_global, w = 1) {
  rr <- r / n
  rs <- sort(rr)
  d <- diff(rs)
  (1 - sum(d) * 10^(-sum((d - 1 / length(d))^2)))^eps_global * w
}

globalweights_fix <- function(n, r, w) {
  w
}
