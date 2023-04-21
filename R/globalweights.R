# Function that calculates the global weights based on the difference
# between response rates
globalweights_diff <- function(n, r, eps_global, w = 1) {
  rr <- r / n
  rs <- sort(rr)
  d <- diff(rs)
  (1 - sum(d) * 10^(-sum((d - 1 / length(d))^2)))^eps_global * w
}

globalweights_fix <- function(n, r, w) {
  w
}
