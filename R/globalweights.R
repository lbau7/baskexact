# Function that calculates the global weights based on the difference
# between response rates
globalweights_diff <- function(n, r, eps_global) {
  rr <- r / n
  rs <- sort(rr)
  d <- diff(rs)
  (1 - sum(d) * 10^(-sum((d - 1 / length(d))^2)))^eps_global
}
