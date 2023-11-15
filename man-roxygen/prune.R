#' @param prune Whether baskets with a number of responses below the
#'   critical pooled value should be pruned before the final analysis.
#'   If this is \code{TRUE} then \code{lambda} is also required and
#'   if \code{globalweight_fun} is not \code{NULL} then
#'   \code{globalweight_fun} and \code{globalweight_params} are also used.
