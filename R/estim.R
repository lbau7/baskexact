#' Posterior Mean and Mean Squared Error
#'
#' Computes the posterior mean and the mean squared error of a basket trial
#' design.
#'
#' @template design
#' @template dotdotdot
#'
#' @return A list containing means of the posterior distributions and
#' the mean squared errors for all baskets.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, p0 = 0.2)
#' estim(design = design, p1 = c(0.2, 0.2, 0.5), n = 15,
#'   weight_fun = weights_fujikawa)
setGeneric("estim",
  function(design, ...) standardGeneric("estim")
)

#' @describeIn estim Posterior mean and mean squared error for a single-stage
#'   basket design.
#'
#' @template design
#' @template p1_pow
#' @template n
#' @param lambda Probability threshold. This is only required when
#'   pruning is applied.
#' @template weights
#' @template globalweights
#' @template dotdotdot
setMethod("estim", "OneStageBasket",
  function(design, p1, n, lambda = NULL, weight_fun, weight_params = list(),
           globalweight_fun = NULL, globalweight_params = list(), ...) {
    check_p1(design = design, p1 = p1, type = "estim")
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda, globalweight_fun = globalweight_fun,
      globalweight_params = list(globalweight_params)))

    events <- arrangements::permutations(0:n, k = design@k, replace = TRUE)
    prob_events <- apply(events, 1, function(x) get_prob(n = n, r = x,
      p = p1))
    post_shapes <- apply(events, 1, function(x) beta_borrow(weight_mat =
        weight_mat, globalweight_fun = globalweight_fun,
      globalweight_params = globalweight_params, design = design, n = n,
      r = x), simplify = FALSE)
    post_means <- t(sapply(post_shapes, mean_beta))

    list(
      Mean = colSums(post_means * prob_events),
      MSE = colSums(t((t(post_means) - p1)^2) * prob_events)
    )
  })
