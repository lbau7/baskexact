#' @include class.R
NULL

#' Test for the Results of a Basket Trial
#'
#' \code{basket_test} evaluates the results of a basket trial and calculates
#' the posterior distributions with and without borrowing.
#'
#' @template design
#' @template dotdotdot
#'
#' @return If \code{details = TRUE}: A list, including matrices of the weights
#' that are used for borrowing, posterior distribution parameters for all
#' baskets without and with borrowing, as well as the posterior probabilities
#' for all baskets without and with borrowing. If \code{details = FALSE}:
#' The posterior probabilities for all baskets with borrowing.
#' @export
#'
#' @examples
#' design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, p0 = 0.2)
#' basket_test(design = design, n = 24, r = c(5, 9, 10), lambda = 0.99,
#'   weight_fun = weights_fujikawa)
setGeneric("basket_test",
  function(design, ...)
    standardGeneric("basket_test")
)

#' @template design
#' @template n
#' @param r The vector of observed responses.
#' @template lambda
#' @template weights
#' @template globalweights
#' @template dotdotdot
#' @params details Whether a detailed list of results or only the vector
#'   of posterior probabilities is returned.
#' @describeIn basket_test Testing for a single-stage basket design.
setMethod("basket_test", "OneStageBasket",
  function(design, n, r, lambda, weight_fun, weight_params = list(),
           globalweight_fun = NULL, globalweight_params = list(),
           details = TRUE, ...) {
    check_params(n = n, lambda = lambda)
    if (any(r > n) | any(r < 0)) stop("responses must be between 0 and n")
    if (length(r) != design@k) stop("r must have length k")
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n, lambda = lambda, globalweight_fun = globalweight_fun,
      globalweight_params = list(globalweight_params)))

    all_combs <- arrangements::combinations(r, 2) + 1
    weights_vec <- weight_mat[all_combs]

    if (!is.null(globalweight_fun)) {
      w <- do.call(globalweight_fun, args = c(n = n, list(r = r),
        globalweight_params))
      weights_vec <- weights_vec * w
    }

    weights <- matrix(0, nrow = design@k, ncol = design@k)
    weights[lower.tri(weights, diag = FALSE)] <- weights_vec
    weights <- weights + t(weights)
    diag(weights) <- 1
    dimnames(weights) <- list(
      sapply(1:design@k, function(x) paste("Basket", x)),
      sapply(1:design@k, function(x) paste("Basket", x))
    )

    shape_post <- matrix(c(design@shape1 + r, design@shape2 + n - r),
      byrow = TRUE, ncol = design@k)
    shape_borrow <- beta_borrow(weight_mat = weight_mat, globalweight_fun =
        globalweight_fun, globalweight_params = globalweight_params,
      design = design, n = n, r = r)

    if (details) {
      dimnames(shape_post) <- list(
        c("shape1", "shape2"),
        sapply(1:design@k, function(x) paste("Basket", x))
      )
      dimnames(shape_borrow) <- list(
        c("shape1", "shape2"),
        sapply(1:design@k, function(x) paste("Basket", x))
      )

      postprob <- post_beta(shape_post, p0 = design@p0)
      postprob_borrow <- post_beta(shape_borrow, p0 = design@p0)

      list(
        weights = weights,
        post_dist_noborrow = shape_post,
        post_dist_borrow = shape_borrow,
        post_prob_noborrow = postprob,
        post_prob_borrow = postprob_borrow
      )
    } else {
      post_beta(shape_borrow, p0 = design@p0)
    }
  })

