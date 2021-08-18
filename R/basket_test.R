#' @include generics.R
NULL

#' @describeIn basket_test Testing for a single-stage basket design.
setMethod("basket_test", "OneStageBasket",
  function(design, n, r, lambda, epsilon, tau, logbase = 2, prune, ...) {
    check_params(n = n, lambda = lambda)
    check_tuning(epsilon = epsilon, tau = tau, logbase = logbase)
    if (any(r > n) | any(r < 0)) stop("responses must be between 0 and n")

    weight_mat <- get_weights(design = design, n = n, epsilon = epsilon,
      tau = tau, logbase = logbase)
    if (prune) {
      crit_pool <- get_crit_pool(design = design, n = n, lambda = lambda)
      weight_mat <- prune_weights(weight_mat = weight_mat, cut = crit_pool)
    }

    all_combs <- arrangements::combinations(r, 2) + 1
    weights_vec <- weight_mat[all_combs]
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
    dimnames(shape_post) <- list(
      c("shape1", "shape2"),
      sapply(1:design@k, function(x) paste("Basket", x))
    )
    shape_borrow <- beta_borrow(k = design@k, r = r, weight_mat = weight_mat,
      shape = shape_post)
    dimnames(shape_borrow) <- list(
      c("shape1", "shape2"),
      sapply(1:design@k, function(x) paste("Basket", x))
    )

    postprob <- post_beta(shape_post, theta0 = design@theta0)
    postprob_borrow <- post_beta(shape_borrow, theta0 = design@theta0)

    list(
      weights = weights,
      post_dist_noborrow = shape_post,
      post_dist_borrow = shape_borrow,
      post_prob_noborrow = postprob,
      post_prob_borrow = postprob_borrow
    )
  })

