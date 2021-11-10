# Computes the posterior distribution with borrowing for a single-stage design
# or an interim analysis of a two-stage design
beta_borrow <- function(weight_mat, ...) {
  UseMethod("beta_borrow", weight_mat)
}

# Borrowing method for Fujikawa's design, where the prior information is also
# shared
#' @export
beta_borrow.fujikawa <- function(weight_mat, design, n, r) {
  all_combs <- arrangements::combinations(r, 2) + 1
  weights_vec <- weight_mat[all_combs]
  shape <- matrix(c(design@shape1 + r, design@shape2 + n - r),
   byrow = TRUE, ncol = design@k)
  weight_beta(k = design@k, weights = weights_vec, shape = shape)
}

# Borrowing method for Power Prior design, where only the observed information
# is shared
#' @export
beta_borrow.pp <- function(weight_mat, design, n, r) {
  all_combs <- arrangements::combinations(r, 2) + 1
  weights_vec <- weight_mat[all_combs]
  shape <- matrix(c(r, n - r), byrow = TRUE, ncol = design@k)
  shape_post <- weight_beta(k = design@k, weights = weights_vec, shape = shape)
  shape_post[1, ] <- shape_post[1, ] + design@shape1
  shape_post[2, ] <- shape_post[2, ] + design@shape2
  shape_post
}

# Computes the posterior distribution with borrowing for the final analysis
# of a two-stage design
beta_borrow_int <- function(weight_mat, ...) {
  UseMethod("beta_borrow_int", weight_mat)
}

# Borrowing method for Fujikawa's design, where the prior information is also
# shared
#' @export
beta_borrow_int.fujikawa <- function(weight_mat, design, n1, r, res_int) {
  r_temp <- get_r_temp(n1 = n1, r = r, res_int = res_int)
  all_combs <- arrangements::combinations(r_temp, 2)
  weights_vec <- weight_mat[all_combs]
  shape <- matrix(c(design@shape1 + r, design@shape2 + n - r),
    byrow = TRUE, ncol = design@k)
  weight_beta(k = design@k, weights = weights_vec, shape = shape)
}

# Borrowing method for Power Prior design, where only the observed information
# is shared
#' @export
beta_borrow_int.pp <- function(weight_mat, design, n1, r, res_int) {
  r_temp <- get_r_temp(n1 = n1, r = r, res_int = res_int)
  all_combs <- arrangements::combinations(r_temp, 2)
  weights_vec <- weight_mat[all_combs]
  shape <- matrix(c(r, n - r), byrow = TRUE, ncol = design@k)
  shape_post <- weight_beta(k = design@k, weights = weights_vec, shape = shape)
  shape_post[1, ] <- shape_post[1, ] + design@shape1
  shape_post[2, ] <- shape_post[2, ] + design@shape2
  shape_post
}
