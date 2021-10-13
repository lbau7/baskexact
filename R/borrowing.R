# Computes the posterior distribution with borrowing for a single-stage design
# or an interim analysis of a two-stage design
beta_borrow <- function(k, r, weight_mat, shape) {
  all_combs <- arrangements::combinations(r, 2) + 1
  weights_vec <- weight_mat[all_combs]
  weight_beta(k = k, weights = weights_vec, shape = shape)
}

# Computes the posterior distribution with borrowing for the final analysis
# of a two-stage design
beta_borrow_int <- function(k, n1, r, res_int, weight_mat, shape) {
  r_temp <- get_r_temp(n1 = n1, r = r, res_int = res_int)
  all_combs <- arrangements::combinations(r_temp, 2)
  weights_vec <- weight_mat[all_combs]
  weight_beta(k = k, weights = weights_vec, shape = shape)
}
