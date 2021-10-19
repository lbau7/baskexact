# Function that conducts the final test in a one-stage design
bskt_final <- function(design, n, lambda, r, weight_mat) {
  shape_post <- matrix(c(design@shape1 + r, design@shape2 + n - r),
    byrow = TRUE, ncol = design@k)
  shape_borrow <- beta_borrow(k = design@k, r = r, weight_mat = weight_mat,
    shape = shape_post)
  post_prob <- post_beta(shape = shape_borrow, theta0 = design@theta0)
  ifelse(post_prob >= lambda, 1, 0)
}

bskt_final_int <- function(design, n, n1, r, res_int, lambda, weight_mat) {
  n_vec <- get_n_vec(n = n, n1 = n1, res_int = res_int)
  shape_post <- matrix(c(design@shape1 + r, design@shape2 + n_vec - r),
    byrow = TRUE, ncol = design@k)
  shape_borrow <- beta_borrow_int(k = design@k, n1 = n1, r = r,
    shape = shape_post, res_int = res_int, weight_mat = weight_mat)
  post_prob <- post_beta(shape = shape_borrow, theta0 = design@theta0)
  get_res_fin(post_prob, res_int, lambda)
}
