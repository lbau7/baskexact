bskt_final <- function(design, n, lambda, r, weight_mat) {
  shape_post <- matrix(c(design@shape1 + r, design@shape2 + n - r),
    byrow = TRUE, ncol = design@k)
  shape_borrow <- beta_borrow(k = design@k, r = r, weight_mat = weight_mat,
    shape = shape_post)
  post_prob <- post_beta(shape = shape_borrow, theta0 = design@theta0)
  ifelse(post_prob >= lambda, 1, 0)
}
