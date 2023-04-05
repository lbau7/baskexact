# Function that conducts the final test in a one-stage design
bskt_final <- function(design, n, lambda, r, weight_mat,
                       globalweight_fun  = NULL, globalweight_params) {
  shape_borrow <- beta_borrow(weight_mat = weight_mat, design = design, n = n,
    r = r, globalweight_fun = globalweight_fun,
    globalweight_params = globalweight_params)
  post_prob <- post_beta(shape = shape_borrow, theta0 = design@theta0)
  ifelse(post_prob >= lambda, 1, 0)
}

# Function that conducts the final test in a two-stage design
bskt_final_int <- function(design, n, n1, r, res_int, lambda, weight_mat) {
  shape_borrow <- beta_borrow_int(weight_mat = weight_mat, design = design,
    n = n, n1 = n1, r = r, res_int = res_int)
  post_prob <- post_beta(shape = shape_borrow, theta0 = design@theta0)
  get_res_fin(post_prob, res_int, lambda)
}
