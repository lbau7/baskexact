# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

get_res_fin <- function(prob, res_int, lambda) {
    .Call(`_baskexact_get_res_fin`, prob, res_int, lambda)
}

get_n_vec <- function(n, n1, res_int) {
    .Call(`_baskexact_get_n_vec`, n, n1, res_int)
}

get_r_temp <- function(n1, r, res_int) {
    .Call(`_baskexact_get_r_temp`, n1, r, res_int)
}

weight_beta <- function(k, weights, shape) {
    .Call(`_baskexact_weight_beta`, k, weights, shape)
}

weight_mat_validate <- function(k, weights) {
    .Call(`_baskexact_weight_mat_validate`, k, weights)
}

