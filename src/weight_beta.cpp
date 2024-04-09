#include "RcppArmadillo.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(rng = false)]]
arma::mat weight_beta(const int& k, const arma::vec& weights,
                      const arma::mat& shape) {
  arma::mat help_mat;
  help_mat.ones(k, k);
  arma::mat weight_mat = trimatl(help_mat, -1);
  weight_mat.elem(find(weight_mat)) = weights;
  weight_mat = weight_mat + trans(weight_mat);
  weight_mat.diag().ones();
  arma::mat shape1 = trans(weight_mat * trans(shape.row(0)));
  arma::mat shape2 = trans(weight_mat * trans(shape.row(1)));
  arma::mat shapepost = join_cols(shape1, shape2);
  return shapepost;
}

// [[Rcpp::export(rng = false)]]
arma::mat weight_mat_validate(const int& k, const arma::vec& weights) {
  arma::mat help_mat;
  help_mat.ones(k, k);
  arma::mat weight_mat = trimatl(help_mat, -1);
  weight_mat.elem(find(weight_mat)) = weights;
  weight_mat = weight_mat + trans(weight_mat);
  weight_mat.diag().ones();
  return weight_mat;
}