#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_res_fin(const NumericVector& prob, const NumericVector& res_int,
                          const double& lambda) {
  int n = prob.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if ((prob[i] >= lambda) && (res_int[i] == 0)) {
      out[i] = 1;
    } else {
      out[i] = -1;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector get_n_vec(const int& n, const int& n1,
                        const NumericVector& res_int) {
  int n0 = res_int.size();
  NumericVector n_vec(n0);
  for (int i = 0; i < n0; ++i) {
    if (res_int[i] == 0) {
      n_vec[i] = n;
    } else {
      n_vec[i] = n1;
    }
  }
  return n_vec;
}

// [[Rcpp::export]]
NumericVector get_r_temp(const int& n1, const NumericVector& r,
                         const NumericVector& res_int) {
  int n = res_int.size();
  NumericVector r_temp(n);
  for (int i = 0; i < n; ++i) {
    if (res_int[i] == 0) {
      r_temp[i] = r[i] + n1 + 2;
    } else {
      r_temp[i] = r[i] + 1;
    }
  }
  return r_temp;
}
