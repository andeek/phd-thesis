#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]
arma::vec binary_two_param_reg_sampler(List data, List params) {
  
  List sums = data["sums"];
  List nums = data["nums"];
  vec u = data["u"];
  
  vec sums_vec_1 = sums[0];
  vec nums_vec_1 = nums[0];
  vec sums_vec_2 = sums[1];
  vec nums_vec_2 = nums[1];
  
  vec res(sums_vec_1.n_elem);
  
  RNGScope scope;
  
  double beta_0 = params["beta_0"];
  double beta_1 = params["beta_1"];
  double eta_1 = params["eta_1"];
  double eta_2 = params["eta_2"];
  
  vec reg = beta_0 + beta_1 * u;
  vec kappa = exp(reg) / (1 + exp(reg));
  
  vec mean_structure = reg + eta_1 * (sums_vec_1 - nums_vec_1 % kappa) + eta_2 * (sums_vec_2 - nums_vec_2 % kappa);
  
  int i;
  
  for(i = 0; i < res.n_elem; ++i) {
    res(i) = as<double>(rbinom(1, 1, exp(mean_structure(i))/(1 + exp(mean_structure(i)))));
  }
  return(res);
}
