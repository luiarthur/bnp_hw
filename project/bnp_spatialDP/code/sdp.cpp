//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <../../../02_dpmm/code/Rcpp/dp/dp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector test (NumericVector Y) {
  //tt <- array(rnorm(18),dim=c(3,3,2)); test(tt) # in R
  IntegerVector arrayDims = Y.attr("dim");
  cube cubeArray(Y.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
  cubeArray(0,0,0) = 518;  
  cout << "slices: " << cubeArray.n_slices <<endl;
  cout << "rows: " << cubeArray.n_rows <<endl;
  cout << "cols: " << cubeArray.n_cols <<endl;
  return(wrap(cubeArray.slice(0))); 
}



//[[Rcpp::export]]
List dppois(NumericVector y, vec s_new, double beta_mu, double beta_s2,
    double tau2_a, double tau2_b, double alpha_a, double alpha_b,
    double sig2_a, double sig2_b, double phi_a, double phi_b,
    double cs_s2, double cs_phi, int B) {

  IntegerVector arrayDims = y.attr("dim");
  cube Y(y.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

  int T = Y.n_rows;
  int n = Y.n_slices;
  vec beta;
  vec tau2;
  vec alpha;
  vec sig2;
  vec phi;
  NumericVector theta; //cube
  int acc_s2 = 0;
  int acc_phi = 0;
  List ret;

  for (int b=1; b<B; b++) {
    for (int t=0; t<T; t++) {
      //update beta
      //update tau2
      //update alpha
      //update sig2
      //update phi
      //update theta
    }
  }

  ret["beta"] = beta;
  ret["tau2"] = tau2;
  ret["alpha"] = alpha;
  ret["sig2"] = sig2;
  ret["phi"] = phi;
  ret["acc_s2"] = acc_s2;
  ret["acc_phi"] = acc_phi;
  ret["theta"] = wrap(theta);

  return ret;
}
