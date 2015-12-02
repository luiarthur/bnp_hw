//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <algorithm> // unique
#include <../../../02_dpmm/code/Rcpp/dp/dp.h> //wsample

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


mat rowSort(mat X) {
  int n = X.n_rows;
  int k = X.n_cols;
  vec pows = zeros<vec>(k);
  vec rsums = zeros<vec>(n);
  uvec ind;

  for (int j=0; j<k; j++) pows[j] = pow(10,j);
  for (int i=0; i<n; i++) {
    rsums[i] = sum(X.row(i) * pows);
  }

  ind = sort_index(rsums);

  return X.rows(ind);
}

mat uniqueRows(mat X) {
  mat sortedX = rowSort(X);
  int n = X.n_rows;
  vec ind = zeros<vec>(n);
  ind[0] = 1;
  for (int i=1; i<n; i++) {
    if ( any(sortedX.row(i) != sortedX.row(i-1))  ) { // if not identical
      ind[i] = 1;
    }
  }
  if ( any(sortedX.row(n-1) != sortedX.row(n-2))  ) { // if not identical
    ind[n-1] = 1;
  }
  return sortedX.rows( find(ind==1) );
}


double update_alpha (double alpha, mat theta, double a, double b, int n) {
  double eta, c, d;
  int ns, ind;
  double alpha_new;
  vec ut = uniqueRows(theta);

  eta = rbeta(1, alpha + 1, n )[0];
  ns = ut.n_rows;
  c = a + ns;
  d = b - log(eta);
  
  ind = wsample({0.0, 1.0}, {c-1, n*d});
  if (ind == 0) {
    alpha_new = rgamma(1, c, 1/d )[0]; // shape, rate
  } else {
    alpha_new = rgamma(1, c-1, 1/d )[0];
  }

  return alpha_new;
}


//[[Rcpp::export]]
List sdp(NumericVector y, vec s_new, double beta_mu, double beta_s2,
    double tau2_a, double tau2_b, double alpha_a, double alpha_b,
    double sig2_a, double sig2_b, double phi_a, double phi_b,
    double cs_beta, double cs_tau2, double cs_sig2, double cs_phi, int B) {

  IntegerVector arrayDims = y.attr("dim");
  cube Y(y.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

  int T = Y.n_rows;
  int n = Y.n_slices;
  vec beta = zeros<vec>(B);
  vec tau2 = zeros<vec>(B);
  vec alpha = zeros<vec>(B);
  vec sig2 = zeros<vec>(B);
  vec phi = zeros<vec>(B);
  cube theta = zeros<cube>(B,T,n);
  int acc_beta = 0;
  int acc_tau2 = 0;
  int acc_sig2 = 0;
  int acc_phi = 0;
  List ret;
  mat tb;

  for (int b=1; b<B; b++) {
    for (int t=0; t<T; t++) {
      //update theta
      tb = theta.slice(b);
      //update beta
      //update tau2
      //update alpha
      update_alpha (alpha[b-1], tb, alpha_a, alpha_b, T); // gibbs 
      //update sig2
      //update phi
    }

    cout <<"\r Progress: " << b << "/" << B << endl;
  }

  ret["beta"] = beta;
  ret["tau2"] = tau2;
  ret["alpha"] = alpha;
  ret["sig2"] = sig2;
  ret["phi"] = phi;
  ret["acc_sig2"] = acc_sig2;
  ret["acc_phi"] = acc_phi;
  ret["theta"] = wrap(theta);

  return ret;
}
