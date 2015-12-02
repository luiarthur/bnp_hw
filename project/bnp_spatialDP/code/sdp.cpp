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

mat mvrnorm(mat M, mat S) {
  int n = M.n_rows;
  mat e = randn(n);
  return M + chol(S).t()*e;
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


double update_beta (double sum_y, mat theta, double tau2, double s2, int T, int n) {
  double denom = tau2 + T*n*s2;
  double m_new = s2 * ( sum_y - sum(sum(theta)) ) / denom;
  double s_new = sqrt( s2*tau2 / denom );

  return rnorm(1, m_new, s_new)[0];
}

double update_tau2 (double a, double b, int n, int T, mat y, mat theta, double beta) {
  double a_new = a + n*T/2;
  mat one_Tn(T,n);
  mat M = y - theta - one_Tn;
  double sum_mm = 0;
  mat m;

  for (int t=0; t<T; t++) {
    m = M.row(t).t() * M.row(t);
    sum_mm += m(0,0);
  }

  double rate_new = b + sum_mm/2;
  double scale_new = 1 / rate_new;

  return rgamma(1, a_new, scale_new)[0];
}

//[[Rcpp::export]]
List sdp(mat Y, vec s_new, mat D, double beta_mu, double beta_s2,
    double tau2_a, double tau2_b, double alpha_a, double alpha_b,
    double sig2_a, double sig2_b, double phi_a, double phi_b,
    double cs_sig2, double cs_phi, int B) {

  int T = Y.n_rows;
  int n = Y.n_cols;
  double sum_y = sum( sum(Y) );
  vec beta = zeros<vec>(B);
  vec tau2 = zeros<vec>(B);
  vec alpha = zeros<vec>(B);
  vec sig2 = zeros<vec>(B);
  vec phi = zeros<vec>(B);
  cube theta = zeros<cube>(B,T,n);
  int acc_sig2 = 0;
  int acc_phi = 0;
  List ret;
  mat tb;

  for (int b=1; b<B; b++) {
    for (int t=0; t<T; t++) {
      //update theta
      tb = theta.slice(b);
      beta[b] = update_beta(sum_y, tb, tau2[b-1], beta_s2, T, n); // gibbs
      tau2[b] = update_tau2(tau2_a, tau2_b, n, T, Y, tb, beta[b]); //gibbs
      alpha[b] = update_alpha(alpha[b-1], tb, alpha_a, alpha_b, T); // gibbs 
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
