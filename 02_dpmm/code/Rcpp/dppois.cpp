//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <dp/dp.h>

using namespace std;
using namespace arma;
using namespace Rcpp ;

//R::rgamma(a,b); // sh,sc

//[[Rcpp::export]]
List dppois(double a_mu, double b_mu, double a_tau, double b_tau,
    double m_beta, double s_beta, int B) { // also need cand_sigs

  mat theta;
  vec mu, tau, beta;
  int acc_mu, acc_beta, acc_tau;
  vec acc_theta;
  List ret;
  
  // Update thetas
  // Update beta
  // Update mu
  // Update tau

  ret["mu"] = mu;
  ret["tau"] = tau;
  ret["beta"] = beta;
  ret["theta"] = theta;
  ret["acc_mu"] = acc_mu;
  ret["acc_tau"] = acc_tau;
  ret["acc_beta"] = acc_beta;
  ret["acc_theta"] = acc_theta;

  return ret;
}
