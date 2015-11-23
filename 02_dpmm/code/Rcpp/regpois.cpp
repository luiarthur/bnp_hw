//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace std;
using namespace arma;
using namespace Rcpp;

double update_theta(double beta, double sum_y, vec x, double a, double b) {
  return rgamma( 1, a + sum_y, 1 / ( b + sum(exp(beta*x)) ) )[0]; // rgamma is shape and scale. I want shape and rate.
}

double update_beta(double beta, double theta, double sum_xy, vec x, double m, double s2, 
    double cs, int* acc) {
  double lg1, lg2, lg, cand, beta_new;

  beta_new = beta;
  cand = rnorm(1,beta,cs)[0];

  lg1 = -pow(cand-m,2) / (2*s2) + cand * sum_xy - theta * sum( exp(cand*x) );
  lg2 = -pow(beta-m,2) / (2*s2) + beta * sum_xy - theta * sum( exp(beta*x) );
  lg = lg1 - lg2;

  if (lg > log(randu()) ) {
    beta_new = cand;
    *acc = *acc + 1;
  }

  return beta_new;
}

//[[Rcpp::export]]
List regpois(vec y, vec x, double zeta, double mu, double m, double s2, double cs_beta, int B) {
  vec theta = ones<vec>(B);
  vec beta = zeros<vec>(B);
  double sum_y = sum(y);
  vec temp_vec = y.t() * x;
  double sum_xy = temp_vec[0];
  int acc_beta = 0;
  List ret;

  for (int b=1; b<B; b++) {
    theta[b] = update_theta(beta[b-1], sum_y, x, zeta, zeta/mu);
    beta[b] = update_beta(beta[b-1], theta[b], sum_xy, x, m, s2, cs_beta, &acc_beta);
    cout << "\r Progress: " << b*100/B << "%";
  }

  ret["theta"] = theta;
  ret["beta"] = beta;
  ret["acc_beta"] = acc_beta*1.0 / B;

  return ret;
}
