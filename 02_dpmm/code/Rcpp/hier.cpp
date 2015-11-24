//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <dp/dp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

// http://adv-r.had.co.nz/Rcpp.html#rcpp-sugar
// Enable C++11 in R: Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
// This is used for randg(shape,scale), vec v = {1,2,3};


mat update_theta (vec theta, double zeta, double mu, double beta, vec x, vec y) {
  int n = theta.size();
  vec theta_new = zeros<vec>(n);
  double shape, scale;
  
  for (int i=0; i<n; i++) {
    shape = y[i] + zeta;
    scale = 1 / ( zeta/mu + exp(x[i]*beta) );
    theta_new[i] = 1 / rgamma(1, shape, scale)[0];
  }

  return reshape(theta_new,1,n);
}

double update_beta (double beta, vec theta, double s_xy, double m, double s2, vec y, vec x, double cs, int* acc) {

  double lg1, lg2, lg, c, beta_new;

  beta_new = beta;
  c = randn() * cs + beta;
  lg1 = -pow(c-m,2)/(2*s2)    +    c*s_xy - sum(theta % exp(c*x));
  lg2 = -pow(beta-m,2)/(2*s2) + beta*s_xy - sum(theta % exp(beta*x));
  lg = lg1 - lg2;

  if (lg > log(randu()) ) {
    beta_new = c;
    *acc = *acc + 1;
  }

  return beta_new;
}

double update_zeta (double zeta, double mu, vec theta, double a, double b, double cs, int* acc) {
  int n = theta.size();
  double lg, lg1, lg2, c, zeta_new;

  zeta_new = zeta;
  c = randn() * cs + zeta; // c=candidate, cs = candidiate sigma

  if (c > 0) {
    lg1 = (a-1)*log(c)    + n*c*log(c/mu)       - n*lgamma(c)     -    c*(b+sum(theta)/mu) +    c*sum(log(theta));
    lg2 = (a-1)*log(zeta) + n*zeta*log(zeta/mu) - n*lgamma(zeta)  - zeta*(b+sum(theta)/mu) + zeta*sum(log(theta));
    lg = lg1 - lg2;

    if ( lg > log(randu()) ) {
      zeta_new = c;
      *acc = *acc + 1;
    }
  }

  return zeta_new;
}


double update_mu (double zeta, vec theta, double a, double b) {
  int n = theta.size();
  return 1/rgamma( 1, a + n*zeta, 1 / (b + sum(theta)*zeta) )[0]; // shape, scale
}


//[[Rcpp::export]]
List hier(vec y, vec x, double a_zeta, double b_zeta, double a_mu, double b_mu,
    double m_beta, double s2_beta, double cs_zeta, double cs_beta, int B) { 

  double s_xy = sum(x%y);
  int acc_beta=0, acc_zeta=0, n=y.size();
  mat theta = ones<mat>(B,n);
  vec zeta, mu, beta, tb;
  mu = zeta = ones<vec>(B);
  beta = zeros<vec>(B);
  List ret;

  for (int b=1; b<B; b++) {
    theta.row(b) = update_theta(vectorise(theta.row(b-1)), zeta[b-1], mu[b-1], beta[b-1], x, y);
    tb = vectorise( theta.row(b) );
    beta[b] = update_beta(beta[b-1], tb, s_xy, m_beta, s2_beta, y, x, cs_beta, &acc_beta); // mh
    zeta[b] = update_zeta (zeta[b-1], mu[b-1], tb, a_zeta, b_zeta, cs_zeta, &acc_zeta); // mh
    mu[b] = update_mu(zeta[b], tb, a_mu, b_mu); // Gibbs
    cout << "\r Progress: " << b*100/B << "%";
  }

  ret["mu"] = mu;
  ret["zeta"] = zeta;
  ret["beta"] = beta;
  ret["theta"] = theta;
  ret["acc_zeta"] = acc_zeta * 1.0 / B;
  ret["acc_beta"] = acc_beta * 1.0 / B;

  return ret;
}

