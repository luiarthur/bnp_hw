//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <dp/dp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
vec test (int n) {
  //uvec inds = find( x == 3 );
  //x.elem(inds) = ones(inds.size()) + 3.0;
  //return x;
  return rbeta(n,2,3);
} 

// http://adv-r.had.co.nz/Rcpp.html#rcpp-sugar
// Enable C++11 in R: Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
// This is used for randg(shape,scale), vec v = {1,2,3};

// Functions for this assignment
double update_zeta (double zeta, double mu, vec theta, double a, double b, double cs, int* acc) {
  vec ut = unique( theta );
  int n = ut.size();
  double lg, lg1, lg2, c, zeta_new;

  zeta_new = zeta;
  c = randn() * cs + zeta; // c=candidate, cs = candidiate sigma

  if (c > 0) {
    lg1 = (a-1)*log(c)    + n*c*log(c/mu)       - n*lgamma(c)     - c*(b+sum(ut)/mu)    + c*sum(log(ut));
    lg2 = (a-1)*log(zeta) + n*zeta*log(zeta/mu) - n*lgamma(zeta)  - zeta*(b+sum(ut)/mu) + zeta*sum(log(ut));
    lg = lg1 - lg2;

    if ( lg > log(randu()) ) {
      zeta_new = c;
      *acc = *acc + 1;
    }
  }

  return zeta_new;
}

double update_mu (double zeta, vec theta, double a, double b) {
  vec ut = unique( theta );
  int n = ut.size();
  return 1/rgamma( 1, a + n*zeta, 1 / (b + sum(ut)*zeta) )[0]; // shape, sc
}

double update_alpha (double alpha, vec theta, double a, double b, int n) {
  double eta, c, d;
  int ns, ind;
  double alpha_new;
  vec ut = unique(theta);

  eta = rbeta(1, alpha + 1, n )[0];
  ns = ut.size();
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

double update_beta (double beta, vec theta, double s_xy, double m, double s2, vec y, vec x, double cs, int* acc) {
    // [beta[ [prod_{i=1:n} (k(y_i|t_i,beta) ]

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

double lpois(double x, double lam) { // without normalizing constant!
  return x*log(lam) - lam -lgamma(x+1);
}

mat update_theta (vec theta, double zeta, double mu, double beta, double alpha, vec x, vec y) {
  vec theta_new = theta;
  vec t_xi, ut_xi;
  int J, ind, n = y.size(), ns;
  vec ut, probs, ys, xs, lprobs;
  uvec uv, inds;
  double lq0, nj, lq;
  double a,b;

  // Update theta_i | theta_{-i}
  for (int i=0; i<n; i++) {
    lq0 = lgamma(y[i]+zeta) - lgamma(y[i]+1) - lgamma(zeta) - y[i]*log( 1 + zeta/mu*exp(-beta*x[i]) ) - zeta*log( 1 + exp(beta*x[i]) * mu/zeta );

    t_xi = theta_new( find(linspace(0,n-1,n) != i) );
    ut_xi = unique(t_xi);
    J = ut_xi.size();
    lprobs = zeros(J+1);
    lprobs[J] = log(alpha) + lq0; // prob of drawing new theta

    for (int j=0; j<J; j++) {
      uv = find(ut_xi[j] == t_xi);
      nj = uv.size();
      lq = lpois( y[i], ut_xi[j] * exp(x[i]*beta) );
      lprobs[j] = log(nj) + lq; // probability of drawing new theta
    }

    probs = exp(lprobs - max(lprobs));
    ind = wsample(linspace(0,J,J+1), probs);
    if (ind == J) {
      a = y[i] + zeta;
      b = zeta/mu + exp(x[i]*beta); // = rate
      theta_new[i] = rgamma(1, a, 1/b)[0]; //shape, scale
    } else {
      theta_new[i] = ut_xi[ind];
    }
  }

  // Update theta | y
  ut = unique(theta_new);
  J = ut.size();
  for (int j=0; j<J; j++) {
    inds = find( theta_new == ut[j] );
    ys = y(inds);
    xs = x(inds);
    ns = inds.size();

    a = sum(ys) + zeta;
    b = sum(exp(xs*beta)) + zeta/mu; // = rate
    theta_new.elem(inds) = zeros(ns) + rgamma(1,a,1/b)[0]; // shape, scale
  }

  return reshape(theta_new,1,32);
}

//[[Rcpp::export]]
List dppois(vec y, vec x, double a_zeta, double b_zeta, double a_mu, double b_mu,
    double a_alpha, double b_alpha, double m_beta, double s2_beta, 
    double cs_zeta, double cs_beta, int B) { 
  // all b_* are rates in the gamma distribution

  double s_xy = sum(x%y);
  int acc_beta=0, acc_zeta=0, n=y.size();
  mat theta = ones<mat>(B,n);
  vec zeta, mu, beta, alpha, tb;
  mu = zeta = alpha = ones<vec>(B);
  beta = zeros<vec>(B);
  List ret;

  // Start MCMC
  for (int b=1; b<B; b++) {
    theta.row(b) = update_theta(vectorise(theta.row(b-1)), zeta[b-1], mu[b-1], beta[b-1], alpha[b-1], x, y); // Gibbs. Done.
    tb = vectorise( theta.row(b) );

    beta[b] = update_beta(beta[b-1], tb, s_xy, m_beta, s2_beta, y, x, cs_beta, &acc_beta); // mh
    zeta[b] = update_zeta (zeta[b-1], mu[b-1], tb, a_zeta, b_zeta, cs_zeta, &acc_zeta); // mh
    mu[b] = update_mu(zeta[b], tb, a_mu, b_mu); // Gibbs
    alpha[b] = update_alpha(alpha[b], tb, a_alpha, b_alpha, n); // Gibbs

    cout << "\r" << b * 100 / B <<"%";
  } // End of MCMC

  ret["mu"] = mu;
  ret["zeta"] = zeta;
  ret["beta"] = beta;
  ret["alpha"] = alpha;
  ret["theta"] = theta;
  ret["acc_zeta"] = acc_zeta * 1.0 / B;
  ret["acc_beta"] = acc_beta * 1.0 / B;

  return ret;
}

// Need a posterior predictive function
