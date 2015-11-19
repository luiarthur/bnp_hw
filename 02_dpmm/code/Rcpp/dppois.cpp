//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <dp/dp.h>

using namespace std;
using namespace arma;
using namespace Rcpp ;

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
double update_mu (double mu, double tau, vec theta, double a, double b, double cs,
    int* acc) {
  // [mu] [prod_{t=uniqueThetas} g_0(t|mu)]
  vec ut = unique( theta );
  int n = ut.size();
  double lg, lg1, lg2, c, mu_new;

  mu_new = mu;
  c = randn() * cs + mu; // c=candidate, cs = candidiate sigma

  if (c > 0) {
    lg1 = (a-1) * log(c)  - b*c  + n*c*log(tau)  - n*lgamma(c)  + c*sum(log(ut));
    lg2 = (a-1) * log(mu) - b*mu + n*mu*log(tau) - n*lgamma(mu) + mu*sum(log(ut));
    lg = lg1 - lg2;

    if ( lg > log(randu()) ) {
      mu_new = c;
      *acc = *acc + 1;
    }
  }

  return mu_new;
  //return 1.0;
}

double rand_beta (double a, double b) { // Only for one draw. Write a function for vectors.
  //randg is in Armadillo, requires c++11
  //the Rcpp rgamma, dbeta are faster for vectors and scalars
  //Julia Distribution is faster than R for draws
  double z1 = randg(1, distr_param(a,1.0))[0];
  double z2 = randg(1, distr_param(b,1.0))[0];
  return z1 / (z1+z2);
}

double update_tau (double mu, vec theta, double a, double b) {
  // [tau] [prod_{t=uniqueThetas} g_0(t|tau)]
  vec ut = unique( theta );
  int n = ut.size();
  return rgamma( 1, a + mu * n, 1 / (b + sum(ut)) )[0]; // shape, sc
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
  
  ind = wsample({0,1}, {c-1, n*d});
  if (ind == 0) {
    alpha_new = rgamma(1, c, 1/d )[0]; // shape, rate
  } else {
    alpha_new = rgamma(1, c-1, 1/d )[0];
  }

  return alpha_new;
}

double update_beta (double beta, vec theta, double m, double s, vec y, vec x, 
    double cs, int* acc) {
    // [beta[ [prod_{i=1:n} (k(y_i|t_i,beta) ]

  double lg1, lg2, lg, c, beta_new;

  beta_new = beta;
  c = randn() * cs + beta;

  lg1 = -pow((c-m)/s,2) / 2    + sum(c * y % log(theta)    % x - theta % exp(c*x));
  lg2 = -pow((beta-m)/s,2) / 2 + sum(beta * y % log(theta) % x - theta % exp(beta*x));
  lg = lg1 - lg2;

  if (lg > log(randu()) ) {
    beta_new = c;
    *acc = *acc + 1;
  }

  return beta_new;
}

double lpois(double x, double lam) { // without normalizing constant!
  return -lam * x*log(lam) - lgamma(x+1);
}

mat update_theta (vec theta, double mu, double tau, double beta, double alpha, vec x, vec y) {
  vec theta_new = theta;
  vec t_xi, ut_xi;
  int J, ind, n = y.size(), ns;
  vec ut, probs, ys, xs, lprobs;
  uvec uv, inds;
  double q0, lq0, n_xi, q;
  double a,b;

  // Update theta_i | theta_{-i}
  for (int i=0; i<n; i++) {
    lq0 = mu * (-x[i]*beta+tau) +lgamma(y[i]+mu) - lgamma(y[i]+1) - lgamma(mu);

    t_xi = theta_new( find(linspace(0,n-1,n) != i) );
    ut_xi = unique(t_xi);
    J = ut_xi.size();
    lprobs = zeros(J+1);
    lprobs[J] = log(alpha) + lq0; // prob of drawing new theta

    for (int j=0; j<J; j++) {
      uv = find(ut_xi[j] == t_xi);
      n_xi = uv.size();
      q = lpois( y[i], ut_xi[j] * exp(x[i]*beta) );
      lprobs[j] = log(n_xi) * q; // probability of drawing new theta
    }

    probs = exp(lprobs - max(lprobs));
    //cout << "probs: "<<probs <<endl;
    ind = wsample(linspace(0,J,J+1), probs);
    if (ind == J) {
      a = y[i]+mu;
      b = 1/( tau+exp(x[i]*beta) );
      theta_new[i] = rgamma(1, a, b)[0]; //shape, scale
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

    a = sum(ys)+mu;
    b = 1 / (sum(exp(xs*beta)) + tau);
    theta_new.elem(inds) = zeros(ns) + rgamma(1,a,b)[0]; // shape, scale
  }

  return reshape(theta_new,1,32);
}

//[[Rcpp::export]]
List dppois(vec y, vec x, double a_mu, double b_mu, double a_tau, double b_tau,
    double a_alpha, double b_alpha, double m_beta, double s_beta, 
    double cs_mu, double cs_beta, int B) { 
  // also need cand_sigs
  // all b_* are rates in the gamma distribution

  int acc_beta=0, acc_mu=0, n=y.size();
  mat theta = ones<mat>(B,n);
  vec mu, tau, beta, alpha, tb;
  mu = tau = alpha = ones<vec>(B);
  beta = zeros<vec>(B);
  List ret;

  // Start MCMC
  for (int b=1; b<B; b++) {
    // Update thetas
    theta.row(b) = update_theta(vectorise(theta.row(b-1)), mu[b-1], tau[b-1],
        beta[b-1], alpha[b-1], x, y);
    tb = vectorise( theta.row(b) );

    beta[b] = update_beta(beta[b-1], tb, m_beta, s_beta, y, x, cs_beta, &acc_beta);
    mu[b] = update_mu(mu[b-1], tau[b-1], tb, a_mu, b_mu, cs_mu, &acc_mu);
    tau[b] = update_tau (mu[b], tb, a_tau, b_tau);
    alpha[b] = update_alpha(alpha[b], tb, a_alpha, b_alpha, n);

    cout << "\r" << b * 100 / B <<"%";
  } // End of MCMC

  ret["mu"] = mu;
  ret["tau"] = tau;
  ret["beta"] = beta;
  ret["alpha"] = alpha;
  ret["theta"] = theta;
  ret["acc_mu"] = acc_mu * 1.0 / B;
  ret["acc_beta"] = acc_beta * 1.0 / B;

  return ret;
}

// Need a posterior predictive function
