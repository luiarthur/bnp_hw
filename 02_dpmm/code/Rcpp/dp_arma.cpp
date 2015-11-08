//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>

using namespace std;
using namespace arma;
using namespace Rcpp ;

const double pi = 3.141592653589793238462;

//[[Rcpp::export]]
mat rdir(int N, vec a) {
  int K = a.size(); 
  mat x(N,K);
  double rowsum_i;
  
  for (int i=0; i<N; i++) {
    rowsum_i = 0;
    for (int j=0; j<K; j++) {
      x(i,j) = R::rgamma(a[j],1.0); // shape & scale
      rowsum_i += x(i,j);
    }
    x.row(i) = x.row(i) / rowsum_i;
  }

  return x;
}


//[[Rcpp::export]]
double ldnorm(double x, double m, double s) {
  return -(log(s) + .5*log(2*pi)) - pow( (x - m) / s , 2) / 2;
}

//[[Rcpp::export]]
double lf(double x, double t, double s) {
  return ldnorm(x,t,s);
}

//[[Rcpp::export]]
double f(double x, double t, double s) {
  return exp(ldnorm(x,t,s));
}


// [[Rcpp::export]]
double wsample( vec x, vec prob ) {
  NumericVector ret = Rcpp::RcppArmadillo::sample(
      as<NumericVector>(wrap(x)), 1, true, 
      as<NumericVector>(wrap(prob))) ;
  return ret[0] ;
}


// [[Rcpp::export]]
int test(vec x, vec y) {
  return sum(x==y);
}

//[[Rcpp::export]]
mat gibbs (vec y, double a, double s, double cs, int burn, int B) {
  int n = y.size();
  int acc_t = 0;
  mat theta = ones(B,n);
  vec G0 = randn(B);

  int k;
  double p, ti, utj;
  vec ut_wo_ti, tb, probs, ts;

  for (int b=1; b<B; b++) {
    // Update theta_i | theta_{-i}
    for (int i=0; i<n; i++) {
      tb = vectorise( theta.row(b) );
      ti = tb[i];
      ut_wo_ti = unique( tb( find(linspace(0,n-1,n) != i) ) );
      k = ut_wo_ti.size();

      if ( any(ut_wo_ti == ti) ) {
        p = G0[b];
      }

      probs = zeros(k+1);
      probs[k] = a * f(y[i], p, s);

      //ts = [ ut_wo_ti, p ]
      ts = zeros(k+1);
      for (int d=0; d<k+1; d++) {
        ts[d] = ut_wo_ti[d];
      }
      ts[k] = p;

      for (int j=0; j<k; j++) {
        utj = ut_wo_ti[j];
        probs[j] = sum( utj == tb ) * f(y[i], utj, s);
      }

      theta(b,i) = wsample( ts , probs );

      if ( b % (B/20) == 0 ) {
        cout << "\r" << round(100*b/B) << "%";
      }
    }
    // Update theta | y
  }
  
  return theta;
}


/*
vec x(n);

for (int i=0; i<n; i++) {
  x[i] = R::rgamma(1.0,1.0);
}
*/
