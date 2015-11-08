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

//[[Rcpp::export]]
double lg0(double x) {
  //return -(log(2*pi) + x*x) / 2;
  return ldnorm(x,0,3);
}


// [[Rcpp::export]]
double wsample( vec x, vec prob ) {
  NumericVector ret = Rcpp::RcppArmadillo::sample(
      as<NumericVector>(wrap(x)), 1, true, 
      as<NumericVector>(wrap(prob))) ;
  return ret[0] ;
}


//[[Rcpp::export]]
mat auxGibbsCpp (vec y, double a, double s, double cs, int B) {
  cout << "Begin Gibbs..." << endl;

  int n = y.size();
  int acc_t = 0;
  mat theta; theta.set_size(B,n);
  vec G0 = randn(B) * 3;

  int k;
  double p, ti, utj, u, cand;
  double lg1, lg2, lg;
  vec ut_wo_ti, tb, probs, ts, ut;

  for (int b=1; b<B; b++) {
    theta.row(b) = theta.row(b-1);
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
      for (int d=0; d<k; d++) {
        ts[d] = ut_wo_ti[d];
      }
      ts[k] = p;

      for (int j=0; j<k; j++) {
        utj = ut_wo_ti[j];
        probs[j] = sum( utj == tb ) * f(y[i], utj, s);
      }

      theta(b,i) = wsample( ts , probs );
    }

    // Update theta | y
    tb = vectorise( theta.row(b) );
    ut = unique( tb );
    for (int j=0; j<ut.size(); j++) {
      u = ut[j];
      uvec ind = find( tb == u );
      cand = randn() * cs + u;
      lg1 = 0;
      lg2 = 0;
      for (int i=0; i<ind.size(); i++) {
        lg1 += lf(y[ind[i]], cand, s) + lg0(cand);
        lg2 += lf(y[ind[i]],    u, s) + lg0(u);
      }
      lg = lg1 - lg2;

      if ( lg > log(randu()) ) {
        // acc_t += 1;
        for (int i=0; i<ind.size(); i++) {
          theta( b, ind[i] ) = cand;
        }
      }
    }

    // Print Progress
    if ( b % (B/20) == 0 ) { // Progress not pring for some reason...
      cout << "\rProgress: " << round(100*b/B) << "%";
    }
  }
  
  cout << "\nDone.\n" << endl;
  return theta;
}

//[[Rcpp::export]]
uvec test(vec x, double y) {
  uvec uv = find(x==y);
  int m = uv.size();
  return uv;
}
