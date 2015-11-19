//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h> // For wsample
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
      x(i,j) = randg(1,distr_param(a[j],1.0))[0]; // shape & scale
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

// [[Rcpp::export]]
double wsample( vec x, vec prob ) {
  NumericVector ret = Rcpp::RcppArmadillo::sample(
      as<NumericVector>(wrap(x)), 1, true, 
      as<NumericVector>(wrap(prob))) ;
  return ret[0] ;
}
