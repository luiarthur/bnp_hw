//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>

using namespace std;

//[[Rcpp::export]]
arma::mat a1 (arma::mat x) {
  return(x.t() * x) ;
}


//[[Rcpp::export]]
arma::mat rdir(int N, arma::vec a) {
  int K = a.size(); 
  arma::mat x(N,K);
  double rowsum_i;
  
  for (int i=0; i<N; i++) {
    rowsum_i = 0;
    for (int j=0; j<K; j++) {
      x(i,j) = R::rgamma(a[j],2.0); // shape & scale
      rowsum_i += x(i,j);
    }
    x.row(i) = x.row(i) / rowsum_i;
  }

  return x;
}


//[[Rcpp::export]]
arma::vec gibbs (int n) {

  arma::vec x(n);

  for (int i=0; i<n; i++) {
    x[i] = R::rgamma(1.0,1.0);
  }

  return x;
}
