//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <algorithm> // unique
#include <../../../02_dpmm/code/Rcpp/dp/dp.h> //wsample

using namespace std;
using namespace arma;
using namespace Rcpp;

// I want to eventually add Time Series, and FIX uniqueRows()

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

//[[Rcpp::export]]
vec mvrnorm(vec M, mat S) {
  int n = M.n_rows;
  mat e = randn(n);
  return vectorise(M + chol(S).t()*e);
}

// For reference
//void bubbleSort(int arr[], int n) {
//  bool swapped = true;
//  int j = 0;
//  int tmp;
//  while (swapped) {
//    swapped = false;
//    j++;
//    for (int i = 0; i < n - j; i++) {
//      if (arr[i] > arr[i + 1]) {
//        tmp = arr[i];
//        arr[i] = arr[i + 1];
//        arr[i + 1] = tmp;
//        swapped = true;
//      }
//    }
//  }
//}


bool r1GTr2(vec r1, vec r2, int c) { // The user should set c = 0.
  bool gt;
  int n = r1.size();

  if (r1[c] > r2[c]) {
    gt = true;
  } else if (r1[c] < r2[c]) {
    gt = false;
  } else if (c == n) { // if last col and still same, return r1 > r2 is true.
    gt = true;
  } else {
    gt = r1GTr2(r1, r2, c+1);
  } 

  return gt;
}

mat sortRows(mat X) { // CHECK OUT: Lexiographical sorts
  int n = X.n_rows;
  bool swapped = true;
  int j = 0;
  mat tmp;

  while (swapped) {
    swapped = false;
    j++;
    for (int i = 0; i < n - j; i++) {
      //if (arr[i] > arr[i + 1]) {
      if ( r1GTr2(vectorise(X.row(i)), vectorise(X.row(i+1)), 0) ) {
        tmp = X.row(i);
        X.row(i) = X.row(i+1);
        X.row(i+1) = tmp;
        swapped = true;
      }
    }
  }
  return X;
}

//[[Rcpp::export]]
mat uniqueRows2(mat X) { // Works, but not OPTIMAL!!!
  int n = X.n_rows;
  X = sortRows(X);

  for (int i=1; i<n; i++) {
    if ( all(X.row(i) == X.row(i-1)) ) {
      X.shed_row(i);
      i--;
      n--;
    }
  }
  return X;
}

//[[Rcpp::export]]
mat uniqueRows(mat Y) { // This works and is fast, but not certain of accuracy
  int n = Y.n_rows;
  int k = Y.n_cols;
  vec x = Y * randu(k);
  vec idx = zeros<vec>(n);

  map<double, int> uRow;
  for (int i=0; i<n; i++) {
    if ( uRow.find(x[i]) == uRow.end() ) { // key doesn't exist
      uRow[x[i]] = 1;
      idx[i] = 1;
    } else { // key found
      uRow[x[i]] = uRow[x[i]] + 1;
    }
  }
  
  return Y.rows( find(idx==1) );
}

//[[Rcpp::export]]
mat uniqueRows1(mat X) { // Works, but Not OPTIMAL!!!
  int n = X.n_rows;
  int i = 0;
  
  while(i < n) {
    for (int j=i+1; j<n; j++) {
      if ( all(X.row(i) == X.row(j)) ){
        X.shed_row(j);
        j--;
        n--;
      }
    }
    i++;
  }

  return X;
}

//[[Rcpp::export]]
double update_alpha (double alpha, int Tstar, double a, double b, int T) {
  double eta, c, d;
  int ind;
  double alpha_new;

  eta = rbeta(1, alpha + 1, T )[0];
  c = a + Tstar;
  d = b - log(eta); // rate
  
  ind = wsample({0.0, 1.0}, {c-1, T*d});
  if (ind == 0.0) {
    alpha_new = rgamma(1, c, 1/d )[0]; // shape, scale
  } else {
    alpha_new = rgamma(1, c-1, 1/d )[0];
  }

  return alpha_new;
}


double update_beta (double sum_y, mat theta, double tau2, double m, double s2, int T, int n) {
  double denom = tau2 + T*n*s2;
  double m_new = (m*tau2 + s2 * ( sum_y - sum(sum(theta)) )) / denom;
  double s_new = sqrt( s2*tau2 / denom );

  return rnorm(1, m_new, s_new)[0];
}

//[[Rcpp::export]]
double update_tau2 (double a, double b, int n, int T, mat y, mat theta, double beta) {
  mat one_Tn(T,n);
  mat M = y - theta - beta*one_Tn;
  double sum_mm = 0;
  mat m;

  //for (int t=0; t<T; t++) {
  //  m = M.row(t).t() * M.row(t);
  //  sum_mm += m(0,0);
  //}
  sum_mm = trace(M.t() * M);

  double a_new = a + n*T/2;
  double rate_new = b + sum_mm/2;
  double scale_new = 1 / rate_new;

  return 1 / rgamma(1, a_new, scale_new)[0];
}

//[[Rcpp::export]]
vec testrg(int n, double a, double b) {
  return rgamma(n, a, b);
}

//[[Rcpp::export]]
mat Hn (double phi, mat D) {
  return exp(-phi * D);
}

//[[Rcpp::export]]
double update_sig2 (double a, double b, int Tstar, int n, mat theta, mat H) {
  mat ut = uniqueRows(theta);
  mat mt;
  double sum_mt = 0;
  mat Hi = H.i();

  //for (int j=0; j<Tstar; j++) {
  //  mt = ut.row(j) * Hi *ut.row(j).t();
  //  sum_mt += mt(0,0);
  //}
  sum_mt = trace(ut * Hi * ut.t()); // concise but slower.

  double a_new = a + .5*n*Tstar;
  double rate_new = b + .5*sum_mt;
  double scale_new = 1 / rate_new;

  return 1 / rgamma(1, a_new, scale_new)[0]; 
}

// Broken!!!???
//double update_phi (double a, double b, mat D, int Tstar, mat theta, double sig2, int L, double cs, int* acc, double phi_curr) {
double update_phi (double a, double b, mat D, int Tstar, mat theta, double sig2, int L) {
  double phi_new; 
  vec lgs = zeros<vec>(L);
  vec phi = linspace<vec>(a,b,L);
  double val, sign;
  mat ut = uniqueRows(theta);
  mat Hp;
  mat mt;
  double sum_mt;
  vec probs;
  double ldet;

  for (int i=0; i<L; i++) {
    sum_mt = 0;
    Hp = Hn(phi[i],D); // pow(H, phi[i]);
    log_det(val, sign, Hp);
    ldet = val; // covariance Matrices ALWAYS have positive determinants!!!
    sum_mt = trace(-ut * Hp.i() * ut.t());
    lgs[i] = -.5*Tstar * ldet + sum_mt/(2*sig2);
  }
  
  probs = exp(lgs - max(lgs));
  phi_new = wsample(phi, probs);
  /*
  double lgcand, lgcurr;
  double cand = randn() * cs + phi_curr;
  double ldetCurr, ldetCand;
  phi_new = phi_curr;
  if (cand > 0) {
    Hp = Hn(cand,D);
    log_det(ldetCand, sign, Hp);
    lgcand = -.5*Tstar * ldetCand - trace(-ut * Hp.i() * ut.t())/(2*sig2);

    Hp = Hn(phi_curr,D);
    log_det(ldetCurr, sign, Hp);
    lgcurr = -.5*Tstar * ldetCurr - trace(-ut * Hp.i() * ut.t())/(2*sig2);
    if (lgcand - lgcurr > log(randu())) {
      phi_new = cand;
      *acc = *acc + 1; 
    }
  }
  */

  return phi_new;
}


//[[Rcpp::export]]
double ldmvnorm(vec x, vec m, mat S) {
  int n = x.size();
  double ldet, val, sign;
  vec out;
  vec b = x - m;
  log_det(val,sign,S);
  ldet = val;
  out = -n/2*log(2*pi) -.5*ldet -.5*b.t()*S.i()*b;

  return out[0];
}


//[[Rcpp::export]]
uvec matchRows (mat X, vec v) {
  uvec uv;
  vec vv = zeros<vec>(X.n_rows);
  mat out;

  for (int i=0; i<X.n_rows; i++) {
    vv[i] = 1 * all(vectorise(X.row(i)) == v);
  }
  
  uv = find( vv == 1);
  return uv;
}

//[[Rcpp::export]]
mat retMatched (mat X, vec v) {
  uvec uv;
  uv = matchRows(X,v);
  return X.rows(uv);
}

mat update_theta(double alpha, double sig2, double phi, double tau2, double beta, mat y,
    mat theta, int T, int n, mat H) {

  double ldet, val, sign;
  double ldetH, val2, sign2;
  mat theta_new = theta;
  double lq, lq0;
  mat t_xt, ut_xt, ut;
  int J, Tj, ind;
  mat tau2I = tau2 * eye<mat>(n,n);
  vec Tind = linspace(0,T-1,T);
  vec lprobs, probs;
  uvec uv, inds; 
  mat S, Sj;
  vec mu, muj;
  mat Mu;
  mat In = eye<mat>(n,n);
  mat Hi = H.i();
  mat tmp;
  vec sum_mt;

  S = (In/tau2 + Hi /sig2).i();
  log_det(val,sign,S);
  ldet = val;

  log_det(val2,sign2,H);
  ldetH = val2;

  for (int t=0; t<T; t++) {
    mu = vectorise(y.row(t)) - ones(n)*beta;
    tmp = .5*ldet + .5*(mu.t()*(In-S/tau2)*mu)/tau2 - n/2 * log(2*pi*tau2*sig2) - .5*ldetH;
    lq0 = tmp(0,0);
    t_xt = theta_new.rows( find(Tind != t) );
    ut_xt = uniqueRows(t_xt);
    J = ut_xt.n_rows;
    lprobs = zeros(J+1);
    probs = zeros(J+1);
    lprobs[J] = log(alpha) + lq0;

    for (int j=0; j<J; j++) {
      uv = matchRows(t_xt, vectorise(ut_xt.row(j)) );
      Tj = uv.size();
      lq = ldmvnorm(vectorise(y.row(t)), vectorise(theta_new.row(j)) + ones(n)*beta, tau2I); // changed this
      lprobs[j] = log(Tj) + lq;
    }

    probs = exp(lprobs - max(lprobs));
    ind = wsample( linspace(0,J,J+1), probs);
    if (ind == J) {
      theta_new.row(t) = reshape( mvrnorm(S*mu/tau2,S), 1, n); // changed this
    } else {
      theta_new.row(t) = ut_xt.row(ind);
    }
  }

  // update theta | y
  ut = uniqueRows(theta_new);
  J = ut.n_rows;
  for (int j=0; j<J; j++) {
    inds = matchRows(theta_new, vectorise( ut.row(j) ));
    Tj = inds.size();
    Sj = (Tj/tau2 * In + Hi/sig2).i();
    Mu = y.rows(inds) - ones(Tj,n)*beta;
    muj = Sj/tau2 * vectorise( sum(Mu,0) ); // sum each column

    theta_new.each_row(inds) = reshape(mvrnorm(muj, Sj),1,n); // Check this !!!
  }

  return theta_new;
}

//[[Rcpp::export]]
List sdp(mat Y, mat D, double beta_m, double beta_s2,
    double tau2_a, double tau2_b, double alpha_a, double alpha_b,
    double sig2_a, double sig2_b, double phi_a, double phi_b, int L, int B) {
  //phi_a = 0, phi_b = 3/(d*max(||si-sj||)), where d is small ~ .01. => a,b about 0,25

  int T = Y.n_rows; // 20 years
  int n = Y.n_cols; // 100 locations
  double sum_y = sum( sum(Y) );
  vec beta = zeros<vec>(B);
  vec tau2 = ones<vec>(B);
  vec alpha = ones<vec>(B);
  vec sig2 = ones<vec>(B);
  vec phi = ones<vec>(B);
  cube theta = zeros<cube>(T,n,B);
  int Tstar;
  mat utheta;
  int acc_sig2 = 0;
  int acc_phi = 0;
  List ret;
  mat tb;

  for (int b=1; b<B; b++) {
    theta.slice(b) = update_theta(alpha[b-1], sig2[b-1], phi[b-1], tau2[b-1], 
        beta[b-1], Y, theta.slice(b-1), T, n, Hn(phi[b-1],D));// ERROR.
    tb = theta.slice(b);
    utheta = uniqueRows(tb);
    Tstar = utheta.n_rows;

    beta[b] = update_beta(sum_y, tb, tau2[b-1], beta_m, beta_s2, T, n); // check
    tau2[b] = 13;//update_tau2(tau2_a, tau2_b, n, T, Y, tb, beta[b]); // check
    alpha[b] = update_alpha(alpha[b-1], Tstar, alpha_a, alpha_b, T);
    sig2[b] = 100;//update_sig2 (sig2_a, sig2_b, Tstar, n, tb, Hn(phi[b-1],D)); // check
    phi[b] = 1;//update_phi(phi_a, phi_b, D, Tstar, tb, sig2[b], L); // check. ERROR.
    //phi[b] = update_phi(phi_a, phi_b, D, Tstar, tb, sig2[b], L,phi_cs,&acc_phi,phi[b-1]); // check. ERROR.

    Rcout << "\rProgress: " << b << "/" << B;
  }
  cout << endl;

  ret["beta"] = beta;
  ret["tau2"] = tau2;
  ret["alpha"] = alpha;
  ret["sig2"] = sig2;
  ret["phi"] = phi;
  ret["theta"] = wrap(theta);
  ret["acc_phi"] = acc_phi*1.0 / B;

  return ret;
}
