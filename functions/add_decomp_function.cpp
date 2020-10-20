#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//Sum_Mat.cpp
//[[Rcpp::export]]
NumericMatrix Sum_Mat(List V){
  Function f("Reduce");
  return f(Named("f") = "+", Named("x") = V);
}

//prod_AB.cpp
//[[Rcpp::export]]
arma::mat prod_AB(arma::mat A, arma::mat B, bool A_t, bool B_t) {
  if(A_t == true) {
    if(B_t == true) return A.t()*B.t();
    if(B_t == false) return A.t()*B;
  }
  if(A_t == false) {
    if(B_t == true) return A*B.t();
    if(B_t == false) return A*B;
  }
}

//update_eta.cpp
//[[Rcpp::export]]
arma::mat update_eta(double delta, arma::mat sum_VV, arma::mat W, arma::mat Theta, 
                     arma::mat sum_VU, arma::mat sum_V, arma::mat Y, arma::mat Z, arma::mat sum_VE) {
  double dim = sum_VV.n_rows;
  arma::mat diag(dim, dim, arma::fill::eye);
  arma::mat First = sum_VV + diag;
  arma::mat Second = W + delta*Theta + sum_VU + delta*sum_V.t()*(Y - Z) - delta * sum_VE;
  arma::mat output = inv(First) * Second;
  return output*(1/delta);
}

//update_theta.cpp
//[[Rcpp::export]]
NumericMatrix update_theta(double delta, double lambda_2, NumericMatrix Eta, NumericMatrix W) {
  double n = W.nrow();
  double m = W.ncol();
  NumericMatrix Theta(n, m);
  for(int g=0; g<m; g++) {
    NumericVector r_g = Eta(_, g) - W(_, g)/delta;
    NumericVector value = 1 - (lambda_2/(delta * abs(r_g)));
    Theta(_, g) = ifelse(value >0, value * r_g, 0);
  }
  return Theta;
}

//update_alpha.cpp
//[[Rcpp::export]]
List update_alpha(double delta, double lambda_1, arma::mat Y, arma::mat X, arma::mat H, 
                  List V, List E, List U) {
  Function f("Reduce");
  double b = V.size();
  arma::mat Sum_V = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = V)));
  arma::mat Sum_E = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = E)));
  arma::mat Sum_U = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = U)));
  arma::mat obj = Y - (1/b)*Sum_V*H - (1/b)*Sum_E - (1/(delta*b))*Sum_U;
  double n = obj.n_rows;
  double p = obj.n_cols;
  arma::mat P;
  arma::mat Q;
  arma::vec d;
  arma::svd(P, d, Q, obj);
  if(p > n) {Q = Q.cols(0, n-1);}
  if(n > p) {P = P.cols(0, p-1);}
  arma::vec d_new = d - lambda_1/(delta*b);
  NumericVector diag_entry = ifelse(as<NumericVector>(wrap(d_new)) > 0, as<NumericVector>(wrap(d_new)), 0);
  arma::mat D_new = diagmat(as<arma::vec>(wrap(diag_entry)));
  arma::mat Z_new = P*D_new*Q.t();
  arma::mat alpha_new = inv(X.t()*X)*X.t()*Z_new;
  return List::create(Named("Z_new") = Z_new, Named("alpha_new") = alpha_new);
}