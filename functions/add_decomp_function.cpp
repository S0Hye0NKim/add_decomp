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
arma::mat update_theta(double delta, double lambda_2, double K, arma::mat Eta, arma::mat W) {
  double n = W.n_rows;
  double m = W.n_cols;
  double p = n/K - 1;
  arma::mat Theta(n, m);
  
  for(int g=0; g<m; g++) {
    arma::vec r_j = Eta.col(g) - (W.col(g)/delta);
    for(int j=0; j<(p+1); j++) {
      arma::vec sub_vec = r_j.subvec(K*j, K*(j+1)-1);
      double accum = 0;
      for(int i =0; i<K; i++) {
        accum += sub_vec[i]*sub_vec[i];
      }
      double norm = sqrt(accum);
      double value = 1 - (lambda_2/(delta*norm));
      if(value >= 0) {
        Theta.submat(K*j, g, K*(j+1)-1, g) = value * sub_vec;
      } else {
        Theta.submat(K*j, g, K*(j+1)-1, g) = arma::vec(K);
      }
    }
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

//update_e.cpp
//[[Rcpp::export]]
List update_e(double delta, List U, List V, arma::mat Y, arma::mat X, arma::mat A, arma::mat H, 
              NumericVector tau_seq) {
  double n = Y.n_rows;
  double b = U.size();
  double m = Y.n_cols;
  arma::mat Z = X*A;
  List e_new = List(b);
  for(int l = 0; l < b; l++) {
    arma::mat mat_E = arma::mat(n, m);
    arma::mat VH = as<arma::mat>(wrap(V[l]))*H;
    arma::mat Error = Y - Z - VH;
    arma::mat mat_U = U[l];
    for(int g = 0; g < m; g++) {
      for(int i = 0; i < n; i++) {
        if(Error(i, g) + mat_U(i, g)/delta < (tau_seq(l)-1)/(n*delta)) {
          mat_E(i, g) = Error(i, g) + mat_U(i, g)/delta - (tau_seq(l)-1)/(n*delta);
        } else if(Error(i, g) + mat_U(i, g)/delta > tau_seq(l)/(n*delta)) {
          mat_E(i, g) = Error(i, g) + mat_U(i, g)/delta - tau_seq(l)/(n*delta);
        } else {mat_E(i, g) = 0;}
      }
    }
    e_new[l] = mat_E;
  }
  return e_new;
}

//update_multipliers.cpp
//[[Rcpp::export]]
List update_multi(List U_old, arma::mat W_old, double delta, arma::mat Y, arma::mat X, arma::mat A, arma::mat H,
                  arma::mat Theta, List V, List E, bool is_u) {
  if(is_u == true) {
    double b = U_old.size();
    List U_new = List(b);
    for(int l = 0; l < b; l++) {
      arma::mat VH = as<arma::mat>(wrap(V[l]))*H;
      arma::mat mat_E = E[l];
      arma::mat mat_U_old = U_old[l];
      arma::mat mat_U = mat_U_old + delta * (Y - X*A - VH - mat_E);
      U_new[l] = mat_U;
    }
    return U_new;
  } else {
    arma::mat W_new = W_old + delta *(Theta - H);
    return List::create(W_new);
  }
}

//calculate_iteration_error.cpp
//[[Rcpp::export]]
double calc_err(List Old, List New) {
  Function f("Reduce");
  double b = Old.size();
  List FR_list = List(b);
  for(int l = 0; l < b; l++) {
    arma::mat Old_mat = Old[l];
    arma::mat New_mat = New[l];
    FR_list[l] = norm(Old_mat - New_mat, "fro");
  }
  double FR_sum = as<double>(wrap(f(Named("f") = "+", Named("x") = FR_list)));
  return FR_sum;
}

//calculate_kronecker_product_V.cpp
//[[Rcpp::export]]
List calc_V(arma::mat X, arma::mat Phi) {
  double n = X.n_rows;
  double p = X.n_cols;
  double K = Phi.n_cols;
  double b = Phi.n_rows;
  Function f("kronecker");
  List V = List(b);
  for(int l = 0; l<b; l++) {
    arma::mat mat_V = arma::mat(n, p*K);
    arma::rowvec Phi_vec = Phi.row(l);
    for(int i = 0; i<n; i++) {
      arma::rowvec X_vec = X.row(i);
      arma::rowvec v_i_l = as<arma::rowvec>(wrap(f(Named("X") = X_vec, Named("Y") = Phi_vec)));
      mat_V.row(i) = v_i_l;
    }
    V[l] = mat_V;
  }
  return V;
}