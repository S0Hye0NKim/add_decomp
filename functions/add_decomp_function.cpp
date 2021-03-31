#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
List add_decomp(double delta, double lambda_1, double lambda_2, double tol_error, int max_iter, arma::mat X,
                arma::mat Y, List V, arma::mat Phi, arma::mat theta_0, arma::mat Z_0, NumericVector tau_seq,
                bool weight) {
  Function f("Reduce");
  //delta = step size
  //lambda_1 = low rank penalty
  //lambda_2 = sparse penalty
  
  //initial value
  double n = X.n_rows;
  double p = X.n_cols - 1;
  double K = Phi.n_cols;
  double b = Phi.n_rows;
  double m = Y.n_cols;
  arma::mat eta_old = theta_0;
  arma::mat theta_old = theta_0;
  arma::mat Z_old = Z_0;
  List e_old(b);
  for(int l =0; l<b; l++) {
    arma::mat e_old_mat = Y  - Z_old - as<arma::mat>(wrap(V[l])) *eta_old;
    e_old[l] = e_old_mat;
  }
  List u_old(b);
  for(int l=0; l<b; l++) {
    arma::mat u_old_mat(n, m);
    u_old_mat.zeros(n, m);
    u_old[l] = u_old_mat;
  }
  arma::mat w_old;
  w_old.zeros((p+1)*K, m);
  arma::mat iter_error(max_iter, 6);
  arma::mat sum_V = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = V)));
  List VV_prod(b);
  for(int l=0; l<b; l++) {
    arma::mat V_mat = V[l];
    VV_prod[l] = V_mat.t() * V_mat;
  }
  arma::mat sum_VV = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = VV_prod)));
  
  //iteration
  for(int iter=0; iter<max_iter; iter++) {
    
    //Process for eta
    List Vu_prod(b);
    List Ve_prod(b);
    for(int l=0; l<b; l++) {
      arma::mat V_mat = V[l];
      arma::mat u_old_mat = u_old[l];
      arma::mat e_old_mat = e_old[l];
      Vu_prod[l] = V_mat.t() * u_old_mat;
      Ve_prod[l] = V_mat.t() * e_old_mat;
    }
    arma::mat sum_Vu = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = Vu_prod)));
    arma::mat sum_Ve = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = Ve_prod)));
    double dim = sum_VV.n_rows;
    arma::mat diag(dim, dim, arma::fill::eye);
    arma::mat First = sum_VV + diag;
    arma::mat Second = w_old + delta * theta_old + sum_Vu + delta*sum_V.t()*(Y-Z_old) - delta*sum_Ve;
    arma::mat output = inv(First) * Second;
    arma::mat eta_new = output*(1/delta);
    
    //Process for theta
    arma::mat theta_new((p+1)*K, m);
    for(int g = 0; g<m; g++) {
      arma::vec r_j = eta_new.col(g) - (w_old.col(g)/delta);
      arma::vec theta_tilde_j = theta_0.col(g);
      for(int j=0; j<(p+1); j++) {
        arma::vec r_j_g = r_j.subvec(K*j, K*(j+1)-1);
        arma::vec theta_tilde_j_g = theta_tilde_j.subvec(K*j, K*(j+1)-1);
        double accum_r_j_g = 0;
        double accum_theta_tilde_j_g = 0;
        for(int i=0; i<K; i++) {
          accum_r_j_g += r_j_g[i]*r_j_g[i];
          accum_theta_tilde_j_g += theta_tilde_j_g[i]*theta_tilde_j_g[i];
        }
        double norm_r_j_g = sqrt(accum_r_j_g);
        double norm_theta_tilde_j_g = sqrt(accum_theta_tilde_j_g);
        double value; 
        if(weight == true) {
          value = 1 - (lambda_2/(delta*norm_r_j_g*norm_theta_tilde_j_g));
        } else {
          value = 1 - (lambda_2/(delta*norm_r_j_g));
        }
        if(value >= 0) {
          theta_new.submat(K*j, g, K*(j+1)-1, g) = value * r_j_g;
        } else {
          theta_new.submat(K*j, g, K*(j+1)-1, g) = arma::vec(K);
        }
      }
    }
    
    //Process for Z = XA
    arma::mat sum_E = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = e_old)));
    arma::mat sum_U = as<arma::mat>(wrap(f(Named("f") = "+", Named("x") = u_old)));
    arma::mat obj = Y - (1/b)*sum_V*eta_new - (1/b)*sum_E + (1/(delta*b))*sum_U;
    arma::mat P;
    arma::mat Q;
    arma::vec d;
    arma::svd(P, d, Q, obj);
    if(m > n) {Q = Q.cols(0, n-1);}
    if(n > m) {P = P.cols(0, m-1);}
    arma::mat Z_0_P;
    arma::mat Z_0_Q;
    arma::vec Z_0_d;
    arma::svd(Z_0_P, Z_0_d, Z_0_Q, Z_0);
    arma::vec d_new;
    
    if(weight == true) {
      d_new = d - lambda_1/(delta*b*Z_0_d);
    } else {
      d_new = d - lambda_1/(delta*b);
    }
    NumericVector diag_entry = ifelse(as<NumericVector>(wrap(d_new)) > 0, 
                                      as<NumericVector>(wrap(d_new)), 0);
    arma::mat D_new = diagmat(as<arma::vec>(wrap(diag_entry)));
    arma::mat Z_new = P*D_new*Q.t();
    
    //Process for e
    List e_new = List(b);
    for(int l = 0; l < b; l++) {
      arma::mat mat_E = arma::mat(n, m);
      arma::mat VH = as<arma::mat>(wrap(V[l]))*eta_new;
      arma::mat Error = Y - Z_new - VH;
      arma::mat mat_U = u_old[l];
      for(int g=0; g<m; g++) {
        for(int i=0; i<n; i++) {
          if(Error(i, g) + mat_U(i, g)/delta < (tau_seq(l)-1)/(n*b*delta)) {
            mat_E(i, g) = Error(i, g) + mat_U(i, g)/delta - (tau_seq(l)-1)/(n*b*delta);
          } else if(Error(i, g) + mat_U(i, g)/delta > tau_seq(l)/(n*b*delta)) {
            mat_E(i, g) = Error(i, g) + mat_U(i, g)/delta - tau_seq(l)/(n*b*delta);
          } else {mat_E(i, g) = 0;}
        }
      }
      e_new[l] = mat_E;
    }
    
    //Process for multiplier u
    List u_new = List(b);
    for(int l=0; l<b; l++) {
      arma::mat VH = as<arma::mat>(wrap(V[l]))*eta_new;
      arma::mat mat_E = e_new[l];
      arma::mat mat_U_old = u_old[l];
      arma::mat mat_U = mat_U_old + delta * (Y - Z_new - VH - mat_E);
      u_new[l] = mat_U;
    }
    
    //Process for multiplier w
    arma::mat w_new = w_old + delta * (theta_new - eta_new);
    
    //Update iteration error
    iter_error(iter, 0) = norm(eta_old - eta_new, "fro"); 
    iter_error(iter, 1) = norm(theta_old - theta_new, "fro");
    iter_error(iter, 2) = norm(Z_old - Z_new, "fro");
    List diff_e = List(b);
    List diff_u = List(b);
    for(int l=0; l<b; l++) {
      double diff_e_norm = norm(as<arma::mat>(wrap(e_old[l])) - as<arma::mat>(wrap(e_new[l])), "fro");
      diff_e[l] = diff_e_norm;
      double diff_u_norm = norm(as<arma::mat>(wrap(u_old[l])) - as<arma::mat>(wrap(u_new[l])), "fro");
      diff_u[l] = diff_u_norm;
    }
    iter_error(iter, 3) = as<double>(wrap(f(Named("f") = "+", Named("x") = diff_e)));
    iter_error(iter, 4) = as<double>(wrap(f(Named("f") = "+", Named("x") = diff_u)));
    iter_error(iter, 5) = norm(w_old - w_new, "fro");
    
    //End the loop
    arma::vec threshold = as<arma::vec>(wrap(iter_error.row(iter)));
    double count = 0;
    for(int i=0; i<6; i++) {
      count += threshold[i];
    }
    if(count < tol_error) break;
    
    //replace old parameter
    eta_old = eta_new;
    theta_old = theta_new;
    Z_old = Z_new;
    e_old = e_new;
    u_old = u_new;
    w_old = w_new;
    
  }
  
  DataFrame params = DataFrame::create(Named("lambda_1")=lambda_1, Named("lambda_2")=lambda_2);
  
  //result
  List estimate = List::create(Named("eta") = eta_old, 
                               Named("theta") = theta_old, 
                               Named("Z") = Z_old, 
                               Named("e") = e_old, 
                               Named("u") = u_old, 
                               Named("w") = w_old, 
                               Named("iter_error") = iter_error, 
                               Named("params") = params);
  
  
  
  
  
  return estimate;
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