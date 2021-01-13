add_decomp <- function(delta, lambda_1, lambda_2, tol_error, max_iter, X, Y, V, Phi, theta_0, alpha_0) {
  # delta = step size
  # lambda_1 = low rank penalty
  # lambda_2 = sparse penalty
  
  # initial value
  eta_old <- theta_0
  theta_old <- eta_old
  alpha_old <- alpha_0
  Z_old <- prod_AB(A = X, B = alpha_old, A_t = FALSE, B_t = FALSE)
  e_old <- list()
  for(l in 1:b) {e_old[[l]] <- Y - Z_old - prod_AB(A = V[[l]], B = eta_old, A_t = FALSE, B_t = FALSE)}
  u_old <- list()
  for(l in 1:b) {u_old[[l]] <- matrix(0, nrow = n, ncol = m)}
  w_old <- matrix(0, nrow = (p+1)*K, ncol = m)
  
  iter_error <- matrix(ncol = 6, nrow = max_iter) %>%
    `colnames<-`(value = c("eta", "theta", "alpha", "e", "u", "w"))
  
  sum_V <- Sum_Mat(V)
  VV_prod <- lapply(V, FUN = function(x) prod_AB(x, x, A_t = TRUE, B_t = FALSE))   # V^TV
  sum_VV <- Sum_Mat(VV_prod)
  
  for(iter in 1:max_iter) {
    # Process for eta
    eta_new <- matrix(nrow = (p+1)*K, ncol = m)
    Vu_prod <- mapply(function(x,y) prod_AB(x, y, A_t = TRUE, B_t = FALSE), V, u_old, SIMPLIFY = FALSE)
    Ve_prod <- mapply(function(x,y) prod_AB(x, y, A_t = TRUE, B_t = FALSE), V, e_old, SIMPLIFY = FALSE)
    eta_new <- update_eta(delta, sum_VV, W = w_old, Theta = theta_old, 
                          sum_VU = Sum_Mat(Vu_prod), sum_V, Y, Z = Z_old, sum_VE = Sum_Mat(Ve_prod))
    
    # Process for theta
    theta_new <- update_theta(delta, lambda_2, K = K, Eta = eta_new, W = w_old)
    
    # Process for Z=XA
    sub_process <- update_alpha(delta, lambda_1, Y = Y, X = X, H = eta_new, V = V, E = e_old, U = u_old)
    Z_new <- sub_process$Z_new
    alpha_new <- sub_process$alpha_new
    
    # Process for e
    e_new <- update_e(delta, U = u_old, V = V, Y = Y, X = X, A = alpha_new, H = eta_new, tau_seq = tau_seq)
    
    # Process for multiplier u
    u_new <- update_multi(U_old = u_old, W_old = w_old, delta, Y, X, A = alpha_new, H = eta_new, 
                          Theta = theta_new, V = V, E = e_new, is_u = TRUE)
    # Process for multiplier w
    w_new <- update_multi(U_old = u_old, W_old = w_old, delta, Y, X, A = alpha_new, H = eta_new, 
                          Theta = theta_new, V = V, E = e_new, is_u = FALSE) %>% .[[1]]
    
    # Update iteration error
    iter_error[iter, "eta"] <- calc_err(list(eta_old), list(eta_new))
    iter_error[iter, "theta"] <- calc_err(list(theta_old), list(theta_new))
    iter_error[iter, "alpha"] <- calc_err(list(alpha_old), list(alpha_new))
    iter_error[iter, "e"] <- calc_err(e_old, e_new)
    iter_error[iter, "u"] <- calc_err(u_old, u_new)
    iter_error[iter, "w"] <- calc_err(list(w_old), list(w_new))
    
    if(sum(iter_error[iter, ]) < tol_error) break
    
    eta_old <- eta_new
    theta_old <- theta_new
    Z_old <- Z_new
    alpha_old <- alpha_new
    e_old <- e_new
    u_old <- u_new
    w_old <- w_new
  }
  
  return(list(eta = eta_new, 
              theta = theta_new, 
              alpha = alpha_new, 
              e = e_new, 
              u = u_new, 
              w = w_new, 
              iter_error = iter_error))
}

# Calculate TP, TN, FP, FN of sparse matrix
check_sp_table <- function(true, est, tol = 0.1^5, table = FALSE) {
  # check sparsity pattern of true and est matrix
  zero_idx_true <- which(abs(true) < tol, arr.ind = TRUE) %>% as_tibble
  zero_idx_est <- lapply(est, FUN = function(x) which(abs(x) < tol, arr.ind = TRUE) %>% as_tibble) %>%
    bind_rows (.id = "tau") %>%
    group_by(row, col) %>%
    summarise(zero = n()) %>%
    filter(zero == b)
  
  num_zero <- which(true==0, arr.ind = TRUE) %>% nrow
  num_nz <- length(true) - num_zero
  
  if(nrow(zero_idx_est) == 0) {
    result <- data.frame(Positive = c(num_nz, num_zero), Negative = c(0, 0))
  } else {
    TN <- semi_join(zero_idx_true, zero_idx_est, by = c("row", "col"))
    FP <- anti_join(zero_idx_true, zero_idx_est, by = c("row", "col"))
    FN <- anti_join(zero_idx_est, zero_idx_true, by = c("row", "col"))
    result <- data.frame(Positive = c(num_nz - nrow(FN), nrow(FP)), Negative = c(nrow(FN), nrow(TN))) %>%
      `rownames<-`(value = c("Positive", "Negative")) %>%
      `colnames<-`(value = c("Est_Positive", "Est_Negative"))
  }
  
  if(table == FALSE) {return(data.frame(FPR = result[2, 1]/num_zero, TPR = result[1, 1]/num_nz))
  } else {return(result)}
}


# Estimate gamma from theta
est_gamma <- function(Phi, theta) {
  l <- nrow(Phi)
  K <- ncol(Phi)
  m <- ncol(theta)
  p <- nrow(theta)/K - 1
  gamma_tau_hat <- list()
  for(l in 1:b) {
    phi_tau <- Phi[l, ]
    gamma_tau_hat[[l]] <- matrix(nrow = (p+1), ncol = m)
    for(i in 1:(p+1)) {
      for(j in 1:m) {
        theta_j <- theta[(1+(i-1)*K):(K*i), j]
        gamma_tau_hat[[l]][i, j] <- theta_j %*% phi_tau
      }
    }
  }
  return(gamma_tau_hat)
}

# Check loss function
check_ft <- function(x, tau) {
  z <- ifelse(x<0, (tau-1)*x, tau*x)
  return(z)
}

# Calculate check loss sum
cal_cl_sum <- function(e, tau_seq) {
  check_loss <- list()
  n <- e[[1]] %>% nrow
  for(i in 1:length(tau_seq)) {
    check_loss[[i]] <- check_ft(e[[i]], tau = tau_seq[i])
  }
  
  result <- lapply(check_loss, FUN = function(x) data.frame(cl_sum = sum(x))) %>%
    unlist %>% sum
  return(result/n)
}