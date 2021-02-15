add_decomp_r <- function(delta, lambda_1, lambda_2, tol_error, max_iter, X, Y, V, Phi, 
                         theta_0, alpha_0, tau_seq) {
  # delta = step size
  # lambda_1 = low rank penalty
  # lambda_2 = sparse penalty
  
  # initial value
  eta_old <- theta_0
  theta_old <- eta_old
  alpha_old <- alpha_0
  Z_old <- X %*% alpha_old
  e_old <- list()
  for(l in 1:b) {e_old[[l]] <- Y - Z_old - V[[l]] %*% eta_old}
  u_old <- list()
  for(l in 1:b) {u_old[[l]] <- matrix(0, nrow = n, ncol = m)}
  w_old <- matrix(0, nrow = (p+1)*K, ncol = m)
  
  iter_error <- matrix(ncol = 6, nrow = max_iter) %>%
    `colnames<-`(value = c("eta", "theta", "alpha", "e", "u", "w"))
  
  sum_V <- Reduce("+", V)
  VV_prod <- lapply(V, FUN = function(x) t(x) %*% x)   # V^TV
  sum_VV <- Reduce("+", VV_prod)
  
  for(iter in 1:max_iter) {
    # Process for eta
    eta_new <- matrix(nrow = (p+1)*K, ncol = m)
    Vu_prod <- mapply(function(x,y) t(x) %*% y, V, u_old, SIMPLIFY = FALSE)
    Ve_prod <- mapply(function(x,y) t(x) %*% y, V, e_old, SIMPLIFY = FALSE)
    eta_new <- (solve(sum_VV+diag(1, (p+1)*K))/delta) %*% (w_old + delta * theta_old + Reduce("+", Vu_prod)
                                                           + delta * t(sum_V) %*% (Y - Z_old)
                                                           - delta * Reduce("+", Ve_prod))
    
    # Process for theta
    theta_new <- matrix(nrow = (p+1)*K, ncol = m)
    for (g in 1:m) {
      for(j in 1:(p+1)) {
        theta_tilde <- theta_0[(K*(j-1) +1):(j*K), g]
        norm_theta_tilde <- (theta_tilde^2) %>% sum %>% sqrt   # weight = 1/norm_theta_tilde
        eta_j_g <- eta_new[(K*(j-1) +1):(j*K), g]
        w_j_g <- w_old[(K*(j-1) +1):(j*K), g]
        r_j_g <- eta_j_g - (w_j_g/delta)
        norm_r_j_g <- (r_j_g^2) %>% sum %>% sqrt
        value <- 1 - (lambda_2/(delta *norm_r_j_g*norm_theta_tilde))
        if(value >= 0) {
          theta_new[(K*(j-1) +1):(j*K), g] <- value * r_j_g
        } else {theta_new[(K*(j-1) +1):(j*K), g] <- 0}
      }
    }
    
    # Process for Z=XA
    Y_list <- list()
    for(i in 1:l) {Y_list[[i]] <- Y}
    VH_list <- lapply(V, FUN = function(x) x %*% eta_new)
    obj_list <- mapply(function(Y, VH, E, U) Y - VH - E - U/delta, Y_list, VH_list, e_old, u_old, SIMPLIFY = FALSE)
    obj <- Reduce("+", obj_list)/b 
    SVD <- svd(obj)
    sing_val_alpha_0 <- svd(alpha_0) %>% .$d  # weight = 1/sing_val_alpha_0
    new_singular <- sapply(SVD$d - lambda_1/(delta*b*sing_val_alpha_0), FUN = function(x) max(x, 0))
    Z_new <- SVD$u %*% diag(new_singular) %*% t(SVD$v)
    alpha_new <- solve(t(X) %*% X) %*% t(X) %*% Z_new
    
    # Process for e
    e_new <- list()
    for(l in 1:b){
      e_new[[l]] <- matrix(nrow = n, ncol = m)
      for(g in 1:m) {
        error <- Y[, g] - Z_new[, g] - V[[l]] %*% eta_new[, g]   #error = Y - XA - VH
        value <- error + u_old[[l]][, g]/delta
        e_new[[l]][, g] <- case_when(value > tau_seq[l]/(n*delta) ~ value - tau_seq[l]/(n*delta), 
                                     value < (tau_seq[l]-1)/(n*delta) ~ value - (tau_seq[l]-1)/(n*delta), 
                                     value >=(tau_seq[l]-1)/(n*delta) & value <= tau_seq[l]/(n*delta) ~ 0)
      }
    }
    
    # Process for multiplier u
    u_new <- list()
    for(l in 1:b) {
      u_new[[l]] <- u_old[[l]] + delta * (Y - Z_new - V[[l]] %*% eta_new - e_new[[l]])
    }
    
    # Process for multiplier w
    w_new <- w_old + delta * (theta_new - eta_new)
    
    # Update iteration error
    iter_error[iter, "eta"] <- Matrix::norm(eta_old - eta_new, type = "F")
    iter_error[iter, "theta"] <- Matrix::norm(theta_old - theta_new, type = "F")
    iter_error[iter, "alpha"] <- Matrix::norm(alpha_old - alpha_new, type = "F")
    e_diff <- mapply(FUN = function(old, new) old - new, e_old, e_new, SIMPLIFY = FALSE)  # sum of frobenius norm
    iter_error[iter, "e"] <- lapply(e_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    u_diff <- mapply(FUN = function(old, new) old - new, u_old, u_new, SIMPLIFY = FALSE)
    iter_error[iter, "u"] <- lapply(u_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    iter_error[iter, "w"] <- Matrix::norm(w_old - w_new, type = "F")
    
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
check_sp_table <- function(true, est, tol = 0.1^5, table = FALSE, tau_seq) {
  # check sparsity pattern of true and est matrix
  b <- length(tau_seq)
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
  b <- nrow(Phi)
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