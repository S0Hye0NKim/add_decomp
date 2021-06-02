scad_deriv <- function(x, lambda, gamma = 3.7) {
  output <- ifelse(abs(x) <= lambda, lambda, 
                   ifelse(abs(x) >= gamma*lambda, 0, (gamma*lambda - abs(x))/(gamma - 1)))
  return(output)
}

add_decomp_r <- function(delta, lambda_1, lambda_2, tol_error, max_iter, X, Y, V, Phi, 
                         theta_0, Z_0, tau_seq, weight = TRUE) {
  # delta = step size
  # lambda_1 = low rank penalty
  # lambda_2 = sparse penalty
  
  n <- nrow(X)
  p <- ncol(X) - 1
  K <- ncol(Phi)
  b <- nrow(Phi)
  m <- ncol(Y)
  
  # initial value
  eta_old <- theta_0
  theta_old <- eta_old
  Z_old <- Z_0
  e_old <- list()
  for(l in 1:b) {e_old[[l]] <- Y - Z_old - V[[l]] %*% eta_old}
  u_old <- list()
  for(l in 1:b) {u_old[[l]] <- matrix(0, nrow = n, ncol = m)}
  w_old <- matrix(0, nrow = (p+1)*K, ncol = m)
  
  iter_error <- matrix(ncol = 6, nrow = max_iter) %>%
    `colnames<-`(value = c("eta", "theta", "Z", "e", "u", "w"))
  
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
        norm_theta_tilde <- ifelse(weight == TRUE, (theta_tilde^2) %>% sum %>% sqrt, 1)  # weight = 1/norm_theta_tilde
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
    for(i in 1:b) {Y_list[[i]] <- Y}
    VH_list <- lapply(V, FUN = function(x) x %*% eta_new)
    obj_list <- mapply(function(Y, VH, E, U) Y - VH - E + U/delta, Y_list, VH_list, e_old, u_old, SIMPLIFY = FALSE)
    obj <- Reduce("+", obj_list)/b 
    SVD <- svd(obj)
    if(weight == TRUE) {
      sing_val_Z_0 <- svd(Z_0) %>% .$d
    } else {sing_val_Z_0 <- rep(1, length(svd(Z_0) %>% .$d))}  # weight = 1/sing_val_Z_0
    new_singular <- sapply(SVD$d - lambda_1/(delta*b*sing_val_Z_0), FUN = function(x) max(x, 0))
    Z_new <- SVD$u %*% diag(new_singular) %*% t(SVD$v)
    
    # Process for e
    e_new <- list()
    for(l in 1:b){
      e_new[[l]] <- matrix(nrow = n, ncol = m)
      for(g in 1:m) {
        error <- Y[, g] - Z_new[, g] - V[[l]] %*% eta_new[, g]   #error = Y - XA - VH
        value <- error + u_old[[l]][, g]/delta
        e_new[[l]][, g] <- case_when(value > tau_seq[l]/(n*b*delta) ~ value - tau_seq[l]/(n*b*delta), 
                                     value < (tau_seq[l]-1)/(n*b*delta) ~ value - (tau_seq[l]-1)/(n*b*delta), 
                                     value >=(tau_seq[l]-1)/(n*b*delta) & value <= tau_seq[l]/(n*b*delta) ~ 0)
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
    iter_error[iter, "Z"] <- Matrix::norm(Z_old - Z_new, type = "F")
    e_diff <- mapply(FUN = function(old, new) old - new, e_old, e_new, SIMPLIFY = FALSE)  # sum of frobenius norm
    iter_error[iter, "e"] <- lapply(e_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    u_diff <- mapply(FUN = function(old, new) old - new, u_old, u_new, SIMPLIFY = FALSE)
    iter_error[iter, "u"] <- lapply(u_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    iter_error[iter, "w"] <- Matrix::norm(w_old - w_new, type = "F")
    
    if(sum(iter_error[iter, ]) < tol_error) break
    
    eta_old <- eta_new
    theta_old <- theta_new
    Z_old <- Z_new
    e_old <- e_new
    u_old <- u_new
    w_old <- w_new
  }
  
  return(list(eta = eta_new, 
              theta = theta_new, 
              Z = Z_new, 
              e = e_new, 
              u = u_new, 
              w = w_new, 
              iter_error = iter_error, 
              params = c(lambda_1, lambda_2) %>% `names<-`(value = c("lambda_1", "lambda_2"))))
}

# Calculate TP, TN, FP, FN of sparse matrix
check_sp_table <- function(true, est, tol = 0.1^5, table = FALSE, tau_seq) {
  # check sparsity pattern of true and est matrix
  b <- length(tau_seq)
  zero_idx_true <- which(abs(true) < tol, arr.ind = TRUE) %>% as_tibble
  zero_idx_est <- lapply(est, FUN = function(x) which(abs(x) < tol, arr.ind = TRUE) %>% as_tibble) %>%
    bind_rows (.id = "tau") %>%
    group_by(row, col) %>%
    summarise(zero = n(), .groups = "keep") %>%
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


# parameter selection via BIC
add_decomp_BIC <- function(X, Y, V, Phi, theta_0, Z_0, tau_seq, tau_seq_real, delta, lamb1_seq, lamb2_seq, 
                           max_iter) {
  m <- ncol(Y)
  p <- ncol(X) - 1
  K <- ncol(V[[1]])/(p+1)
  n <- nrow(X)
  idx_tau <- tau_seq %in% tau_seq_real
  
  # iteration for lamb1_seq and lamb2_seq
  simulation <- list()
  for(lamb1_idx in 1:length(lamb1_seq)) {
    simulation[[lamb1_idx]] <- list()
    for(lamb2_idx in 1:length(lamb2_seq)) {
      simulation[[lamb1_idx]][[lamb2_idx]] <- add_decomp_r(delta = delta, lambda_1 = lamb1_seq[lamb1_idx], 
                                                           lambda_2 = lamb2_seq[lamb2_idx], tol_error = 0.1^5, 
                                                           max_iter = max_iter, X = X, Y = Y, V = V, Phi = Phi, 
                                                           theta_0, Z_0, tau_seq = tau_seq, weight = TRUE)
    }
  }
  
  simulation <- lapply(simulation, FUN = function(x) `names<-`(x, value = paste0("lambda_2=", lamb2_seq)))
  names(simulation) <- paste0("lambda_1=", lamb1_seq)
  
  BIC <- list()
  for(i in 1:length(lamb1_seq)) {
    BIC[[i]] <- list()
    for(j in 1:length(lamb2_seq)) {
      result <- simulation[[i]][[j]]
      est_error <- lapply(V[idx_tau], FUN = function(x) (Y - result$Z - x %*% result$theta)
                          %>% as.vector())
      check_loss_err <- mapply(FUN = function(x, tau) check_ft(x, tau), x = est_error, 
                               tau = as.list(tau_seq_real), SIMPLIFY = FALSE) %>%
        lapply(FUN = function(x) sum(x)) %>% unlist %>% sum
      gamma_tau_hat <- est_gamma(Phi[idx_tau, ], result$theta)
      zero_idx_est <- lapply(gamma_tau_hat, FUN = function(x) which(abs(x) < 0.1^5, arr.ind = TRUE) %>% as_tibble) %>%
        bind_rows(.id = "tau") %>%
        group_by(row, col) %>%
        summarise(zero = n(), .groups = "keep") %>%
        filter(zero == sum(idx_tau))
      num_nz <- nrow(zero_idx_est)
      S_hat <- (gamma_tau_hat[[1]] %>% dim %>% prod) - num_nz 
      num_nz_intercept <- m - (zero_idx_est %>% filter(row == 1) %>% nrow)
      
      BIC[[i]][[j]] <- data.frame(log_Q = log(check_loss_err), 
                                  r_hat = rankMatrix(result$Z)[1], 
                                  S_hat = S_hat, 
                                  num_nz_intercept = num_nz_intercept)
    }
  }
  
  r_X <- rankMatrix(X)
  params_table <- 
    lapply(BIC, FUN = function(x) `names<-`(x, value = lamb2_seq) %>% 
             bind_rows(.id = "lambda_2")) %>%
    `names<-`(value = lamb1_seq) %>%
    bind_rows(.id = "lambda_1") 
  
  output <- list(table = params_table, 
                 simulation = simulation)
  
  
  return(output)
}


# parameter selection via BIC
LR_model_BIC <- function(X, Y, Z_0, tau_seq, tau_seq_real, delta, lamb_seq, max_iter, r_X) {
  m <- ncol(Y)
  p <- ncol(X) - 1
  n <- nrow(X)
  idx_tau <- tau_seq %in% tau_seq_real
  
  # iteration for lamb_seq
  simulation <- list()
  for(lamb_idx in 1:length(lamb_seq)) {
    simulation[[lamb_idx]] <- LR_model_r(delta = delta, lambda = lamb_seq[lamb_idx], tol_error = 0.1^5, 
                                       max_iter = max_iter, X = X, Y = Y, Z_0 = Z_0, tau_seq = tau_seq, weight = TRUE)
  }
  
  names(simulation) <- paste0("lambda=", lamb_seq)
  
  BIC <- list()
  for(i in 1:length(lamb_seq)) {
    result <- simulation[[i]]
    est_error <- Y - result$Z
    check_loss_err <- lapply(as.list(tau_seq_real), FUN = function(x) check_ft(est_error, x)) %>%
      Reduce("+", .) %>% as.vector %>% sum
    BIC[[i]] <- data.frame(log_Q = log(check_loss_err), r_hat = rankMatrix(result$Z)[1])
  }

  names(BIC) <- lamb_seq
  BIC_data <- BIC %>% bind_rows(.id = "lambda") %>%
    mutate(term = (r_hat * max(r_X, m))/(2*n*m), 
           BIC_log_sum = log_Q + log(p+m)*term, 
           BIC_log_p = log_Q + log(p)*term, 
           BIC_log_n = log_Q + log(n)*term, 
           BIC_llog_p = log_Q + log(log(p))*term,
           BIC_llog_n = log_Q + log(log(n))*term) %>%
    group_by(lambda) %>%
    select_at(vars(starts_with("BIC"))) %>%
    ungroup()
  
  BIC_val_min <- apply(select_at(BIC_data, vars(starts_with("BIC"))), 2, min) %>%
    `names<-`(value = c("log_sum", "log_p", "log_n", "llog_p", "llog_n"))
  
  min_BIC <- filter(BIC_data, BIC_log_sum == BIC_val_min["log_sum"]| BIC_log_p == BIC_val_min["log_p"] |
                      BIC_log_n == BIC_val_min["log_n"] | BIC_llog_p == BIC_val_min["llog_p"] | 
                      BIC_llog_n == BIC_val_min["llog_n"])
  
  output <- list(min_BIC = min_BIC, 
                 BIC_data = BIC_data, 
                 simulation = simulation)
  
  return(output)
}

# parameter selection via BIC
SP_model_BIC <- function(X, Y, V, Phi, theta_0, tau_seq, tau_seq_real, delta, lamb_seq, max_iter) {
  m <- ncol(Y)
  p <- ncol(X) - 1
  K <- ncol(V[[1]])/(p+1)
  n <- nrow(X)
  idx_tau <- tau_seq %in% tau_seq_real
  
  # iteration for lamb1_seq and lamb2_seq
  simulation <- list()
  for(lamb_idx in 1:length(lamb_seq)) {
    simulation[[lamb_idx]] <- SP_model_r(delta = delta, lambda = lamb_seq[lamb_idx], tol_error = 0.1^5, 
                                       max_iter = max_iter, X = X, Y = Y, V = V, Phi = Phi, theta_0 = theta_0, 
                                       tau_seq = tau_seq, weight = TRUE)
  }
  
  names(simulation) <- paste0("lambda=", lamb_seq)
  
  BIC <- list()
  for(i in 1:length(lamb_seq)) {
    result <- simulation[[i]]
    est_error <- lapply(V[idx_tau], FUN = function(x) (Y - x %*% result$theta)
                        %>% as.vector())
    check_loss_err <- mapply(FUN = function(x, tau) check_ft(x, tau), x = est_error, 
                             tau = as.list(tau_seq_real), SIMPLIFY = FALSE) %>%
      lapply(FUN = function(x) sum(x)) %>% unlist %>% sum
    gamma_tau_hat <- est_gamma(Phi[idx_tau, ], result$theta)
    S_hat <- check_sp_table(true = matrix(0, nrow = (p+1), ncol = m), 
                            est = gamma_tau_hat, table = TRUE, tol = 0.1^5, tau_seq = tau_seq_real) %>%
      .$Est_Positive %>% sum
    BIC[[i]] <- data.frame(log_Q = log(check_loss_err), S_hat = S_hat)
  }
  
  names(BIC) <- lamb_seq
  BIC_data <- BIC %>% bind_rows(.id = "lambda") %>%
    mutate(term = (K * S_hat)/(2*n*m), 
           BIC_log_sum = log_Q + log(p+m)*term, 
           BIC_log_p = log_Q + log(p)*term, 
           BIC_log_n = log_Q + log(n)*term, 
           BIC_llog_p = log_Q + log(log(p))*term,
           BIC_llog_n = log_Q + log(log(n))*term) %>%
    group_by(lambda) %>%
    select_at(vars(starts_with("BIC"))) %>%
    ungroup()
  
  BIC_val_min <- apply(select_at(BIC_data, vars(starts_with("BIC"))), 2, min) %>%
    `names<-`(value = c("log_sum", "log_p", "log_n", "llog_p", "llog_n"))
  
  min_BIC <- filter(BIC_data, BIC_log_sum == BIC_val_min["log_sum"]| BIC_log_p == BIC_val_min["log_p"] |
                      BIC_log_n == BIC_val_min["log_n"] | BIC_llog_p == BIC_val_min["llog_p"] | 
                      BIC_llog_n == BIC_val_min["llog_n"])
  
  output <- list(min_BIC = min_BIC, 
                 BIC_data = BIC_data, 
                 simulation = simulation)
  
  return(output)
}

# Low rank model 
LR_model_r <- function(delta, lambda, tol_error, max_iter, X, Y, Z_0, tau_seq, weight) {
  n <- nrow(X)
  p <- ncol(X) - 1
  b <- length(tau_seq)
  m <- ncol(Y)
  Z_old <- Z_0
  e_old <- list()
  u_old <- list()
  for(l in 1:b) {
    e_old[[l]] <- Y- Z_old
    u_old[[l]] <- matrix(0, nrow = n, ncol = m)
  }
  iter_error <- matrix(ncol = 3, nrow = max_iter) %>%
    `colnames<-`(value = c("Z", "e", "u"))
  
  for(iter in 1:max_iter) {
    # Process for Z
    sum_E <- Reduce("+", e_old)
    sum_U <- Reduce("+", u_old)
    obj <- Y - sum_E/b + (1/(delta*b))*sum_U
    SVD_obj <- svd(obj)
    SVD_Z_0 <- svd(Z_0)
    if(weight == TRUE) {
      d_new <- SVD_obj$d - lambda/(delta*b*SVD_Z_0$d)
    } else {
      d_new <- SVD_obj$d - lambda/(delta*b)
    }
    diag_entry <- ifelse(d_new > 0, d_new, 0)
    Z_new <- SVD_obj$u %*% diag(diag_entry) %*% t(SVD_obj$v)
    
    # Process for e
    e_new <- list()
    for(l in 1:b){
      e_new[[l]] <- matrix(nrow = n, ncol = m)
      for(g in 1:m) {
        error <- Y[, g] - Z_new[, g]   #error = Y - XA
        value <- error + u_old[[l]][, g]/delta
        e_new[[l]][, g] <- case_when(value > tau_seq[l]/(n*b*delta) ~ value - tau_seq[l]/(n*b*delta), 
                                     value < (tau_seq[l]-1)/(n*b*delta) ~ value - (tau_seq[l]-1)/(n*b*delta), 
                                     value >=(tau_seq[l]-1)/(n*b*delta) & value <= tau_seq[l]/(n*b*delta) ~ 0)
      }
    }
    
    # Process for u
    u_new <- list()
    for(l in 1:b) {
      u_new[[l]] <- u_old[[l]] + delta * (Y - Z_new - e_new[[l]])
    }
    
    # update iteration error
    iter_error[iter, "Z"] <- Matrix::norm(Z_old - Z_new, type = "F")
    e_diff <- mapply(FUN = function(old, new) old - new, e_old, e_new, SIMPLIFY = FALSE)  
    iter_error[iter, "e"] <- lapply(e_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    u_diff <- mapply(FUN = function(old, new) old - new, u_old, u_new, SIMPLIFY = FALSE)
    iter_error[iter, "u"] <- lapply(u_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    
    if(sum(iter_error[iter, ]) < tol_error) break
    
    Z_old <- Z_new
    e_old <- e_new
    u_old <- u_new
  }
  
  return(list(Z = Z_new, 
              e = e_new, 
              u = u_new, 
              iter_error = iter_error, 
              params = lambda))
  
}

# Sparse model
SP_model_r <- function(delta, lambda, tol_error, max_iter, X, Y, V, Phi, 
                       theta_0, tau_seq, weight = TRUE) {
  n <- nrow(X)
  m <- ncol(Y)
  p <- ncol(X) - 1
  b <- length(tau_seq)
  
  # initial value
  eta_old <- theta_0
  theta_old <- eta_old
  e_old <- list()
  for(l in 1:b) {e_old[[l]] <- Y - V[[l]] %*% eta_old}
  u_old <- list()
  for(l in 1:b) {u_old[[l]] <- matrix(0, nrow = n, ncol = m)}
  w_old <- matrix(0, nrow = (p+1)*K, ncol = m)
  
  iter_error <- matrix(ncol = 5, nrow = max_iter) %>%
    `colnames<-`(value = c("eta", "theta", "e", "u", "w"))
  
  sum_V <- Reduce("+", V)
  VV_prod <- lapply(V, FUN = function(x) t(x) %*% x)   # V^TV
  sum_VV <- Reduce("+", VV_prod)
  
  for(iter in 1:max_iter) {
    # Process for eta
    eta_new <- matrix(nrow = (p+1)*K, ncol = m)
    Vu_prod <- mapply(function(x,y) t(x) %*% y, V, u_old, SIMPLIFY = FALSE)
    Ve_prod <- mapply(function(x,y) t(x) %*% y, V, e_old, SIMPLIFY = FALSE)
    eta_new <- (solve(sum_VV+diag(1, (p+1)*K))/delta) %*% (w_old + delta * theta_old + Reduce("+", Vu_prod)
                                                           + delta * t(sum_V) %*% (Y)
                                                           - delta * Reduce("+", Ve_prod))
    
    # Process for theta
    theta_new <- matrix(nrow = (p+1)*K, ncol = m)
    for (g in 1:m) {
      for(j in 1:(p+1)) {
        theta_tilde <- theta_0[(K*(j-1) +1):(j*K), g]
        norm_theta_tilde <- ifelse(weight == TRUE, (theta_tilde^2) %>% sum %>% sqrt, 1)  # weight = 1/norm_theta_tilde
        eta_j_g <- eta_new[(K*(j-1) +1):(j*K), g]
        w_j_g <- w_old[(K*(j-1) +1):(j*K), g]
        r_j_g <- eta_j_g - (w_j_g/delta)
        norm_r_j_g <- (r_j_g^2) %>% sum %>% sqrt
        value <- 1 - (lambda/(delta *norm_r_j_g*norm_theta_tilde))
        if(value >= 0) {
          theta_new[(K*(j-1) +1):(j*K), g] <- value * r_j_g
        } else {theta_new[(K*(j-1) +1):(j*K), g] <- 0}
      }
    }
    
    # Process for e
    e_new <- list()
    for(l in 1:b){
      e_new[[l]] <- matrix(nrow = n, ncol = m)
      for(g in 1:m) {
        error <- Y[, g]  - V[[l]] %*% eta_new[, g]   #error = Y - VH
        value <- error + u_old[[l]][, g]/delta
        e_new[[l]][, g] <- case_when(value > tau_seq[l]/(n*b*delta) ~ value - tau_seq[l]/(n*b*delta), 
                                     value < (tau_seq[l]-1)/(n*b*delta) ~ value - (tau_seq[l]-1)/(n*b*delta), 
                                     value >=(tau_seq[l]-1)/(n*b*delta) & value <= tau_seq[l]/(n*b*delta) ~ 0)
      }
    }
    
    # Process for multiplier u
    u_new <- list()
    for(l in 1:b) {
      u_new[[l]] <- u_old[[l]] + delta * (Y - V[[l]] %*% eta_new - e_new[[l]])
    }
    
    # Process for multiplier w
    w_new <- w_old + delta * (theta_new - eta_new)
    
    # Update iteration error
    iter_error[iter, "eta"] <- Matrix::norm(eta_old - eta_new, type = "F")
    iter_error[iter, "theta"] <- Matrix::norm(theta_old - theta_new, type = "F")
    e_diff <- mapply(FUN = function(old, new) old - new, e_old, e_new, SIMPLIFY = FALSE)  # sum of frobenius norm
    iter_error[iter, "e"] <- lapply(e_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    u_diff <- mapply(FUN = function(old, new) old - new, u_old, u_new, SIMPLIFY = FALSE)
    iter_error[iter, "u"] <- lapply(u_diff, FUN = function(x) Matrix::norm(x, type = "F")) %>% Reduce("+", .)
    iter_error[iter, "w"] <- Matrix::norm(w_old - w_new, type = "F")
    
    if(sum(iter_error[iter, ]) < tol_error) break
    
    eta_old <- eta_new
    theta_old <- theta_new
    e_old <- e_new
    u_old <- u_new
    w_old <- w_new
  }
  
  return(list(eta = eta_new, 
              theta = theta_new, 
              e = e_new, 
              u = u_new, 
              w = w_new, 
              iter_error = iter_error, 
              params = lambda))
}
