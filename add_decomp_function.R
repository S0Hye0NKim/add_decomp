add_decomp <- function(delta, lambda_1, lambda_2, tol_error, max_iter) {
  # delta = step size
  # lambda_1 = low rank penalty
  # lambda_2 = sparse penalty
  
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
      r_g <- eta_new[, g] - w_old[, g]/delta
      value <- 1 - (lambda_2/(delta * abs(r_g)))
      theta_new[, g] <- ifelse(value > 0, value * r_g, 0)
    }
    
    # Process for Z=XA
    Y_list <- list()
    for(i in 1:l) {Y_list[[i]] <- Y}
    VH_list <- lapply(V, FUN = function(x) x %*% eta_new)
    obj_list <- mapply(function(Y, VH, E, U) Y - VH - E - U/delta, Y_list, VH_list, e_old, u_old, SIMPLIFY = FALSE)
    obj <- Reduce("+", obj_list)/b 
    SVD <- svd(obj)
    new_singular <- sapply(SVD$d - lambda_1/(delta*b), FUN = function(x) max(x, 0))
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