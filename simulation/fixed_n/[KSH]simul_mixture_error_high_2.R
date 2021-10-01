# Simulation with normal error
rm(list = ls())

library(dplyr)
library(splines)
library(Matrix)
library(foreach)
library(doParallel)
library(Rcpp)
library(glmnet)
library(fda)
library(expm)

sourceCpp("[KSH]add_decomp_function.cpp")
source("https://raw.githubusercontent.com/S0Hye0NKim/add_decomp/master/functions/add_decomp_function.R")


# 1. low dim : (n,p) = (400, 100)

#################
## 1-0. Set up ##
#################

## Generate data
set.seed(3)
n <- 400
m <- 10
p <- 800
b <- 15
num_rank <- 3
num_rank_X <- 300
simul_times <- 33

sigma_mat <- matrix(nrow = p, ncol = p)
for(j in 1:p) {
  for(k in 1:p) {
    sigma_mat[j, k] <- 0.5^(abs(j-k))
  }
}

# sparse matrix
col_ind <- sample(1:m, size = m, replace = FALSE)
row_ind <- sample(2:(p+1), size = m)
sp_mat <- matrix(0, nrow = p+1, ncol = m)
for(i in 1:m) {
  sp_mat[row_ind[i], col_ind[i]] <- rnorm(1, mean = 5, sd = 0.1)
}
sp_mat[1, ] <- rnorm(m, mean = 5, sd = 0.1)

# low rank matrix
L1 <- matrix(rnorm(p*num_rank, mean = 0, sd = 0.4), nrow = p)
L2 <- matrix(rnorm(m*num_rank, mean = 0, sd = 0.4), nrow = m)
LR_mat <- L1 %*% t(L2)

X_list <- list()
eps_list <- list()

set.seed(2)
for(simul in 1:simul_times) {
  X1 <- matrix(rnorm(n*num_rank_X, mean = 0, sd = 1), nrow = n)
  X2 <- matrix(rnorm(p*num_rank_X, mean = 0, sd = 1), nrow = p)
  X <- X1 %*% t(X2) 
  X_list[[simul]] <- X %*% expm::sqrtm(sigma_mat) %>% cbind(rep(1, n), .)
  
  components <- sample(1:2, size = n*m, replace = TRUE)
  mus <- c(0, 1)
  
  eps_entry <- rnorm(n = n*m, mean = mus[components], sd = 1)
  eps_list[[simul]] <- matrix(eps_entry, nrow = n, ncol = m)
}

Y_list <- mapply(FUN = function(X, LR, SP, eps) X[, -1] %*% LR_mat + X %*% SP + eps, 
                 X_list, list(LR_mat), list(sp_mat), eps_list, SIMPLIFY = FALSE)


### Calculate kronecker product
K <- 5
tau_seq <- seq(from = 0.35, to = 0.65, length.out = b)
tau_seq_real <- tau_seq[tau_seq >= 0.4 & tau_seq  <= 0.6]
idx_tau <- (tau_seq >= "0.4" & tau_seq <= "0.6")

knots_seq <- seq(min(tau_seq)- 0.02, max(tau_seq) + 0.02, length.out = K)
Phi <- fda::bsplineS(tau_seq, breaks= knots_seq, norder=2, nderiv=0, returnMatrix=FALSE)

V_list <- list()
for(simul in 1:simul_times) {
  V_list[[simul]] <- calc_V(X_list[[simul]], Phi)
}

##################################
## 1-1. Simulation - add_decomp ##
##################################

simul_high_add_decomp_2 <- vector("list", length = simul_times)
for(simul in 1:simul_times) {
  
  tau_seq <- seq(from = 0.35, to = 0.65, length.out = b)
  X <- X_list[[simul]]
  Y <- Y_list[[simul]]
  V <- V_list[[simul]]
  
  lasso_coef <- matrix(nrow = p+1, ncol = m)
  for(g in 1:m) {
    cv.lasso <- cv.glmnet(x = X[, -1], y = Y[, g], 
                          alpha = 1, type.measure = "mae")
    lasso_model <- glmnet(X[, -1], Y[, g], 
                          family = "gaussian", alpha = 1, lambda = 0.01)
    lasso_coef[, g] <- c(lasso_model$a0, as.vector(lasso_model$beta))
  }
  
  theta_init <- matrix(nrow = (p+1)*K, ncol = m)
  for(g in 1:m) {
    for(j in 0:p) {
      theta_init[((j*K)+1):((j+1)*K), g] <- lasso_coef[j+1, g]
    }
  }
  
  Y_modified <- Y - X%*%lasso_coef
  ridge_coef <- matrix(nrow = p+1, ncol = m)
  for(g in 1:m) {
    cv.ridge <- cv.glmnet(x = X[, -1], y = Y_modified[, g], alpha = 0, type.measure = "mae")
    ridge_model <- glmnet(X[, -1], Y_modified[, g], family = "gaussian", alpha = 0, lambda = cv.ridge$lambda.min)
    ridge_coef[, g] <- c(ridge_model$a0, as.vector(ridge_model$beta))
  }
  alpha_init <- ridge_coef
  
  init_val <- add_decomp(delta = 1, lambda_1 = 0.01, lambda_2 = 0.001, tol_error = 0.1^5, max_iter = 50,
                         X = X, Y = Y, V = V, Phi = Phi, 
                         theta_0 = theta_init, Z_0 = X%*%alpha_init, tau_seq = tau_seq, weight = FALSE)
  
  log_lamb1 <- c( seq(-0.9, 0.1, length.out = 20))
  lamb1_seq <- exp(log_lamb1)
  log_lamb2 <- c(seq(4, 6.3, length.out = 20))
  lamb2_seq <- exp(log_lamb2)
  
  BIC_table <- list()
  
  for(idx in 1:length(lamb1_seq)) {
    lamb1 <- lamb1_seq[idx]
    cl <- makeCluster(20) #not to overload your computer
    registerDoParallel(cl) # Ready to parallel
    
    temp_BIC <- foreach(lamb2 = lamb2_seq, .noexport = "add_decomp") %dopar% {
      library(dplyr)
      library(splines)
      library(Matrix)
      library(glmnet)
      library(fda)
      library(Rcpp)
      sourceCpp("[KSH]add_decomp_function.cpp")
      
      BIC_simul <- add_decomp_BIC(X, Y, V, Phi, theta_0 = init_val$theta, Z_0 = init_val$Z, tau_seq, tau_seq_real,
                                  lamb1_seq = lamb1, lamb2_seq = lamb2, max_iter = 50, delta = 1, fun_type = "cpp")
      BIC_simul$table
    }
    stopCluster(cl)
    BIC_table[[idx]] <- temp_BIC
  }
  
  r_X <- rankMatrix(X[, -1])[1]
  BIC_params <- BIC_table %>% lapply(FUN = function(x) bind_rows(x)) %>%
    bind_rows() %>%
    mutate(LR_part = r_hat * max(r_X, m) / (2*n*m),
           S_hat_net = S_hat - num_nz_intercept,
           LR = log(p) * log(log(n)) * LR_part, 
           SP = log(p) * log(log(n)) * K * S_hat_net / (2*n*m),
           BIC = log_Q + LR + SP) %>%
    mutate_all(as.numeric) %>%
    filter(S_hat_net != 0) %>%
    arrange(BIC) %>%  
    head(1)
  
  result <- add_decomp(delta = 1, lambda_1 = BIC_params$lambda_1, lambda_2 = BIC_params$lambda_2, 
                       tol_error = 0.1^5, max_iter = 50, X, Y, V, Phi, 
                       theta_0 = init_val$theta, Z_0 = init_val$Z, tau_seq = tau_seq, weight = TRUE)
  
  simul_high_add_decomp_2[[simul]] <- result
}

################################
## 1-2. Simulation - LR_model ##
################################

cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
simul_high_LR_model_2 <- foreach(simul = 1:simul_times, .noexport = "add_decomp") %dopar% {
  library(dplyr)
  library(splines)
  library(Matrix)
  library(glmnet)
  library(fda)
  
  tau_seq <- seq(from = 0.35, to = 0.65, length.out = b)
  X <- X_list[[simul]]
  Y <- Y_list[[simul]]
  
  ridge_coef <- matrix(nrow = p+1, ncol = m)
  for(g in 1:m) {
    cv.ridge <- cv.glmnet(x = X[, -1], y = Y[, g], alpha = 0, type.measure = "mae")
    ridge_model <- glmnet(X[, -1], Y[, g], family = "gaussian", alpha = 0, lambda = cv.ridge$lambda.min)
    ridge_coef[, g] <- c(ridge_model$a0, as.vector(ridge_model$beta))
  }
  first_init_LR <- ridge_coef
  
  init_val_LR <- LR_model_r(delta = 1, lambda = 100, tol_error = 0.1^5, max_iter = 50, 
                            X = X, Y = Y, Z_0 = X %*% first_init_LR, tau_seq = tau_seq, weight = FALSE)
  
  lamb_seq <- seq(0.1, 1, length.out = 25)
  r_X <- rankMatrix(X[, -1])
  BIC_simul <- LR_model_BIC(X, Y, Z_0 = init_val_LR$Z, tau_seq, tau_seq_real, lamb_seq, max_iter = 50, delta = 1, r_X = rankMatrix(X[, -1]))
  
  BIC_params <- BIC_simul$min_BIC %>%
    arrange(BIC_log_p) %>%
    head(1)
  
  result <- BIC_simul$simulation[[which(lamb_seq == BIC_params$lambda)]]
  
  result
  
}
stopCluster(cl)


################################
## 1-3. Simulation - SP_model ##
################################

simul_high_SP_model_2 <- vector("list", length = simul_times)
for(simul in 1:simul_times) {
  
  tau_seq <- seq(from = 0.35, to = 0.65, length.out = b)
  X <- X_list[[simul]]
  Y <- Y_list[[simul]]
  V <- V_list[[simul]]
  
  lasso_coef <- matrix(nrow = p+1, ncol = m)
  for(g in 1:m) {
    cv.lasso <- cv.glmnet(x = X[, -1], y = Y[, g], alpha = 1, type.measure = "mae")
    lasso_model <- glmnet(X[, -1], Y[, g], family = "gaussian", alpha = 1, lambda = cv.lasso$lambda.min)
    lasso_coef[, g] <- c(lasso_model$a0, as.vector(lasso_model$beta))
  }
  
  first_init_SP <- matrix(nrow = (p+1)*K, ncol = m)
  for(g in 1:m) {
    for(j in 0:p) {
      first_init_SP[((j*K)+1):((j+1)*K), g] <- lasso_coef[j+1, g]
    }
  }
  
  init_val_SP <- SP_model(delta = 1, lambda = 0.05, tol_error = 0.1^5, max_iter = 50, 
                          X = X, Y = Y, V = V, Phi = Phi, theta_0 = first_init_SP, tau_seq = tau_seq, weight = FALSE)
  
  log_lamb <- c(seq(5, 8, length.out = 20))
  lamb_seq <- exp(log_lamb)
  
  BIC_table <- list()
  cl <- makeCluster(20) #not to overload your computer
  registerDoParallel(cl) # Ready to parallel
  
  BIC_table <- foreach(lambda = lamb_seq, .noexport = "SP_model") %dopar% {
    library(dplyr)
    library(splines)
    library(Matrix)
    library(glmnet)
    library(fda)
    library(Rcpp)
    sourceCpp("[KSH]add_decomp_function.cpp")
    
    BIC_simul <- SP_model_BIC(X, Y, V, Phi, theta_0 = init_val_SP$theta, 
                              tau_seq, tau_seq_real, lamb_seq = lambda, max_iter = 50, delta = 1, fun_type = "cpp")
    BIC_simul$BIC_data
  }
  stopCluster(cl)
  
  BIC_params <- BIC_table %>%
    bind_rows() %>%
    arrange(BIC_log_p) %>%
    mutate_all(as.numeric) %>%
    head(1)
  
  result <- SP_model(delta = 1, lambda = BIC_params$lambda, tol_error = 0.1^5, max_iter = 50, 
                     X = X, Y = Y, V = V, Phi = Phi, theta_0 = init_val_SP$theta, tau_seq = tau_seq, weight = TRUE)
  simul_high_SP_model_2[[simul]] <- result
}


###############
## Save Data ##
###############

save(simul_high_add_decomp_2, simul_high_LR_model_2, simul_high_SP_model_2, 
     LR_mat, sp_mat, Phi, tau_seq, tau_seq_real, X_list, file = "ksh_simul_mixture_error_high_2.RData")


