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
set.seed(1)
n <- 400
m <- 10
p <- 100
b <- 15
num_rank <- 5
num_rank_X <- 30
simul_times <- 100

sigma_mat <- matrix(nrow = p, ncol = p)
for(j in 1:p) {
  for(k in 1:p) {
    sigma_mat[j, k] <- 0.5^(abs(j-k))
  }
}

# sparse matrix
col_ind <- sample(1:m, size = m, replace = FALSE)
row_ind <- sample(1:(p+1), size = m)
sp_mat <- matrix(0, nrow = p+1, ncol = m)
for(i in 1:m) {
  sp_mat[row_ind[i], col_ind[i]] <- rnorm(1, mean = 5, sd = 0.1)
}

# low rank matrix
L1 <- matrix(rnorm((p+1)*num_rank, mean = 0, sd = 0.3), nrow = p+1)
L2 <- matrix(rnorm(m*num_rank, mean = 0, sd = 0.3), nrow = m)
LR_mat <- L1 %*% t(L2)

X_list <- list()
eps_list <- list()

for(simul in 1:simul_times) {
  X1 <- matrix(rnorm(n*num_rank_X, mean = 0, sd = 1), nrow = n)
  X2 <- matrix(rnorm(p*num_rank_X, mean = 0, sd = 1), nrow = p)
  X <- X1 %*% t(X2) 
  X_list[[simul]] <- X %*% expm::sqrtm(sigma_mat) %>% cbind(rep(1, n), .)
  
  sigma <- runif(n,0.3,0.7)
  eps_list[[simul]] <- MASS::mvrnorm(n = m, mu = rep(0,n), Sigma = diag(sigma^2,n)) %>% t()
}

Y_list <- mapply(FUN = function(X, LR, SP, eps) X %*% (LR + SP) + eps, 
                 X_list, list(LR_mat), list(sp_mat), eps_list, SIMPLIFY = FALSE)

### Calculate kronecker product
K <- 10
tau_seq <- seq(from = 0.35, to = 0.65, length.out = b)
tau_seq_real <- tau_seq[tau_seq >= 0.4 & tau_seq  <= 0.6]

knots_seq <- seq(min(tau_seq)- 0.02, max(tau_seq) + 0.02, length.out = K)
Phi <- fda::bsplineS(tau_seq, breaks= knots_seq, norder=2, nderiv=0, returnMatrix=FALSE)

V_list <- list()
for(simul in 1:simul_times) {
  V_list[[simul]] <- calc_V(X_list[[simul]], Phi)
}

##################################
## 1-1. Simulation - add_decomp ##
##################################

cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
simul_low_add_decomp <- foreach(simul = 1:simul_times, .noexport = "add_decomp") %dopar% {
  library(dplyr)
  library(splines)
  library(Matrix)
  library(glmnet)
  library(fda)
  
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
  
  init_val <- add_decomp_r(delta = 1, lambda_1 = 0.03, lambda_2 = 0.01, tol_error = 0.1^5, max_iter = 50,
                           X = X, Y = Y, V = V, Phi = Phi, 
                           theta_0 = theta_init, Z_0 = X%*%alpha_init, tau_seq = tau_seq, weight = FALSE)
  
  lamb1_seq <- c( seq(0.001, 0.15, length.out = 10))
  lamb2_seq <- c(seq(100, 300, length.out = 10))
  BIC_simul <- add_decomp_BIC(X, Y, V, Phi, theta_0 = init_val$theta, Z_0 = init_val$Z, tau_seq, tau_seq_real, 
                              lamb1_seq = lamb1_seq, lamb2_seq = lamb2_seq, max_iter = 50)
  
  
  BIC_params <- BIC_simul$table %>%
    mutate(LR = log(log(log(p))) * LR_part, 
           SP = log(n)*log(p) * SP_part, 
           BIC = log_Q + LR + SP) %>%
    arrange(BIC) %>% 
    head(1)
  
  result <- BIC_simul$simulation[[which(lamb1_seq == BIC_params$lambda_1)]][[which(lamb2_seq == BIC_params$lambda_2)]]
  
  result
  
}
stopCluster(cl)

################################
## 1-2. Simulation - LR_model ##
################################

cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
simul_low_LR_model <- foreach(simul = 1:simul_times, .noexport = "add_decomp") %dopar% {
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
  
  lamb_seq <- seq(0.1^20, 0.1, length.out = 25)
  BIC_simul <- LR_model_BIC(X, Y, Z_0 = init_val_LR$Z, tau_seq, tau_seq_real, lamb_seq, max_iter = 50)
  
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

cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
simul_low_SP_model <- foreach(simul = 1:simul_times, .noexport = "add_decomp") %dopar% {
  library(dplyr)
  library(splines)
  library(Matrix)
  library(glmnet)
  library(fda)
  
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
  
  init_val_SP <- SP_model_r(delta = 1, lambda = 0.05, tol_error = 0.1^5, max_iter = 50, 
                            X = X, Y = Y, V = V, Phi = Phi, theta_0 = first_init_SP, tau_seq = tau_seq, weight = FALSE)
  
  lamb_seq <- seq(0.05, 5, length.out = 20)
  BIC_simul <- SP_model_BIC(X, Y, V, Phi, theta_0 = init_val_SP$theta, 
                            tau_seq, tau_seq_real, lamb_seq, max_iter = 50)
  
  BIC_params <- BIC_simul$min_BIC %>%
    arrange(BIC_log_p) %>%
    head(1)
  
  result <- BIC_simul$simulation[[which(lamb_seq == BIC_params$lambda)]]
  
  result
  
}
stopCluster(cl)


###############
## Save Data ##
###############

save(simul_low_add_decomp, simul_low_LR_model, simul_low_SP_model, 
     LR_mat, sp_mat, Phi, tau_seq, tau_seq_real, X_list, file = "ksh_simul_normal_error_low.RData")


