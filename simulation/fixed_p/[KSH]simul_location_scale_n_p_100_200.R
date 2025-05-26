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
library(rqPen)

sourceCpp("[KSH]add_decomp_function.cpp")
source("https://raw.githubusercontent.com/S0Hye0NKim/add_decomp/master/functions/add_decomp_function.R")


#################
## 1-0. Set up ##
#################

## Generate data
set.seed(3)
n <- 100
m <- 10
p <- 200
b <- 15
num_rank <- 3
num_rank_X <- 50
simul_times <- 50

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
  sp_mat[row_ind[i], col_ind[i]] <- rnorm(1, mean = 8, sd = 0.1)
}
sp_mat[1, ] <- rnorm(m, mean = 8, sd = 0.1)

inc_func_idx <- row_ind[which(col_ind == 1)] - 1


# low rank matrix
L1 <- matrix(rnorm(p*num_rank, mean = 0, sd = 0.4), nrow = p)
L2 <- matrix(rnorm(m*num_rank, mean = 0, sd = 0.4), nrow = m)
LR_mat <- L1 %*% t(L2)

X_list <- list()
eps_list <- list()



for(simul in 1:simul_times) {
  X1 <- matrix(rnorm(n*num_rank_X, mean = 0, sd = 1), nrow = n)
  X2 <- matrix(rnorm(p*num_rank_X, mean = 0, sd = 1), nrow = p)
  X <- X1 %*% t(X2) %*% expm::sqrtm(sigma_mat)
  X[, inc_func_idx] <- pnorm(X[, inc_func_idx]) * sqrt(12)
  X_list[[simul]] <- X %>% cbind(rep(1, n), .)
  
  eps_list[[simul]] <- matrix(rnorm(n*m, mean = 0, sd = 1), nrow = n, ncol = m)
}


Y_list <- mapply(FUN = function(X, LR, SP, eps) X[, -1] %*% LR + X %*% SP + matrix(c(3 * X[, (inc_func_idx + 1)], rep(1, n*(m-1))), nrow = n, ncol = m) * eps , 
                 X_list, list(LR_mat), list(sp_mat), eps_list, SIMPLIFY = FALSE)


### Calculate kronecker product
K <- (2 * n^(1/5)) %>% round
tau_seq <- seq(from = 0.35, to = 0.65, length.out = b)
tau_seq_real <- tau_seq[tau_seq >= 0.4 & tau_seq  <= 0.6]
idx_tau <- (tau_seq >= "0.4" & tau_seq <= "0.6")

knots_seq <- seq(min(tau_seq)- 0.02, max(tau_seq) + 0.02, length.out = K)
Phi <- fda::bsplineS(tau_seq, breaks= knots_seq, norder=2, nderiv=0, returnMatrix=FALSE)

V_list <- list()
for(simul in 1:simul_times) {
  V_list[[simul]] <- calc_V(X_list[[simul]], Phi)
}


### quantile regression for initial value
cl <- makeCluster(20) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
QR_Lasso <- foreach(simul = 1:simul_times) %:%
  foreach(g = 1:m) %dopar% {
    library(dplyr)
    library(splines)
    library(rqPen)

    X <- X_list[[simul]]
    Y <- Y_list[[simul]]

    QR_model <- rqPen::rq.lasso.fit.mult(x = X[, -1], y = Y[, g], tau_seq = tau_seq, lambda = 0.0001)
    QR_Lasso_coef <- lapply(QR_model, FUN = function(x) x$coefficients %>% data.frame(value = .) %>%
                            tibble::rownames_to_column(var = "variable")) %>%
      `names<-`(tau_seq) %>%
      bind_rows(.id = "tau") %>%
      mutate(variable = ifelse(variable == "intercept", "x0", variable))
    QR_Lasso_coef
  }
stopCluster(cl)

QR_Lasso_data <- lapply(QR_Lasso, FUN = function(x) x %>% `names<-`(value = 1:m) %>% bind_rows(.id = "group")) %>%
  `names<-`(value = 1:simul_times) %>%
  bind_rows(.id = "simulation")


#########################
## First Initial value ##
#########################
init_val_AD <- list()
init_val_LR <- list()
init_val_SP <- list()

for(simul in 1:simul_times) {
    X <- X_list[[simul]]
    Y <- Y_list[[simul]]
    V <- V_list[[simul]]
  
    first_init_SP <- matrix(nrow = (p+1)*K, ncol = m)
    for(g in 1:m) {  
        for(j in 0:p) {
            QR_coef <- QR_Lasso_data %>% filter(simulation == simul, group == g, variable == paste0("x", j)) %>% .$value
            idx <- seq(2, b-1, length.out = K)
      
            est_theta <- solve(Phi[idx, ], QR_coef[idx])

            first_init_SP[((j*K)+1):((j+1)*K), g] <- est_theta
        }
    }
  
    init_val_SP[[simul]] <- SP_model(delta = 1, lambda = 0.01, tol_error = 0.1^5, max_iter = 50, 
                            X = X, Y = Y, V = V, Phi = Phi, theta_0 = first_init_SP, tau_seq = tau_seq, weight = FALSE)

    lasso_coef <- matrix(nrow = p+1, ncol = m)
    for(g in 1:m) {
        cv.lasso <- cv.glmnet(x = X[, -1], y = Y[, g], 
                              alpha = 1, type.measure = "mae")
        lasso_model <- glmnet(X[, -1], Y[, g], 
                              family = "gaussian", alpha = 1, lambda = cv.lasso$lambda.min)
        lasso_coef[, g] <- c(lasso_model$a0, as.vector(lasso_model$beta))
    }

    Y_modified <- Y - X%*%lasso_coef
    ridge_coef_AD <- matrix(nrow = p+1, ncol = m)
    for(g in 1:m) {
        cv.ridge <- cv.glmnet(x = X[, -1], y = Y_modified[, g], alpha = 0, type.measure = "mae")
        ridge_model <- glmnet(X[, -1], Y_modified[, g], family = "gaussian", alpha = 0, lambda = cv.ridge$lambda.min)
        ridge_coef_AD[, g] <- c(ridge_model$a0, as.vector(ridge_model$beta))
    }
    alpha_init <- ridge_coef_AD


    init_val_AD[[simul]] <- add_decomp(delta = 1, lambda_1 = 0.01, lambda_2 = 0.01, tol_error = 0.1^5, max_iter = 50,
                                       X = X, Y = Y, V = V, Phi = Phi, 
                                       theta_0 = init_val_SP[[simul]]$theta, Z_0 = X%*%alpha_init, tau_seq = tau_seq, weight = FALSE)
    
    ridge_coef_LR <- matrix(nrow = p+1, ncol = m)
    for(g in 1:m) {
        cv.ridge <- cv.glmnet(x = X[, -1], y = Y[, g], alpha = 0, type.measure = "mae")
        ridge_model <- glmnet(X[, -1], Y[, g], family = "gaussian", alpha = 0, lambda = cv.ridge$lambda.min)
        ridge_coef_LR[, g] <- c(ridge_model$a0, as.vector(ridge_model$beta))
    }
  first_init_LR <- ridge_coef_LR
  
  init_val_LR[[simul]] <- LR_model(delta = 1, lambda = 0.01, tol_error = 0.1^5, max_iter = 50, 
                            X = X, Y = Y, Z_0 = X %*% first_init_LR, tau_seq = tau_seq, weight = FALSE)
}



##################################
## 1. Additive decomposed model ##
##################################
simul_add_decomp <- vector("list", length = simul_times)
for(simul in 1:simul_times) {

  X <- X_list[[simul]]
  Y <- Y_list[[simul]]
  V <- V_list[[simul]]
  init_val <- init_val_AD[[simul]]
  
  log_lamb1 <- c( seq(1.5, 2.3, length.out = 20))
  lamb1_seq <- exp(log_lamb1)
  log_lamb2 <- c(seq(4, 4.3, length.out = 20))
  lamb2_seq <- exp(log_lamb2)
  
  BIC_table <- list()
  
  for(idx in 1:length(lamb1_seq)) {
    lamb1 <- lamb1_seq[idx]
    cl <- makeCluster(20) #not to overload your computer
    registerDoParallel(cl) # Ready to parallel
    
    temp_BIC <- foreach(lamb2 = lamb2_seq, .noexport = "add_decomp") %dopar% {
      library(dplyr)
      library(Matrix)
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
  
  simul_add_decomp[[simul]] <- result
}


#######################
## 2. Low-rank Model ##
#######################

cl <- makeCluster(20) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
simul_LR_model <- foreach(simul = 1:simul_times, .noexport = "add_decomp") %dopar% {
  library(dplyr)
  library(Matrix)
  
  X <- X_list[[simul]]
  Y <- Y_list[[simul]]
  init_val <- init_val_LR[[simul]]
  
  lamb_seq <- seq(0.1, 5, length.out = 20)
  r_X <- rankMatrix(X[, -1])
  BIC_simul <- LR_model_BIC(X, Y, Z_0 = init_val$Z, tau_seq, tau_seq_real, lamb_seq, max_iter = 50, delta = 1, r_X = rankMatrix(X[, -1]))
  
  r_X <- rankMatrix(X[, -1])[1]
  BIC_params <- BIC_simul$BIC_data %>%
                  mutate(LR_part = r_hat * max(r_X, m) / (2*n*m), 
                         LR = log(p) * log(log(n)) * LR_part, 
                         BIC = log_Q + LR) %>%
                  mutate_all(as.numeric) %>%
                  arrange(BIC) %>% 
                  head(1)

  result <- BIC_simul$simulation[[which(lamb_seq == BIC_params$lambda)]]
  
  result
  
}
stopCluster(cl)


#####################
## 3. Sparse Model ##
#####################
simul_SP_model <- vector("list", length = simul_times)
for(simul in 1:simul_times) {
  
  X <- X_list[[simul]]
  Y <- Y_list[[simul]]
  V <- V_list[[simul]]
  init_val <- init_val_SP[[simul]]
  
  log_lamb <- c(seq(3, 4.3, length.out = 20))
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
    
    BIC_simul <- SP_model_BIC(X, Y, V, Phi, theta_0 = init_val$theta, 
                              tau_seq, tau_seq_real, lamb_seq = lambda, max_iter = 50, delta = 1, fun_type = "cpp")
    BIC_simul$BIC_data
  }
  stopCluster(cl)
  
  BIC_params <- BIC_table %>%
    bind_rows() %>%
    mutate(S_hat_net = S_hat - num_nz_intercept,
           SP = log(p) * log(log(n)) * K * S_hat_net / (2*n*m),
           BIC = log_Q + SP) %>%
    mutate_all(as.numeric) %>%
    filter(S_hat_net != 0) %>%
    arrange(BIC) %>%  
    head(1)

  result <- SP_model(delta = 1, lambda = BIC_params$lambda, tol_error = 0.1^5, max_iter = 50, 
                     X = X, Y = Y, V = V, Phi = Phi, theta_0 = init_val$theta, tau_seq = tau_seq, weight = TRUE)
  simul_SP_model[[simul]] <- result
}


###############
## Save Data ##
###############

save(simul_add_decomp, simul_LR_model, simul_SP_model, est_gamma, check_sp_table, 
     LR_mat, sp_mat, Phi, tau_seq, tau_seq_real, X_list, file = "ksh_simul_location_scale_n_p_100_200.RData")

