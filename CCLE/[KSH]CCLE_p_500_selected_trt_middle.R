# Packages
rm(list = ls())
library(dplyr)
library(stringr)
library(Matrix)
library(data.table)
library(Rcpp)
library(splines)
library(fda)
library(glmnet)
library(foreach)
library(doParallel)
library(tibble)
library(corpcor)
library(rqPen)

source("https://raw.githubusercontent.com/S0Hye0NKim/add_decomp/master/functions/add_decomp_function.R")
sourceCpp("[KSH]add_decomp_function.cpp")

#########################
## 0. Data & Screening ##
#########################

## 0-(1) Screening Y
Y <- fread("CCLE_Y_new.csv")

m <- ncol(Y)

Y_sd <- Y[, lapply(.SD, sd)] %>% as.matrix %>% t()
data.frame(Y_sd) %>%
  mutate(group = rownames(.)) %>%
  dplyr::select(group, Y_sd) %>%
  arrange(Y_sd) %>% head(5)

Y <- Y[, -"Panobinostat"] %>%
  column_to_rownames(var = "CCLE_Name") %>%
  as.matrix

## Select treatment which has wide boxplot
quant_diff <- apply(Y, 2, quantile, probs = c(0.25,0.75)) %>%
  apply(2, diff)

drug_idx <- colnames(Y)[quant_diff != 0]

Y <- Y[, drug_idx]

m <- ncol(Y)
n <- nrow(Y)

## 0-(2) Screening X

### 0-(2)-(1) Standard deviation of X
X <- fread("CCLE_X_new.csv")
X <- X[, -"Broad_ID"]

X_sd <- X[, -"CCLE_Name"][, lapply(.SD, FUN = sd)] %>%
  as.matrix() %>% t()
idx_scr_sd <-  data.frame(X_sd) %>%
  mutate(variable = rownames(.)) %>%
  dplyr::select(variable, X_sd) %>%
  arrange(-X_sd) %>% 
  .[1:10000, ] %>% .$variable

X <- dplyr::select(X, c(CCLE_Name,all_of(idx_scr_sd))) %>%
  column_to_rownames(var = "CCLE_Name")

p <- 500

### 0-(2)-(2) Correlation of X and Y
idx_scr_cor <- cor(X, Y) %>% abs %>%
  apply(MARGIN = 1, FUN = mean) %>%
  data.frame(variable = names(.), cor = .) %>%
  arrange(-cor) %>% first(p) %>% .$variable

X <- dplyr::select(X, all_of(idx_scr_cor)) %>% 
  apply(MARGIN = 2, FUN = function(x) (x - mean(x))/(sd(x)*sqrt((n-1)/n))) %>%
  cbind(1, .)
colnames(X)[1] <- "intercept"
gene_ex <- colnames(X)

### 0-(3) spline
K <- (2 * n^(1/5)) %>% round
b <- 30
tau_seq <- seq(from = 0.1, to = 0.9, length.out = b)
idx_tau <- (tau_seq >= "0.25" & tau_seq <= "0.75")
tau_seq_real <- tau_seq[idx_tau]

knots_seq <- seq(min(tau_seq) - 0.02, max(tau_seq) + 0.02, length.out = K)
Phi <- fda::bsplineS(tau_seq, breaks= knots_seq, norder=2, nderiv=0, returnMatrix=FALSE)

V <- calc_V(as.matrix(X), Phi)


### Treatment / gene expression data
trt_data <- data.frame(col = 1:ncol(Y), trt_nm = colnames(Y))
gene_ex_data <- data.frame(col = 1:ncol(X), gene_ex_nm = colnames(X))

#########################
## First Initial value ##
#########################


### quantile regression for initial value
set.seed(3)
X <- as.matrix(X)
Y <- as.matrix(Y)

cl <- makeCluster(20) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
QR_Lasso <- foreach(g = 1:m) %dopar% {
  library(dplyr)
  library(splines)
  library(rqPen)

  median <- median(Y[, g])

  init_lamb <- ifelse(median == 8, 0.0001, 0.001)

  QR_model <- rqPen::rq.lasso.fit.mult(x = X[, -1], y = Y[, g], tau_seq = tau_seq, lambda = init_lamb)
  QR_Lasso_coef <- lapply(QR_model, FUN = function(x) x$coefficients %>% data.frame(value = .) %>%
                            tibble::rownames_to_column(var = "variable")) %>%
      `names<-`(tau_seq) %>%
      bind_rows(.id = "tau")

  QR_Lasso_coef
}
stopCluster(cl)

QR_Lasso_data <- QR_Lasso %>% `names<-`(value = 1:m) %>% bind_rows(.id = "group") %>%
  left_join(gene_ex_data, by = c("variable" = "gene_ex_nm")) %>%
  rename(gene_ex_nm = variable) %>%
  mutate(variable = paste0("x", col - 1))

# First initial value

first_init_SP <- matrix(nrow = (p+1)*K, ncol = m)
for(g in 1:m) {  
  for(j in 0:p) {
    QR_coef <- QR_Lasso_data %>% filter(group == g, variable == paste0("x", j)) %>% .$value
    idx <- seq(2, b-1, length.out = K)
      
    est_theta <- solve(Phi[idx, ], QR_coef[idx])

    first_init_SP[((j*K)+1):((j+1)*K), g] <- est_theta
  }
}
  
init_val_SP <- SP_model(delta = 1, lambda = 0.01, tol_error = 0.1^5, max_iter = 50, 
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

init_val_AD <- add_decomp(delta = 1, lambda_1 = 0.01, lambda_2 = 0.01, tol_error = 0.1^5, max_iter = 50,
                                       X = X, Y = Y, V = V, Phi = Phi, 
                                       theta_0 = init_val_SP$theta, Z_0 = X%*%alpha_init, tau_seq = tau_seq, weight = FALSE)


##################################
## 1. Additive decomposed model ##
##################################
log_lamb1 <- c( seq(3.3, 3.7, length.out = 20))
lamb1_seq <- exp(log_lamb1)
log_lamb2 <- c(seq(1, 2.2, length.out = 20))
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
      
    BIC_simul <- add_decomp_BIC(X, Y, V, Phi, theta_0 = init_val_AD$theta, Z_0 = init_val_AD$Z, tau_seq, tau_seq_real,
                                lamb1_seq = lamb1, lamb2_seq = lamb2, max_iter = 50, delta = 1, fun_type = "cpp")
    BIC_simul$table
  }
  
  stopCluster(cl)
  BIC_table[[idx]] <- temp_BIC
}
  
r_X <- rankMatrix(X[, -1])[1]

GIC_table <- BIC_table %>% lapply(FUN = function(x) bind_rows(x)) %>%
  bind_rows() %>%
  mutate(LR_part = r_hat * max(1, m) / (2*n*m),
         S_hat_net = S_hat - num_nz_intercept,
         LR = log(p) * log(log(n)) * LR_part, 
         SP = log(p) * log(log(n)) * K * S_hat_net / (2*n*m),
         BIC = log_Q + LR + SP) %>%
  mutate_all(as.numeric) %>%
  filter(S_hat_net != 0) %>%
  arrange(BIC)
  
BIC_params <- GIC_table %>% head(1)
  
result_AD <- add_decomp(delta = 1, lambda_1 = BIC_params$lambda_1, lambda_2 = BIC_params$lambda_2, 
                       tol_error = 0.1^5, max_iter = 50, X, Y, V, Phi, 
                       theta_0 = init_val_AD$theta, Z_0 = init_val_AD$Z, tau_seq = tau_seq, weight = TRUE)


###############
## Save data ##
###############

save.image(file = "ksh_CCLE_p_500_selected_trt_middle.RData")