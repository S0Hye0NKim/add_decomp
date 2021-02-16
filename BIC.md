BIC
================
Sohyeon Kim
2021 2 16



``` r
library(tidyverse)
library(ggplot2)
library(splines)
library(Matrix)
library(pander)
library(gridExtra)
library(foreach)
library(doParallel)
library(Rcpp)
library(glmnet)
library(fda)

sourceCpp("functions/add_decomp_function.cpp")
source("https://raw.githubusercontent.com/S0Hye0NKim/add_decomp/master/functions/add_decomp_function.R")
```

# Toy Example

``` r
set.seed(0)
n <- 400
m <- 10
p <- 30
#num_nz <- 15   #round((p+1)*m*(1/5)) 
b <- 31
num_rank <- 5

X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) %>% cbind(1, .)   #add intercept term in X

#sp_mat <- rsparsematrix(nrow = p+1, ncol = m, nnz = num_nz)             # make sparse matrix
#for(i in 1:(p+1)) {
  #for(j in 1:m) {
    #if(abs(sp_mat[i, j]) != 0) {sp_mat[i, j] <- rnorm(1, mean = 5, sd = 0.1)}
  #}
#}
col_ind <- sample(1:m, size = m, replace = FALSE)
row_ind <- sample(1:(p+1), size = m)
sp_mat <- matrix(0, nrow = p+1, ncol = m)
for(i in 1:m) {
  sp_mat[row_ind[i], col_ind[i]] <- rnorm(1, mean = 5, sd = 0.1)
}
num_zero <- which(sp_mat==0, arr.ind = TRUE) %>% nrow
num_nz <- (p+1)*m - num_zero

LR_mat <- matrix(rnorm((p+1)*m, mean = 1, sd = 0.1), ncol = m) # make low rank matrix using SVD
SVD <- svd(LR_mat)
D_mat <- diag(SVD$d)
idx <- (num_rank+1):m
D_mat[idx, idx] <- 0
LR_mat <- SVD$u %*% D_mat %*% t(SVD$v)          

true_B <- sp_mat + LR_mat            # B(tau) = sparse matrix + low rank matrix

eps <- matrix(rnorm(n*m, mean = 0, sd = 0.1), nrow = n)

Y <- X%*%true_B + eps
```

``` r
K <- 10
tau_seq <- seq(from = 0.35, to = 0.65, length.out = b)
tau_seq_real <- tau_seq[tau_seq > 0.39 & tau_seq < 0.61]

knots_seq <- seq(min(tau_seq)- 0.02, max(tau_seq) + 0.02, length.out = K)
Phi <- fda::bsplineS(tau_seq, breaks= knots_seq, norder=2, nderiv=0, returnMatrix=FALSE)
V <- calc_V(X, Phi)

# Initial value
lasso_coef <- matrix(nrow = p+1, ncol = m)
for(g in 1:m) {
  cv.lasso <- cv.glmnet(x = X[, -1], y = Y[, g], alpha = 1, type.measure = "mae")
  lasso_model <- glmnet(X[, -1], Y[, g], family = "gaussian", alpha = 1, lambda = cv.lasso$lambda.min)
  lasso_coef[, g] <- c(lasso_model$a0, as.vector(lasso_model$beta))
}


theta_init <- matrix(nrow = (p+1)*K, ncol = m)
for(g in 1:m) {
  for(j in 0:p) {
    theta_init[((j*K)+1):((j+1)*K), g] <- lasso_coef[j+1, g]
  }
}

Y_modified <- Y - X%*%lasso_coef

lin_model <- lm(Y_modified~., data = data.frame(X[, -1]))
alpha_init <- lin_model$coefficients
```

## add\_decomp

``` r
lamb1_seq <- seq(0.001, 0.01, by = 0.001)
lamb2_seq <- seq(1500, 1700, by = 10)


start <- Sys.time()
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl) # Ready to parallel

simulation <- foreach(lambda_1 = lamb1_seq, .noexport = "add_decomp") %:%
  foreach(lambda_2 = lamb2_seq, .noexport = "add_decomp") %dopar% {
    library(tidyverse)
    library(splines)
    library(Matrix)
    library(Rcpp)
    library(glmnet)
    library(fda)
    sourceCpp("functions/add_decomp_function.cpp")
    
    result <- add_decomp(delta = 1, lambda_1 = lambda_1, lambda_2 = lambda_2, tol_error = 0.001, 
                         max_iter = 50, X = X, Y = Y, V = V, Phi = Phi, 
                         theta_0 = theta_init, alpha_0 = alpha_init, tau_seq = tau_seq)
    
    result
  }



stopCluster(cl)
end <- Sys.time()
```

## BIC

  
![\\begin{aligned}BIC=&\\;log\\bigg(\\sum\_{g=1}^m\\sum\_{\\ell=1}^b\\sum\_{i=1}^n\\frac{1}{n}\\rho\_{\\tau\_\\ell}(Y\_i^{(g)}-X\\hat{\\alpha}^{(g)}-V^{(\\ell)}\\hat{\\theta}^{(g)})\\bigg)\\\\&+\\Big(Penalty\\Big)\\Big(\\frac{\\hat{r}(m+p-\\hat{r})+K\\hat{S}}{2n}\\Big)\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7DBIC%3D%26%5C%3Blog%5Cbigg%28%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%5Cfrac%7B1%7D%7Bn%7D%5Crho_%7B%5Ctau_%5Cell%7D%28Y_i%5E%7B%28g%29%7D-X%5Chat%7B%5Calpha%7D%5E%7B%28g%29%7D-V%5E%7B%28%5Cell%29%7D%5Chat%7B%5Ctheta%7D%5E%7B%28g%29%7D%29%5Cbigg%29%5C%5C%26%2B%5CBig%28Penalty%5CBig%29%5CBig%28%5Cfrac%7B%5Chat%7Br%7D%28m%2Bp-%5Chat%7Br%7D%29%2BK%5Chat%7BS%7D%7D%7B2n%7D%5CBig%29%5Cend%7Baligned%7D
"\\begin{aligned}BIC=&\\;log\\bigg(\\sum_{g=1}^m\\sum_{\\ell=1}^b\\sum_{i=1}^n\\frac{1}{n}\\rho_{\\tau_\\ell}(Y_i^{(g)}-X\\hat{\\alpha}^{(g)}-V^{(\\ell)}\\hat{\\theta}^{(g)})\\bigg)\\\\&+\\Big(Penalty\\Big)\\Big(\\frac{\\hat{r}(m+p-\\hat{r})+K\\hat{S}}{2n}\\Big)\\end{aligned}")  

  - n : \# of samples.

  - p : \# of gene.

  - m : \# of response.

  - ![\\hat{r}](https://latex.codecogs.com/png.latex?%5Chat%7Br%7D
    "\\hat{r}") :
    ![\\text{rank}(\\hat{A})](https://latex.codecogs.com/png.latex?%5Ctext%7Brank%7D%28%5Chat%7BA%7D%29
    "\\text{rank}(\\hat{A})").

  - K : \# of spline function.

  - ![\\hat{S}](https://latex.codecogs.com/png.latex?%5Chat%7BS%7D
    "\\hat{S}") : \# of non-zero entries in
    ![\\hat{\\Gamma}(\\tau)](https://latex.codecogs.com/png.latex?%5Chat%7B%5CGamma%7D%28%5Ctau%29
    "\\hat{\\Gamma}(\\tau)") regardless of quantiles.

<!-- end list -->

``` r
idx_tau <- (tau_seq > 0.39 & tau_seq <= 0.61)

BIC <- list()
for(i in 1:length(lamb1_seq)) {
  BIC[[i]] <- list()
  for(j in 1:length(lamb2_seq)) {
    result <- simulation[[i]][[j]]
    est_error <- lapply(V[idx_tau], FUN = function(x) (Y - X %*% result$alpha - x %*% result$theta)
                    %>% as.vector())
    check_loss_err <- mapply(FUN = function(x, tau) check_ft(x, tau), x = est_error, 
                           tau = as.list(tau_seq_real), SIMPLIFY = FALSE) %>%
      lapply(FUN = function(x) sum(x)/n) %>% unlist %>% sum
    gamma_tau_hat <- est_gamma(Phi[idx_tau, ], result$theta)
    S_hat <- check_sp_table(true = matrix(0, nrow = (p+1), ncol = m), 
               est = gamma_tau_hat, table = TRUE, tol = 0.1^5, tau_seq = tau_seq_real) %>%
      .$Est_Positive %>% sum
    BIC[[i]][[j]] <- data.frame(log_Q = log(check_loss_err), 
                                r_hat = rankMatrix(result$alpha)[1], 
                                S_hat = S_hat)
  }
}

BIC_data <- lapply(BIC, FUN = function(x) `names<-`(x, value = lamb2_seq) %>% bind_rows(.id = "lambda_2")) %>%
  `names<-`(value = lamb1_seq) %>%
  bind_rows(.id = "lambda_1") %>%
  mutate(term = (r_hat * (m + p - r_hat) + K * S_hat)/(2*n), 
         BIC_log_sum = log_Q + log(p+m)*term, 
         BIC_log_p = log_Q + log(p)*term, 
         BIC_log_n = log_Q + log(n)*term, 
         BIC_llog_p = log_Q + log(log(p))*term,
         BIC_llog_n = log_Q + log(log(n))*term) %>%
  group_by(lambda_1, lambda_2) %>%
  select_at(vars(starts_with("BIC"))) %>%
  ungroup()
```

``` r
BIC_val_min <- apply(select_at(BIC_data, vars(starts_with("BIC"))), 2, min) %>%
  `names<-`(value = c("log_sum", "log_p", "log_n", "llog_p", "llog_n"))
```

``` r
filter(BIC_data, BIC_log_sum == BIC_val_min["log_sum"]| BIC_log_p == BIC_val_min["log_p"] |
         BIC_log_n == BIC_val_min["log_n"] | BIC_llog_p == BIC_val_min["llog_p"] | 
         BIC_llog_n == BIC_val_min["llog_n"]) %>%
  pander()
```

| lambda\_1 | lambda\_2 | BIC\_log\_sum | BIC\_log\_p | BIC\_log\_n | BIC\_llog\_p | BIC\_llog\_n |
| :-------: | :-------: | :-----------: | :---------: | :---------: | :----------: | :----------: |
|   0.002   |   1560    |     5.626     |    5.576    |    6.026    |    5.198     |    5.296     |
