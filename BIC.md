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

``` r
delta = 1
eta_old <- theta_init
  theta_old <- eta_old
  alpha_old <- alpha_init
  Z_old <- X %*% alpha_old
  e_old <- list()
  for(l in 1:b) {e_old[[l]] <- Y - Z_old - V[[l]] %*% eta_old}
  u_old <- list()
  for(l in 1:b) {u_old[[l]] <- matrix(0, nrow = n, ncol = m)}
  w_old <- matrix(0, nrow = (p+1)*K, ncol = m)

Y_list <- list()
    for(l in 1:b) {Y_list[[l]] <- Y}
    VH_list <- lapply(V, FUN = function(x) x %*% theta_init)
    obj_list <- mapply(function(Y, VH, E, U) Y - VH - E - U/delta, Y_list, VH_list, e_old, u_old, SIMPLIFY = FALSE)
    obj <- Reduce("+", obj_list)/b 
    SVD <- svd(obj)
    sing_val_alpha_0 <- svd(alpha_init) %>% .$d  # weight = 1/sing_val_alpha_0
delta*b*SVD$d[1]*sing_val_alpha_0[1] 
```

    ## [1] 323.7015

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
idx_tau <- tau_seq %in% tau_seq_real

BIC_func(X, Y, V, Phi, theta_0 = theta_init, alpha_0 = alpha_init, tau_seq, tau_seq_real = tau_seq[idx_tau],
         lamb1_seq = lamb1_seq, lamb2_seq = lamb2_seq, table = FALSE, max_iter = 50)
```

# boundary for lambda2

  
![Q=\\sum\_{g=1}^m\\sum\_{\\ell=1}^b\\sum\_{i=1}^n\\frac{1}{n}\\rho\_{\\tau\_\\ell}\\bigg(Y\_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum\_{j=0}^px\_{ij}\\Phi(\\tau\_\\ell)^T\\boldsymbol{\\theta}\_j^{(g)}\\bigg)+\\lambda\_1||\\textbf{XA}||\_{\*\\psi}+\\lambda\_2\\sum\_{g=1}^m\\sum\_{j=1}^p\\phi\_j^{(g)}||\\boldsymbol{\\theta}^{(g)}\_j||\_2](https://latex.codecogs.com/png.latex?Q%3D%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%5Cfrac%7B1%7D%7Bn%7D%5Crho_%7B%5Ctau_%5Cell%7D%5Cbigg%28Y_i%5E%7B%28g%29%7D-%5Ctextbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Csum_%7Bj%3D0%7D%5Epx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%5Cbigg%29%2B%5Clambda_1%7C%7C%5Ctextbf%7BXA%7D%7C%7C_%7B%2A%5Cpsi%7D%2B%5Clambda_2%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7Bj%3D1%7D%5Ep%5Cphi_j%5E%7B%28g%29%7D%7C%7C%5Cboldsymbol%7B%5Ctheta%7D%5E%7B%28g%29%7D_j%7C%7C_2
"Q=\\sum_{g=1}^m\\sum_{\\ell=1}^b\\sum_{i=1}^n\\frac{1}{n}\\rho_{\\tau_\\ell}\\bigg(Y_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum_{j=0}^px_{ij}\\Phi(\\tau_\\ell)^T\\boldsymbol{\\theta}_j^{(g)}\\bigg)+\\lambda_1||\\textbf{XA}||_{*\\psi}+\\lambda_2\\sum_{g=1}^m\\sum_{j=1}^p\\phi_j^{(g)}||\\boldsymbol{\\theta}^{(g)}_j||_2")  

  
![\\frac{\\partial
Q}{\\partial\\boldsymbol{\\theta}\_j^{(g)}}=\\sum\_{\\ell=1}^b\\sum\_{i=1}^n-x\_{ij}\\Phi(\\tau\_\\ell)\\frac{t\_i^{(\\ell)(g)}}{n}+\\lambda\_2\\phi\_j^{(g)}s\_j^{(g)}=0\\;\\text{
by kkt
condition}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%7D%3D%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En-x_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5Cfrac%7Bt_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7Bn%7D%2B%5Clambda_2%5Cphi_j%5E%7B%28g%29%7Ds_j%5E%7B%28g%29%7D%3D0%5C%3B%5Ctext%7B%20by%20kkt%20condition%7D
"\\frac{\\partial Q}{\\partial\\boldsymbol{\\theta}_j^{(g)}}=\\sum_{\\ell=1}^b\\sum_{i=1}^n-x_{ij}\\Phi(\\tau_\\ell)\\frac{t_i^{(\\ell)(g)}}{n}+\\lambda_2\\phi_j^{(g)}s_j^{(g)}=0\\;\\text{ by kkt condition}")  
where

  
![\\begin{aligned}t\_i^{(\\ell)(g)}&=\\frac{\\partial
\\rho\_{\\tau\_{\\ell}}(e\_i^{(\\ell)(g)})}{\\partial
e\_i^{(\\ell)(g)}}\\\\&=\\begin{cases}\\tau\_\\ell-1\\quad&\\text{if
}e\_i^{(\\ell)(g)}\<0\\\\\\{c\\in\\mathbb{R}:\\tau\_\\ell-1\\le
c\\le\\tau\_\\ell\\}\\quad&\\text{if
}e\_i^{(\\ell)(g)}=0\\\\\\tau\_\\ell\\quad&\\text{if
}e\_i^{(\\ell)(g)}\>0\\end{cases}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7Dt_i%5E%7B%28%5Cell%29%28g%29%7D%26%3D%5Cfrac%7B%5Cpartial%20%5Crho_%7B%5Ctau_%7B%5Cell%7D%7D%28e_i%5E%7B%28%5Cell%29%28g%29%7D%29%7D%7B%5Cpartial%20e_i%5E%7B%28%5Cell%29%28g%29%7D%7D%5C%5C%26%3D%5Cbegin%7Bcases%7D%5Ctau_%5Cell-1%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3C0%5C%5C%5C%7Bc%5Cin%5Cmathbb%7BR%7D%3A%5Ctau_%5Cell-1%5Cle%20c%5Cle%5Ctau_%5Cell%5C%7D%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3D0%5C%5C%5Ctau_%5Cell%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3E0%5Cend%7Bcases%7D%5Cend%7Baligned%7D
"\\begin{aligned}t_i^{(\\ell)(g)}&=\\frac{\\partial \\rho_{\\tau_{\\ell}}(e_i^{(\\ell)(g)})}{\\partial e_i^{(\\ell)(g)}}\\\\&=\\begin{cases}\\tau_\\ell-1\\quad&\\text{if }e_i^{(\\ell)(g)}\<0\\\\\\{c\\in\\mathbb{R}:\\tau_\\ell-1\\le c\\le\\tau_\\ell\\}\\quad&\\text{if }e_i^{(\\ell)(g)}=0\\\\\\tau_\\ell\\quad&\\text{if }e_i^{(\\ell)(g)}\>0\\end{cases}\\end{aligned}")  

where
![e\_i^{(\\ell)(g)}=Y\_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum\_{j=0}^px\_{ij}\\Phi(\\tau\_\\ell)^T\\boldsymbol{\\theta}\_j^{(g)}](https://latex.codecogs.com/png.latex?e_i%5E%7B%28%5Cell%29%28g%29%7D%3DY_i%5E%7B%28g%29%7D-%5Ctextbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Csum_%7Bj%3D0%7D%5Epx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D
"e_i^{(\\ell)(g)}=Y_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum_{j=0}^px_{ij}\\Phi(\\tau_\\ell)^T\\boldsymbol{\\theta}_j^{(g)}")

  
![\\rightarrow
|t\_i^{(\\ell)(g)}|\\le\\text{max}(\\tau\_\\ell, 1-\\tau\_\\ell)](https://latex.codecogs.com/png.latex?%5Crightarrow%20%7Ct_i%5E%7B%28%5Cell%29%28g%29%7D%7C%5Cle%5Ctext%7Bmax%7D%28%5Ctau_%5Cell%2C%201-%5Ctau_%5Cell%29
"\\rightarrow |t_i^{(\\ell)(g)}|\\le\\text{max}(\\tau_\\ell, 1-\\tau_\\ell)")  

  
![\\begin{aligned}s\_j^{(g)}&=\\frac{\\partial||\\boldsymbol{\\theta}\_j^{(g)}||\_2}{\\partial\\boldsymbol{\\theta}\_j^{(g)}}\\\\&=\\begin{cases}\\frac{1}{||\\theta\_j^{(g)}||\_2}\\theta\_j^{(g)}&\\text{if
}\\boldsymbol{\\theta}\_j^{(g)}\\ne0\\\\k:||k||\_2\\le1&\\text{if
}\\boldsymbol{\\theta}\_j^{(g)}=0\\end{cases}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7Ds_j%5E%7B%28g%29%7D%26%3D%5Cfrac%7B%5Cpartial%7C%7C%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%7C%7C_2%7D%7B%5Cpartial%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%7D%5C%5C%26%3D%5Cbegin%7Bcases%7D%5Cfrac%7B1%7D%7B%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%7D%5Ctheta_j%5E%7B%28g%29%7D%26%5Ctext%7Bif%20%7D%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%5Cne0%5C%5Ck%3A%7C%7Ck%7C%7C_2%5Cle1%26%5Ctext%7Bif%20%7D%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%3D0%5Cend%7Bcases%7D%5Cend%7Baligned%7D
"\\begin{aligned}s_j^{(g)}&=\\frac{\\partial||\\boldsymbol{\\theta}_j^{(g)}||_2}{\\partial\\boldsymbol{\\theta}_j^{(g)}}\\\\&=\\begin{cases}\\frac{1}{||\\theta_j^{(g)}||_2}\\theta_j^{(g)}&\\text{if }\\boldsymbol{\\theta}_j^{(g)}\\ne0\\\\k:||k||_2\\le1&\\text{if }\\boldsymbol{\\theta}_j^{(g)}=0\\end{cases}\\end{aligned}")  

  
![\\begin{aligned}\\boldsymbol{\\theta}\_j^{(g)}=0&\\iff\\lambda\_2\\phi\_j^{(g)}k=\\sum\_{\\ell=1}^b\\sum\_{i=1}^nx\_{ij}\\Phi(\\tau\_\\ell)\\frac{t\_i^{(\\ell)(g)}}{n}\\\\&\\iff\\lambda\_2|\\phi\_j^{(g)}|\\cdot||k||\_2=||\\sum\_{\\ell=1}^b\\sum\_{i=1}^nx\_{ij}\\Phi(\\tau\_\\ell)\\frac{t\_i^{(\\ell)(g)}}{n}||\_2\\le\\lambda\_2|\\phi\_j^{(g)}|\\\\&\\iff||\\sum\_{\\ell=1}^b\\sum\_{i=1}^nx\_{ij}\\Phi(\\tau\_\\ell)\\frac{t\_i^{(\\ell)(g)}}{n|\\phi\_j^{(g)}|}||\\le\\lambda\_2\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%3D0%26%5Ciff%5Clambda_2%5Cphi_j%5E%7B%28g%29%7Dk%3D%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5Enx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5Cfrac%7Bt_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7Bn%7D%5C%5C%26%5Ciff%5Clambda_2%7C%5Cphi_j%5E%7B%28g%29%7D%7C%5Ccdot%7C%7Ck%7C%7C_2%3D%7C%7C%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5Enx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5Cfrac%7Bt_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7Bn%7D%7C%7C_2%5Cle%5Clambda_2%7C%5Cphi_j%5E%7B%28g%29%7D%7C%5C%5C%26%5Ciff%7C%7C%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5Enx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5Cfrac%7Bt_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7Bn%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%7C%7C%5Cle%5Clambda_2%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{\\theta}_j^{(g)}=0&\\iff\\lambda_2\\phi_j^{(g)}k=\\sum_{\\ell=1}^b\\sum_{i=1}^nx_{ij}\\Phi(\\tau_\\ell)\\frac{t_i^{(\\ell)(g)}}{n}\\\\&\\iff\\lambda_2|\\phi_j^{(g)}|\\cdot||k||_2=||\\sum_{\\ell=1}^b\\sum_{i=1}^nx_{ij}\\Phi(\\tau_\\ell)\\frac{t_i^{(\\ell)(g)}}{n}||_2\\le\\lambda_2|\\phi_j^{(g)}|\\\\&\\iff||\\sum_{\\ell=1}^b\\sum_{i=1}^nx_{ij}\\Phi(\\tau_\\ell)\\frac{t_i^{(\\ell)(g)}}{n|\\phi_j^{(g)}|}||\\le\\lambda_2\\end{aligned}")  

  
![\\begin{aligned}||\\sum\_{\\ell=1}^b\\sum\_{i=1}^nx\_{ij}\\Phi(\\tau\_\\ell)\\frac{t\_i^{(\\ell)(g)}}{n|\\phi\_j^{(g)}|}||\_2&\\le\\frac{1}{n|\\phi\_j^{(g)}|}||\\sum\_{\\ell=1}^b\\sum\_{i=1}^nx\_{ij}\\Phi(\\tau\_\\ell)||\_2\\\\&\\qquad(\\because
t\_i^{(\\ell)(g)}\\le\\text{max}(\\tau\_\\ell-1,\\tau\_\\ell)\\le1)\\\\&\\le\\frac{1}{n|\\phi\_j^{(g)}|}\\sum\_{\\ell=1}^b\\sum\_{i=1}^n|x\_{ij}|\\cdot||\\Phi(\\tau\_\\ell)||\_2\\\\&\\qquad(\\because\\text{
triangle
inequality)}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%7C%7C%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5Enx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5Cfrac%7Bt_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7Bn%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%7C%7C_2%26%5Cle%5Cfrac%7B1%7D%7Bn%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%7C%7C%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5Enx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%7C%7C_2%5C%5C%26%5Cqquad%28%5Cbecause%20t_i%5E%7B%28%5Cell%29%28g%29%7D%5Cle%5Ctext%7Bmax%7D%28%5Ctau_%5Cell-1%2C%5Ctau_%5Cell%29%5Cle1%29%5C%5C%26%5Cle%5Cfrac%7B1%7D%7Bn%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%7Cx_%7Bij%7D%7C%5Ccdot%7C%7C%5CPhi%28%5Ctau_%5Cell%29%7C%7C_2%5C%5C%26%5Cqquad%28%5Cbecause%5Ctext%7B%20triangle%20inequality%29%7D%5Cend%7Baligned%7D
"\\begin{aligned}||\\sum_{\\ell=1}^b\\sum_{i=1}^nx_{ij}\\Phi(\\tau_\\ell)\\frac{t_i^{(\\ell)(g)}}{n|\\phi_j^{(g)}|}||_2&\\le\\frac{1}{n|\\phi_j^{(g)}|}||\\sum_{\\ell=1}^b\\sum_{i=1}^nx_{ij}\\Phi(\\tau_\\ell)||_2\\\\&\\qquad(\\because t_i^{(\\ell)(g)}\\le\\text{max}(\\tau_\\ell-1,\\tau_\\ell)\\le1)\\\\&\\le\\frac{1}{n|\\phi_j^{(g)}|}\\sum_{\\ell=1}^b\\sum_{i=1}^n|x_{ij}|\\cdot||\\Phi(\\tau_\\ell)||_2\\\\&\\qquad(\\because\\text{ triangle inequality)}\\end{aligned}")  

  
![\\sum\_{g=1}^m\\sum\_{\\ell=1}^b\\sum\_{i=1}^n\\frac{1}{n}\\rho\_{\\tau\_\\ell}\\bigg(Y\_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum\_{j=0}^px\_{ij}\\Phi(\\tau\_\\ell)^T\\boldsymbol{\\theta}\_j^{(g)}\\bigg)=\\sum\_{\\ell=1}^b\\frac{1}{n}||\\sqrt{\\rho\_{\\tau\_\\ell}(\\textbf{E})}||\_F^2\\\\\\text{where
}\\sqrt{\\rho\_{\\tau\_\\ell}(\\textbf{E})}=\\bigg\\{\\sqrt{\\rho\_{\\tau\_\\ell}(Y\_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum\_{j=0}^px\_{ij}\\Phi(\\tau\_\\ell)^T\\boldsymbol{\\theta}\_j^{(g)})}\\bigg\\}\_{i,g}](https://latex.codecogs.com/png.latex?%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%5Cfrac%7B1%7D%7Bn%7D%5Crho_%7B%5Ctau_%5Cell%7D%5Cbigg%28Y_i%5E%7B%28g%29%7D-%5Ctextbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Csum_%7Bj%3D0%7D%5Epx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%5Cbigg%29%3D%5Csum_%7B%5Cell%3D1%7D%5Eb%5Cfrac%7B1%7D%7Bn%7D%7C%7C%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28%5Ctextbf%7BE%7D%29%7D%7C%7C_F%5E2%5C%5C%5Ctext%7Bwhere%20%7D%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28%5Ctextbf%7BE%7D%29%7D%3D%5Cbigg%5C%7B%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28Y_i%5E%7B%28g%29%7D-%5Ctextbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Csum_%7Bj%3D0%7D%5Epx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%29%7D%5Cbigg%5C%7D_%7Bi%2Cg%7D
"\\sum_{g=1}^m\\sum_{\\ell=1}^b\\sum_{i=1}^n\\frac{1}{n}\\rho_{\\tau_\\ell}\\bigg(Y_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum_{j=0}^px_{ij}\\Phi(\\tau_\\ell)^T\\boldsymbol{\\theta}_j^{(g)}\\bigg)=\\sum_{\\ell=1}^b\\frac{1}{n}||\\sqrt{\\rho_{\\tau_\\ell}(\\textbf{E})}||_F^2\\\\\\text{where }\\sqrt{\\rho_{\\tau_\\ell}(\\textbf{E})}=\\bigg\\{\\sqrt{\\rho_{\\tau_\\ell}(Y_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum_{j=0}^px_{ij}\\Phi(\\tau_\\ell)^T\\boldsymbol{\\theta}_j^{(g)})}\\bigg\\}_{i,g}")  

  
![\\frac{\\partial Q}{\\partial
\\textbf{Z}}=-\\frac{1}{n}\\sum\_{\\ell=1}^n\\frac{\\partial\\,||\\sqrt{\\rho\_{\\tau\_\\ell}(\\textbf{E})}||\_F^2}{\\partial\\textbf{E}}+\\lambda\_1\\frac{\\partial||\\textbf{Z}||\_{\*,\\psi}}{\\partial\\,\\textbf{Z}}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20%5Ctextbf%7BZ%7D%7D%3D-%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7B%5Cell%3D1%7D%5En%5Cfrac%7B%5Cpartial%5C%2C%7C%7C%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28%5Ctextbf%7BE%7D%29%7D%7C%7C_F%5E2%7D%7B%5Cpartial%5Ctextbf%7BE%7D%7D%2B%5Clambda_1%5Cfrac%7B%5Cpartial%7C%7C%5Ctextbf%7BZ%7D%7C%7C_%7B%2A%2C%5Cpsi%7D%7D%7B%5Cpartial%5C%2C%5Ctextbf%7BZ%7D%7D
"\\frac{\\partial Q}{\\partial \\textbf{Z}}=-\\frac{1}{n}\\sum_{\\ell=1}^n\\frac{\\partial\\,||\\sqrt{\\rho_{\\tau_\\ell}(\\textbf{E})}||_F^2}{\\partial\\textbf{E}}+\\lambda_1\\frac{\\partial||\\textbf{Z}||_{*,\\psi}}{\\partial\\,\\textbf{Z}}")
