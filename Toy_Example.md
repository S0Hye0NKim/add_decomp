Toy Example
================
Sohyeon Kim



``` r
library(tidyverse)
library(ggplot2)
library(splines)
library(Matrix)
library(pander)
```

# 0\. Initial Setting

  
![\\begin{aligned}Y &\\in \\mathbb{R}^{n\\times g}\\\\Y &= XB(\\tau) +
\\epsilon(\\tau) , \\; X\\sim N(0, 1)\\\\ X &\\in\\mathbb{R}^{n \\times
(p+1)}, \\; \\;B(\\tau) \\in \\mathbb{R}^{(p+1)\\times g} , \\; and \\;
Q\_\\tau(\\epsilon(\\tau)|X)
= 0\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7DY%20%26%5Cin%20%5Cmathbb%7BR%7D%5E%7Bn%5Ctimes%20g%7D%5C%5CY%20%26%3D%20XB%28%5Ctau%29%20%2B%20%5Cepsilon%28%5Ctau%29%20%2C%20%5C%3B%20X%5Csim%20N%280%2C%201%29%5C%5C%20X%20%26%5Cin%5Cmathbb%7BR%7D%5E%7Bn%20%5Ctimes%20%28p%2B1%29%7D%2C%20%5C%3B%20%5C%3BB%28%5Ctau%29%20%5Cin%20%5Cmathbb%7BR%7D%5E%7B%28p%2B1%29%5Ctimes%20g%7D%20%2C%20%5C%3B%20and%20%5C%3B%20Q_%5Ctau%28%5Cepsilon%28%5Ctau%29%7CX%29%20%3D%200%5Cend%7Baligned%7D
"\\begin{aligned}Y &\\in \\mathbb{R}^{n\\times g}\\\\Y &= XB(\\tau) + \\epsilon(\\tau) , \\; X\\sim N(0, 1)\\\\ X &\\in\\mathbb{R}^{n \\times (p+1)}, \\; \\;B(\\tau) \\in \\mathbb{R}^{(p+1)\\times g} , \\; and \\; Q_\\tau(\\epsilon(\\tau)|X) = 0\\end{aligned}")  

  
![\\begin{aligned}\\beta\_0(\\tau) &= \\beta\_0 + Q\_\\tau(\\epsilon)
\\\\ \\epsilon(\\tau) &= \\epsilon- Q\_\\tau(\\epsilon), \\; \\epsilon
\\sim
N(0, 0.1^2)\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cbeta_0%28%5Ctau%29%20%26%3D%20%5Cbeta_0%20%2B%20Q_%5Ctau%28%5Cepsilon%29%20%5C%5C%20%5Cepsilon%28%5Ctau%29%20%26%3D%20%5Cepsilon-%20Q_%5Ctau%28%5Cepsilon%29%2C%20%5C%3B%20%5Cepsilon%20%5Csim%20N%280%2C%200.1%5E2%29%5Cend%7Baligned%7D
"\\begin{aligned}\\beta_0(\\tau) &= \\beta_0 + Q_\\tau(\\epsilon) \\\\ \\epsilon(\\tau) &= \\epsilon- Q_\\tau(\\epsilon), \\; \\epsilon \\sim N(0, 0.1^2)\\end{aligned}")  

  - n : \# of observation.
  - m : \# of group.
  - p : \# of covariate.
  - r : nonzero entry in true
    ![\\Gamma(\\tau)](https://latex.codecogs.com/png.latex?%5CGamma%28%5Ctau%29
    "\\Gamma(\\tau)")
  - K : \# of basis function.
  - b : \# of regional quantiel we consider.

For toy example, I will use n = 300, p = 50, g = 20.

``` r
set.seed(0)
n <- 300
m <- 15
p <- 20
num_nz <- round((p+1)*m*(2/3)) 
b <- 6
num_rank <- 5

X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) %>% cbind(1, .)   #add intercept term in X

row_idx <- sample(1:(p+1), size = num_nz, replace = TRUE)     # make sparse matrix following N(4, 1)
col_idx <- sample(1:m, size = num_nz, replace = TRUE)
rand.x <- rnorm(num_nz, mean = 4, sd = 1)
sp_mat <- sparseMatrix(row_idx, col_idx, x = rand.x)
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

![B(\\tau)](https://latex.codecogs.com/png.latex?B%28%5Ctau%29
"B(\\tau)") = sparse matrix + low rank matrix

# 1\. Preliminary

### 1-1. Check function

``` r
check_ft <- function(x, tau) {
  z <- ifelse(x<0, (tau-1)*x, tau*x)
  return(z)
}
```

  
![\\rho\_\\tau(x) = (\\tau -
\\mathbb{I}(x\<0))x](https://latex.codecogs.com/png.latex?%5Crho_%5Ctau%28x%29%20%3D%20%28%5Ctau%20-%20%5Cmathbb%7BI%7D%28x%3C0%29%29x
"\\rho_\\tau(x) = (\\tau - \\mathbb{I}(x\<0))x")  

### 1-2. Basis function

``` r
K <- 10

tau_seq <- seq(from = 0.4, to = 0.6, length.out = b)
knots_seq <- seq(from = tau_seq[1] - 0.02, to = tau_seq[b] + 0.02, length.out = K - 2) %>% .[-c(1, (K - 2))]
boundary <- c(tau_seq[1] - 0.02, tau_seq[b] + 0.02)
Phi <- bs(tau_seq, knots = knots_seq, degree = 3, intercept = TRUE, Boundary.knots = boundary)

attr(Phi, "knots")
```

    ## [1] 0.4142857 0.4485714 0.4828571 0.5171429 0.5514286 0.5857143

  
![\\begin{aligned}\\boldsymbol{v}\_i^{(\\ell)} &= (\\boldsymbol{x}\_i
\\otimes \\Phi(\\tau\_{\\ell}))^T\\\\ &=(\\Phi(\\tau\_\\ell)^T,
x\_{i1}\\Phi(\\tau\_\\ell)^T, \\dots, x\_{ip}\\Phi(\\tau\_\\ell)^T) \\\\
&\\in
\\mathbb{R}^{(p+1)K}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7Bv%7D_i%5E%7B%28%5Cell%29%7D%20%26%3D%20%28%5Cboldsymbol%7Bx%7D_i%20%5Cotimes%20%5CPhi%28%5Ctau_%7B%5Cell%7D%29%29%5ET%5C%5C%20%26%3D%28%5CPhi%28%5Ctau_%5Cell%29%5ET%2C%20x_%7Bi1%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%2C%20%5Cdots%2C%20x_%7Bip%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%29%20%5C%5C%20%26%5Cin%20%5Cmathbb%7BR%7D%5E%7B%28p%2B1%29K%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{v}_i^{(\\ell)} &= (\\boldsymbol{x}_i \\otimes \\Phi(\\tau_{\\ell}))^T\\\\ &=(\\Phi(\\tau_\\ell)^T, x_{i1}\\Phi(\\tau_\\ell)^T, \\dots, x_{ip}\\Phi(\\tau_\\ell)^T) \\\\ &\\in \\mathbb{R}^{(p+1)K}\\end{aligned}")  

### 1-3. New design matrix for regional quantile

``` r
V <- list()
for(l in 1:b){
  V_l <- matrix(nrow = n, ncol = (p+1)*K)
  for(i in 1:n) {
    v_i_l <- kronecker(X[i, ], Phi[l, ])
    V_l[i, ] <- v_i_l 
  }
  V[[l]] <- V_l
}
```

  - Dimension for V : ![b \\times n \\times
    (p+1)K](https://latex.codecogs.com/png.latex?b%20%5Ctimes%20n%20%5Ctimes%20%28p%2B1%29K
    "b \\times n \\times (p+1)K")

  
![\\boldsymbol{V}^{(\\ell)} = \\begin{bmatrix} v\_1^{(\\ell)}
\\\\v\_2^{(\\ell)} \\\\ \\vdots \\\\ v\_n^{(\\ell)}\\end{bmatrix} \\in
\\mathbb{R}^{n\\times
(p+1)K}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%20%3D%20%5Cbegin%7Bbmatrix%7D%20v_1%5E%7B%28%5Cell%29%7D%20%5C%5Cv_2%5E%7B%28%5Cell%29%7D%20%5C%5C%20%5Cvdots%20%5C%5C%20v_n%5E%7B%28%5Cell%29%7D%5Cend%7Bbmatrix%7D%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bn%5Ctimes%20%28p%2B1%29K%7D
"\\boldsymbol{V}^{(\\ell)} = \\begin{bmatrix} v_1^{(\\ell)} \\\\v_2^{(\\ell)} \\\\ \\vdots \\\\ v_n^{(\\ell)}\\end{bmatrix} \\in \\mathbb{R}^{n\\times (p+1)K}")  

### 1-4. Initial value for estimator

  
![\\boldsymbol{\\eta}^{(g)}, \\boldsymbol{\\theta}^{(g)},
\\boldsymbol{w}^{(g)}\\in \\mathbb{R}^{(p+1)\\times K}\\\\
\\boldsymbol{\\alpha}^{(g)} \\in \\mathbb{R}^{p+1} \\\\
\\boldsymbol{e}^{(\\ell)(g)}, \\boldsymbol{u}^{(\\ell)(g)}, \\in
\\mathbb{R}^{n}\\\\subject\\; to \\; \\boldsymbol{\\eta}^{(g)} -
\\boldsymbol{\\theta}^{(g)} = 0, \\; and \\;\\boldsymbol{Y}^{(g)} -
\\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)} -
\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)} -
\\boldsymbol{e}^{(\\ell)(g)}
= 0](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%7D%2C%20%5Cboldsymbol%7B%5Ctheta%7D%5E%7B%28g%29%7D%2C%20%5Cboldsymbol%7Bw%7D%5E%7B%28g%29%7D%5Cin%20%5Cmathbb%7BR%7D%5E%7B%28p%2B1%29%5Ctimes%20K%7D%5C%5C%20%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bp%2B1%7D%20%5C%5C%20%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%28g%29%7D%2C%20%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28g%29%7D%2C%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bn%7D%5C%5Csubject%5C%3B%20to%20%5C%3B%20%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%7D%20-%20%5Cboldsymbol%7B%5Ctheta%7D%5E%7B%28g%29%7D%20%3D%200%2C%20%5C%3B%20and%20%5C%3B%5Cboldsymbol%7BY%7D%5E%7B%28g%29%7D%20-%20%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D%20-%20%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%7D%20-%20%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%28g%29%7D%20%3D%200
"\\boldsymbol{\\eta}^{(g)}, \\boldsymbol{\\theta}^{(g)}, \\boldsymbol{w}^{(g)}\\in \\mathbb{R}^{(p+1)\\times K}\\\\ \\boldsymbol{\\alpha}^{(g)} \\in \\mathbb{R}^{p+1} \\\\ \\boldsymbol{e}^{(\\ell)(g)}, \\boldsymbol{u}^{(\\ell)(g)}, \\in \\mathbb{R}^{n}\\\\subject\\; to \\; \\boldsymbol{\\eta}^{(g)} - \\boldsymbol{\\theta}^{(g)} = 0, \\; and \\;\\boldsymbol{Y}^{(g)} - \\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)} - \\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)} - \\boldsymbol{e}^{(\\ell)(g)} = 0")  

# 2\. Algorithm

  
![\\begin{aligned}\\boldsymbol{Y}&=\\boldsymbol{XA}+\\boldsymbol{X\\Gamma}(\\tau)\\\\&=\\boldsymbol{Z}+\\boldsymbol{V}\_{\\tau\_\\ell}\\boldsymbol{\\Theta}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7BY%7D%26%3D%5Cboldsymbol%7BXA%7D%2B%5Cboldsymbol%7BX%5CGamma%7D%28%5Ctau%29%5C%5C%26%3D%5Cboldsymbol%7BZ%7D%2B%5Cboldsymbol%7BV%7D_%7B%5Ctau_%5Cell%7D%5Cboldsymbol%7B%5CTheta%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{Y}&=\\boldsymbol{XA}+\\boldsymbol{X\\Gamma}(\\tau)\\\\&=\\boldsymbol{Z}+\\boldsymbol{V}_{\\tau_\\ell}\\boldsymbol{\\Theta}\\end{aligned}")  

``` r
add_decomp <- function(delta, lambda_1, lambda_2, tol_error, max_iter) {
  # delta = step size
  # lambda_1 = low rank penalty
  # lambda_2 = sparse penalty
  
  # initial value
  eta_old <- matrix(1, nrow = (p+1)*K, ncol = m) 
  theta_old <- eta_old
  alpha_old <- matrix(0, nrow = p+1, ncol = m)
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
```

  
![\\hat{\\Gamma}\_{11}(\\tau\_1)=\\cdots=\\hat{\\Gamma}\_{11}(\\tau\_\\ell)=0=\\Gamma\_{11}(\\tau)](https://latex.codecogs.com/png.latex?%5Chat%7B%5CGamma%7D_%7B11%7D%28%5Ctau_1%29%3D%5Ccdots%3D%5Chat%7B%5CGamma%7D_%7B11%7D%28%5Ctau_%5Cell%29%3D0%3D%5CGamma_%7B11%7D%28%5Ctau%29
"\\hat{\\Gamma}_{11}(\\tau_1)=\\cdots=\\hat{\\Gamma}_{11}(\\tau_\\ell)=0=\\Gamma_{11}(\\tau)")  

``` r
check_sp_table <- function(true, est, tol = 0.1^5, table = FALSE, nnz = num_nz, nz = num_zero) {
  # check sparsity pattern of true and est matrix
  zero_idx_true <- which(abs(true) < tol, arr.ind = TRUE) %>% as_tibble
  zero_idx_est <- lapply(est, FUN = function(x) which(abs(x) < tol, arr.ind = TRUE) %>% as_tibble) %>%
    bind_rows %>%
    split(., .$row) %>%
    lapply(FUN = function(x) x$col %>% unique %>% tibble(col = .)) %>%
    bind_rows(.id = "row") %>% mutate_all(as.numeric)
  
  if(nrow(zero_idx_est) == 0) {
    result <- data.frame(Positive = c(nnz, nz), Negative = c(0, 0))
  } else {
    TN <- semi_join(zero_idx_true, zero_idx_est, by = c("row", "col"))
    FP <- anti_join(zero_idx_true, zero_idx_est, by = c("row", "col"))
    FN <- anti_join(zero_idx_est, zero_idx_true, by = c("row", "col"))
    result <- data.frame(Positive = c(nnz - nrow(FN), nrow(FP)), Negative = c(nrow(FN), nrow(TN))) %>%
      `rownames<-`(value = c("Positive", "Negative"))
  }
  
  if(table == FALSE) {return(data.frame(FPR = result[2, 1]/nz, TPR = result[1, 1]/nnz))
    } else {return(result)}
}
```

``` r
max_iter <- 50
delta <- 0.5
lambda_1 <- 0.15
lambda_2 <- 2.5
tol_error <- 0.1

result <- add_decomp(delta, lambda_1, lambda_2, tol_error, max_iter)
```

![\\lambda\_1, \\lambda\_2,
\\delta](https://latex.codecogs.com/png.latex?%5Clambda_1%2C%20%5Clambda_2%2C%20%5Cdelta
"\\lambda_1, \\lambda_2, \\delta") 어떻게?

``` r
result$iter_error %>% 
  na.omit() %>%
  data.frame %>%
  mutate(iter = 1:nrow(.)) %>%
  filter(iter > 15) %>% 
  gather(key = "estimator", value = "value", -iter) %>%
  ggplot() +
  geom_line(aes(x = iter, y = value, group = estimator, color = estimator)) +
  facet_wrap(~estimator, scales = "free_y")
```

![](Toy_Example_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

  - Stopping criteria

  
![\\begin{aligned}&||\\eta^{k+1}-\\eta^{k}||\_2^2+||\\theta^{k+1}-\\theta^k||\_2^2+||\\alpha^{k+1}-\\alpha^k||\_2^2\\\\&+||e^{k+1}-e^k||\_2^2+||u^{k+1}-u^k||\_2^2+||w^{k+1}-w^k||\_2^2
\\le
\\epsilon\_{tol}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%7C%7C%5Ceta%5E%7Bk%2B1%7D-%5Ceta%5E%7Bk%7D%7C%7C_2%5E2%2B%7C%7C%5Ctheta%5E%7Bk%2B1%7D-%5Ctheta%5Ek%7C%7C_2%5E2%2B%7C%7C%5Calpha%5E%7Bk%2B1%7D-%5Calpha%5Ek%7C%7C_2%5E2%5C%5C%26%2B%7C%7Ce%5E%7Bk%2B1%7D-e%5Ek%7C%7C_2%5E2%2B%7C%7Cu%5E%7Bk%2B1%7D-u%5Ek%7C%7C_2%5E2%2B%7C%7Cw%5E%7Bk%2B1%7D-w%5Ek%7C%7C_2%5E2%20%5Cle%20%5Cepsilon_%7Btol%7D%5Cend%7Baligned%7D
"\\begin{aligned}&||\\eta^{k+1}-\\eta^{k}||_2^2+||\\theta^{k+1}-\\theta^k||_2^2+||\\alpha^{k+1}-\\alpha^k||_2^2\\\\&+||e^{k+1}-e^k||_2^2+||u^{k+1}-u^k||_2^2+||w^{k+1}-w^k||_2^2 \\le \\epsilon_{tol}\\end{aligned}")  

# 3\. Evaluation

![\\Gamma(\\tau)](https://latex.codecogs.com/png.latex?%5CGamma%28%5Ctau%29
"\\Gamma(\\tau)")는 ![\\tau](https://latex.codecogs.com/png.latex?%5Ctau
"\\tau")에 상관 없이 동일하므로, quantile 상관 없이
![Y\\approx\\hat{Y}\_{\\tau\_\\ell}](https://latex.codecogs.com/png.latex?Y%5Capprox%5Chat%7BY%7D_%7B%5Ctau_%5Cell%7D
"Y\\approx\\hat{Y}_{\\tau_\\ell}") 동일.

``` r
Y_hat_tau <- lapply(V, FUN = function(x) x %*% result$theta + X %*% result$alpha) %>%
  `names<-`(value = tau_seq)
```

### Low rank matrix

  
![\\frac{\\sigma\_1+\\cdots+\\sigma\_5}{\\sigma\_1+\\cdots+\\sigma\_{20}}\\quad\\text{where
}\\sigma\_i\\text{ is a singular
value}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Csigma_1%2B%5Ccdots%2B%5Csigma_5%7D%7B%5Csigma_1%2B%5Ccdots%2B%5Csigma_%7B20%7D%7D%5Cquad%5Ctext%7Bwhere%20%7D%5Csigma_i%5Ctext%7B%20is%20a%20singular%20value%7D
"\\frac{\\sigma_1+\\cdots+\\sigma_5}{\\sigma_1+\\cdots+\\sigma_{20}}\\quad\\text{where }\\sigma_i\\text{ is a singular value}")  

``` r
sing_val <- svd(result$alpha)$d
sum(sing_val[1:5])/sum(sing_val)
```

    ## [1] 0.9859056

``` r
rankMatrix(result$alpha)
```

    ## [1] 11
    ## attr(,"method")
    ## [1] "tolNorm2"
    ## attr(,"useGrad")
    ## [1] FALSE
    ## attr(,"tol")
    ## [1] 4.662937e-15

  
![\\frac{||\\hat{A}-A||\_F}{||A||\_F}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%7C%7C%5Chat%7BA%7D-A%7C%7C_F%7D%7B%7C%7CA%7C%7C_F%7D
"\\frac{||\\hat{A}-A||_F}{||A||_F}")  

``` r
norm(result$alpha - LR_mat, "F")/norm(LR_mat, "F")
```

    ## [1] 0.2098985

### Sparse matrix

![\\hat{\\gamma}\_j^{(g)}(\\tau\_\\ell)=\\hat{\\theta}\_{j}^{(g)T}\\phi\_s(\\tau)](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cgamma%7D_j%5E%7B%28g%29%7D%28%5Ctau_%5Cell%29%3D%5Chat%7B%5Ctheta%7D_%7Bj%7D%5E%7B%28g%29T%7D%5Cphi_s%28%5Ctau%29
"\\hat{\\gamma}_j^{(g)}(\\tau_\\ell)=\\hat{\\theta}_{j}^{(g)T}\\phi_s(\\tau)")

``` r
gamma_tau_hat <- list()
for(l in 1:b) {
  phi_tau <- Phi[l, ]
  gamma_tau_hat[[l]] <- matrix(nrow = (p+1), ncol = m)
  for(i in 1:(p+1)) {
    for(j in 1:m) {
      theta_j <- result$theta[(1+(i-1)*K):(K*i), j]
      gamma_tau_hat[[l]][i, j] <- theta_j %*% phi_tau
    }
  }
}
```

### sparsity pattern check

``` r
data.frame(c("TP", "FP"), c("FN", "TN")) %>% 
  `rownames<-`(value = c("True:non_zero", "True:zero")) %>%
  `colnames<-`(value = c("Est:non_zero", "Est:zero")) %>% pander
```

|                    | Est:non\_zero | Est:zero |
| :----------------: | :-----------: | :------: |
| **True:non\_zero** |      TP       |    FN    |
|   **True:zero**    |      FP       |    TN    |

``` r
check_sp_table(true = sp_mat, est = gamma_tau_hat, table = FALSE)
```

    ##        FPR       TPR
    ## 1 0.695122 0.7417219

### Frobenius norm ratio

  
![\\frac{||\\hat{A}-A||\_F}{||A||\_F}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%7C%7C%5Chat%7BA%7D-A%7C%7C_F%7D%7B%7C%7CA%7C%7C_F%7D
"\\frac{||\\hat{A}-A||_F}{||A||_F}")  

``` r
lamb_1_seq <- seq(from = 0.05, to = 0.3, by = 0.05)
simul_lamb1 <- list()

for(idx in 1:length(lamb_1_seq)) {
  simul_lamb1[[idx]] <- add_decomp(delta = 0.5, lambda_1 = lamb_1_seq[idx], lambda_2 = 2.5,
                                   tol_error = 0.1, max_iter = 50)
  if((idx %% 2) == 1) print(paste0("iter = ", idx))
}
```

``` r
alpha_hat <- simul_lamb1 %>% lapply(FUN = function(x) x$alpha)

alpha_hat %>%
  lapply(FUN = function(x) norm(LR_mat - x, type = "F")/norm(LR_mat, type = "F") %>% data.frame(value = .)) %>%
  `names<-`(value = lamb_1_seq) %>%
  bind_rows(.id = "lambda_1") %>%
  ggplot() +
  geom_line(aes(x = as.numeric(lambda_1), y = value)) +
  labs(title = expression("Frobenius norm ratio when "~delta~"=0.1, "~lambda[2]~"=0.5, "~lambda[1]~"in (0.05,0.3)"), 
       x = expression(lambda[1]))
```

``` r
lapply(alpha_hat, FUN = function(x) svd(x) %>% .$d) %>%
  lapply(FUN = function(x) data.frame(ratio = (sum(x[1:num_rank])/sum(x)), rank = sum(abs(x) > 0.1^5))) %>%
  `names<-`(value = lamb_1_seq) %>%
  bind_rows(.id = "lambda1")
```

### ROC curve

``` r
lamb_2_seq <- seq(from = 1, to = 5, by = 0.25)
simul <- list()

for(idx in 1:length(lamb_2_seq)) {
  simul[[idx]] <- add_decomp(delta = 0.5, lambda_1 = 0.2, lambda_2 = lamb_2_seq[idx], tol_error = 0.1, max_iter = 50)
  if((idx %% 2) == 1) print(paste0("iter = ", idx))
}
```

``` r
theta_hat <- simul %>%
  lapply(FUN = function(x) x$theta)

gamma_tau_hat <- list()
for(idx in 1:length(lamb_2_seq)) {
  gamma_tau_hat[[idx]] <- list()
  for(l in 1:b) {
    phi_tau <- Phi[l, ]
    gamma_tau_hat[[idx]][[l]] <- matrix(nrow = (p+1), ncol = m)
    for(i in 1:(p+1)) {
      for(j in 1:m) {
        theta_j_g <- theta_hat[[idx]][(1+(i-1)*K):(K*i), j]
        gamma_tau_hat[[idx]][[l]][i, j] <- theta_j_g %*% phi_tau
      }
    }
  }
}

ROC_table <- list()
for(idx in 1:length(lamb_2_seq)) {
  ROC_table[[idx]] <- check_sp_table(true = sp_mat, est = gamma_tau_hat[[idx]]) 
} 

ROC_table %>%
  `names<-`(value = lamb_2_seq) %>%
  bind_rows(.id = "lamb_2") %>%
  arrange(FPR, TPR) %>%
  ggplot() +
  geom_line(aes(x = FPR, y = TPR)) +
  geom_point(aes(x = FPR, y = TPR)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed(ratio = 1) +
  ggtitle(label = expression(delta ~ " =0.25, "~lambda[1] ~ "=0.2, "~lambda[2] ~ "in (0.5, 1.5)")) +
  scale_y_continuous(limits = c(0, 1))
```
