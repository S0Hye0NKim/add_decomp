Toy Example
================
Sohyeon Kim
7/21/2020



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
    ![B(\\tau)](https://latex.codecogs.com/png.latex?B%28%5Ctau%29
    "B(\\tau)")
  - K : \# of basis function.
  - b : \# of regional quantiel we consider.

For toy example, I will use n = 300, p = 50, g = 20.

``` r
set.seed(0)
n <- 300
m <- 20
p <- 50
r <- 500
b <- 21
num_rank <- 5

X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) %>% cbind(1, .)   #add intercept term in X
sp_mat <- rsparsematrix(nrow = p+1, ncol = m, nnz = r)             # make sparse matrix
LR_mat <- matrix(rnorm(m*(p+1), mean = 0, sd = 1), ncol = m)      # make low rank matrix using SVD
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

#### Question

1.  Y를 generate할 때, 사용하는
    ![B(\\tau)](https://latex.codecogs.com/png.latex?B%28%5Ctau%29
    "B(\\tau)") 값에 tau depend?
2.  우선 ![\\tau
    = 0.5](https://latex.codecogs.com/png.latex?%5Ctau%20%3D%200.5
    "\\tau = 0.5") 로 generate 시켰는데, 이렇게 하는게 맞나?

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
K <- 15

tau_seq <- seq(from = 0.4, to = 0.6, length.out = b) 
Phi <- bs(tau_seq, df = K, degree = 3, intercept = TRUE)
```

  
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

``` r
eta_old <- matrix(0, nrow = (p+1)*K, ncol = m) 
theta_old <- eta_old
alpha_old <- matrix(0, nrow = p+1, ncol = m)
Z_old <- X %*% alpha_old
e_old <- list()
for(l in 1:b) {e_old[[l]] <- Y - Z_old - V[[l]] %*% eta_old}
u_old <- list()
for(l in 1:b) {u_old[[l]] <- matrix(0, nrow = n, ncol = m)}
w_old <- matrix(0, nrow = (p+1)*K, ncol = m)
```

# 2\. Algorithm

  
![\\begin{aligned}\\boldsymbol{Y}&=\\boldsymbol{XA}+\\boldsymbol{X\\Gamma}(\\tau)\\\\&=\\boldsymbol{Z}+\\boldsymbol{V}\_{\\tau\_\\ell}\\boldsymbol{\\Theta}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7BY%7D%26%3D%5Cboldsymbol%7BXA%7D%2B%5Cboldsymbol%7BX%5CGamma%7D%28%5Ctau%29%5C%5C%26%3D%5Cboldsymbol%7BZ%7D%2B%5Cboldsymbol%7BV%7D_%7B%5Ctau_%5Cell%7D%5Cboldsymbol%7B%5CTheta%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{Y}&=\\boldsymbol{XA}+\\boldsymbol{X\\Gamma}(\\tau)\\\\&=\\boldsymbol{Z}+\\boldsymbol{V}_{\\tau_\\ell}\\boldsymbol{\\Theta}\\end{aligned}")  

``` r
max_iter <- 50
delta <- 5
lambda_1 <- 15          #low rank penalty
lambda_2 <- 28           #sparse penalty
tol_error <- 10
iter_error <- matrix(ncol = 6, nrow = max_iter) %>%
  `colnames<-`(value = c("eta", "theta", "alpha", "e", "u", "w"))

sum_V <- Reduce("+", V)
VV_prod <- lapply(V, FUN = function(x) t(x) %*% x)   # V^TV
sum_VV <- Reduce("+", VV_prod)

for(iter in 1:max_iter){
  # Process for eta
  eta_new <- matrix(nrow = (p+1)*K, ncol = m)
  Vu_prod <- mapply(function(x,y) t(x) %*% y, V, u_old, SIMPLIFY = FALSE)
  Ve_prod <- mapply(function(x,y) t(x) %*% y, V, e_old, SIMPLIFY = FALSE)
  eta_new <- (solve(sum_VV+diag(1, (p+1)*K))/delta) %*% (w_old + delta * theta_old + Reduce("+", Vu_prod)
                                                         + delta * t(sum_V) %*% (Y - Z_old)
                                                         - delta * Reduce("+", Ve_prod))
  
  # Process for theta
  theta_new <- matrix(nrow = (p+1)*K, ncol = m)
  threshold <- lambda_2/delta
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
     e_new[[l]][, g] <- case_when(value > tau_seq[l]/delta ~ value - tau_seq[l]/delta, 
                                  value < (tau_seq[l]-1)/delta ~ value - (tau_seq[l]-1)/delta, 
                                  value >=(tau_seq[l]-1)/delta & value <= tau_seq[l]/delta ~ 0)
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
  
  #if(iter_error[iter, ] < tol_error) break
  
  eta_old <- eta_new
  theta_old <- theta_new
  Z_old <- Z_new
  alpha_old <- alpha_new
  e_old <- e_new
  u_old <- u_new
  w_old <- w_new
  
  if(iter%%10 == 0) {print(paste0("iter = ", iter))}
}
```

![\\lambda\_1, \\lambda\_2,
\\delta](https://latex.codecogs.com/png.latex?%5Clambda_1%2C%20%5Clambda_2%2C%20%5Cdelta
"\\lambda_1, \\lambda_2, \\delta") 어떻게?

``` r
iter_error %>% data.frame %>%
  mutate(iter = 1:nrow(.)) %>%
  filter(iter %in% 15:50) %>% 
  gather(key = "estimator", value = "value", -iter) %>%
  ggplot() +
  geom_line(aes(x = iter, y = value, group = estimator, color = estimator)) +
  facet_wrap(~estimator, scales = "free_y")
```

![](Toy_Example_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

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
Y_hat_tau <- lapply(V, FUN = function(x) x %*% theta_new + Z_new) %>%
  `names<-`(value = tau_seq)
```

### Low rank matrix

  
![\\frac{\\sigma\_1+\\cdots+\\sigma\_5}{\\sigma\_1+\\cdots+\\sigma\_{20}}\\quad\\text{where
}\\sigma\_i\\text{ is a singular
value}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Csigma_1%2B%5Ccdots%2B%5Csigma_5%7D%7B%5Csigma_1%2B%5Ccdots%2B%5Csigma_%7B20%7D%7D%5Cquad%5Ctext%7Bwhere%20%7D%5Csigma_i%5Ctext%7B%20is%20a%20singular%20value%7D
"\\frac{\\sigma_1+\\cdots+\\sigma_5}{\\sigma_1+\\cdots+\\sigma_{20}}\\quad\\text{where }\\sigma_i\\text{ is a singular value}")  

``` r
sing_val <- svd(alpha_new)$d
sum(sing_val[1:5])/sum(sing_val)
```

    ## [1] 0.4385151

  
![||A-\\hat{A}||\_F](https://latex.codecogs.com/png.latex?%7C%7CA-%5Chat%7BA%7D%7C%7C_F
"||A-\\hat{A}||_F")  

``` r
(LR_mat - alpha_new) %>%
  Matrix::norm(type = "F")
```

    ## [1] 44.15708

### Sparse matrix

  
![||V\_{\\tau\_\\ell}\\hat{\\Theta}-X\\Gamma(\\tau)||\_F](https://latex.codecogs.com/png.latex?%7C%7CV_%7B%5Ctau_%5Cell%7D%5Chat%7B%5CTheta%7D-X%5CGamma%28%5Ctau%29%7C%7C_F
"||V_{\\tau_\\ell}\\hat{\\Theta}-X\\Gamma(\\tau)||_F")  

``` r
lapply(V, FUN = function(x) ((x %*% theta_new) - X %*% sp_mat) %>% Matrix::norm(type = "F")) %>%
  `names<-`(value = tau_seq) %>%
  bind_cols %>% t() 
```

    ##          [,1]
    ## 0.4  325.6633
    ## 0.41 343.5448
    ## 0.42 309.8238
    ## 0.43 299.9283
    ## 0.44 300.1387
    ## 0.45 302.5406
    ## 0.46 302.0807
    ## 0.47 301.2192
    ## 0.48 301.8234
    ## 0.49 301.8836
    ## 0.5  301.6634
    ## 0.51 302.4676
    ## 0.52 302.9132
    ## 0.53 302.7039
    ## 0.54 303.9290
    ## 0.55 304.7996
    ## 0.56 302.9462
    ## 0.57 303.3842
    ## 0.58 313.3453
    ## 0.59 345.5079
    ## 0.6  329.1110

![X\\hat{\\Gamma}(\\tau)=V\_{\\tau\_\\ell}\\hat{\\Theta}](https://latex.codecogs.com/png.latex?X%5Chat%7B%5CGamma%7D%28%5Ctau%29%3DV_%7B%5Ctau_%5Cell%7D%5Chat%7B%5CTheta%7D
"X\\hat{\\Gamma}(\\tau)=V_{\\tau_\\ell}\\hat{\\Theta}")

![\\hat{\\Gamma}(\\tau)=(X^TX)^{-1}X^TV\_{\\tau\_\\ell}\\hat{\\Theta}](https://latex.codecogs.com/png.latex?%5Chat%7B%5CGamma%7D%28%5Ctau%29%3D%28X%5ETX%29%5E%7B-1%7DX%5ETV_%7B%5Ctau_%5Cell%7D%5Chat%7B%5CTheta%7D
"\\hat{\\Gamma}(\\tau)=(X^TX)^{-1}X^TV_{\\tau_\\ell}\\hat{\\Theta}")

``` r
gamma_tau_hat <- lapply(V, FUN = function(x) solve(t(X) %*% X) %*% t(X) %*% x %*% theta_new) %>%
  `names<-`(value = tau_seq)
lapply(gamma_tau_hat, FUN = function(x) (x - sp_mat) %>% Matrix::norm(type = "F")) %>%
  bind_cols %>% t()
```

    ##          [,1]
    ## 0.4  19.36359
    ## 0.41 20.31675
    ## 0.42 18.36523
    ## 0.43 17.67847
    ## 0.44 17.66829
    ## 0.45 17.83623
    ## 0.46 17.79435
    ## 0.47 17.72865
    ## 0.48 17.77671
    ## 0.49 17.77540
    ## 0.5  17.74744
    ## 0.51 17.79991
    ## 0.52 17.82164
    ## 0.53 17.78757
    ## 0.54 17.86407
    ## 0.55 17.91877
    ## 0.56 17.77074
    ## 0.57 17.80109
    ## 0.58 18.48764
    ## 0.59 20.38757
    ## 0.6  19.47473

![\\\#\\text{ of non-zero entry in
}\\Gamma(\\tau)=500](https://latex.codecogs.com/png.latex?%5C%23%5Ctext%7B%20of%20non-zero%20entry%20in%20%7D%5CGamma%28%5Ctau%29%3D500
"\\#\\text{ of non-zero entry in }\\Gamma(\\tau)=500")

![x\<0.1^5](https://latex.codecogs.com/png.latex?x%3C0.1%5E5
"x\<0.1^5")이면 0으로 간주.

### sparsity pattern check

  
![\\begin{aligned}\\text{common}&:=\\text{zero index in
}\\Gamma(\\tau)=\\text{zero index in
}\\hat{\\Gamma}(\\tau)\\\\\\text{diff}&:=\\text{zero index in
}\\Gamma(\\tau)\\ne\\text{zero index in
}\\hat{\\Gamma}(\\tau)\\\\\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Ctext%7Bcommon%7D%26%3A%3D%5Ctext%7Bzero%20index%20in%20%7D%5CGamma%28%5Ctau%29%3D%5Ctext%7Bzero%20index%20in%20%7D%5Chat%7B%5CGamma%7D%28%5Ctau%29%5C%5C%5Ctext%7Bdiff%7D%26%3A%3D%5Ctext%7Bzero%20index%20in%20%7D%5CGamma%28%5Ctau%29%5Cne%5Ctext%7Bzero%20index%20in%20%7D%5Chat%7B%5CGamma%7D%28%5Ctau%29%5C%5C%5Cend%7Baligned%7D
"\\begin{aligned}\\text{common}&:=\\text{zero index in }\\Gamma(\\tau)=\\text{zero index in }\\hat{\\Gamma}(\\tau)\\\\\\text{diff}&:=\\text{zero index in }\\Gamma(\\tau)\\ne\\text{zero index in }\\hat{\\Gamma}(\\tau)\\\\\\end{aligned}")  

``` r
num_zero <- which(sp_mat==0, arr.ind = TRUE) %>% nrow
print(paste0("The number of zero entry in sp_mat is ", num_zero))
```

    ## [1] "The number of zero entry in sp_mat is 520"

``` r
check_sp_pattern <- function(true, est, tol = 0.1^5, idx = FALSE) {
  # diff1 means that true = 0 but est != 0
  # diff2 means that true != 0 but est = 0
  zero_idx_true <- which(abs(true) < tol, arr.ind = TRUE) %>% as_tibble
  zero_idx_est <- which(abs(est) < tol, arr.ind = TRUE) %>% as_tibble
  common <- semi_join(zero_idx_true, zero_idx_est, by = c("row", "col"))
  diff_1 <- anti_join(zero_idx_true, zero_idx_est, by = c("row", "col"))
  diff_2 <- anti_join(zero_idx_est, zero_idx_true, by = c("row", "col"))
  if(idx == FALSE) {return(data.frame(diff1 = nrow(diff_1), common = nrow(common), diff2 = nrow(diff_2)))
    } else {return(list(common = common, diff = diff))}
}
```

``` r
gamma_tau_hat %>%
  lapply(FUN = function(x) check_sp_pattern(true = sp_mat, est = x)) %>%
  bind_rows(.id = "tau") %>%
  mutate(accuracy = common/num_zero)
```

    ##     tau diff1 common diff2  accuracy
    ## 1   0.4   133    387   262 0.7442308
    ## 2  0.41   226    294   188 0.5653846
    ## 3  0.42   225    295   188 0.5673077
    ## 4  0.43   226    294   188 0.5653846
    ## 5  0.44   226    294   187 0.5653846
    ## 6  0.45   226    294   188 0.5653846
    ## 7  0.46   227    293   188 0.5634615
    ## 8  0.47   205    315   195 0.6057692
    ## 9  0.48   205    315   196 0.6057692
    ## 10 0.49   205    315   195 0.6057692
    ## 11  0.5   203    317   198 0.6096154
    ## 12 0.51   205    315   196 0.6057692
    ## 13 0.52   207    313   196 0.6019231
    ## 14 0.53   207    313   196 0.6019231
    ## 15 0.54   221    299   189 0.5750000
    ## 16 0.55   220    300   189 0.5769231
    ## 17 0.56   221    299   189 0.5750000
    ## 18 0.57   221    299   189 0.5750000
    ## 19 0.58   220    300   189 0.5769231
    ## 20 0.59   221    299   189 0.5750000
    ## 21  0.6   137    383   265 0.7365385
