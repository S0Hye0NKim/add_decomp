Toy Example
================
Sohyeon Kim
7/21/2020



``` r
library(tidyverse)
library(ggplot2)
library(splines)
library(Matrix)
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
b <- 20

X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) %>% cbind(1, .)   #add intercept term in X
sp_mat <- rsparsematrix(nrow = p+1, ncol = m, nnz = r)             # make sparse matrix
LR_mat <- matrix(rnorm(m*(p+1), mean = 0, sd = 1), ncol = m)      # make low rank matrix using SVD
SVD <- svd(LR_mat)
D_mat <- diag(SVD$d)
idx <- sample(1:m, 15)
D_mat[idx, idx] <- 0
LR_mat <- SVD$u %*% D_mat %*% t(SVD$v)          

true_B <- sp_mat + LR_mat            # B(tau) = sparse matrix + low rank matrix

eps <- matrix(rnorm(n*m, mean = 0, sd = 0.1), nrow = n)
#tau <- 0.5
#true_B_tau <- true_B
#true_B_tau[1, ] <- true_B[1, ] + qnorm(tau, mean = 0, sd = 0.1)
#eps_tau <- eps - qnorm(tau, mean = 0, sd = 0.1)

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
eta_old <- matrix(rpois(m*(p+1)*K, lambda = 5), nrow = (p+1)*K, ncol = m) 
theta_old <- matrix(rpois(m*(p+1)*K, lambda = 3), nrow = (p+1)*K, ncol = m) 
alpha_old <- matrix(1, nrow = p+1, ncol = m)
Z_old <- X %*% alpha_old
e_old <- list()
for(l in 1:b) {e_old[[l]] <- Y - X %*% alpha_old - V[[l]] %*% eta_old}
u_old <- list()
for(l in 1:b) {u_old[[l]] <- matrix(1, nrow = n, ncol = m)}
w_old <- matrix(1, nrow = (p+1)*K, ncol = m)
```

# 2\. Algorithm

  
![\\begin{aligned}\\boldsymbol{Y}&=\\boldsymbol{XA}+\\boldsymbol{X\\Gamma}(\\tau)\\\\&=\\boldsymbol{Z}+\\boldsymbol{V}\_{\\tau\_\\ell}\\boldsymbol{\\Theta}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7BY%7D%26%3D%5Cboldsymbol%7BXA%7D%2B%5Cboldsymbol%7BX%5CGamma%7D%28%5Ctau%29%5C%5C%26%3D%5Cboldsymbol%7BZ%7D%2B%5Cboldsymbol%7BV%7D_%7B%5Ctau_%5Cell%7D%5Cboldsymbol%7B%5CTheta%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{Y}&=\\boldsymbol{XA}+\\boldsymbol{X\\Gamma}(\\tau)\\\\&=\\boldsymbol{Z}+\\boldsymbol{V}_{\\tau_\\ell}\\boldsymbol{\\Theta}\\end{aligned}")  

``` r
max_iter <- 50
delta <- 0.05
lambda_1 <- 0.5
lambda_2 <- 0.5
tol_error <- 0.01

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
    value <- eta_new[, g] - w_old[, g]/delta 
    theta_new[, g] <- case_when(value < -threshold ~ value + threshold, 
                                abs(value) < threshold ~ 0, 
                                value > threshold ~ value - threshold)
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
  
  # Process for e
  e_new <- list()
  for(l in 1:b){
    e_new[[l]] <- matrix(nrow = n, ncol = m)
    common_val <- Y - Z_new - V[[l]] %*% eta_new
    for(g in 1:m) {
      e_new[[l]][, g] <- case_when(e_old[[l]][, g] < 0 ~ common_val[, g] - (u_old[[l]][, g] + tau_seq[[l]] -1)/delta,
                                   e_old[[l]][, g] == 0 ~ common_val[, g] - u_old[[l]][, g]/delta, 
                                   e_old[[l]][, g] > 0 ~ common_val[, g] - (u_old[[l]][, g] + tau_seq[[l]])/delta)
    }
  }
  
  # Process for multiplier u
  u_new <- list()
  for(l in 1:b) {
    u_new[[l]] <- u_old[[l]] + delta * (Y - Z_new - V[[l]] %*% eta_new - e_new[[l]])
  }
  # Process for multiplier w
  w_new <- w_old + delta * (theta_new - eta_new)
  
  #if(e1 < tol_error) {break}
  
  eta_old <- eta_new
  theta_old <- theta_new
  Z_old <- Z_new
  e_old <- e_new
  u_old <- u_new
  w_old <- w_new
}
```

![\\lambda\_1,
\\lambda\_2](https://latex.codecogs.com/png.latex?%5Clambda_1%2C%20%5Clambda_2
"\\lambda_1, \\lambda_2")는 임의로 0.5로 설정.

## Question

  - Stopping criteria

  
![\\begin{aligned}&||\\eta^{k+1}-\\eta^{k}||\_2^2+||\\theta^{k+1}-\\theta^k||\_2^2+||\\alpha^{k+1}-\\alpha^k||\_2^2\\\\&+||e^{k+1}-e^k||\_2^2+||u^{k+1}-u^k||\_2^2+||w^{k+1}-w^k||\_2^2
\\le
\\epsilon\_{tol}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%7C%7C%5Ceta%5E%7Bk%2B1%7D-%5Ceta%5E%7Bk%7D%7C%7C_2%5E2%2B%7C%7C%5Ctheta%5E%7Bk%2B1%7D-%5Ctheta%5Ek%7C%7C_2%5E2%2B%7C%7C%5Calpha%5E%7Bk%2B1%7D-%5Calpha%5Ek%7C%7C_2%5E2%5C%5C%26%2B%7C%7Ce%5E%7Bk%2B1%7D-e%5Ek%7C%7C_2%5E2%2B%7C%7Cu%5E%7Bk%2B1%7D-u%5Ek%7C%7C_2%5E2%2B%7C%7Cw%5E%7Bk%2B1%7D-w%5Ek%7C%7C_2%5E2%20%5Cle%20%5Cepsilon_%7Btol%7D%5Cend%7Baligned%7D
"\\begin{aligned}&||\\eta^{k+1}-\\eta^{k}||_2^2+||\\theta^{k+1}-\\theta^k||_2^2+||\\alpha^{k+1}-\\alpha^k||_2^2\\\\&+||e^{k+1}-e^k||_2^2+||u^{k+1}-u^k||_2^2+||w^{k+1}-w^k||_2^2 \\le \\epsilon_{tol}\\end{aligned}")
