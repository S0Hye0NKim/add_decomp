Toy Example
================
Sohyeon Kim
7/21/2020



``` r
library(tidyverse)
library(ggplot2)
library(splines)
```

# Initial Setting

  
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
r <- 30
b <- 20

X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) %>% cbind(1, .)     #add intercept term in X
true_B <- matrix(0, nrow = p+1, ncol = m) 
nonzero_row <- sample(1:(p+1), size = r, replace = TRUE)
nonzero_col <- c(1:m, sample(1:m, size = r - m, replace = TRUE))
for(i in 1:r) {true_B[nonzero_row[i], nonzero_col[i]] = rnorm(n = 1, mean = 3, sd = 1)}

eps <- matrix(rnorm(n*m, mean = 0, sd = 0.1), nrow = n)
tau <- 0.5
true_B_tau <- true_B
true_B_tau[1, ] <- true_B[1, ] + qnorm(tau, mean = 0, sd = 0.1)
eps_tau <- eps - qnorm(tau, mean = 0, sd = 0.1)

Y <- X%*%true_B_tau + eps_tau
```

![B(\\tau)](https://latex.codecogs.com/png.latex?B%28%5Ctau%29
"B(\\tau)")는 각 열들을 기준으로 최소 1개의 non-zero entry, 하지만 특정 행에 대해서는 모두 0인 값을
가진다.

## Question

1.  Y를 generate할 때, 사용하는
    ![B(\\tau)](https://latex.codecogs.com/png.latex?B%28%5Ctau%29
    "B(\\tau)") 값에 tau depend?
2.  우선 ![\\tau
    = 0.5](https://latex.codecogs.com/png.latex?%5Ctau%20%3D%200.5
    "\\tau = 0.5") 로 generate 시켰는데, 이렇게 하는게 맞나?

# Check function

``` r
check_ft <- function(x, tau) {
  z <- ifelse(x<0, (tau-1)*x, tau*x)
  return(z)
}
```

  
![\\rho\_\\tau(x) = (\\tau -
\\mathbb{I}(x\<0))x](https://latex.codecogs.com/png.latex?%5Crho_%5Ctau%28x%29%20%3D%20%28%5Ctau%20-%20%5Cmathbb%7BI%7D%28x%3C0%29%29x
"\\rho_\\tau(x) = (\\tau - \\mathbb{I}(x\<0))x")  

# Basis function

``` r
knot <- 10
degree <- 3
K <- knot + degree + 1 # the number of basis function

tau_seq <- seq(from = 0.25, to = 0.75, length.out = b + 1) 
Phi <- bs(tau_seq, df = K, intercept = TRUE)
```

## Questions

1.  tau seq 이렇게 설정하는게 맞나?
2.  degree of freedom 설정

  
![\\begin{aligned}\\boldsymbol{v}\_i^{(\\ell)} &= (\\boldsymbol{x}\_i
\\otimes \\Phi(\\tau\_{\\ell}))^T\\\\ &=(\\Phi(\\tau\_\\ell)^T,
x\_{i1}\\Phi(\\tau\_\\ell)^T, \\dots, x\_{ip}\\Phi(\\tau\_\\ell)^T) \\\\
&\\in
\\mathbb{R}^{(p+1)K}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7Bv%7D_i%5E%7B%28%5Cell%29%7D%20%26%3D%20%28%5Cboldsymbol%7Bx%7D_i%20%5Cotimes%20%5CPhi%28%5Ctau_%7B%5Cell%7D%29%29%5ET%5C%5C%20%26%3D%28%5CPhi%28%5Ctau_%5Cell%29%5ET%2C%20x_%7Bi1%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%2C%20%5Cdots%2C%20x_%7Bip%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%29%20%5C%5C%20%26%5Cin%20%5Cmathbb%7BR%7D%5E%7B%28p%2B1%29K%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{v}_i^{(\\ell)} &= (\\boldsymbol{x}_i \\otimes \\Phi(\\tau_{\\ell}))^T\\\\ &=(\\Phi(\\tau_\\ell)^T, x_{i1}\\Phi(\\tau_\\ell)^T, \\dots, x_{ip}\\Phi(\\tau_\\ell)^T) \\\\ &\\in \\mathbb{R}^{(p+1)K}\\end{aligned}")  

# New design matrix for regional quantile

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

# Initial value for estimator

  
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
eta_init <- matrix(rpois(m*(p+1)*K, lambda = 5), nrow = (p+1)*K, ncol = m) # theta와 동일하게?
theta_init <- matrix(rpois(m*(p+1)*K, lambda = 3), nrow = (p+1)*K, ncol = m) # Quant reg값?
alpha_init <- matrix(1, nrow = p+1, ncol = m)
e_init <- list()
for(l in 1:b) {e_init[[l]] <- Y - X %*% alpha_init - V[[l]] %*% eta_init}
u_init <- matrix(1, nrow = n, ncol = 1)
w_init <- matrix(1, nrow = (p+1)*K, ncol = 1)
```

## Question

1.  dual problem 전개할 때, ![\\theta\_0 = \\eta\_0, \\;e\_0 =
    Y-X\\alpha=V-\\eta](https://latex.codecogs.com/png.latex?%5Ctheta_0%20%3D%20%5Ceta_0%2C%20%5C%3Be_0%20%3D%20Y-X%5Calpha%3DV-%5Ceta
    "\\theta_0 = \\eta_0, \\;e_0 = Y-X\\alpha=V-\\eta") ??? 초기값 설정 어떻게?
    -\> 우선 ![\\theta,
    \\eta](https://latex.codecogs.com/png.latex?%5Ctheta%2C%20%5Ceta
    "\\theta, \\eta")는 각각 설정하고, e는 equality condition과 동일하게 설정.
2.  Objective function이 convex한지? convex 안하면 초기값에 따라 local에 빠질 수도…
    global을 찾자\!

# closed form

  
![ \\frac{\\partial Q}{\\partial \\boldsymbol{\\eta}^{(g)T}} =
-\\sum\_{\\ell}^b \\boldsymbol{V}^{(\\ell)T} \\boldsymbol{u}^{(\\ell)} -
\\delta\\sum\_{\\ell=1}^bV^{(\\ell)T}(\\boldsymbol{Y} -
\\boldsymbol{X}\\boldsymbol{\\alpha} -
\\boldsymbol{e}^{(\\ell)})+\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta^{(g)}}
- \\boldsymbol{w} + \\delta\\boldsymbol{\\eta}^{(g)} -
\\delta\\boldsymbol{\\theta} = 0
](https://latex.codecogs.com/png.latex?%20%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29T%7D%7D%20%3D%20-%5Csum_%7B%5Cell%7D%5Eb%20%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%20%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%7D%20-%20%5Cdelta%5Csum_%7B%5Cell%3D1%7D%5EbV%5E%7B%28%5Cell%29T%7D%28%5Cboldsymbol%7BY%7D%20-%20%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%20-%20%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%7D%29%2B%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%5E%7B%28g%29%7D%7D%20-%20%5Cboldsymbol%7Bw%7D%20%2B%20%5Cdelta%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%7D%20-%20%5Cdelta%5Cboldsymbol%7B%5Ctheta%7D%20%20%3D%200%20
" \\frac{\\partial Q}{\\partial \\boldsymbol{\\eta}^{(g)T}} = -\\sum_{\\ell}^b \\boldsymbol{V}^{(\\ell)T} \\boldsymbol{u}^{(\\ell)} - \\delta\\sum_{\\ell=1}^bV^{(\\ell)T}(\\boldsymbol{Y} - \\boldsymbol{X}\\boldsymbol{\\alpha} - \\boldsymbol{e}^{(\\ell)})+\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta^{(g)}} - \\boldsymbol{w} + \\delta\\boldsymbol{\\eta}^{(g)} - \\delta\\boldsymbol{\\theta}  = 0 ")  

Closed form for ![\\eta](https://latex.codecogs.com/png.latex?%5Ceta
"\\eta")

  
![\\hat{\\boldsymbol{\\eta}}^{(g)} =
\\frac{1}{\\delta}(\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}
+ \\boldsymbol{I})^{-1}(\\boldsymbol{w} +
\\delta\\boldsymbol{\\theta}+\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}
+ \\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y} -
\\boldsymbol{X}\\boldsymbol{\\alpha} -
\\boldsymbol{e}))](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cboldsymbol%7B%5Ceta%7D%7D%5E%7B%28g%29%7D%20%3D%20%5Cfrac%7B1%7D%7B%5Cdelta%7D%28%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%20%2B%20%5Cboldsymbol%7BI%7D%29%5E%7B-1%7D%28%5Cboldsymbol%7Bw%7D%20%2B%20%5Cdelta%5Cboldsymbol%7B%5Ctheta%7D%2B%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7Bu%7D%20%2B%20%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%28%5Cboldsymbol%7BY%7D%20-%20%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%20-%20%5Cboldsymbol%7Be%7D%29%29
"\\hat{\\boldsymbol{\\eta}}^{(g)} = \\frac{1}{\\delta}(\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)} + \\boldsymbol{I})^{-1}(\\boldsymbol{w} + \\delta\\boldsymbol{\\theta}+\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u} + \\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y} - \\boldsymbol{X}\\boldsymbol{\\alpha} - \\boldsymbol{e}))")  

Closed form for ![\\eta](https://latex.codecogs.com/png.latex?%5Ceta
"\\eta")

  
![\\frac{\\partial Q}{\\partial \\boldsymbol{\\theta}^{(g)T}} =
\\boldsymbol{w} + \\delta\\boldsymbol{\\theta} -
\\delta\\boldsymbol{\\eta} + \\frac{\\partial
\\lambda\_2||\\boldsymbol{\\theta}||\_2}{\\partial
\\boldsymbol{\\theta}}
= 0](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20%5Cboldsymbol%7B%5Ctheta%7D%5E%7B%28g%29T%7D%7D%20%3D%20%5Cboldsymbol%7Bw%7D%20%2B%20%5Cdelta%5Cboldsymbol%7B%5Ctheta%7D%20-%20%5Cdelta%5Cboldsymbol%7B%5Ceta%7D%20%2B%20%5Cfrac%7B%5Cpartial%20%5Clambda_2%7C%7C%5Cboldsymbol%7B%5Ctheta%7D%7C%7C_2%7D%7B%5Cpartial%20%5Cboldsymbol%7B%5Ctheta%7D%7D%20%3D%200
"\\frac{\\partial Q}{\\partial \\boldsymbol{\\theta}^{(g)T}} = \\boldsymbol{w} + \\delta\\boldsymbol{\\theta} - \\delta\\boldsymbol{\\eta} + \\frac{\\partial \\lambda_2||\\boldsymbol{\\theta}||_2}{\\partial \\boldsymbol{\\theta}} = 0")  
Closed form for ![\\theta](https://latex.codecogs.com/png.latex?%5Ctheta
"\\theta") -\> penalty term 미분…?

  
![\\begin{aligned} \\frac{\\partial Q}{\\partial \\alpha^{(g)T}} &=
-\\sum\_{\\ell}^b\\boldsymbol{X}\\boldsymbol{u}^{(\\ell)}
-\\delta\\boldsymbol{X}^T\\sum\_{\\ell}^b(\\boldsymbol{Y} -
\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta} -
\\boldsymbol{e}^{(\\ell)}) +
\\delta\\boldsymbol{X}^T\\boldsymbol{X}\\boldsymbol{\\alpha} +
\\frac{\\partial \\lambda\_1||\\boldsymbol{A}||\_\*}{\\partial
\\boldsymbol{\\alpha}} \\\\
&= 0\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%20%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20%5Calpha%5E%7B%28g%29T%7D%7D%20%26%3D%20-%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BX%7D%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%7D%20-%5Cdelta%5Cboldsymbol%7BX%7D%5ET%5Csum_%7B%5Cell%7D%5Eb%28%5Cboldsymbol%7BY%7D%20-%20%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%7D%20-%20%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%7D%29%20%2B%20%5Cdelta%5Cboldsymbol%7BX%7D%5ET%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%20%2B%20%5Cfrac%7B%5Cpartial%20%5Clambda_1%7C%7C%5Cboldsymbol%7BA%7D%7C%7C_%2A%7D%7B%5Cpartial%20%5Cboldsymbol%7B%5Calpha%7D%7D%20%5C%5C%20%26%3D%200%5Cend%7Baligned%7D
"\\begin{aligned} \\frac{\\partial Q}{\\partial \\alpha^{(g)T}} &= -\\sum_{\\ell}^b\\boldsymbol{X}\\boldsymbol{u}^{(\\ell)} -\\delta\\boldsymbol{X}^T\\sum_{\\ell}^b(\\boldsymbol{Y} - \\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta} - \\boldsymbol{e}^{(\\ell)}) + \\delta\\boldsymbol{X}^T\\boldsymbol{X}\\boldsymbol{\\alpha} + \\frac{\\partial \\lambda_1||\\boldsymbol{A}||_*}{\\partial \\boldsymbol{\\alpha}} \\\\ &= 0\\end{aligned}")  

Closed form for ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha
"\\alpha") -\> penalty term 미분…?

## Questions

1.  Penalty term 미분 어떻게??
2.  check function 미분 어떻게??
