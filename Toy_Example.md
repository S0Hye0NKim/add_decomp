Toy Example
================
Sohyeon Kim
7/21/2020



``` r
library(tidyverse)
library(ggplot2)
```

\#Initial Setting

  
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
  - g : \# of group.
  - p : \# of covariate.
  - m : nonzero entry in true
    ![B(\\tau)](https://latex.codecogs.com/png.latex?B%28%5Ctau%29
    "B(\\tau)")

For toy example, I will use n = 300, p = 50, g = 20.

``` r
set.seed(0)
n <- 300
g <- 20
p <- 50
m <- 30

X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) %>% cbind(1, .)     #add intercept term in X
true_B <- matrix(0, nrow = p+1, ncol = g) 
nonzero_row <- sample(1:(p+1), size = m, replace = TRUE)
nonzero_col <- c(1:g, sample(1:g, size = m - g, replace = TRUE))
for(i in 1:m) {true_B[nonzero_row[i], nonzero_col[i]] = rnorm(n = 1, mean = 3, sd = 1)}

eps <- matrix(rnorm(n*g, mean = 0, sd = 0.1), nrow = n)
tau <- 0.5
true_B_tau <- true_B
true_B_tau[1, ] <- true_B[1, ] + qnorm(tau, mean = 0, sd = 0.1)
eps_tau <- eps - qnorm(tau, mean = 0, sd = 0.1)

Y <- X%*%true_B_tau + eps_tau
```

![B(\\tau)](https://latex.codecogs.com/png.latex?B%28%5Ctau%29
"B(\\tau)")는 각 열들을 기준으로 최소 1개의 non-zero entry, 하지만 특정 행에 대해서는 모두 0인 값을
가진다.