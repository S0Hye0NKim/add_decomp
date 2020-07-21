

    library(tidyverse)
    library(ggplot2)

\#Initial Setting

$$\\begin{aligned}Y &\\in \\mathbb{R}^{n\\times g}\\\\Y &= XB(\\tau) + \\epsilon(\\tau) , \\\\ X &\\in\\mathbb{R}^{n \\times (p+1)}, \\; \\;B(\\tau) \\in \\mathbb{R}^{(p+1)\\times g} , \\; and \\; Q\_\\tau(\\epsilon(\\tau)|X) = 0\\end{aligned}$$

-   n : \# of observation.
-   g : \# of group.
-   p : \# of covariate.
-   m : nonzero entry in true *B*(*τ*)

For toy example, I will use n = 300, p = 50, g = 20.

    set.seed(0)
    n <- 300
    g <- 20
    p <- 50
    m <- 30

    X <- matrix(rnorm(n*p, mean = 0, sd = 1), nrow = n) %>% cbind(1, .)     #add intercept term in X
    true_B_tau <- matrix(0, nrow = p+1, ncol = g) 
    nonzero_row <- sample(1:(p+1), size = m, replace = TRUE)
    nonzero_col <- c(1:g, sample(1:g, size = m - g, replace = TRUE))
    for(i in 1:m) {true_B_tau[nonzero_row[i], nonzero_col[i]] = rnorm(n = 1, mean = 3, sd = 1)}

*B*(*τ*)는 각 열들을 기준으로 최소 1개의 non-zero entry, 하지만 특정
행에 대해서는 모두 0인 값을 가진다.
