# add_decomp : additive decomposition into low-rank and sparse matrices in multivariate regional quantile regression with a large number of covariates.

## Overview
We propose multiple response regional quantile regression by imposing structural condition on the underlying coefficient matrix. In the CCLE data analyses, we assume that only a few genes are relevant to the effect of drug resistance responses and some genes could have similar effects on multiple responses. We propose a penalized multivariate quantile regression by decomposing the quantile coefficient function into the low-rank and sparse matrix. Low-rank part is a constant function of quantile levels, which represents the global pattern of the coefficient function, whereas the sparse matrix can be smoothly varying by quantile levels, which represents a local and specific pattern of the coefficient function. 

---

## Main functions

[calc_V](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.cpp) : Calculate kronecker product between X and spline function. (C++)

[add_decomp](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.cpp) : Main ADMM algorithm estimating low-rank and sparse coefficient matrix in regional quatnile regression. (C++)

[LR_model](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.cpp) : ADMM algorithm estimating only low-rank coefficient matrix. (C++)

[SP_model](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.cpp) : ADMM algorithm estimating only sparse coefficient matrix. (C++)

[add_decmop_BIC](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.R) : Calculate the GIC for low-rank and sparse penalties. (R)

[LR_model_BIC](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.R) : Calculate the GIC for only low-rank penalty. (R)

[SP_model_BIC](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.R) : Calculate the GIC for only sparse penalty. (R)

[check_sp_table](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.R) : Calculate TP, TN, FP, FN of sparse matrix. (R)

[est_gamma](https://github.com/S0Hye0NKim/add_decomp/blob/master/functions/add_decomp_function.R) : Estimate smooth function of sparse matrix. (R)

---

## Directory
All the functions used in the proposed algorithm are located in the directory "[**functions**](https://github.com/S0Hye0NKim/add_decomp/tree/master/functions)"

All the codes used in simulation are located in the directory "[**simulation**](https://github.com/S0Hye0NKim/add_decomp/tree/master/simulation)"

All the codes used in CCLE analysis are located in the directory "[**CCLE**](https://github.com/S0Hye0NKim/add_decomp/tree/master/CCLE)"

---
## Version information

### Program
- R : 3.6.3
  
### Package
- dplyr : 1.0.5
- stringr : 1.4.0
- Matrix : 1.3.2
- data.table : 1.13.6
- splines : 3.6.3
- fda : 5.19
- glmnet : 4.1.1
- foreach : 1.5.1
- doParallel : 1.0.16
- gplots : 3.1.1
- Rcpp : 1.0.6
- RcpArmadillo : 0.10.2.2.0