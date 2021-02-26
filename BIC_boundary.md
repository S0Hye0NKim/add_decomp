BIC parameter boundary
================
Sohyeon Kim

  
![Q=\\sum\_{g=1}^m\\sum\_{\\ell=1}^b\\sum\_{i=1}^n\\frac{1}{n}\\rho\_{\\tau\_\\ell}\\bigg(Y\_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum\_{j=0}^px\_{ij}\\Phi(\\tau\_\\ell)^T\\boldsymbol{\\theta}\_j^{(g)}\\bigg)+\\lambda\_1||\\textbf{XA}||\_{\*\\psi}+\\lambda\_2\\sum\_{g=1}^m\\sum\_{j=1}^p\\phi\_j^{(g)}||\\boldsymbol{\\theta}^{(g)}\_j||\_2](https://latex.codecogs.com/png.latex?Q%3D%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%5Cfrac%7B1%7D%7Bn%7D%5Crho_%7B%5Ctau_%5Cell%7D%5Cbigg%28Y_i%5E%7B%28g%29%7D-%5Ctextbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Csum_%7Bj%3D0%7D%5Epx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%5Cbigg%29%2B%5Clambda_1%7C%7C%5Ctextbf%7BXA%7D%7C%7C_%7B%2A%5Cpsi%7D%2B%5Clambda_2%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7Bj%3D1%7D%5Ep%5Cphi_j%5E%7B%28g%29%7D%7C%7C%5Cboldsymbol%7B%5Ctheta%7D%5E%7B%28g%29%7D_j%7C%7C_2
"Q=\\sum_{g=1}^m\\sum_{\\ell=1}^b\\sum_{i=1}^n\\frac{1}{n}\\rho_{\\tau_\\ell}\\bigg(Y_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum_{j=0}^px_{ij}\\Phi(\\tau_\\ell)^T\\boldsymbol{\\theta}_j^{(g)}\\bigg)+\\lambda_1||\\textbf{XA}||_{*\\psi}+\\lambda_2\\sum_{g=1}^m\\sum_{j=1}^p\\phi_j^{(g)}||\\boldsymbol{\\theta}^{(g)}_j||_2")  

# Boundary for lambda2

  
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

  
![\\begin{aligned}||\\sum\_{\\ell=1}^b\\sum\_{i=1}^nx\_{ij}\\Phi(\\tau\_\\ell)\\frac{t\_i^{(\\ell)(g)}}{n|\\phi\_j^{(g)}|}||\_2&\\le\\frac{1}{n|\\phi\_j^{(g)}|}||\\sum\_{\\ell=1}^b\\sum\_{i=1}^n|x\_{ij}|\\cdot\\Phi(\\tau\_\\ell)||\_2\\\\&\\qquad(\\because
|t\_i^{(\\ell)(g)}|\\le\\text{max}(|\\tau\_\\ell-1|,\\tau\_\\ell)\\le1)\\\\&\\le\\frac{1}{n|\\phi\_j^{(g)}|}\\sum\_{i=1}^n|x\_{ij}|\\cdot\\sum\_{\\ell=1}^b||\\Phi(\\tau\_\\ell)||\_2\\\\&\\qquad(\\because\\text{
triangle
inequality)}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%7C%7C%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5Enx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5Cfrac%7Bt_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7Bn%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%7C%7C_2%26%5Cle%5Cfrac%7B1%7D%7Bn%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%7C%7C%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%7Cx_%7Bij%7D%7C%5Ccdot%5CPhi%28%5Ctau_%5Cell%29%7C%7C_2%5C%5C%26%5Cqquad%28%5Cbecause%20%7Ct_i%5E%7B%28%5Cell%29%28g%29%7D%7C%5Cle%5Ctext%7Bmax%7D%28%7C%5Ctau_%5Cell-1%7C%2C%5Ctau_%5Cell%29%5Cle1%29%5C%5C%26%5Cle%5Cfrac%7B1%7D%7Bn%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%5Csum_%7Bi%3D1%7D%5En%7Cx_%7Bij%7D%7C%5Ccdot%5Csum_%7B%5Cell%3D1%7D%5Eb%7C%7C%5CPhi%28%5Ctau_%5Cell%29%7C%7C_2%5C%5C%26%5Cqquad%28%5Cbecause%5Ctext%7B%20triangle%20inequality%29%7D%5Cend%7Baligned%7D
"\\begin{aligned}||\\sum_{\\ell=1}^b\\sum_{i=1}^nx_{ij}\\Phi(\\tau_\\ell)\\frac{t_i^{(\\ell)(g)}}{n|\\phi_j^{(g)}|}||_2&\\le\\frac{1}{n|\\phi_j^{(g)}|}||\\sum_{\\ell=1}^b\\sum_{i=1}^n|x_{ij}|\\cdot\\Phi(\\tau_\\ell)||_2\\\\&\\qquad(\\because |t_i^{(\\ell)(g)}|\\le\\text{max}(|\\tau_\\ell-1|,\\tau_\\ell)\\le1)\\\\&\\le\\frac{1}{n|\\phi_j^{(g)}|}\\sum_{i=1}^n|x_{ij}|\\cdot\\sum_{\\ell=1}^b||\\Phi(\\tau_\\ell)||_2\\\\&\\qquad(\\because\\text{ triangle inequality)}\\end{aligned}")  

  
![\\therefore\\lambda\_2\\ge\\text{max}\_{g,j}\\frac{1}{n\\,|\\phi\_j^{(g)}|}\\sum\_{\\ell=1}^b\\sum\_{i=1}^n|x\_{ij}|\\cdot||\\Phi(\\tau\_\\ell)||\_2\\rightarrow\\boldsymbol{\\Theta}=0](https://latex.codecogs.com/png.latex?%5Ctherefore%5Clambda_2%5Cge%5Ctext%7Bmax%7D_%7Bg%2Cj%7D%5Cfrac%7B1%7D%7Bn%5C%2C%7C%5Cphi_j%5E%7B%28g%29%7D%7C%7D%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%7Cx_%7Bij%7D%7C%5Ccdot%7C%7C%5CPhi%28%5Ctau_%5Cell%29%7C%7C_2%5Crightarrow%5Cboldsymbol%7B%5CTheta%7D%3D0
"\\therefore\\lambda_2\\ge\\text{max}_{g,j}\\frac{1}{n\\,|\\phi_j^{(g)}|}\\sum_{\\ell=1}^b\\sum_{i=1}^n|x_{ij}|\\cdot||\\Phi(\\tau_\\ell)||_2\\rightarrow\\boldsymbol{\\Theta}=0")  

-----

# Boundary for lambda1

  
![\\sum\_{g=1}^m\\sum\_{\\ell=1}^b\\sum\_{i=1}^n\\frac{1}{n}\\rho\_{\\tau\_\\ell}\\bigg(Y\_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum\_{j=0}^px\_{ij}\\Phi(\\tau\_\\ell)^T\\boldsymbol{\\theta}\_j^{(g)}\\bigg)=\\sum\_{\\ell=1}^b\\frac{1}{n}||\\sqrt{\\rho\_{\\tau\_\\ell}(\\textbf{E})}||\_F^2](https://latex.codecogs.com/png.latex?%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7B%5Cell%3D1%7D%5Eb%5Csum_%7Bi%3D1%7D%5En%5Cfrac%7B1%7D%7Bn%7D%5Crho_%7B%5Ctau_%5Cell%7D%5Cbigg%28Y_i%5E%7B%28g%29%7D-%5Ctextbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Csum_%7Bj%3D0%7D%5Epx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%5Cbigg%29%3D%5Csum_%7B%5Cell%3D1%7D%5Eb%5Cfrac%7B1%7D%7Bn%7D%7C%7C%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28%5Ctextbf%7BE%7D%29%7D%7C%7C_F%5E2
"\\sum_{g=1}^m\\sum_{\\ell=1}^b\\sum_{i=1}^n\\frac{1}{n}\\rho_{\\tau_\\ell}\\bigg(Y_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum_{j=0}^px_{ij}\\Phi(\\tau_\\ell)^T\\boldsymbol{\\theta}_j^{(g)}\\bigg)=\\sum_{\\ell=1}^b\\frac{1}{n}||\\sqrt{\\rho_{\\tau_\\ell}(\\textbf{E})}||_F^2")  

  
![\\text{where
}\\sqrt{\\rho\_{\\tau\_\\ell}(\\textbf{E})}=\\bigg\\{\\sqrt{\\rho\_{\\tau\_\\ell}(Y\_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum\_{j=0}^px\_{ij}\\Phi(\\tau\_\\ell)^T\\boldsymbol{\\theta}\_j^{(g)})}\\bigg\\}\_{i,g}](https://latex.codecogs.com/png.latex?%5Ctext%7Bwhere%20%7D%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28%5Ctextbf%7BE%7D%29%7D%3D%5Cbigg%5C%7B%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28Y_i%5E%7B%28g%29%7D-%5Ctextbf%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Csum_%7Bj%3D0%7D%5Epx_%7Bij%7D%5CPhi%28%5Ctau_%5Cell%29%5ET%5Cboldsymbol%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%29%7D%5Cbigg%5C%7D_%7Bi%2Cg%7D
"\\text{where }\\sqrt{\\rho_{\\tau_\\ell}(\\textbf{E})}=\\bigg\\{\\sqrt{\\rho_{\\tau_\\ell}(Y_i^{(g)}-\\textbf{X}\\boldsymbol{\\alpha}^{(g)}-\\sum_{j=0}^px_{ij}\\Phi(\\tau_\\ell)^T\\boldsymbol{\\theta}_j^{(g)})}\\bigg\\}_{i,g}")  

  
![\\frac{\\partial
Q}{\\partial\\textbf{Z}}=-\\frac{1}{n}\\sum\_{\\ell=1}^n\\frac{\\partial\\,||\\sqrt{\\rho\_{\\tau\_\\ell}(\\textbf{E})}||\_F^2}{\\partial\\textbf{E}}+\\lambda\_1\\frac{\\partial||\\textbf{Z}||\_{\*,\\psi}}{\\partial\\,\\textbf{Z}}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%5Ctextbf%7BZ%7D%7D%3D-%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7B%5Cell%3D1%7D%5En%5Cfrac%7B%5Cpartial%5C%2C%7C%7C%5Csqrt%7B%5Crho_%7B%5Ctau_%5Cell%7D%28%5Ctextbf%7BE%7D%29%7D%7C%7C_F%5E2%7D%7B%5Cpartial%5Ctextbf%7BE%7D%7D%2B%5Clambda_1%5Cfrac%7B%5Cpartial%7C%7C%5Ctextbf%7BZ%7D%7C%7C_%7B%2A%2C%5Cpsi%7D%7D%7B%5Cpartial%5C%2C%5Ctextbf%7BZ%7D%7D
"\\frac{\\partial Q}{\\partial\\textbf{Z}}=-\\frac{1}{n}\\sum_{\\ell=1}^n\\frac{\\partial\\,||\\sqrt{\\rho_{\\tau_\\ell}(\\textbf{E})}||_F^2}{\\partial\\textbf{E}}+\\lambda_1\\frac{\\partial||\\textbf{Z}||_{*,\\psi}}{\\partial\\,\\textbf{Z}}")
