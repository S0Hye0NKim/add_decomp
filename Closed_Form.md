Closed Form
================
Sohyeon Kim
7/28/2020

## 1\. Closed form for eta.

  
![\\begin{aligned}\\frac{\\partial
Q}{\\partial\\boldsymbol{\\eta}^{(g)T}}&=-\\sum\_{\\ell}^b
\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}^{(\\ell)}-\\delta\\sum\_{\\ell=1}^bV^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}^{(\\ell)})+\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta^{(g)}}-\\boldsymbol{w}+\\delta\\boldsymbol{\\eta}^{(g)}-\\delta\\boldsymbol{\\theta}\\\\&=0\\\\\\hat{\\boldsymbol{\\eta}}^{(g)}&=\\frac{1}{\\delta}(\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\boldsymbol{I})^{-1}(\\boldsymbol{w}+\\delta\\boldsymbol{\\theta}+\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}+\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}))\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29T%7D%7D%26%3D-%5Csum_%7B%5Cell%7D%5Eb%20%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%7D-%5Cdelta%5Csum_%7B%5Cell%3D1%7D%5EbV%5E%7B%28%5Cell%29T%7D%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D-%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%7D%29%2B%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%5E%7B%28g%29%7D%7D-%5Cboldsymbol%7Bw%7D%2B%5Cdelta%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%7D-%5Cdelta%5Cboldsymbol%7B%5Ctheta%7D%5C%5C%26%3D0%5C%5C%5Chat%7B%5Cboldsymbol%7B%5Ceta%7D%7D%5E%7B%28g%29%7D%26%3D%5Cfrac%7B1%7D%7B%5Cdelta%7D%28%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%2B%5Cboldsymbol%7BI%7D%29%5E%7B-1%7D%28%5Cboldsymbol%7Bw%7D%2B%5Cdelta%5Cboldsymbol%7B%5Ctheta%7D%2B%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7Bu%7D%2B%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D-%5Cboldsymbol%7Be%7D%29%29%5Cend%7Baligned%7D
"\\begin{aligned}\\frac{\\partial Q}{\\partial\\boldsymbol{\\eta}^{(g)T}}&=-\\sum_{\\ell}^b \\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}^{(\\ell)}-\\delta\\sum_{\\ell=1}^bV^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}^{(\\ell)})+\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta^{(g)}}-\\boldsymbol{w}+\\delta\\boldsymbol{\\eta}^{(g)}-\\delta\\boldsymbol{\\theta}\\\\&=0\\\\\\hat{\\boldsymbol{\\eta}}^{(g)}&=\\frac{1}{\\delta}(\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\boldsymbol{I})^{-1}(\\boldsymbol{w}+\\delta\\boldsymbol{\\theta}+\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}+\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}))\\end{aligned}")  

  
![\\boldsymbol{H}^{k+1}(=\[\\boldsymbol{\\eta}^{(1)}\\cdots\\boldsymbol{\\eta}^{(m)}\])\\\\:=\\frac{1}{\\delta}(\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\boldsymbol{I})^{-1}(\\boldsymbol{W}+\\delta\\boldsymbol{\\Theta}+\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{U}^{(\\ell)}+\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{Z}-\\boldsymbol{E}^{(\\ell)}))](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D%28%3D%5B%5Cboldsymbol%7B%5Ceta%7D%5E%7B%281%29%7D%5Ccdots%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28m%29%7D%5D%29%5C%5C%3A%3D%5Cfrac%7B1%7D%7B%5Cdelta%7D%28%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%2B%5Cboldsymbol%7BI%7D%29%5E%7B-1%7D%28%5Cboldsymbol%7BW%7D%2B%5Cdelta%5Cboldsymbol%7B%5CTheta%7D%2B%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%2B%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BZ%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%7D%29%29
"\\boldsymbol{H}^{k+1}(=[\\boldsymbol{\\eta}^{(1)}\\cdots\\boldsymbol{\\eta}^{(m)}])\\\\:=\\frac{1}{\\delta}(\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\boldsymbol{I})^{-1}(\\boldsymbol{W}+\\delta\\boldsymbol{\\Theta}+\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{U}^{(\\ell)}+\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{Z}-\\boldsymbol{E}^{(\\ell)}))")  

  
![\\frac{\\partial^2
Q}{\\partial\\boldsymbol{\\eta}^2}=\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\delta\\boldsymbol{I}\\succeq0](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%5E2%20Q%7D%7B%5Cpartial%5Cboldsymbol%7B%5Ceta%7D%5E2%7D%3D%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%2B%5Cdelta%5Cboldsymbol%7BI%7D%5Csucceq0
"\\frac{\\partial^2 Q}{\\partial\\boldsymbol{\\eta}^2}=\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\delta\\boldsymbol{I}\\succeq0")  

  - Second derivative function is semi-postiive definite?

## 2\. Closed form for theta

  
![\\begin{aligned}\\theta\_j^{(g)^{k+1}}:=\&argmin\_{\\theta\_j^{(g)}\\in\\mathbb{R}^{K}}\\sum\_{g=1}^m\\boldsymbol{w}^{(g)^kT}\\theta\_j^{(g)}+
\\frac{\\delta}{2}\\sum\_{g=1}^m||\\theta\_j^{(g)} -
\\eta\_j^{(g)^{k+1}}||\_2^2 +
\\lambda\_2\\sum\_{g=1}^m\\phi\_j^{(g)}||\\theta\_j^{(g)}||\_{2}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Ctheta_j%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%3A%3D%26argmin_%7B%5Ctheta_j%5E%7B%28g%29%7D%5Cin%5Cmathbb%7BR%7D%5E%7BK%7D%7D%5Csum_%7Bg%3D1%7D%5Em%5Cboldsymbol%7Bw%7D%5E%7B%28g%29%5EkT%7D%5Ctheta_j%5E%7B%28g%29%7D%2B%20%5Cfrac%7B%5Cdelta%7D%7B2%7D%5Csum_%7Bg%3D1%7D%5Em%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%20-%20%5Ceta_j%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%7C%7C_2%5E2%20%2B%20%5Clambda_2%5Csum_%7Bg%3D1%7D%5Em%5Cphi_j%5E%7B%28g%29%7D%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_%7B2%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\theta_j^{(g)^{k+1}}:=&argmin_{\\theta_j^{(g)}\\in\\mathbb{R}^{K}}\\sum_{g=1}^m\\boldsymbol{w}^{(g)^kT}\\theta_j^{(g)}+ \\frac{\\delta}{2}\\sum_{g=1}^m||\\theta_j^{(g)} - \\eta_j^{(g)^{k+1}}||_2^2 + \\lambda_2\\sum_{g=1}^m\\phi_j^{(g)}||\\theta_j^{(g)}||_{2}\\end{aligned}")  

where
![\\phi\_j^{(g)}=\\frac{1}{||\\tilde{\\theta}\_j^{(g)}||\_2}](https://latex.codecogs.com/png.latex?%5Cphi_j%5E%7B%28g%29%7D%3D%5Cfrac%7B1%7D%7B%7C%7C%5Ctilde%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%7C%7C_2%7D
"\\phi_j^{(g)}=\\frac{1}{||\\tilde{\\theta}_j^{(g)}||_2}")

  
![\\begin{aligned}\\frac{\\partial Q}{\\partial
\\theta\_j^{(g)T}}&=w\_j+\\delta\\theta\_j^{(g)}-\\delta\\eta\_j^{(g)}+\\lambda\_2\\cdot
\\phi\_j^{(g)}\\frac{\\partial
||\\theta\_j^{(g)}||\_2}{\\partial\\theta\_j^{(g)}}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20%5Ctheta_j%5E%7B%28g%29T%7D%7D%26%3Dw_j%2B%5Cdelta%5Ctheta_j%5E%7B%28g%29%7D-%5Cdelta%5Ceta_j%5E%7B%28g%29%7D%2B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%5Cfrac%7B%5Cpartial%20%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%7D%7B%5Cpartial%5Ctheta_j%5E%7B%28g%29%7D%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\frac{\\partial Q}{\\partial \\theta_j^{(g)T}}&=w_j+\\delta\\theta_j^{(g)}-\\delta\\eta_j^{(g)}+\\lambda_2\\cdot \\phi_j^{(g)}\\frac{\\partial ||\\theta_j^{(g)}||_2}{\\partial\\theta_j^{(g)}}\\end{aligned}")  

  
![\\frac{\\partial||\\theta\_j^{(g)}||\_2}{\\partial\\theta\_j^{(g)}}
=\\begin{cases}\\theta\_j^{(g)}/||\\theta\_j^{(g)}||\_2&\\quad \\text{if
} \\theta\_j^{(g)} \\ne 0\\\\\\{u:||u||\_2\\le 1\\}&\\quad\\text{if }
\\theta\_j^{(g)}=0\\end{cases}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%7D%7B%5Cpartial%5Ctheta_j%5E%7B%28g%29%7D%7D%20%3D%5Cbegin%7Bcases%7D%5Ctheta_j%5E%7B%28g%29%7D%2F%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%26%5Cquad%20%5Ctext%7Bif%20%7D%20%5Ctheta_j%5E%7B%28g%29%7D%20%5Cne%200%5C%5C%5C%7Bu%3A%7C%7Cu%7C%7C_2%5Cle%201%5C%7D%26%5Cquad%5Ctext%7Bif%20%7D%20%5Ctheta_j%5E%7B%28g%29%7D%3D0%5Cend%7Bcases%7D
"\\frac{\\partial||\\theta_j^{(g)}||_2}{\\partial\\theta_j^{(g)}} =\\begin{cases}\\theta_j^{(g)}/||\\theta_j^{(g)}||_2&\\quad \\text{if } \\theta_j^{(g)} \\ne 0\\\\\\{u:||u||_2\\le 1\\}&\\quad\\text{if } \\theta_j^{(g)}=0\\end{cases}")  

-----

if
![\\theta\_j^{(g)}\\ne0](https://latex.codecogs.com/png.latex?%5Ctheta_j%5E%7B%28g%29%7D%5Cne0
"\\theta_j^{(g)}\\ne0")

  
![\\begin{aligned}&(1+\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}/\\delta}{||\\theta\_j^{(g)}||\_2})\\theta\_j^{(g)}=\\underbrace{\\eta\_j^{(g)}-\\frac{w\_j^{(g)}}{\\delta}}\_{r\_j^{(g)}}\\\\&(1+\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}/\\delta}{||\\theta\_j^{(g)}||\_2})||\\theta\_j^{(g)}||\_2=||r\_j^{(g)}||\_2\\\\&\\rightarrow||\\theta\_j^{(g)}||\_2=||r\_j^{(g)}||\_2-\\lambda\_2\\cdot
\\phi\_j^{(g)}/\\delta\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%281%2B%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%7D%7B%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%7D%29%5Ctheta_j%5E%7B%28g%29%7D%3D%5Cunderbrace%7B%5Ceta_j%5E%7B%28g%29%7D-%5Cfrac%7Bw_j%5E%7B%28g%29%7D%7D%7B%5Cdelta%7D%7D_%7Br_j%5E%7B%28g%29%7D%7D%5C%5C%26%281%2B%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%7D%7B%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%7D%29%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%3D%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2%5C%5C%26%5Crightarrow%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%3D%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2-%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%5Cend%7Baligned%7D
"\\begin{aligned}&(1+\\frac{\\lambda_2\\cdot \\phi_j^{(g)}/\\delta}{||\\theta_j^{(g)}||_2})\\theta_j^{(g)}=\\underbrace{\\eta_j^{(g)}-\\frac{w_j^{(g)}}{\\delta}}_{r_j^{(g)}}\\\\&(1+\\frac{\\lambda_2\\cdot \\phi_j^{(g)}/\\delta}{||\\theta_j^{(g)}||_2})||\\theta_j^{(g)}||_2=||r_j^{(g)}||_2\\\\&\\rightarrow||\\theta_j^{(g)}||_2=||r_j^{(g)}||_2-\\lambda_2\\cdot \\phi_j^{(g)}/\\delta\\end{aligned}")  

  
![\\begin{aligned}&||\\theta\_j^{(g)}||\_2=||r\_j^{(g)}||\_2-\\lambda\_2\\cdot
\\phi\_j^{(g)}/\\delta\>0\\\\&\\rightarrow||r\_j^{(g)}||\_2\>\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}}{\\delta}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%3D%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2-%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%3E0%5C%5C%26%5Crightarrow%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2%3E%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%7D%7B%5Cdelta%7D%5Cend%7Baligned%7D
"\\begin{aligned}&||\\theta_j^{(g)}||_2=||r_j^{(g)}||_2-\\lambda_2\\cdot \\phi_j^{(g)}/\\delta\>0\\\\&\\rightarrow||r_j^{(g)}||_2\>\\frac{\\lambda_2\\cdot \\phi_j^{(g)}}{\\delta}\\end{aligned}")  

  
![\\begin{aligned}\\theta\_j^{(g)}&=r\_j^{(g)}(1+\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}/\\delta}{||\\theta\_j^{(g)}||\_2})^{-1}\\\\&=(1-\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}/\\delta}{||r\_j^{(g)}||\_2})r\_j^{(g)}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Ctheta_j%5E%7B%28g%29%7D%26%3Dr_j%5E%7B%28g%29%7D%281%2B%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%7D%7B%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%7D%29%5E%7B-1%7D%5C%5C%26%3D%281-%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%7D%7B%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2%7D%29r_j%5E%7B%28g%29%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\theta_j^{(g)}&=r_j^{(g)}(1+\\frac{\\lambda_2\\cdot \\phi_j^{(g)}/\\delta}{||\\theta_j^{(g)}||_2})^{-1}\\\\&=(1-\\frac{\\lambda_2\\cdot \\phi_j^{(g)}/\\delta}{||r_j^{(g)}||_2})r_j^{(g)}\\end{aligned}")  

-----

if ![\\theta=0](https://latex.codecogs.com/png.latex?%5Ctheta%3D0
"\\theta=0")

  
![w\_j^{(g)}-\\delta\\eta\_j^{(g)}+\\lambda\_2\\phi\_j^{(g)}\\cdot
u=0\\quad\\text{where
}||u||\_2\\le1](https://latex.codecogs.com/png.latex?w_j%5E%7B%28g%29%7D-%5Cdelta%5Ceta_j%5E%7B%28g%29%7D%2B%5Clambda_2%5Cphi_j%5E%7B%28g%29%7D%5Ccdot%20u%3D0%5Cquad%5Ctext%7Bwhere%20%7D%7C%7Cu%7C%7C_2%5Cle1
"w_j^{(g)}-\\delta\\eta_j^{(g)}+\\lambda_2\\phi_j^{(g)}\\cdot u=0\\quad\\text{where }||u||_2\\le1")  
  
![\\begin{aligned}\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}}{\\delta}u&=\\eta\_j^{(g)}-\\frac{w\_j^{(g)}}{\\delta}=r\_j^{(g)}\\\\\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}}{\\delta}||u||\_2&=||\\eta\_j^{(g)}-\\frac{w\_j^{(g)}}{\\delta}||\_2=||r\_j^{(g)}||\_2\\le\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}}{\\delta}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%7D%7B%5Cdelta%7Du%26%3D%5Ceta_j%5E%7B%28g%29%7D-%5Cfrac%7Bw_j%5E%7B%28g%29%7D%7D%7B%5Cdelta%7D%3Dr_j%5E%7B%28g%29%7D%5C%5C%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%7D%7B%5Cdelta%7D%7C%7Cu%7C%7C_2%26%3D%7C%7C%5Ceta_j%5E%7B%28g%29%7D-%5Cfrac%7Bw_j%5E%7B%28g%29%7D%7D%7B%5Cdelta%7D%7C%7C_2%3D%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2%5Cle%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%7D%7B%5Cdelta%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\frac{\\lambda_2\\cdot \\phi_j^{(g)}}{\\delta}u&=\\eta_j^{(g)}-\\frac{w_j^{(g)}}{\\delta}=r_j^{(g)}\\\\\\frac{\\lambda_2\\cdot \\phi_j^{(g)}}{\\delta}||u||_2&=||\\eta_j^{(g)}-\\frac{w_j^{(g)}}{\\delta}||_2=||r_j^{(g)}||_2\\le\\frac{\\lambda_2\\cdot \\phi_j^{(g)}}{\\delta}\\end{aligned}")  

-----

  
![\\therefore\\hat{\\theta}^{(g)}\_j=(1-\\frac{\\lambda\_2\\cdot
\\phi\_j^{(g)}/\\delta}{||r\_j^{(g)}||\_2})\_+r\_j^{(g)}\\quad\\text{where
}r\_j^{(g)}=\\eta\_j^{(g)}-w\_j^{(g)}/\\delta](https://latex.codecogs.com/png.latex?%5Ctherefore%5Chat%7B%5Ctheta%7D%5E%7B%28g%29%7D_j%3D%281-%5Cfrac%7B%5Clambda_2%5Ccdot%20%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%7D%7B%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2%7D%29_%2Br_j%5E%7B%28g%29%7D%5Cquad%5Ctext%7Bwhere%20%7Dr_j%5E%7B%28g%29%7D%3D%5Ceta_j%5E%7B%28g%29%7D-w_j%5E%7B%28g%29%7D%2F%5Cdelta
"\\therefore\\hat{\\theta}^{(g)}_j=(1-\\frac{\\lambda_2\\cdot \\phi_j^{(g)}/\\delta}{||r_j^{(g)}||_2})_+r_j^{(g)}\\quad\\text{where }r_j^{(g)}=\\eta_j^{(g)}-w_j^{(g)}/\\delta")  

  
![\\begin{aligned}&\\rightarrow\\text{if
}\\;(1-\\frac{\\lambda\_2\\cdot\\phi\_j^{(g)}/\\delta}{||r\_j^{(g)}||\_2})\>0,\\quad\\hat{\\theta}\_j^{(g)}\\ne0\\\\&\\rightarrow\\text{if
}\\;\\lambda\_2\<\\delta\\cdot||r\_j^{(g)}||\_2\\cdot||\\tilde{\\theta}\_j^{(g)}||\_2,\\quad\\hat{\\theta}\_j^{(g)}\\ne0\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%5Crightarrow%5Ctext%7Bif%20%7D%5C%3B%281-%5Cfrac%7B%5Clambda_2%5Ccdot%5Cphi_j%5E%7B%28g%29%7D%2F%5Cdelta%7D%7B%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2%7D%29%3E0%2C%5Cquad%5Chat%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%5Cne0%5C%5C%26%5Crightarrow%5Ctext%7Bif%20%7D%5C%3B%5Clambda_2%3C%5Cdelta%5Ccdot%7C%7Cr_j%5E%7B%28g%29%7D%7C%7C_2%5Ccdot%7C%7C%5Ctilde%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%7C%7C_2%2C%5Cquad%5Chat%7B%5Ctheta%7D_j%5E%7B%28g%29%7D%5Cne0%5Cend%7Baligned%7D
"\\begin{aligned}&\\rightarrow\\text{if }\\;(1-\\frac{\\lambda_2\\cdot\\phi_j^{(g)}/\\delta}{||r_j^{(g)}||_2})\>0,\\quad\\hat{\\theta}_j^{(g)}\\ne0\\\\&\\rightarrow\\text{if }\\;\\lambda_2\<\\delta\\cdot||r_j^{(g)}||_2\\cdot||\\tilde{\\theta}_j^{(g)}||_2,\\quad\\hat{\\theta}_j^{(g)}\\ne0\\end{aligned}")  

  
![\\begin{aligned}\\frac{\\partial^2 Q}{\\partial
\\theta\_j^{2(g)T}}&=\\underbrace{\\delta+\\begin{cases}1/||\\theta\_j||\_2&\\quad\\text{if
} \\theta\_j \\ne 0\\\\0&\\quad\\text{if }
\\theta\_j=0\\end{cases}}\_{\>0}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cfrac%7B%5Cpartial%5E2%20Q%7D%7B%5Cpartial%20%5Ctheta_j%5E%7B2%28g%29T%7D%7D%26%3D%5Cunderbrace%7B%5Cdelta%2B%5Cbegin%7Bcases%7D1%2F%7C%7C%5Ctheta_j%7C%7C_2%26%5Cquad%5Ctext%7Bif%20%7D%20%5Ctheta_j%20%5Cne%200%5C%5C0%26%5Cquad%5Ctext%7Bif%20%7D%20%5Ctheta_j%3D0%5Cend%7Bcases%7D%7D_%7B%3E0%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\frac{\\partial^2 Q}{\\partial \\theta_j^{2(g)T}}&=\\underbrace{\\delta+\\begin{cases}1/||\\theta_j||_2&\\quad\\text{if } \\theta_j \\ne 0\\\\0&\\quad\\text{if } \\theta_j=0\\end{cases}}_{\>0}\\end{aligned}")  

## 3\. Closed form for alpha

  - Vector form

  
![\\begin{aligned}\\boldsymbol{X\\alpha}^{(g)^{k+1}}:=argmin\_{\\alpha\\in\\mathbb{R}^{p+1}}&-\\sum\_{g=1}^m\\sum\_{\\ell=1}^b\\boldsymbol{u}^{(\\ell)(g)^{k}T}(\\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)})+\\lambda\_1||\\boldsymbol{XA}||\_{\*\\boldsymbol{\\psi}}\\\\&+\\frac{\\delta}{2}\\sum\_{g=1}^m\\sum\_{\\ell=1}^b||\\boldsymbol{Y}^{(g)}-\\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)^{k+1}}-\\boldsymbol{e}^{(\\ell)(g)^k}||\_2^2\\\\&\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7BX%5Calpha%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%3A%3Dargmin_%7B%5Calpha%5Cin%5Cmathbb%7BR%7D%5E%7Bp%2B1%7D%7D%26-%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7B%5Cell%3D1%7D%5Eb%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28g%29%5E%7Bk%7DT%7D%28%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D%29%2B%5Clambda_1%7C%7C%5Cboldsymbol%7BXA%7D%7C%7C_%7B%2A%5Cboldsymbol%7B%5Cpsi%7D%7D%5C%5C%26%2B%5Cfrac%7B%5Cdelta%7D%7B2%7D%5Csum_%7Bg%3D1%7D%5Em%5Csum_%7B%5Cell%3D1%7D%5Eb%7C%7C%5Cboldsymbol%7BY%7D%5E%7B%28g%29%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D-%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%28g%29%5Ek%7D%7C%7C_2%5E2%5C%5C%26%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{X\\alpha}^{(g)^{k+1}}:=argmin_{\\alpha\\in\\mathbb{R}^{p+1}}&-\\sum_{g=1}^m\\sum_{\\ell=1}^b\\boldsymbol{u}^{(\\ell)(g)^{k}T}(\\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)})+\\lambda_1||\\boldsymbol{XA}||_{*\\boldsymbol{\\psi}}\\\\&+\\frac{\\delta}{2}\\sum_{g=1}^m\\sum_{\\ell=1}^b||\\boldsymbol{Y}^{(g)}-\\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)^{k+1}}-\\boldsymbol{e}^{(\\ell)(g)^k}||_2^2\\\\&\\end{aligned}")  

where
![\\lambda\_1||\\cdot||\_{\*\\boldsymbol{\\psi}}](https://latex.codecogs.com/png.latex?%5Clambda_1%7C%7C%5Ccdot%7C%7C_%7B%2A%5Cboldsymbol%7B%5Cpsi%7D%7D
"\\lambda_1||\\cdot||_{*\\boldsymbol{\\psi}}") is adaptive nuclear norm.

  - Matrix form

  
![\\begin{aligned}&(\\boldsymbol{XA})^{^{k+1}}:=argmin\_{\\boldsymbol{XA}}-\\sum\_\\ell^b\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})+\\frac{\\delta}{2}\\sum\_{\\ell=1}^b||\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{A}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}||\_F^2+\\lambda\_1||\\boldsymbol{XA}||\_{\*\\boldsymbol{\\psi}}\\\\&\\text{where}\\;\\boldsymbol{U}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{u}^{(\\ell)(1)}\\dots\\boldsymbol{u}^{(\\ell)(m)}\\end{bmatrix},\\quad\\boldsymbol{Y}=\\begin{bmatrix}\\boldsymbol{Y}^{(1)}\\dots\\boldsymbol{Y}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{H}=\\begin{bmatrix}\\boldsymbol{\\eta}^{(1)}\\dots\\boldsymbol{\\eta}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{E}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{e}^{(\\ell)(1)}\\dots\\boldsymbol{e}^{(\\ell)(m)}\\end{bmatrix}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%28%5Cboldsymbol%7BXA%7D%29%5E%7B%5E%7Bk%2B1%7D%7D%3A%3Dargmin_%7B%5Cboldsymbol%7BXA%7D%7D-%5Csum_%5Cell%5Eb%5Ctext%7Btr%7D%28%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%5ET%7D%5Cboldsymbol%7BXA%7D%29%2B%5Cfrac%7B%5Cdelta%7D%7B2%7D%5Csum_%7B%5Cell%3D1%7D%5Eb%7C%7C%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7BA%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D%7C%7C_F%5E2%2B%5Clambda_1%7C%7C%5Cboldsymbol%7BXA%7D%7C%7C_%7B%2A%5Cboldsymbol%7B%5Cpsi%7D%7D%5C%5C%26%5Ctext%7Bwhere%7D%5C%3B%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%281%29%7D%5Cdots%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28m%29%7D%5Cend%7Bbmatrix%7D%2C%5Cquad%5Cboldsymbol%7BY%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7BY%7D%5E%7B%281%29%7D%5Cdots%5Cboldsymbol%7BY%7D%5E%7B%28m%29%7D%5Cend%7Bbmatrix%7D%2C%5Cquad%5Cboldsymbol%7BH%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7B%5Ceta%7D%5E%7B%281%29%7D%5Cdots%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28m%29%7D%5Cend%7Bbmatrix%7D%2C%5Cquad%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%281%29%7D%5Cdots%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%28m%29%7D%5Cend%7Bbmatrix%7D%5Cend%7Baligned%7D
"\\begin{aligned}&(\\boldsymbol{XA})^{^{k+1}}:=argmin_{\\boldsymbol{XA}}-\\sum_\\ell^b\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})+\\frac{\\delta}{2}\\sum_{\\ell=1}^b||\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{A}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}||_F^2+\\lambda_1||\\boldsymbol{XA}||_{*\\boldsymbol{\\psi}}\\\\&\\text{where}\\;\\boldsymbol{U}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{u}^{(\\ell)(1)}\\dots\\boldsymbol{u}^{(\\ell)(m)}\\end{bmatrix},\\quad\\boldsymbol{Y}=\\begin{bmatrix}\\boldsymbol{Y}^{(1)}\\dots\\boldsymbol{Y}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{H}=\\begin{bmatrix}\\boldsymbol{\\eta}^{(1)}\\dots\\boldsymbol{\\eta}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{E}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{e}^{(\\ell)(1)}\\dots\\boldsymbol{e}^{(\\ell)(m)}\\end{bmatrix}\\end{aligned}")  

-----

  
![\\begin{aligned}&\\sum\_{\\ell=1}^b\\{\\frac{\\delta}{2}||\\boldsymbol{Y}-\\underbrace{\\boldsymbol{XA}}\_{\\boldsymbol{Z}}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}||\_F^2-\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})\\}\\\\&=\\sum\_{\\ell}^b\\text{tr}\\bigg\\{\\frac{\\delta}{2}(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})^T(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})-(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{Z})\\bigg\\}\\\\&=\\text{tr}\\Bigg\\{\\frac{\\delta}{2}\\sum\_{\\ell}^b\\bigg\[(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})^T(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})-\\frac{2}{\\delta}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{Z})\\bigg\]\\Bigg\\}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%5Csum_%7B%5Cell%3D1%7D%5Eb%5C%7B%5Cfrac%7B%5Cdelta%7D%7B2%7D%7C%7C%5Cboldsymbol%7BY%7D-%5Cunderbrace%7B%5Cboldsymbol%7BXA%7D%7D_%7B%5Cboldsymbol%7BZ%7D%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D%7C%7C_F%5E2-%5Ctext%7Btr%7D%28%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%5ET%7D%5Cboldsymbol%7BXA%7D%29%5C%7D%5C%5C%26%3D%5Csum_%7B%5Cell%7D%5Eb%5Ctext%7Btr%7D%5Cbigg%5C%7B%5Cfrac%7B%5Cdelta%7D%7B2%7D%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cboldsymbol%7BZ%7D%29%5ET%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cboldsymbol%7BZ%7D%29-%28%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%5ET%7D%5Cboldsymbol%7BZ%7D%29%5Cbigg%5C%7D%5C%5C%26%3D%5Ctext%7Btr%7D%5CBigg%5C%7B%5Cfrac%7B%5Cdelta%7D%7B2%7D%5Csum_%7B%5Cell%7D%5Eb%5Cbigg%5B%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cboldsymbol%7BZ%7D%29%5ET%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cboldsymbol%7BZ%7D%29-%5Cfrac%7B2%7D%7B%5Cdelta%7D%28%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%5ET%7D%5Cboldsymbol%7BZ%7D%29%5Cbigg%5D%5CBigg%5C%7D%5Cend%7Baligned%7D
"\\begin{aligned}&\\sum_{\\ell=1}^b\\{\\frac{\\delta}{2}||\\boldsymbol{Y}-\\underbrace{\\boldsymbol{XA}}_{\\boldsymbol{Z}}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}||_F^2-\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})\\}\\\\&=\\sum_{\\ell}^b\\text{tr}\\bigg\\{\\frac{\\delta}{2}(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})^T(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})-(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{Z})\\bigg\\}\\\\&=\\text{tr}\\Bigg\\{\\frac{\\delta}{2}\\sum_{\\ell}^b\\bigg[(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})^T(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\boldsymbol{Z})-\\frac{2}{\\delta}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{Z})\\bigg]\\Bigg\\}\\end{aligned}")  

  
![\\begin{aligned}&\\propto\\text{tr}\\Bigg\\{\\frac{\\delta\\cdot
b}{2}\\bigg\[\\boldsymbol{Z}^T\\boldsymbol{Z}-\\frac{2}{b}\\sum\_{\\ell}^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})^T\\boldsymbol{Z}\\bigg\]\\Bigg\\}\\\\&\\propto\\frac{\\delta\\cdot
b}{2}\\text{tr}\\Bigg\\{\\bigg\[\\frac{1}{b}\\sum\_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}\\bigg\]^T\\bigg\[\\frac{1}{b}\\sum\_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}\\bigg\]\\Bigg\\}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%5Cpropto%5Ctext%7Btr%7D%5CBigg%5C%7B%5Cfrac%7B%5Cdelta%5Ccdot%20b%7D%7B2%7D%5Cbigg%5B%5Cboldsymbol%7BZ%7D%5ET%5Cboldsymbol%7BZ%7D-%5Cfrac%7B2%7D%7Bb%7D%5Csum_%7B%5Cell%7D%5Eb%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%29%5ET%5Cboldsymbol%7BZ%7D%5Cbigg%5D%5CBigg%5C%7D%5C%5C%26%5Cpropto%5Cfrac%7B%5Cdelta%5Ccdot%20b%7D%7B2%7D%5Ctext%7Btr%7D%5CBigg%5C%7B%5Cbigg%5B%5Cfrac%7B1%7D%7Bb%7D%5Csum_%5Cell%5Eb%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%29-%5Cboldsymbol%7BZ%7D%5Cbigg%5D%5ET%5Cbigg%5B%5Cfrac%7B1%7D%7Bb%7D%5Csum_%5Cell%5Eb%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%29-%5Cboldsymbol%7BZ%7D%5Cbigg%5D%5CBigg%5C%7D%5Cend%7Baligned%7D
"\\begin{aligned}&\\propto\\text{tr}\\Bigg\\{\\frac{\\delta\\cdot b}{2}\\bigg[\\boldsymbol{Z}^T\\boldsymbol{Z}-\\frac{2}{b}\\sum_{\\ell}^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})^T\\boldsymbol{Z}\\bigg]\\Bigg\\}\\\\&\\propto\\frac{\\delta\\cdot b}{2}\\text{tr}\\Bigg\\{\\bigg[\\frac{1}{b}\\sum_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}\\bigg]^T\\bigg[\\frac{1}{b}\\sum_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}\\bigg]\\Bigg\\}\\end{aligned}")  

  
![\\begin{aligned}&=\\frac{\\delta\\cdot
b}{2}||\\frac{1}{b}\\sum\_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}||\_F^2\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%3D%5Cfrac%7B%5Cdelta%5Ccdot%20b%7D%7B2%7D%7C%7C%5Cfrac%7B1%7D%7Bb%7D%5Csum_%5Cell%5Eb%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%29-%5Cboldsymbol%7BZ%7D%7C%7C_F%5E2%5Cend%7Baligned%7D
"\\begin{aligned}&=\\frac{\\delta\\cdot b}{2}||\\frac{1}{b}\\sum_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}||_F^2\\end{aligned}")  

-----

  
![\\begin{aligned}\\therefore\\quad&(\\boldsymbol{XA})^{k+1}(=\\boldsymbol{Z}^{k+1})\\\\&\\quad\\quad:=\\text{argmin}\_{Z}\\frac{\\delta\\cdot
b}{2}||\\frac{1}{b}\\sum\_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}||\_F^2+\\lambda\_1||\\boldsymbol{Z}||\_{\*\\boldsymbol{\\psi}}\\\\&\\quad\\quad=\\text{argmin}\_{Z}\\frac{1}{2}||\\frac{1}{b}\\sum\_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}||\_F^2+\\frac{\\lambda\_1}{\\delta\\cdot
b}||\\boldsymbol{Z}||\_{\*\\boldsymbol{\\psi}}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Ctherefore%5Cquad%26%28%5Cboldsymbol%7BXA%7D%29%5E%7Bk%2B1%7D%28%3D%5Cboldsymbol%7BZ%7D%5E%7Bk%2B1%7D%29%5C%5C%26%5Cquad%5Cquad%3A%3D%5Ctext%7Bargmin%7D_%7BZ%7D%5Cfrac%7B%5Cdelta%5Ccdot%20b%7D%7B2%7D%7C%7C%5Cfrac%7B1%7D%7Bb%7D%5Csum_%5Cell%5Eb%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%29-%5Cboldsymbol%7BZ%7D%7C%7C_F%5E2%2B%5Clambda_1%7C%7C%5Cboldsymbol%7BZ%7D%7C%7C_%7B%2A%5Cboldsymbol%7B%5Cpsi%7D%7D%5C%5C%26%5Cquad%5Cquad%3D%5Ctext%7Bargmin%7D_%7BZ%7D%5Cfrac%7B1%7D%7B2%7D%7C%7C%5Cfrac%7B1%7D%7Bb%7D%5Csum_%5Cell%5Eb%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%29-%5Cboldsymbol%7BZ%7D%7C%7C_F%5E2%2B%5Cfrac%7B%5Clambda_1%7D%7B%5Cdelta%5Ccdot%20b%7D%7C%7C%5Cboldsymbol%7BZ%7D%7C%7C_%7B%2A%5Cboldsymbol%7B%5Cpsi%7D%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\therefore\\quad&(\\boldsymbol{XA})^{k+1}(=\\boldsymbol{Z}^{k+1})\\\\&\\quad\\quad:=\\text{argmin}_{Z}\\frac{\\delta\\cdot b}{2}||\\frac{1}{b}\\sum_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}||_F^2+\\lambda_1||\\boldsymbol{Z}||_{*\\boldsymbol{\\psi}}\\\\&\\quad\\quad=\\text{argmin}_{Z}\\frac{1}{2}||\\frac{1}{b}\\sum_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})-\\boldsymbol{Z}||_F^2+\\frac{\\lambda_1}{\\delta\\cdot b}||\\boldsymbol{Z}||_{*\\boldsymbol{\\psi}}\\end{aligned}")  

  
![\\begin{aligned}&\\quad\\quad=\\boldsymbol{P}\\boldsymbol{D}\_{\\lambda\_1/\\delta\\cdot
b,\\boldsymbol{\\psi}}(\\Sigma)\\boldsymbol{Q}^T\\\\&\\text{where
}\\boldsymbol{P\\Sigma Q}^T\\text{is the SVD of
}\\frac{1}{b}\\sum\_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})\\\\&\\text{and
}\\boldsymbol{D}\_{\\lambda\_1/\\delta\\cdot
b,\\boldsymbol{\\psi}}(\\Sigma)=\\text{diag}(\\{\\sigma\_i-\\frac{\\lambda\_1\\boldsymbol{\\psi}\_i}{(\\delta\\cdot
b)}\\}\_+)\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%5Cquad%5Cquad%3D%5Cboldsymbol%7BP%7D%5Cboldsymbol%7BD%7D_%7B%5Clambda_1%2F%5Cdelta%5Ccdot%20b%2C%5Cboldsymbol%7B%5Cpsi%7D%7D%28%5CSigma%29%5Cboldsymbol%7BQ%7D%5ET%5C%5C%26%5Ctext%7Bwhere%20%7D%5Cboldsymbol%7BP%5CSigma%20Q%7D%5ET%5Ctext%7Bis%20the%20SVD%20of%20%7D%5Cfrac%7B1%7D%7Bb%7D%5Csum_%5Cell%5Eb%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%29%5C%5C%26%5Ctext%7Band%20%7D%5Cboldsymbol%7BD%7D_%7B%5Clambda_1%2F%5Cdelta%5Ccdot%20b%2C%5Cboldsymbol%7B%5Cpsi%7D%7D%28%5CSigma%29%3D%5Ctext%7Bdiag%7D%28%5C%7B%5Csigma_i-%5Cfrac%7B%5Clambda_1%5Cboldsymbol%7B%5Cpsi%7D_i%7D%7B%28%5Cdelta%5Ccdot%20b%29%7D%5C%7D_%2B%29%5Cend%7Baligned%7D
"\\begin{aligned}&\\quad\\quad=\\boldsymbol{P}\\boldsymbol{D}_{\\lambda_1/\\delta\\cdot b,\\boldsymbol{\\psi}}(\\Sigma)\\boldsymbol{Q}^T\\\\&\\text{where }\\boldsymbol{P\\Sigma Q}^T\\text{is the SVD of }\\frac{1}{b}\\sum_\\ell^b(\\boldsymbol{Y}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}-\\frac{1}{\\delta}\\boldsymbol{U}^{(\\ell)})\\\\&\\text{and }\\boldsymbol{D}_{\\lambda_1/\\delta\\cdot b,\\boldsymbol{\\psi}}(\\Sigma)=\\text{diag}(\\{\\sigma_i-\\frac{\\lambda_1\\boldsymbol{\\psi}_i}{(\\delta\\cdot b)}\\}_+)\\end{aligned}")  

where
![\\boldsymbol{\\psi}\_i=\\frac{1}{\\tilde{\\sigma}\_i}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%7B%5Cpsi%7D_i%3D%5Cfrac%7B1%7D%7B%5Ctilde%7B%5Csigma%7D_i%7D
"\\boldsymbol{\\psi}_i=\\frac{1}{\\tilde{\\sigma}_i}")

  
![\\begin{aligned}&\\rightarrow\\text{if
}\\;\\sigma\_i-\\frac{\\lambda\_1\\boldsymbol{\\psi}\_i}{(\\delta\\cdot
b)}\>0,\\quad\\textbf{Z}\\ne0\\\\&\\rightarrow\\text{if
}\\;\\lambda\_1\<\\delta\\cdot
b\\cdot\\sigma\_i\\cdot\\tilde{\\sigma}\_i,\\quad\\textbf{Z}\\ne0\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%5Crightarrow%5Ctext%7Bif%20%7D%5C%3B%5Csigma_i-%5Cfrac%7B%5Clambda_1%5Cboldsymbol%7B%5Cpsi%7D_i%7D%7B%28%5Cdelta%5Ccdot%20b%29%7D%3E0%2C%5Cquad%5Ctextbf%7BZ%7D%5Cne0%5C%5C%26%5Crightarrow%5Ctext%7Bif%20%7D%5C%3B%5Clambda_1%3C%5Cdelta%5Ccdot%20b%5Ccdot%5Csigma_i%5Ccdot%5Ctilde%7B%5Csigma%7D_i%2C%5Cquad%5Ctextbf%7BZ%7D%5Cne0%5Cend%7Baligned%7D
"\\begin{aligned}&\\rightarrow\\text{if }\\;\\sigma_i-\\frac{\\lambda_1\\boldsymbol{\\psi}_i}{(\\delta\\cdot b)}\>0,\\quad\\textbf{Z}\\ne0\\\\&\\rightarrow\\text{if }\\;\\lambda_1\<\\delta\\cdot b\\cdot\\sigma_i\\cdot\\tilde{\\sigma}_i,\\quad\\textbf{Z}\\ne0\\end{aligned}")  

## 4\. Closed form for e

  
![e\_i^{(\\ell)(g)} = \\text{argmin
}\\frac{1}{n}\\rho\_{\\tau\_{\\ell}}(e\_i^{(\\ell)(g)})-u\_i^{(\\ell)(g)}e\_i^{(\\ell)(g)}+\\frac{\\delta}{2}(Y\_i^{(g)}-X\_i\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}-e\_i^{(\\ell)(g)})^2](https://latex.codecogs.com/png.latex?e_i%5E%7B%28%5Cell%29%28g%29%7D%20%3D%20%5Ctext%7Bargmin%20%7D%5Cfrac%7B1%7D%7Bn%7D%5Crho_%7B%5Ctau_%7B%5Cell%7D%7D%28e_i%5E%7B%28%5Cell%29%28g%29%7D%29-u_i%5E%7B%28%5Cell%29%28g%29%7De_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7B%5Cdelta%7D%7B2%7D%28Y_i%5E%7B%28g%29%7D-X_i%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D-e_i%5E%7B%28%5Cell%29%28g%29%7D%29%5E2
"e_i^{(\\ell)(g)} = \\text{argmin }\\frac{1}{n}\\rho_{\\tau_{\\ell}}(e_i^{(\\ell)(g)})-u_i^{(\\ell)(g)}e_i^{(\\ell)(g)}+\\frac{\\delta}{2}(Y_i^{(g)}-X_i\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}-e_i^{(\\ell)(g)})^2")  

  
![\\frac{\\partial Q}{\\partial
e\_i^{(\\ell)(g)}}=-\\delta(Y\_i^{(g)}-X\_i\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}-e\_i^{(\\ell)(g)})-u\_i^{(\\ell)(g)}+\\frac{1}{n}\\frac{\\partial
\\rho\_{\\tau\_{\\ell}}(e\_i^{(\\ell)(g)})}{\\partial
e\_i^{(\\ell)(g)}}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20e_i%5E%7B%28%5Cell%29%28g%29%7D%7D%3D-%5Cdelta%28Y_i%5E%7B%28g%29%7D-X_i%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D-e_i%5E%7B%28%5Cell%29%28g%29%7D%29-u_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7B1%7D%7Bn%7D%5Cfrac%7B%5Cpartial%20%5Crho_%7B%5Ctau_%7B%5Cell%7D%7D%28e_i%5E%7B%28%5Cell%29%28g%29%7D%29%7D%7B%5Cpartial%20e_i%5E%7B%28%5Cell%29%28g%29%7D%7D
"\\frac{\\partial Q}{\\partial e_i^{(\\ell)(g)}}=-\\delta(Y_i^{(g)}-X_i\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}-e_i^{(\\ell)(g)})-u_i^{(\\ell)(g)}+\\frac{1}{n}\\frac{\\partial \\rho_{\\tau_{\\ell}}(e_i^{(\\ell)(g)})}{\\partial e_i^{(\\ell)(g)}}")  

  
![\\frac{\\partial^2 Q}{\\partial
e\_i^{2(\\ell)(g)}}=\\delta\>0](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%5E2%20Q%7D%7B%5Cpartial%20e_i%5E%7B2%28%5Cell%29%28g%29%7D%7D%3D%5Cdelta%3E0
"\\frac{\\partial^2 Q}{\\partial e_i^{2(\\ell)(g)}}=\\delta\>0")  

  
![\\frac{\\partial \\rho\_{\\tau\_{\\ell}}(e\_i^{(\\ell)(g)})}{\\partial
e\_i^{(\\ell)(g)}}=\\begin{cases}\\tau\_\\ell-1\\quad&\\text{if
}e\_i^{(\\ell)(g)}\<0\\\\\\{c\\in\\mathbb{R}:\\tau\_\\ell-1\\le
c\\le\\tau\_\\ell\\}\\quad&\\text{if
}e\_i^{(\\ell)(g)}=0\\\\\\tau\_\\ell\\quad&\\text{if
}e\_i^{(\\ell)(g)}\>0\\end{cases}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20%5Crho_%7B%5Ctau_%7B%5Cell%7D%7D%28e_i%5E%7B%28%5Cell%29%28g%29%7D%29%7D%7B%5Cpartial%20e_i%5E%7B%28%5Cell%29%28g%29%7D%7D%3D%5Cbegin%7Bcases%7D%5Ctau_%5Cell-1%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3C0%5C%5C%5C%7Bc%5Cin%5Cmathbb%7BR%7D%3A%5Ctau_%5Cell-1%5Cle%20c%5Cle%5Ctau_%5Cell%5C%7D%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3D0%5C%5C%5Ctau_%5Cell%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3E0%5Cend%7Bcases%7D
"\\frac{\\partial \\rho_{\\tau_{\\ell}}(e_i^{(\\ell)(g)})}{\\partial e_i^{(\\ell)(g)}}=\\begin{cases}\\tau_\\ell-1\\quad&\\text{if }e_i^{(\\ell)(g)}\<0\\\\\\{c\\in\\mathbb{R}:\\tau_\\ell-1\\le c\\le\\tau_\\ell\\}\\quad&\\text{if }e_i^{(\\ell)(g)}=0\\\\\\tau_\\ell\\quad&\\text{if }e_i^{(\\ell)(g)}\>0\\end{cases}")  

Let
![\\epsilon\_i^{(\\ell)(g)}=Y\_i^{(g)}-X\_i^T\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}](https://latex.codecogs.com/png.latex?%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%3DY_i%5E%7B%28g%29%7D-X_i%5ET%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D
"\\epsilon_i^{(\\ell)(g)}=Y_i^{(g)}-X_i^T\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}")

-----

if
![e\_i^{(\\ell)(g)}\<0](https://latex.codecogs.com/png.latex?e_i%5E%7B%28%5Cell%29%28g%29%7D%3C0
"e_i^{(\\ell)(g)}\<0")

  
![\\begin{aligned}&-(\\epsilon\_i^{(\\ell)(g)}-e\_i^{(\\ell)(g)})-u\_i^{(\\ell)(g)}/\\delta+(\\tau-1)/n\\delta=0\\\\\&e\_i=\\epsilon\_i^{(\\ell)(g)}+\\frac{u\_i^{(\\ell)(g)}}{\\delta}-\\frac{\\tau-1}{n\\delta}\<0\\quad\\text{by
hypo}\\\\&\\leftrightarrow\\;\\epsilon\_i^{(\\ell)(g)}+\\frac{u\_i^{(\\ell)(g)}}{\\delta}\<\\frac{\\tau-1}{n\\delta}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26-%28%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D-e_i%5E%7B%28%5Cell%29%28g%29%7D%29-u_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta%2B%28%5Ctau-1%29%2Fn%5Cdelta%3D0%5C%5C%26e_i%3D%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7Bu_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7B%5Cdelta%7D-%5Cfrac%7B%5Ctau-1%7D%7Bn%5Cdelta%7D%3C0%5Cquad%5Ctext%7Bby%20hypo%7D%5C%5C%26%5Cleftrightarrow%5C%3B%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7Bu_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7B%5Cdelta%7D%3C%5Cfrac%7B%5Ctau-1%7D%7Bn%5Cdelta%7D%5Cend%7Baligned%7D
"\\begin{aligned}&-(\\epsilon_i^{(\\ell)(g)}-e_i^{(\\ell)(g)})-u_i^{(\\ell)(g)}/\\delta+(\\tau-1)/n\\delta=0\\\\&e_i=\\epsilon_i^{(\\ell)(g)}+\\frac{u_i^{(\\ell)(g)}}{\\delta}-\\frac{\\tau-1}{n\\delta}\<0\\quad\\text{by hypo}\\\\&\\leftrightarrow\\;\\epsilon_i^{(\\ell)(g)}+\\frac{u_i^{(\\ell)(g)}}{\\delta}\<\\frac{\\tau-1}{n\\delta}\\end{aligned}")  

-----

if
![e\_i^{(\\ell)(g)}\>0](https://latex.codecogs.com/png.latex?e_i%5E%7B%28%5Cell%29%28g%29%7D%3E0
"e_i^{(\\ell)(g)}\>0")

  
![-(\\epsilon\_i^{(\\ell)(g)}-e\_i^{(\\ell)(g)})-u\_i^{(\\ell)(g)}/\\delta+\\tau/n\\delta=0\\\\e\_i=\\epsilon\_i^{(\\ell)(g)}+\\frac{u\_i^{(\\ell)(g)}}{\\delta}-\\frac{\\tau}{n\\delta}\>0\\quad\\text{by
hypo}\\\\\\leftrightarrow\\;\\epsilon\_i^{(\\ell)(g)}+\\frac{u\_i^{(\\ell)(g)}}{\\delta}\<\\frac{\\tau}{n\\delta}\\\\](https://latex.codecogs.com/png.latex?-%28%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D-e_i%5E%7B%28%5Cell%29%28g%29%7D%29-u_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta%2B%5Ctau%2Fn%5Cdelta%3D0%5C%5Ce_i%3D%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7Bu_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7B%5Cdelta%7D-%5Cfrac%7B%5Ctau%7D%7Bn%5Cdelta%7D%3E0%5Cquad%5Ctext%7Bby%20hypo%7D%5C%5C%5Cleftrightarrow%5C%3B%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7Bu_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7B%5Cdelta%7D%3C%5Cfrac%7B%5Ctau%7D%7Bn%5Cdelta%7D%5C%5C
"-(\\epsilon_i^{(\\ell)(g)}-e_i^{(\\ell)(g)})-u_i^{(\\ell)(g)}/\\delta+\\tau/n\\delta=0\\\\e_i=\\epsilon_i^{(\\ell)(g)}+\\frac{u_i^{(\\ell)(g)}}{\\delta}-\\frac{\\tau}{n\\delta}\>0\\quad\\text{by hypo}\\\\\\leftrightarrow\\;\\epsilon_i^{(\\ell)(g)}+\\frac{u_i^{(\\ell)(g)}}{\\delta}\<\\frac{\\tau}{n\\delta}\\\\")  

-----

if
![e\_i^{(\\ell)(g)}=0](https://latex.codecogs.com/png.latex?e_i%5E%7B%28%5Cell%29%28g%29%7D%3D0
"e_i^{(\\ell)(g)}=0")

  
![-(\\epsilon\_i^{(\\ell)(g)})-u\_i^{(\\ell)(g)}/\\delta+c/n\\delta=0\\quad\\text{where
}-\\tau\_\\ell\\le
c\\le1-\\tau\_\\ell](https://latex.codecogs.com/png.latex?-%28%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%29-u_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta%2Bc%2Fn%5Cdelta%3D0%5Cquad%5Ctext%7Bwhere%20%7D-%5Ctau_%5Cell%5Cle%20c%5Cle1-%5Ctau_%5Cell
"-(\\epsilon_i^{(\\ell)(g)})-u_i^{(\\ell)(g)}/\\delta+c/n\\delta=0\\quad\\text{where }-\\tau_\\ell\\le c\\le1-\\tau_\\ell")  

  
![\\frac{\\tau\_\\ell-1}{n\\delta}\\le-(\\epsilon\_i^{(\\ell)(g)})-\\frac{u\_i^{(\\ell)(g)}}{\\delta}=-\\frac{c}{\\delta}\\le\\frac{\\tau\_\\ell}{n\\delta}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Ctau_%5Cell-1%7D%7Bn%5Cdelta%7D%5Cle-%28%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%29-%5Cfrac%7Bu_i%5E%7B%28%5Cell%29%28g%29%7D%7D%7B%5Cdelta%7D%3D-%5Cfrac%7Bc%7D%7B%5Cdelta%7D%5Cle%5Cfrac%7B%5Ctau_%5Cell%7D%7Bn%5Cdelta%7D
"\\frac{\\tau_\\ell-1}{n\\delta}\\le-(\\epsilon_i^{(\\ell)(g)})-\\frac{u_i^{(\\ell)(g)}}{\\delta}=-\\frac{c}{\\delta}\\le\\frac{\\tau_\\ell}{n\\delta}")  

-----

  
![\\therefore
e\_i^{(\\ell)(g)^{k+1}}=\\begin{cases}\\epsilon\_i^{(\\ell)(g)}+u\_i^{(\\ell)(g)}/\\delta-(\\tau\_\\ell-1)/n\\delta&\\text{if
}\\epsilon\_i^{(\\ell)(g)}+u\_i^{(\\ell)(g)}/\\delta\<(\\tau\_\\ell-1)/n\\delta\\\\0\\quad&\\text{if
}(\\tau\_\\ell-1)/n\\delta\<\\epsilon\_i^{(\\ell)(g)}+u\_i^{(\\ell)(g)}/\\delta\<\\tau\_\\ell/n\\delta\\\\\\epsilon\_i^{(\\ell)(g)}+u\_i^{(\\ell)(g)}/\\delta-\\tau\_\\ell/n\\delta\\quad&\\text{if
}\\epsilon\_i^{(\\ell)(g)}+u\_i^{(\\ell)(g)}/\\delta\>\\tau\_\\ell/n\\delta\\end{cases}](https://latex.codecogs.com/png.latex?%5Ctherefore%20e_i%5E%7B%28%5Cell%29%28g%29%5E%7Bk%2B1%7D%7D%3D%5Cbegin%7Bcases%7D%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2Bu_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta-%28%5Ctau_%5Cell-1%29%2Fn%5Cdelta%26%5Ctext%7Bif%20%7D%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2Bu_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta%3C%28%5Ctau_%5Cell-1%29%2Fn%5Cdelta%5C%5C0%5Cquad%26%5Ctext%7Bif%20%7D%28%5Ctau_%5Cell-1%29%2Fn%5Cdelta%3C%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2Bu_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta%3C%5Ctau_%5Cell%2Fn%5Cdelta%5C%5C%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2Bu_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta-%5Ctau_%5Cell%2Fn%5Cdelta%5Cquad%26%5Ctext%7Bif%20%7D%5Cepsilon_i%5E%7B%28%5Cell%29%28g%29%7D%2Bu_i%5E%7B%28%5Cell%29%28g%29%7D%2F%5Cdelta%3E%5Ctau_%5Cell%2Fn%5Cdelta%5Cend%7Bcases%7D
"\\therefore e_i^{(\\ell)(g)^{k+1}}=\\begin{cases}\\epsilon_i^{(\\ell)(g)}+u_i^{(\\ell)(g)}/\\delta-(\\tau_\\ell-1)/n\\delta&\\text{if }\\epsilon_i^{(\\ell)(g)}+u_i^{(\\ell)(g)}/\\delta\<(\\tau_\\ell-1)/n\\delta\\\\0\\quad&\\text{if }(\\tau_\\ell-1)/n\\delta\<\\epsilon_i^{(\\ell)(g)}+u_i^{(\\ell)(g)}/\\delta\<\\tau_\\ell/n\\delta\\\\\\epsilon_i^{(\\ell)(g)}+u_i^{(\\ell)(g)}/\\delta-\\tau_\\ell/n\\delta\\quad&\\text{if }\\epsilon_i^{(\\ell)(g)}+u_i^{(\\ell)(g)}/\\delta\>\\tau_\\ell/n\\delta\\end{cases}")  

## 5\. Closed form for multipliers

  
![\\begin{aligned}\\boldsymbol{u}^{(\\ell)(g)^{k+1}}:=&\\boldsymbol{u}^{(\\ell)(g)^k}+\\delta(\\boldsymbol{Y}^{(g)}-\\boldsymbol{Z}\[\\quad,g\]-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)^{k+1}}-\\boldsymbol{e}^{(\\ell)(g)^{k+1}})\\\\\\boldsymbol{w}^{(g)^{k+1}}:=&\\boldsymbol{w}^{(g)^{k}}+\\delta(\\boldsymbol{\\theta}^{(g)^{k+1}}-\\boldsymbol{\\eta}^{(g)^{k+1}})\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28g%29%5E%7Bk%2B1%7D%7D%3A%3D%26%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28g%29%5Ek%7D%2B%5Cdelta%28%5Cboldsymbol%7BY%7D%5E%7B%28g%29%7D-%5Cboldsymbol%7BZ%7D%5B%5Cquad%2Cg%5D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D-%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%28g%29%5E%7Bk%2B1%7D%7D%29%5C%5C%5Cboldsymbol%7Bw%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%3A%3D%26%5Cboldsymbol%7Bw%7D%5E%7B%28g%29%5E%7Bk%7D%7D%2B%5Cdelta%28%5Cboldsymbol%7B%5Ctheta%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D-%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%29%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{u}^{(\\ell)(g)^{k+1}}:=&\\boldsymbol{u}^{(\\ell)(g)^k}+\\delta(\\boldsymbol{Y}^{(g)}-\\boldsymbol{Z}[\\quad,g]-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)^{k+1}}-\\boldsymbol{e}^{(\\ell)(g)^{k+1}})\\\\\\boldsymbol{w}^{(g)^{k+1}}:=&\\boldsymbol{w}^{(g)^{k}}+\\delta(\\boldsymbol{\\theta}^{(g)^{k+1}}-\\boldsymbol{\\eta}^{(g)^{k+1}})\\end{aligned}")
