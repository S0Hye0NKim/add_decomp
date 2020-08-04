Closed Form
================
Sohyeon Kim
7/28/2020

## 1\. Closed form for eta.

  
![\\begin{aligned}\\frac{\\partial
Q}{\\partial\\boldsymbol{\\eta}^{(g)T}}&=-\\sum\_{\\ell}^b
\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}^{(\\ell)}-\\delta\\sum\_{\\ell=1}^bV^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}^{(\\ell)})+\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta^{(g)}}-\\boldsymbol{w}+\\delta\\boldsymbol{\\eta}^{(g)}-\\delta\\boldsymbol{\\theta}\\\\&=0\\\\\\hat{\\boldsymbol{\\eta}}^{(g)}&=\\frac{1}{\\delta}(\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\boldsymbol{I})^{-1}(\\boldsymbol{w}+\\delta\\boldsymbol{\\theta}+\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}+\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}))\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29T%7D%7D%26%3D-%5Csum_%7B%5Cell%7D%5Eb%20%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%7D-%5Cdelta%5Csum_%7B%5Cell%3D1%7D%5EbV%5E%7B%28%5Cell%29T%7D%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D-%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%7D%29%2B%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%5E%7B%28g%29%7D%7D-%5Cboldsymbol%7Bw%7D%2B%5Cdelta%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%7D-%5Cdelta%5Cboldsymbol%7B%5Ctheta%7D%5C%5C%26%3D0%5C%5C%5Chat%7B%5Cboldsymbol%7B%5Ceta%7D%7D%5E%7B%28g%29%7D%26%3D%5Cfrac%7B1%7D%7B%5Cdelta%7D%28%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%2B%5Cboldsymbol%7BI%7D%29%5E%7B-1%7D%28%5Cboldsymbol%7Bw%7D%2B%5Cdelta%5Cboldsymbol%7B%5Ctheta%7D%2B%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7Bu%7D%2B%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%28%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D-%5Cboldsymbol%7Be%7D%29%29%5Cend%7Baligned%7D
"\\begin{aligned}\\frac{\\partial Q}{\\partial\\boldsymbol{\\eta}^{(g)T}}&=-\\sum_{\\ell}^b \\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}^{(\\ell)}-\\delta\\sum_{\\ell=1}^bV^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}^{(\\ell)})+\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta^{(g)}}-\\boldsymbol{w}+\\delta\\boldsymbol{\\eta}^{(g)}-\\delta\\boldsymbol{\\theta}\\\\&=0\\\\\\hat{\\boldsymbol{\\eta}}^{(g)}&=\\frac{1}{\\delta}(\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\boldsymbol{I})^{-1}(\\boldsymbol{w}+\\delta\\boldsymbol{\\theta}+\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{u}+\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}(\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{\\alpha}-\\boldsymbol{e}))\\end{aligned}")  

  
![\\frac{\\partial^2
Q}{\\partial\\boldsymbol{\\eta}^2}=\\delta\\sum\_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\delta\\boldsymbol{I}\\succeq0](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%5E2%20Q%7D%7B%5Cpartial%5Cboldsymbol%7B%5Ceta%7D%5E2%7D%3D%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29T%7D%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%2B%5Cdelta%5Cboldsymbol%7BI%7D%5Csucceq0
"\\frac{\\partial^2 Q}{\\partial\\boldsymbol{\\eta}^2}=\\delta\\sum_{\\ell}^b\\boldsymbol{V}^{(\\ell)T}\\boldsymbol{V}^{(\\ell)}+\\delta\\boldsymbol{I}\\succeq0")  

  - Second derivative function is semi-postiive definite?

## 2\. Closed form for theta

  
![\\begin{aligned}\\theta\_j^{(g)^{k+1}}:=\&argmin\_{\\theta^{(g)}}\\sum\_{g=1}^m\\boldsymbol{w}^{(g)^kT}\\theta\_j^{(g)}+
\\frac{\\delta}{2}\\sum\_{g=1}^m||\\theta\_j^{(g)} -
\\eta\_j^{(g)^{k+1}}||\_2^2 +
\\lambda\_2\\sum\_{g=1}^m||\\theta\_j^{(g)}||\_2\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Ctheta_j%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%3A%3D%26argmin_%7B%5Ctheta%5E%7B%28g%29%7D%7D%5Csum_%7Bg%3D1%7D%5Em%5Cboldsymbol%7Bw%7D%5E%7B%28g%29%5EkT%7D%5Ctheta_j%5E%7B%28g%29%7D%2B%20%5Cfrac%7B%5Cdelta%7D%7B2%7D%5Csum_%7Bg%3D1%7D%5Em%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%20-%20%5Ceta_j%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%7C%7C_2%5E2%20%2B%20%5Clambda_2%5Csum_%7Bg%3D1%7D%5Em%7C%7C%5Ctheta_j%5E%7B%28g%29%7D%7C%7C_2%5Cend%7Baligned%7D
"\\begin{aligned}\\theta_j^{(g)^{k+1}}:=&argmin_{\\theta^{(g)}}\\sum_{g=1}^m\\boldsymbol{w}^{(g)^kT}\\theta_j^{(g)}+ \\frac{\\delta}{2}\\sum_{g=1}^m||\\theta_j^{(g)} - \\eta_j^{(g)^{k+1}}||_2^2 + \\lambda_2\\sum_{g=1}^m||\\theta_j^{(g)}||_2\\end{aligned}")  

  
![\\begin{aligned}\\frac{\\partial Q}{\\partial
\\theta\_j^{(g)T}}&=w\_j+\\delta\\theta\_j-\\delta\\eta\_j+\\frac{\\partial
\\lambda\_2||\\theta\_j||\_2}{\\partial\\theta\_j}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20%5Ctheta_j%5E%7B%28g%29T%7D%7D%26%3Dw_j%2B%5Cdelta%5Ctheta_j-%5Cdelta%5Ceta_j%2B%5Cfrac%7B%5Cpartial%20%5Clambda_2%7C%7C%5Ctheta_j%7C%7C_2%7D%7B%5Cpartial%5Ctheta_j%7D%5Cend%7Baligned%7D
"\\begin{aligned}\\frac{\\partial Q}{\\partial \\theta_j^{(g)T}}&=w_j+\\delta\\theta_j-\\delta\\eta_j+\\frac{\\partial \\lambda_2||\\theta_j||_2}{\\partial\\theta_j}\\end{aligned}")  

  
![\\frac{\\partial||\\theta\_j||\_2}{\\partial\\theta\_j}
=\\begin{cases}\\theta\_j/||\\theta\_j||&\\quad \\text{if } \\theta\_j
\\ne 0\\\\\\{u:||u||\\le 1\\}&\\quad\\text{if }
\\theta\_j=0\\end{cases}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%7C%7C%5Ctheta_j%7C%7C_2%7D%7B%5Cpartial%5Ctheta_j%7D%20%3D%5Cbegin%7Bcases%7D%5Ctheta_j%2F%7C%7C%5Ctheta_j%7C%7C%26%5Cquad%20%5Ctext%7Bif%20%7D%20%5Ctheta_j%20%5Cne%200%5C%5C%5C%7Bu%3A%7C%7Cu%7C%7C%5Cle%201%5C%7D%26%5Cquad%5Ctext%7Bif%20%7D%20%5Ctheta_j%3D0%5Cend%7Bcases%7D
"\\frac{\\partial||\\theta_j||_2}{\\partial\\theta_j} =\\begin{cases}\\theta_j/||\\theta_j||&\\quad \\text{if } \\theta_j \\ne 0\\\\\\{u:||u||\\le 1\\}&\\quad\\text{if } \\theta_j=0\\end{cases}")  

  
![\\theta\_j^{(g)^{k+1}}=\\begin{cases}\\eta\_j^{(g)}-w^{(g)}\_j/\\delta+\\lambda\_2/\\delta&\\quad
\\text{if
}w\_j^{(g)}/\\delta-\\eta\_j^{(g)}\\le-\\lambda\_2/\\delta\\\\0&\\quad\\text{if
}||w\_j^{(g)}/\\delta-\\eta\_j^{(g)}||\\le\\lambda\_2/\\delta\\\\\\eta\_j^{(g)}-w\_j^{(g)}/\\delta-\\lambda\_2/\\delta&\\quad
\\text{if
}w\_j^{(g)}/\\delta-\\eta\_j^{(g)}\\ge\\lambda\_2/\\delta\\end{cases}](https://latex.codecogs.com/png.latex?%5Ctheta_j%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%3D%5Cbegin%7Bcases%7D%5Ceta_j%5E%7B%28g%29%7D-w%5E%7B%28g%29%7D_j%2F%5Cdelta%2B%5Clambda_2%2F%5Cdelta%26%5Cquad%20%5Ctext%7Bif%20%7Dw_j%5E%7B%28g%29%7D%2F%5Cdelta-%5Ceta_j%5E%7B%28g%29%7D%5Cle-%5Clambda_2%2F%5Cdelta%5C%5C0%26%5Cquad%5Ctext%7Bif%20%7D%7C%7Cw_j%5E%7B%28g%29%7D%2F%5Cdelta-%5Ceta_j%5E%7B%28g%29%7D%7C%7C%5Cle%5Clambda_2%2F%5Cdelta%5C%5C%5Ceta_j%5E%7B%28g%29%7D-w_j%5E%7B%28g%29%7D%2F%5Cdelta-%5Clambda_2%2F%5Cdelta%26%5Cquad%20%5Ctext%7Bif%20%7Dw_j%5E%7B%28g%29%7D%2F%5Cdelta-%5Ceta_j%5E%7B%28g%29%7D%5Cge%5Clambda_2%2F%5Cdelta%5Cend%7Bcases%7D
"\\theta_j^{(g)^{k+1}}=\\begin{cases}\\eta_j^{(g)}-w^{(g)}_j/\\delta+\\lambda_2/\\delta&\\quad \\text{if }w_j^{(g)}/\\delta-\\eta_j^{(g)}\\le-\\lambda_2/\\delta\\\\0&\\quad\\text{if }||w_j^{(g)}/\\delta-\\eta_j^{(g)}||\\le\\lambda_2/\\delta\\\\\\eta_j^{(g)}-w_j^{(g)}/\\delta-\\lambda_2/\\delta&\\quad \\text{if }w_j^{(g)}/\\delta-\\eta_j^{(g)}\\ge\\lambda_2/\\delta\\end{cases}")  

## 3\. Closed form for alpha

  
![\\begin{aligned}&\\boldsymbol{A}^{^{k+1}}:=argmin\_{\\boldsymbol{A}}-\\sum\_\\ell^b\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})+\\frac{\\delta}{2}\\sum\_{\\ell=1}^b||\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{A}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}||\_2^2+||\\boldsymbol{A}||\_\*\\\\&\\text{where}\\;\\boldsymbol{U}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{u}^{(\\ell)(1)}\\dots\\boldsymbol{u}^{(\\ell)(m)}\\end{bmatrix},\\quad\\boldsymbol{Y}=\\begin{bmatrix}\\boldsymbol{Y}^{(1)}\\dots\\boldsymbol{Y}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{H}=\\begin{bmatrix}\\boldsymbol{\\eta}^{(1)}\\dots\\boldsymbol{\\eta}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{E}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{e}^{(\\ell)(1)}\\dots\\boldsymbol{e}^{(\\ell)(m)}\\end{bmatrix}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%5Cboldsymbol%7BA%7D%5E%7B%5E%7Bk%2B1%7D%7D%3A%3Dargmin_%7B%5Cboldsymbol%7BA%7D%7D-%5Csum_%5Cell%5Eb%5Ctext%7Btr%7D%28%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%5ET%7D%5Cboldsymbol%7BXA%7D%29%2B%5Cfrac%7B%5Cdelta%7D%7B2%7D%5Csum_%7B%5Cell%3D1%7D%5Eb%7C%7C%5Cboldsymbol%7BY%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7BA%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%5E%7Bk%2B1%7D-%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%5Ek%7D%7C%7C_2%5E2%2B%7C%7C%5Cboldsymbol%7BA%7D%7C%7C_%2A%5C%5C%26%5Ctext%7Bwhere%7D%5C%3B%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%281%29%7D%5Cdots%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28m%29%7D%5Cend%7Bbmatrix%7D%2C%5Cquad%5Cboldsymbol%7BY%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7BY%7D%5E%7B%281%29%7D%5Cdots%5Cboldsymbol%7BY%7D%5E%7B%28m%29%7D%5Cend%7Bbmatrix%7D%2C%5Cquad%5Cboldsymbol%7BH%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7B%5Ceta%7D%5E%7B%281%29%7D%5Cdots%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28m%29%7D%5Cend%7Bbmatrix%7D%2C%5Cquad%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%7D%3D%5Cbegin%7Bbmatrix%7D%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%281%29%7D%5Cdots%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%28m%29%7D%5Cend%7Bbmatrix%7D%5Cend%7Baligned%7D
"\\begin{aligned}&\\boldsymbol{A}^{^{k+1}}:=argmin_{\\boldsymbol{A}}-\\sum_\\ell^b\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})+\\frac{\\delta}{2}\\sum_{\\ell=1}^b||\\boldsymbol{Y}-\\boldsymbol{X}\\boldsymbol{A}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H}^{k+1}-\\boldsymbol{E}^{(\\ell)^k}||_2^2+||\\boldsymbol{A}||_*\\\\&\\text{where}\\;\\boldsymbol{U}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{u}^{(\\ell)(1)}\\dots\\boldsymbol{u}^{(\\ell)(m)}\\end{bmatrix},\\quad\\boldsymbol{Y}=\\begin{bmatrix}\\boldsymbol{Y}^{(1)}\\dots\\boldsymbol{Y}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{H}=\\begin{bmatrix}\\boldsymbol{\\eta}^{(1)}\\dots\\boldsymbol{\\eta}^{(m)}\\end{bmatrix},\\quad\\boldsymbol{E}^{(\\ell)}=\\begin{bmatrix}\\boldsymbol{e}^{(\\ell)(1)}\\dots\\boldsymbol{e}^{(\\ell)(m)}\\end{bmatrix}\\end{aligned}")  

  
![\\begin{aligned}&\\frac{\\partial\\,\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})}{\\partial\\boldsymbol{A}}=\\frac{\\partial\\,\\text{tr}(\\boldsymbol{AU}^{(\\ell)^T}\\boldsymbol{X})}{\\partial\\boldsymbol{A}}=\\boldsymbol{X}^T\\boldsymbol{U}^{(\\ell)}\\\\&\\partial||\\boldsymbol{A}||\_\*=\\{\\boldsymbol{B}\\boldsymbol{C}^T+\\boldsymbol{W}:\\boldsymbol{W}\\in\\mathbb{R}^{(p+1)\\times
n},\\,\\boldsymbol{U}^T\\boldsymbol{W}=0,\\,\\boldsymbol{W}\\boldsymbol{V}=0,\\,||\\boldsymbol{W}||\_2\\le1\\}\\\\&\\text{where}\\quad\\boldsymbol{A}=\\boldsymbol{B\\Sigma
C}^T\\;\\text{by SVD},\\quad||\\cdot||\_2:\\text{spectral norm(largest
singular
value)}\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%26%5Cfrac%7B%5Cpartial%5C%2C%5Ctext%7Btr%7D%28%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%5ET%7D%5Cboldsymbol%7BXA%7D%29%7D%7B%5Cpartial%5Cboldsymbol%7BA%7D%7D%3D%5Cfrac%7B%5Cpartial%5C%2C%5Ctext%7Btr%7D%28%5Cboldsymbol%7BAU%7D%5E%7B%28%5Cell%29%5ET%7D%5Cboldsymbol%7BX%7D%29%7D%7B%5Cpartial%5Cboldsymbol%7BA%7D%7D%3D%5Cboldsymbol%7BX%7D%5ET%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%5C%5C%26%5Cpartial%7C%7C%5Cboldsymbol%7BA%7D%7C%7C_%2A%3D%5C%7B%5Cboldsymbol%7BB%7D%5Cboldsymbol%7BC%7D%5ET%2B%5Cboldsymbol%7BW%7D%3A%5Cboldsymbol%7BW%7D%5Cin%5Cmathbb%7BR%7D%5E%7B%28p%2B1%29%5Ctimes%20n%7D%2C%5C%2C%5Cboldsymbol%7BU%7D%5ET%5Cboldsymbol%7BW%7D%3D0%2C%5C%2C%5Cboldsymbol%7BW%7D%5Cboldsymbol%7BV%7D%3D0%2C%5C%2C%7C%7C%5Cboldsymbol%7BW%7D%7C%7C_2%5Cle1%5C%7D%5C%5C%26%5Ctext%7Bwhere%7D%5Cquad%5Cboldsymbol%7BA%7D%3D%5Cboldsymbol%7BB%5CSigma%20C%7D%5ET%5C%3B%5Ctext%7Bby%20%20SVD%7D%2C%5Cquad%7C%7C%5Ccdot%7C%7C_2%3A%5Ctext%7Bspectral%20norm%28largest%20singular%20value%29%7D%5Cend%7Baligned%7D
"\\begin{aligned}&\\frac{\\partial\\,\\text{tr}(\\boldsymbol{U}^{(\\ell)^T}\\boldsymbol{XA})}{\\partial\\boldsymbol{A}}=\\frac{\\partial\\,\\text{tr}(\\boldsymbol{AU}^{(\\ell)^T}\\boldsymbol{X})}{\\partial\\boldsymbol{A}}=\\boldsymbol{X}^T\\boldsymbol{U}^{(\\ell)}\\\\&\\partial||\\boldsymbol{A}||_*=\\{\\boldsymbol{B}\\boldsymbol{C}^T+\\boldsymbol{W}:\\boldsymbol{W}\\in\\mathbb{R}^{(p+1)\\times n},\\,\\boldsymbol{U}^T\\boldsymbol{W}=0,\\,\\boldsymbol{W}\\boldsymbol{V}=0,\\,||\\boldsymbol{W}||_2\\le1\\}\\\\&\\text{where}\\quad\\boldsymbol{A}=\\boldsymbol{B\\Sigma C}^T\\;\\text{by  SVD},\\quad||\\cdot||_2:\\text{spectral norm(largest singular value)}\\end{aligned}")  

  - convert vector form to matrix form\!

  
![\\begin{aligned} \\frac{\\partial Q}{\\partial \\boldsymbol{A}} &=
-\\sum\_{\\ell}^b\\boldsymbol{X}^T\\boldsymbol{U}^{(\\ell)}
-\\delta\\boldsymbol{X}^T\\sum\_{\\ell}^b(\\boldsymbol{Y} -
\\boldsymbol{V}^{(\\ell)}\\boldsymbol{H} - \\boldsymbol{E}^{(\\ell)}) +
\\delta\\sum\_{\\ell}^b\\boldsymbol{X}^T\\boldsymbol{X}\\boldsymbol{A}
+\\boldsymbol{BC}^T+\\boldsymbol{W}\\\\&=0\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%20%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20%5Cboldsymbol%7BA%7D%7D%20%26%3D%20-%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BX%7D%5ET%5Cboldsymbol%7BU%7D%5E%7B%28%5Cell%29%7D%20-%5Cdelta%5Cboldsymbol%7BX%7D%5ET%5Csum_%7B%5Cell%7D%5Eb%28%5Cboldsymbol%7BY%7D%20-%20%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7BH%7D%20-%20%5Cboldsymbol%7BE%7D%5E%7B%28%5Cell%29%7D%29%20%2B%20%5Cdelta%5Csum_%7B%5Cell%7D%5Eb%5Cboldsymbol%7BX%7D%5ET%5Cboldsymbol%7BX%7D%5Cboldsymbol%7BA%7D%20%2B%5Cboldsymbol%7BBC%7D%5ET%2B%5Cboldsymbol%7BW%7D%5C%5C%26%3D0%5Cend%7Baligned%7D
"\\begin{aligned} \\frac{\\partial Q}{\\partial \\boldsymbol{A}} &= -\\sum_{\\ell}^b\\boldsymbol{X}^T\\boldsymbol{U}^{(\\ell)} -\\delta\\boldsymbol{X}^T\\sum_{\\ell}^b(\\boldsymbol{Y} - \\boldsymbol{V}^{(\\ell)}\\boldsymbol{H} - \\boldsymbol{E}^{(\\ell)}) + \\delta\\sum_{\\ell}^b\\boldsymbol{X}^T\\boldsymbol{X}\\boldsymbol{A} +\\boldsymbol{BC}^T+\\boldsymbol{W}\\\\&=0\\end{aligned}")  

  - optimize the equation through singular value thresholding

## 4\. Closed form for e

  
![e\_i^{(\\ell)(g)} = \\text{argmin
}\\rho\_{\\tau\_{\\ell}}(e\_i^{(\\ell)(g)})-u\_i^{(\\ell)(g)}e\_i^{(\\ell)(g)}+\\frac{\\delta}{2}(Y\_i^{(g)}-X\_i\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}-e\_i^{(\\ell)(g)})^2](https://latex.codecogs.com/png.latex?e_i%5E%7B%28%5Cell%29%28g%29%7D%20%3D%20%5Ctext%7Bargmin%20%7D%5Crho_%7B%5Ctau_%7B%5Cell%7D%7D%28e_i%5E%7B%28%5Cell%29%28g%29%7D%29-u_i%5E%7B%28%5Cell%29%28g%29%7De_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7B%5Cdelta%7D%7B2%7D%28Y_i%5E%7B%28g%29%7D-X_i%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D-e_i%5E%7B%28%5Cell%29%28g%29%7D%29%5E2
"e_i^{(\\ell)(g)} = \\text{argmin }\\rho_{\\tau_{\\ell}}(e_i^{(\\ell)(g)})-u_i^{(\\ell)(g)}e_i^{(\\ell)(g)}+\\frac{\\delta}{2}(Y_i^{(g)}-X_i\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}-e_i^{(\\ell)(g)})^2")  

  
![\\frac{\\partial Q}{\\partial
e\_i^{(\\ell)(g)}}=-\\delta(Y\_i^{(g)}-X\_i\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}-e\_i^{(\\ell)(g)})-u\_i^{(\\ell)(g)}+\\frac{\\partial
\\rho\_{\\tau\_{\\ell}}(e\_i^{(\\ell)(g)})}{\\partial
e\_i^{(\\ell)(g)}}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20e_i%5E%7B%28%5Cell%29%28g%29%7D%7D%3D-%5Cdelta%28Y_i%5E%7B%28g%29%7D-X_i%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D-e_i%5E%7B%28%5Cell%29%28g%29%7D%29-u_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Cfrac%7B%5Cpartial%20%5Crho_%7B%5Ctau_%7B%5Cell%7D%7D%28e_i%5E%7B%28%5Cell%29%28g%29%7D%29%7D%7B%5Cpartial%20e_i%5E%7B%28%5Cell%29%28g%29%7D%7D
"\\frac{\\partial Q}{\\partial e_i^{(\\ell)(g)}}=-\\delta(Y_i^{(g)}-X_i\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}-e_i^{(\\ell)(g)})-u_i^{(\\ell)(g)}+\\frac{\\partial \\rho_{\\tau_{\\ell}}(e_i^{(\\ell)(g)})}{\\partial e_i^{(\\ell)(g)}}")  

  
![\\frac{\\partial \\rho\_{\\tau\_{\\ell}}(e\_i^{(\\ell)(g)})}{\\partial
e\_i^{(\\ell)(g)}}=\\begin{cases}\\tau\_\\ell-1\\quad&\\text{if
}e\_i^{(\\ell)(g)}\<0\\\\0\\quad&\\text{if
}e\_i^{(\\ell)(g)}=0\\\\\\tau\_\\ell\\quad&\\text{if
}e\_i^{(\\ell)(g)}\>0\\end{cases}](https://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20%5Crho_%7B%5Ctau_%7B%5Cell%7D%7D%28e_i%5E%7B%28%5Cell%29%28g%29%7D%29%7D%7B%5Cpartial%20e_i%5E%7B%28%5Cell%29%28g%29%7D%7D%3D%5Cbegin%7Bcases%7D%5Ctau_%5Cell-1%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3C0%5C%5C0%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3D0%5C%5C%5Ctau_%5Cell%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%7D%3E0%5Cend%7Bcases%7D
"\\frac{\\partial \\rho_{\\tau_{\\ell}}(e_i^{(\\ell)(g)})}{\\partial e_i^{(\\ell)(g)}}=\\begin{cases}\\tau_\\ell-1\\quad&\\text{if }e_i^{(\\ell)(g)}\<0\\\\0\\quad&\\text{if }e_i^{(\\ell)(g)}=0\\\\\\tau_\\ell\\quad&\\text{if }e_i^{(\\ell)(g)}\>0\\end{cases}")  

  
![e\_i^{(\\ell)(g)^{k+1}}=\\begin{cases}Y\_i^{(g)}-X\_i\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}-\\frac{1}{\\delta}(u\_i^{(\\ell)(g)}+\\tau\_\\ell-1)&\\text{if
}e\_i^{(\\ell)(g)^k}\<0\\\\Y\_i^{(g)}-X\_i\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}-\\frac{1}{\\delta}(u\_i^{(\\ell)(g)})\\quad&\\text{if
}e\_i^{(\\ell)(g)^k}=0\\\\Y\_i^{(g)}-X\_i\\alpha^{(g)}-V\_i^{(\\ell)}\\eta^{(g)}-\\frac{1}{\\delta}(u\_i^{(\\ell)(g)}+\\tau\_\\ell)\\quad&\\text{if
}e\_i^{(\\ell)(g)^k}\>0\\end{cases}](https://latex.codecogs.com/png.latex?e_i%5E%7B%28%5Cell%29%28g%29%5E%7Bk%2B1%7D%7D%3D%5Cbegin%7Bcases%7DY_i%5E%7B%28g%29%7D-X_i%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%28u_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Ctau_%5Cell-1%29%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%5Ek%7D%3C0%5C%5CY_i%5E%7B%28g%29%7D-X_i%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%28u_i%5E%7B%28%5Cell%29%28g%29%7D%29%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%5Ek%7D%3D0%5C%5CY_i%5E%7B%28g%29%7D-X_i%5Calpha%5E%7B%28g%29%7D-V_i%5E%7B%28%5Cell%29%7D%5Ceta%5E%7B%28g%29%7D-%5Cfrac%7B1%7D%7B%5Cdelta%7D%28u_i%5E%7B%28%5Cell%29%28g%29%7D%2B%5Ctau_%5Cell%29%5Cquad%26%5Ctext%7Bif%20%7De_i%5E%7B%28%5Cell%29%28g%29%5Ek%7D%3E0%5Cend%7Bcases%7D
"e_i^{(\\ell)(g)^{k+1}}=\\begin{cases}Y_i^{(g)}-X_i\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}-\\frac{1}{\\delta}(u_i^{(\\ell)(g)}+\\tau_\\ell-1)&\\text{if }e_i^{(\\ell)(g)^k}\<0\\\\Y_i^{(g)}-X_i\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}-\\frac{1}{\\delta}(u_i^{(\\ell)(g)})\\quad&\\text{if }e_i^{(\\ell)(g)^k}=0\\\\Y_i^{(g)}-X_i\\alpha^{(g)}-V_i^{(\\ell)}\\eta^{(g)}-\\frac{1}{\\delta}(u_i^{(\\ell)(g)}+\\tau_\\ell)\\quad&\\text{if }e_i^{(\\ell)(g)^k}\>0\\end{cases}")  

## 5\. Closed form for multipliers

  
![\\begin{aligned}\\boldsymbol{u}^{(\\ell)(g)^{k+1}}:=&\\boldsymbol{u}^{(\\ell)(g)^k}+\\delta(\\boldsymbol{Y}^{(g)}-\\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)^{k+1}}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)^{k+1}}-\\boldsymbol{e}^{(\\ell)(g)^{k+1}})\\\\\\boldsymbol{w}^{(g)^{k+1}}:=&\\boldsymbol{w}^{(g)^{k}}+\\delta(\\boldsymbol{\\theta}^{(g)^{k+1}}-\\boldsymbol{\\eta}^{(g)^{k+1}})\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28g%29%5E%7Bk%2B1%7D%7D%3A%3D%26%5Cboldsymbol%7Bu%7D%5E%7B%28%5Cell%29%28g%29%5Ek%7D%2B%5Cdelta%28%5Cboldsymbol%7BY%7D%5E%7B%28g%29%7D-%5Cboldsymbol%7BX%7D%5Cboldsymbol%7B%5Calpha%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D-%5Cboldsymbol%7BV%7D%5E%7B%28%5Cell%29%7D%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D-%5Cboldsymbol%7Be%7D%5E%7B%28%5Cell%29%28g%29%5E%7Bk%2B1%7D%7D%29%5C%5C%5Cboldsymbol%7Bw%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%3A%3D%26%5Cboldsymbol%7Bw%7D%5E%7B%28g%29%5E%7Bk%7D%7D%2B%5Cdelta%28%5Cboldsymbol%7B%5Ctheta%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D-%5Cboldsymbol%7B%5Ceta%7D%5E%7B%28g%29%5E%7Bk%2B1%7D%7D%29%5Cend%7Baligned%7D
"\\begin{aligned}\\boldsymbol{u}^{(\\ell)(g)^{k+1}}:=&\\boldsymbol{u}^{(\\ell)(g)^k}+\\delta(\\boldsymbol{Y}^{(g)}-\\boldsymbol{X}\\boldsymbol{\\alpha}^{(g)^{k+1}}-\\boldsymbol{V}^{(\\ell)}\\boldsymbol{\\eta}^{(g)^{k+1}}-\\boldsymbol{e}^{(\\ell)(g)^{k+1}})\\\\\\boldsymbol{w}^{(g)^{k+1}}:=&\\boldsymbol{w}^{(g)^{k}}+\\delta(\\boldsymbol{\\theta}^{(g)^{k+1}}-\\boldsymbol{\\eta}^{(g)^{k+1}})\\end{aligned}")
