数値微分の練習

* forward_difference
<img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x)}{h}" />
<img src="https://latex.codecogs.com/gif.latex?h=2^{n},2^{n-1},2^{n-2},2^{n-3},\dots,2^{n-52}" />
<img src="https://latex.codecogs.com/gif.latex?n=\lfloor\log_2|x|\rfloor" />

* backward_difference
<img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x)-f(x-h)}{h}" />
<img src="https://latex.codecogs.com/gif.latex?h=2^{n},2^{n-1},2^{n-2},2^{n-3},\dots,2^{n-52}" />
<img src="https://latex.codecogs.com/gif.latex?n=\lfloor\log_2|x|\rfloor" />

* central_difference
<img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" />
<img src="https://latex.codecogs.com/gif.latex?h=2^{n},2^{n-1},2^{n-2},2^{n-3},\dots,2^{n-52}" />
<img src="https://latex.codecogs.com/gif.latex?n=\lfloor\log_2|x|\rfloor" />

* dualnumber
<img src="https://latex.codecogs.com/gif.latex?f(x&plus;\epsilon)=\begin{pmatrix}&space;f(x)&space;&&space;f'(x)\\&space;0&space;&&space;f(x)&space;\end{pmatrix}" />
<img src="https://latex.codecogs.com/gif.latex?\epsilon=\begin{pmatrix}&space;0&space;&&space;1\\&space;0&space;&&space;0&space;\end{pmatrix}\neq\begin{pmatrix}&space;0&space;&&space;0\\&space;0&space;&&space;0&space;\end{pmatrix}" />
<img src="https://latex.codecogs.com/gif.latex?\epsilon^2=\begin{pmatrix}&space;0&space;&&space;0\\&space;0&space;&&space;0&space;\end{pmatrix}" />
