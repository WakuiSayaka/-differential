数値微分の練習

* forward_difference
<img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x)}{h}" />
<img src="https://latex.codecogs.com/gif.latex?h=2^{-1},2^{-2},2^{-3},2^{-4},\dots,2^{-52}" />

* backward_difference
<img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x)-f(x-h)}{h}" />
<img src="https://latex.codecogs.com/gif.latex?h=2^{-1},2^{-2},2^{-3},2^{-4},\dots,2^{-52}" />

* central_difference
<img src="https://latex.codecogs.com/gif.latex?f'(x)&space;\approx&space;\frac{f(x&plus;h)-f(x-h)}{2h}" />
<img src="https://latex.codecogs.com/gif.latex?h=2^{-1},2^{-2},2^{-3},2^{-4},\dots,2^{-52}" />

* dualnumber
<img src="https://latex.codecogs.com/gif.latex?f(x&plus;\epsilon)=\begin{pmatrix}&space;f(x)&space;&&space;f'(x)\\&space;0&space;&&space;f(x)&space;\end{pmatrix}" />
<img src="https://latex.codecogs.com/gif.latex?\epsilon=\begin{pmatrix}&space;0&space;&&space;1\\&space;0&space;&&space;0&space;\end{pmatrix}\neq\begin{pmatrix}&space;0&space;&&space;0\\&space;0&space;&&space;0&space;\end{pmatrix}" />
<img src="https://latex.codecogs.com/gif.latex?\epsilon^2=\begin{pmatrix}&space;0&space;&&space;0\\&space;0&space;&&space;0&space;\end{pmatrix}" />
