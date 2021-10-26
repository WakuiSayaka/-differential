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


* dualnumber/DualNumber.hpp

Matrix(double)  Matrix型へのキャスト

Matrix型とMatrix型の四則演算

Matrix型とdouble型の四則演算

double型とMatrix型の四則演算

Matrix型への代入(右辺はMatrix型とdouble型)

Matrix型への代入演算子(右辺はMatrix型とdouble型)

Matrix.DualNumber();       Matrix型の変数に行列εを代入

<img src="https://latex.codecogs.com/gif.latex?\epsilon&space;=\begin{pmatrix}&space;0&space;&&space;1\\&space;0&space;&&space;1&space;\end{pmatrix}" />

Matrix.IdentityMatrix();   Matrix型の変数に単位行列Iを代入

<img src="https://latex.codecogs.com/gif.latex?I&space;=&space;\begin{pmatrix}&space;1&space;&&space;0\\&space;0&space;&&space;1&space;\end{pmatrix}" />
