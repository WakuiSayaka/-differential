#include "DualNumber.hpp"


// 確認用の関数
// return exp(x);
// return exp(2.0*x);
// return sin(x);
// return sin(2.0*x);
// return sin(x) * sin(x);
// return cos(x);
// return cos(2.0 * x);
// return cos(x) * cos(x);
// return tan(x);
// return tan(2.0*x);
// return tan(x)*tan(x);
// return log(x);
// return log(2.0*x);
// return log(x)*log(x);
// return sqrt(x);
// return pow(x,-0.5);
// return pow(x,1.25);
// return pow(x,-1.25);
// return sinh(x);
// return cosh(x);
// return tanh(x);



//f(x) 実数
double func(double x) {
  return 4.0*x*x + 2.0*x + 3.0;
}

//f(x) 行列(2x2)
Matrix func(Matrix x) {
  return 4.0*x*x + 2.0*x + 3.0;
}

//f'(x)
double differential_func(double rx) {
  Matrix x;
  x.DualNumber(); //ε
  x += rx;  //ε + x

  return func(x).GetDN();
}

int main(){
  double x;
  //f(x)とf'(x)を出力する為の任意の値xを入力
  x = 2.0;

  cout << "x="     << x                     << '\n';
  cout << "f(x)="  << func(x)               << '\n';
  cout << "f'(x)=" << differential_func(x)  << '\n';
  return 0;
}
