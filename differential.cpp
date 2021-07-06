#include "DualNumber.hpp"

double func(double x) {
  double res;
  res = 4.0*x*x + 2.0*x + 3.0;
  // res = exp(x);
  // res = exp(2.0*x);
  // res = sin(x);
  // res = sin(2.0*x);
  // res = sin(x) * sin(x);
  // res = cos(x);
  // res = cos(2.0 * x);
  // res = cos(x) * cos(x);
  // res = tan(x);
  // res = tan(2.0*x);
  // res = tan(x)*tan(x);
  // res = log(x);
  // res = log(2.0*x);
  // res = log(x)*log(x);
  // res = sqrt(x);
  // res = pow(x,-0.5);
  // res = pow(x,1.25);
  // res = pow(x,-1.25);
  // res = sinh(x);
  // res = cosh(x);
  // res = tanh(x);
  return res;
}

double differential_func(double rx) {
  Matrix res;
  Matrix x;

  x.DualNumber();
  x = x + rx;  //ε + x

  res = 4.0*x*x + 2.0*x + 3.0;
  // res = exp(x);
  // res = exp(2.0*x);
  // res = sin(x);
  // res = sin(2.0*x);
  // res = sin(x) * sin(x);
  // res = cos(x);
  // res = cos(2.0 * x);
  // res = cos(x) * cos(x);
  // res = tan(x);
  // res = tan(2.0*x);
  // res = tan(x)*tan(x);
  // res = log(x);
  // res = log(2.0*x);
  // res = log(x)*log(x);
  // res = sqrt(x);
  // res = pow(x,-0.5);
  // res = pow(x,1.25);
  // res = pow(x,-1.25);
  // res = sinh(x);
  // res = cosh(x);
  // res = tanh(x);
  return res.GetDN();
}



 // main 関数
int main(){
  double x = 2.0;
  double fx,fdx;
  double test;

  fx   = func(x);
  fdx  = differential_func(x);

  cout << "x="     << x    << '\n';
  cout << "f(x)="  << fx   << '\n';
  cout << "f'(x)=" << fdx  << '\n';
  return 0;
}
