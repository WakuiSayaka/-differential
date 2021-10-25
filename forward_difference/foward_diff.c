#include <stdio.h>
#include <math.h>


double func(double x) {
  return x*x;
}

//前進差分
double foward_diff(double h,double x) {
  return (func(x + h) - func(x)) / h ;
}

int main(void) {
  double h,x;

  //任意の値
  x = 1.0;

  h = 1e-1;
  for (int i = 0; i < 15; i++) {
    //f'(x)を出力
    printf("[%2d]%.15lf\th=%e\n",i+1,foward_diff(h,x),h);
    h = h/10.0;
  }

  return 0;
}
