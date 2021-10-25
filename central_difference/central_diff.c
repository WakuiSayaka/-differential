#include <stdio.h>
#include <math.h>


double func(double x) {
  return x*x;
}

//中心差分
double central_diff(double h,double x) {
  return (func(x + h) - func(x - h)) * 0.5 / h ;
}

int main(void) {
  double h,x;

  //任意の値
  x = 1.0;

  h = 0.5;
  for (int i = 0; i < 52; i++) {
    //f'(x)を出力
    printf("[%2d]%.15lf\th=%.15e\n",i+1,central_diff(h,x),h);
    h *= 0.5;
  }

  return 0;
}
