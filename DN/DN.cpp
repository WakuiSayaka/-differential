#include <iostream>
#include <cmath>
using namespace std;


//再帰の形で書き直せそう
double factrial(int n) {
  if (n < 2) {
    return 1;
  } else {
    int temp=1;
    for (int i = 2; i <= n; i++) {
      temp *= i;
    }
    return temp;
  }
}


//DN(2x2)のみ,四則演算(+-*/)のみ
class Matrix{
    double Mat[2][2];
public:
  ~Matrix(void) {};

  Matrix(void){
    Mat[0][0] = 0.0; Mat[0][1] = 0.0;
    Mat[1][0] = 0.0; Mat[1][1] = 0.0;
  }

  Matrix(const double obj) {
    Matrix temp;
    temp = obj;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        this->Mat[i][j] = temp.Mat[i][j];
      }
    }
  }


  inline const Matrix operator + (void) const {
    return *this;
  }

  const Matrix operator - (void) const {
    Matrix res;

    res = -1.0 * (*this);

    return res;
  }

  friend const Matrix operator * (const Matrix lhs,const Matrix rhs) {
    Matrix res;
    double temp;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        temp = 0.0;
        for (int k = 0; k < 2; k++) {
          temp += lhs.Mat[i][k]*rhs.Mat[k][j];
        }
        res.Mat[i][j] = temp;
      }
    }
    return res;
  }

  friend const Matrix operator * (const Matrix lhs,const double rhs) {
    Matrix res;
    Matrix temp = rhs;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = lhs.Mat[i][j]*rhs;
      }
    }
    return res;
  }

  friend const Matrix operator * (const double lhs,const Matrix rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = lhs*rhs.Mat[i][j];
      }
    }
    return res;
  }


  friend const Matrix operator + (const Matrix lhs,const Matrix rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = lhs.Mat[i][j] + rhs.Mat[i][j];
      }
    }
    return res;
  }


  friend const Matrix operator - (const Matrix lhs,const Matrix rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = lhs.Mat[i][j] - rhs.Mat[i][j];
      }
    }
    return res;
  }


  friend const Matrix operator / (const Matrix lhs,const Matrix rhs) {
    Matrix res;
    Matrix rhs_inv;
    double temp;

    temp = rhs.Mat[0][0]*rhs.Mat[1][1] - rhs.Mat[0][1]*rhs.Mat[1][0];
    if (temp == 0) {
      cout << "逆行列が存在しない" << '\n';
      return res;
    } else {
      rhs_inv.Mat[0][0] = rhs.Mat[1][1];
      rhs_inv.Mat[0][1] = -rhs.Mat[0][1];
      rhs_inv.Mat[1][0] = -rhs.Mat[1][0];
      rhs_inv.Mat[1][1] = rhs.Mat[0][0];
      rhs_inv = (1.0/temp) * rhs_inv;
    }
    res = lhs * rhs_inv;

    return res;
  }


  Matrix operator = (const Matrix obj) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        this->Mat[i][j] = obj.Mat[i][j];
      }
    }
    return *this;
  }

  Matrix operator = (const double obj) {
    Matrix res;
    this->Mat[0][0] = obj;
    this->Mat[1][1] = obj;
    return *this;
  }

  Matrix Inverse(void){
    Matrix res;
    double temp;

    temp = this->Mat[0][0]*this->Mat[1][1] - this->Mat[0][1]*this->Mat[1][0];
    if (temp == 0) {
      cout << "逆行列が存在しない" << '\n';
      return *this;
    } else {
      res.Mat[0][0] = this->Mat[1][1];
      res.Mat[0][1] = -this->Mat[0][1];
      res.Mat[1][0] = -this->Mat[1][0];
      res.Mat[1][1] = this->Mat[0][0];
      res = (1.0/temp)*res;
    }
    return res;
  }

  Matrix M_pow(int n) {
    Matrix res;
    if (!n) {
      res.IdentityMatrix();
      return res;
    } else {
      res = *this;
      for (int i = 1; i < n; i++) {
        res = res * (*this);
      }
      return res;
    }
  }


  //DNの場合、ε^2=0で収束するのを利用
  //exp(a+ε)の計算は私はわからないのでexp(a+ε)=exp(a)+exp(ε)を使った
  //exp(a)はdouble型の用意された関数を使って計算してからMatrix {{exp(a),0},{0,exp(a)}}にする
  //exp(ε)の行列の計算はε^2=0のおかげですぐ収束
  //三角関数
  //sin(x+ε) = sin(x)*cos(ε) + cos(x)*sin(ε)
  //cos(x+ε) = cos(x)*cos(ε) - sin(x)*sin(ε)
  //sin(x),cos(x)はdouble型の用意された関数を使って計算してからMatrix {{sin(x),0},{0,sin(x)}},{{cos(x),0},{0,cos(x)}}
  //sin(ε),cos(ε)は
  //sin(A):=sum_k=0^inf ( ((-1)^k)/((2k+1)!)  ) * A^(2k+1)
  //cos(A):=sum_k=0^inf (  ((-1)^k)/((2k+1)!) ) * A^2k
  //k>=1の項はε^2=0により0になる
  //sin(ε)=ε,cos(ε)=I
  //まとめると
  //sin(x+ε) = sin(x) + cos(x)*ε
  //cos(x+ε) = cos(x) - sin(x)*ε
  //対数関数(自然対数)
  //x+ε=x*(1+ε/x)
  //log(x*(1+ε/x))=log(x) + log(1+ε/x)
  //log(x)はdouble型の用意された関数を使って計算してからMatrix{{log(x),0},{0,log(x)}}
  //log(1+ε/x)は
  //log(1+X) := sum_k=1^inf ((-1)^(k-1) * (x^k))/k
  //k>=2の項はε^2=0により0になる
  //log(1+ε/x) = (1/x)*ε
  //まとめると
  //log(x+ε) = log(x) + (1/x)*ε
  Matrix M_exp(void) {
		Matrix res;
    double rn;
    double dn;
    if (this->Mat[1][0] != 0) {
      cout << "対象外" << '\n';
      return res;
    }

    if((this->Mat[0][0] == this->Mat[1][1]) && (this->Mat[0][0] != 0.0)) {
      rn  = this->Mat[0][0];
      res = Matrix(exp(rn));
      if (this->Mat[0][1]!=0.0) {
        Matrix temp;
        Matrix m_dn;
        dn = this->Mat[0][1];
        m_dn.DualNumber();
        m_dn = dn*m_dn;
        for (int i = 0; i <= 1; i++) {
          temp = temp + (1.0/factrial(i))*m_dn.M_pow(i);
        }
        res = res * temp;
      }
      return res;
    } else {
      for (int i = 0; i <= 1; i++) {
        res = res + (1.0/factrial(i))*this->M_pow(i);
      }
      return res;
    }
	}

  void DualNumber(){
    Mat[0][0]=0.0; Mat[0][1]=1.0;
    Mat[1][0]=0.0; Mat[1][1]=0.0;
  }

  void IdentityMatrix(){
    Mat[0][0]=1.0; Mat[0][1]=0.0;
    Mat[1][0]=0.0; Mat[1][1]=1.0;
  }

  void Zero(){
    Mat[0][0]=0.0; Mat[0][1]=0.0;
    Mat[1][0]=0.0; Mat[1][1]=0.0;
  }



  double GetDN() {
    return Mat[0][1];
  }




  void show(){
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        cout << Mat[i][j];
        if (j!=1) {
          cout << ",";
        }
      }
      cout << '\n';
    }
  }
};

Matrix exp(Matrix obj) {
  Matrix res;
  res = obj.M_exp();
  return res;
}

Matrix pow(Matrix obj,int n) {
  Matrix res;
  res = obj.M_pow(n);
  return res;
}


Matrix Matrix_inv(Matrix obj) {
  Matrix res;
  res = obj.Inverse();
  return res;
}








double func(double x) {
  double res;

  res = 4*x*x + 2*x + 3;
  // res = exp(2.0*x);
  return res;
}

double differential_func(double rx) {
  Matrix res;
  Matrix x;

  x.DualNumber();
  x = x + Matrix(rx);

  res = Matrix(4.0)*x*x + Matrix(2.0)*x + Matrix(3.0);
  // res = exp(2.0*x);

  return res.GetDN();
}

 // main 関数
int main(){
  double x = 1.0;
  double fx,fdx;
  Matrix res;


  fx  = func(x);
  fdx = differential_func(x);

  cout << "f(x)=4x*x + 2x + 3" << '\n';
  cout << "x=" << x << '\n';
  cout << "f(x)=" << fx << '\n';
  cout << "f'(x)=" << fdx << '\n';


  return 0;
}
