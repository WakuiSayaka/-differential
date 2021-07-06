#include <iostream>
#include <cmath>
using namespace std;


//階乗
double factrial(int n) {
  if(n > 0) {
    return n * factrial(n-1);
  } else {
    return 1;
  }
}

// //二項係数 nCr
// double comb(int n,int r) {
//   if(n < 0 || r < 0 || n < r) return 0.0;
//   if(n == r || !r) return 1.0;
//   if (r > n-r) {
//     r = n-r;
//   }
//   return comb(n - 1, r - 1) * (double)n / (double)r;
// }
//
// //ベルヌーイ数（できなかった）
// double Bernoulli(int n){
//   if (n > 0) {
//     double res;
//     for (int i = 0; i < n; i++) {
//       res += comb(n+1,i)*Bernoulli(i);
//     }
//     return -1.0/(double)(n+1.0)*res;
//   } else {
//     return 1;
//   }
// }


//powを実数に対応させる為に必要
bool is_integer( double x ){
  return std::floor(x)==x;
}


//DN(2x2)のみ,1階微分のみ
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

  friend const Matrix operator + (const Matrix lhs,const double rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = lhs.Mat[i][j];
        if (i == j) res.Mat[i][j] += rhs;
      }
    }
    return res;
  }

  friend const Matrix operator + (const double lhs,const Matrix rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = rhs.Mat[i][j];
        if (i == j) res.Mat[i][j] += lhs;
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

  friend const Matrix operator - (const Matrix lhs,const double rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = lhs.Mat[i][j];
        if (i == j) res.Mat[i][j] -= rhs;
      }
    }
    return res;
  }

  friend const Matrix operator - (const double lhs,const Matrix rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = -rhs.Mat[i][j];
        if (i == j) res.Mat[i][j] += lhs;
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

  friend const Matrix operator / (const Matrix lhs,const double rhs) {
    Matrix res;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        res.Mat[i][j] = lhs.Mat[i][j]/rhs;
      }
    }
    return res;
  }


  friend const Matrix operator / (const double lhs,const Matrix rhs) {
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

  Matrix operator += (const Matrix obj) {
    *this = *this + obj;
    return *this;
  }

  Matrix operator += (const double obj) {
    *this = *this + obj;
    return *this;
  }

  Matrix operator -= (const Matrix obj) {
    *this = *this - obj;
    return *this;
  }

  Matrix operator -= (const double obj) {
    *this = *this - obj;
    return *this;
  }

  Matrix operator *= (const Matrix obj) {
    *this = *this * obj;
    return *this;
  }

  Matrix operator *= (const double obj) {
    *this = *this * obj;
    return *this;
  }

  Matrix operator /= (const Matrix obj) {
    *this = *this / obj;
    return *this;
  }

  Matrix operator /= (const double obj) {
    *this = *this / obj;
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

  Matrix M_pow(double n) {
    Matrix res;
    if (is_integer(n)) {
      if (!n) {
        res.IdentityMatrix();
        return res;
      } else if(n > 0) {
        res = *this;
        for (int i = 1; i < n; i++) {
          res = res * (*this);
        }
        return res;
      } else {
        res = *this;
        for (int i = 1; i < -n; i++) {
          res = res * (*this);
        }
        return 1.0/res;
      }
    } else {
      if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
        cout << "対象外" << '\n';
        return res;
      }
      if ((this->Mat[0][0] <= 0.0) || (this->Mat[0][1] < 0.0)) {
        cout << "a^x a>0 を満たしていない。" << '\n';
        return res;
      }
      res = n * this->M_log();
      res = res.M_exp();
      return res;
    }
  }


  //exp(a+ε)=exp(a)*exp(ε)
  //exp(a)はdouble型の用意された関数を使って計算してからMatrix {{exp(a),0},{0,exp(a)}}にする
  //exp(ε)の行列の計算はε^2=0ですぐ終わる
  Matrix M_exp(void) {
		Matrix res;
    double rn;
    double dn;
    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    if(this->Mat[0][0]) {
      rn  = this->Mat[0][0];
      res = Matrix(exp(rn));
      if (this->Mat[0][1]) {
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

  //三角関数(sin,cos)
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
  Matrix Ep_sin(void) {
    Matrix res;
    for (int i = 0; i <= 1; i++) {
      res = res + (pow(-1.0,(double)i) / factrial(2.0*i + 1.0)) * this->M_pow(2*i + 1);
    }
    return res;
  }
  Matrix Ep_cos(void) {
    Matrix res;
    for (int i = 0; i <= 1; i++) {
      res = res + (pow(-1.0,(double)i) / factrial(2.0*i + 1.0)) * this->M_pow(2*i);
    }
    return res;
  }
  Matrix M_sin(void) {
    Matrix res;
    double rn;
    double dn;
    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    if(this->Mat[0][0]) {
      Matrix sinx,cosx;
      rn  = this->Mat[0][0];
      sinx = Matrix(sin(rn));
      cosx = Matrix(cos(rn));
      if (this->Mat[0][1]) {
        Matrix m_dn;
        dn = this->Mat[0][1];
        m_dn.DualNumber();
        m_dn = dn*m_dn;
        //sin(x+ε) = sin(x)*cos(ε) + cos(x)*sin(ε)
        res = sinx*m_dn.Ep_cos() + cosx*m_dn.Ep_sin();
      }
      return res;
    } else {
      res = this->Ep_sin();
      return res;
    }
  }

  Matrix M_cos(void) {
    Matrix res;
    double rn;
    double dn;
    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    if(this->Mat[0][0]) {
      Matrix sinx,cosx;
      rn  = this->Mat[0][0];
      sinx = Matrix(sin(rn));
      cosx = Matrix(cos(rn));
      if (this->Mat[0][1]) {
        Matrix m_dn;
        dn = this->Mat[0][1];
        m_dn.DualNumber();
        m_dn = dn*m_dn;
        //cos(x+ε) = cos(x)*cos(ε) - sin(x)*sin(ε)
        res = cosx*m_dn.Ep_cos() - sinx*m_dn.Ep_sin();
      }
      return res;
    } else {
      res = this->Ep_cos();
      return res;
    }
  }


  //三角関数(tan) （やめた）
  //tan(x+ε) = ( tan(x) + tan(ε) ) / ( 1 - tan(x) * tan(ε) )
  //tan(x)はdouble型の用意された関数を使って計算してからMatrix {{tan(x),0},{0,tan(x)}}
  //tan(ε)は
  //tan(A) := sum_n=1^inf ( ( (-1)^(n-1) * 2^2n * (2^2n -1) * B_2n ) / (2n)! ) * x^(2n-1)
  //(B_nはベルヌーイ数)
  //n>=2の項はε^2=0により0になる
  //tan(ε) = ε
  //まとめると
  //tan(x + ε) = (tan(x) + ε) / (1 - tan(x)ε) = {{tan(x),1},{0,tan(x)}}/{{1,-tan(x)},{0,1}}
  //{{1,-tan(x)},{0,1}}の逆行列は{{1,-tan(x)},{0,1}}なので
  //{{tan(x),1},{0,tan(x)}}*{{1,-tan(x)},{0,1}}={{tan(x),1-tan^2(x)},{0,tan(x)}}
  //cos^2(x)=1/(1-tan^2(x))なので1/cos^2(x)=(1-tan^2(x))
  //{{tan(x),1/cos^2(x)},{0,tan(x)}}=tan(x) + (1/cos^2(x))ε
  //三角関数(tan) （こっちにした）
  //tan(x+ε) = sin(x+ε)/cos(x+ε)
  Matrix M_tan(void) {
    Matrix res;
    double rn;
    double dn;
    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    if(this->Mat[0][0]) {
      Matrix sinx,cosx;
      Matrix sinxe,cosxe;
      rn  = this->Mat[0][0];
      sinx = Matrix(sin(rn));
      cosx = Matrix(cos(rn));
      if (this->Mat[0][1]) {
        Matrix m_dn;
        dn = this->Mat[0][1];
        m_dn.DualNumber();
        m_dn = dn*m_dn;
        //sin(x+ε) = sin(x)*cos(ε) + cos(x)*sin(ε)
        sinxe = sinx*m_dn.Ep_cos() + cosx*m_dn.Ep_sin();
        //cos(x+ε) = cos(x)*cos(ε) - sin(x)*sin(ε)
        cosxe = cosx*m_dn.Ep_cos() - sinx*m_dn.Ep_sin();
        res = sinxe/cosxe;
      }
      return res;
    } else {
      res = this->Ep_sin()/this->Ep_cos();
      return res;
    }
  }

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
  Matrix M_log(void) {
		Matrix res;
    double rn;
    double dn;
    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    if(this->Mat[0][0] <= 0) {
      cout << "log(x) x>0 を満たしていない。" << '\n';
    }


    rn  = this->Mat[0][0];
    res = Matrix(log(rn));
    if (this->Mat[0][1]) {
      Matrix m_dn;
      //ε/x
      dn = this->Mat[0][1]/rn;
      m_dn.DualNumber();
      m_dn = dn*m_dn;
      //log(1+ε/x) の 1+ の部分はここで相殺される
      for (int i = 1; i <= 1; i++) {
        res = res + (pow(-1.0,(double)i-1) * m_dn.M_pow(i) )/(double)i ;
      }
    }
    return res;
	}

  //sqrt(x)=exp(0.5*log(x)) ただし x>0 に注意
  Matrix M_sqrt(void){
    Matrix res;
    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }
    if ((this->Mat[0][0] <= 0.0) || (this->Mat[0][1] < 0.0)) {
      cout << "sqrt(x) x>0 を満たしていない。" << '\n';
      return res;
    }
    res = 0.5 * this->M_log();
    res = res.M_exp();
    return res;
  }

  Matrix M_sinh(void) {
    Matrix res;

    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    //0除算防止
    if (this->Mat[0][0] >= 0) {
      res = 0.5 * (this->M_exp() - 1.0/this->M_exp() );
    } else {
      Matrix negative;
      negative = -*this;
      res = 0.5 * (1.0/negative.M_exp() - negative.M_exp() );
    }
    return res;
  }

  Matrix M_cosh(void) {
    Matrix res;

    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    //0除算防止
    if (this->Mat[0][0] >= 0) {
      res = 0.5 * (this->M_exp() + 1.0/this->M_exp() );
    } else {
      Matrix negative;
      negative = -*this;
      res = 0.5 * (1.0/negative.M_exp() + negative.M_exp() );
    }
    return res;
  }

  Matrix M_tanh(void) {
    Matrix res;

    if ((this->Mat[1][0] != 0) || (this->Mat[0][0] != this->Mat[1][1])) {
      cout << "対象外" << '\n';
      return res;
    }

    //オーバーフロー防止
    if (this->Mat[0][0] > 0.5) {
      Matrix temp = 2.0 * (*this);
      res = 1.0 - 2.0/(1.0 + temp.M_exp());
    } else if (this->Mat[0][0] < -0.5) {
      Matrix temp = -2.0 * (*this);
      res = 2.0/(1.0 + temp.M_exp()) - 1.0;
    } else {
      res = this->M_sinh()/this->M_cosh();
    }
    return res;
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
        if (!j) {
          cout << ",";
        }
      }
      cout << '\n';
    }
  }
};

void show(Matrix obj) {
  obj.show();
}

Matrix exp(Matrix obj) {
  Matrix res;
  res = obj.M_exp();
  return res;
}

Matrix sin(Matrix obj) {
  Matrix res;
  res = obj.M_sin();
  return res;
}

Matrix cos(Matrix obj) {
  Matrix res;
  res = obj.M_cos();
  return res;
}

Matrix tan(Matrix obj) {
  Matrix res;
  res = obj.M_tan();
  return res;
}

Matrix log(Matrix obj) {
  Matrix res;
  res = obj.M_log();
  return res;
}

Matrix sqrt(Matrix obj) {
  Matrix res;
  res = obj.M_sqrt();
  return res;
}

Matrix pow(Matrix obj,double n) {
  Matrix res;
  res = obj.M_pow(n);
  return res;
}

Matrix sinh(Matrix obj) {
  Matrix res;
  res = obj.M_sinh();
  return res;
}

Matrix cosh(Matrix obj) {
  Matrix res;
  res = obj.M_cosh();
  return res;
}

Matrix tanh(Matrix obj) {
  Matrix res;
  res = obj.M_tanh();
  return res;
}

Matrix inv(Matrix obj) {
  Matrix res;
  res = obj.Inverse();
  return res;
}


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

//fdxが合っているかの確認用
double Test_differential_func(double x) {
  double res;
  res = 8.0 * x + 2.0;                                 //(4.0 * x * x + 2.0 * x + 3.0)'
  // res = exp(x);                                        // (exp(x))'
  // res = 2.0 * exp(2.0 * x);                            // (exp(2.0 * x))'
  // res = cos(x);                                        // (sin(x))'
  // res = 2.0*cos(2.0*x);                                // (sin(2.0 * x))'
  // res = sin(2.0*x);                                    // (sin(x) * sin(x))'
  // res = -sin(x);                                       // (cos(x))'
  // res = -2.0 * sin(2.0*x);                             // (cos(2.0 * x))'
  // res = -sin(2.0*x);                                   // (cos(x) * cos(x))'
  // res = 1.0/(cos(x)*cos(x));                           // (tan(x))'
  // res = 2.0/(cos(2.0*x)*cos(2.0*x));                   // (tan(2.0 * x))'
  // res = 2.0 * sin(x) / (cos(x) * cos(x) * cos(x));     // (tan(x) * tan(x))'
  // res = 1.0 / x;                                       // (log(x))'
  // res = 1.0 / x;                                       // (log(2.0 * x))'
  // res = 2.0 * log(x) / x;                              // (log(x) * log(x))'
  // res = 0.5/sqrt(x);                                   // (sqrt(x))'
  // res = -0.5*pow(x,-3.0/2.0);                          // (1.0/sqrt(x))'
  // res = 1.25*pow(x,0.25);                              // (pow(x,1.25))'
  // res = -1.25*pow(x,-2.25);                            // (pow(x,-1.25))'
  // res = cosh(x);                                       //(sinh(x))'
  // res = sinh(x);                                       //(cosh(x))'
  // res = 1.0/(cosh(x)*cosh(x));                         //(tanh(x))'
  return res;
}

 // main 関数
int main(){
  double x = 2.0;
  double fx,fdx;
  double test;

  fx   = func(x);
  fdx  = differential_func(x);
  test = Test_differential_func(x);

  cout << "x="     << x    << '\n';
  cout << "f(x)="  << fx   << '\n';
  cout << "f'(x)=" << fdx  << '\n';
  cout << "f'(x)=" << test << '\n';

  return 0;
}
