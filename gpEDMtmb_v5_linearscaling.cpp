
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_MATRIX(xcpue);
  DATA_MATRIX(xcatch);
  DATA_SCALAR(ysigma2);
  PARAMETER_VECTOR(logW);
  PARAMETER(logtau);
  PARAMETER(logg);
  PARAMETER(logb);
  //
  vector<Type> W = exp(logW);
  Type b = exp(logb);
  //
  Type tau = (ysigma2-0.001)/(1+exp(-logtau))+0.001;
  Type g = (ysigma2-0.001)/(1+exp(-logg))+0.001;
  REPORT(ysigma2);
  //
  matrix<Type> X(xcpue.rows(),xcpue.cols());
  X = xcpue-(xcatch*b);
  //
  vector<Type> v1(X.cols());
  vector<Type> v2(X.cols());
  matrix<Type> K(X.rows(),X.rows());
  //
  for(int i=0; i<K.rows(); i++) {
    for(int j=0; j<K.rows(); j++) {
      v1 = X.row(i)-X.row(j);
      v1 = v1.abs();
      v2 = (W*v1)*(W*v1);
      K(i,j) = tau*exp(-sum(v2));
    }
  }
  matrix<Type> G(K);
  G.fill(0.0);
  vector<Type> gvec(K.rows());
  gvec.fill(g);
  G.diagonal() = gvec;
  K = K+G;
  Type nll = density::MVNORM_t<Type>(K)(y);
  REPORT(K);
  vector<Type> v3 = W*W;
  Type v4 = (-0.5)*sum(v3)/(3.14159/2);
  Type v5 = dbeta(tau/(ysigma2*1.01),Type(1.1),Type(1.1),true);
  Type v6 = dbeta(g/(ysigma2*1.01),Type(1.1),Type(1.1),true);
  nll -= v4;
  nll -= v5;
  nll -= v6;
  return nll;
}
