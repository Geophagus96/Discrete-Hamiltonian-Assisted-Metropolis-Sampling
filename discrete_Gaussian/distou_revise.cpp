#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double distou(const int &i, const vec &cumd){
  
  vec number(1, fill::randu);
  double start, end;
  double u;
  if(i == 1){
    start = 0.0;
    end = cumd(0);
  }else{
    start = cumd((i-2));
    end = cumd((i-1));
  }
  u = start + (end-start)*number(0);
  return u;
  
}

// [[Rcpp::export]]
int utodis (const double&u, const vec &cumd) {
  vec a(1, fill::zeros);
  vec c;
  c = join_cols(a, cumd);
  int low = 0, high = c.n_elem-2;
  while (low <= high) {
    int mid = (low + high) / 2;
    if (c(mid) <= u && c((mid+1)) > u) {
      return (mid+1);
    } else if (c(mid) < u) {
      low = mid + 1;
    } else {
      high = mid - 1;
    }
  }
  return c.n_elem - 1;
}

//[[Rcpp::export]]
double jlless (const double &pi1, const double &pi, const double &pl1, const double &pl, const double &gamma){
  if((-1+pi1+gamma) < (-1+pl1)){
    if((-1+pi+gamma) < (-1+pl1)){
      return 0;
    }
    else if(((-1+pl1)<=(-1+pi+gamma)) && ((-1+pi+gamma)<(-1+pl))){
      return pow((pi+gamma-pl1),2)/2;
    }
    else if(((-1+pl)<=(-1+pi+gamma))&&((-1+pi+gamma)<pl1)){
      return pow((pl-pl1),2)/2+(pl-pl1)*(pi+gamma-pl);
    }
    else if((pl1<=(-1+pi+gamma))&&((-1+pi+gamma)<pl)){
      return pow((pl-pl1),2)/2+(pl-pl1)*(pi+gamma-pl)+pow((-1+pi+gamma-pl1),2)/2;
    }
    else{
      return pow((pl-pl1),2)+(pl-pl1)*(2*(pi+gamma-pl)-1);
    }
  }
  else if(((-1+pl1)<=(-1+pi1+gamma))&&((-1+pi1+gamma)<(-1+pl))){
    if(((-1+pl1)<=(-1+pi+gamma)) && ((-1+pi+gamma)<(-1+pl))){
      return pow((pi+gamma-pl1),2)/2-pow((pi1+gamma-pl1),2)/2;
    }
    else if(((-1+pl)<=(-1+pi+gamma))&&((-1+pi+gamma)<pl1)){
      return pow(pl-pl1,2)/2-pow((pi1+gamma-pl1),2)/2+(pl-pl1)*(pi+gamma-pl);
    }
    else if((pl1<=(-1+pi+gamma))&&((-1+pi+gamma)<pl)){
      return pow(pl-pl1,2)/2-pow((pi1+gamma-pl1),2)/2+(pl-pl1)*(pi+gamma-pl)+pow((-1+pi+gamma-pl1),2)/2;
    }
    else{
      return pow((pl-pl1),2)-pow((pi1+gamma-pl1),2)/2+(pl-pl1)*(2*(pi+gamma-pl)-1);
    }
  }
  else if(((-1+pl)<=(-1+pi1+gamma))&&((-1+pi1+gamma)<pl1)){
    if(((-1+pl)<=(-1+pi+gamma))&&((-1+pi+gamma)<pl1)){
      return (pl-pl1)*(pi-pi1);
    }
    else if((pl1<=(-1+pi+gamma))&&((-1+pi+gamma)<pl)){
      return (pl-pl1)*(pi-pi1)+pow((-1+pi+gamma-pl1),2)/2;
    }
    else{
      return (pl-pl1)*(pi-pi1)+pow((pl-pl1),2)/2+(pl-pl1)*(-1+pi+gamma-pl);
    }
  }
  else if((pl1<=(-1+pi1+gamma))&&((-1+pi1+gamma)<pl)){
    if ((pl1<=(-1+pi+gamma))&&((-1+pi+gamma)<pl)){
      return ((-1+pi1+gamma)+(-1+pi+gamma)-2*(2*pl1-pl))*(pi-pi1)/2;
    }
    else{
      return (pl+(-1+pi+gamma)-2*(2*pl1-pl))*(pl-(-1+pi1+gamma))/2+2*(pl-pl1)*(-1+pi+gamma-pl);
    }
  }
  else{
    return 2*(pl-pl1)*(pi-pi1);
  }
}

//[[Rcpp::export]]
double jlwithin(const double &u, const double &u1, const double &pl1, const double &pl, const double &gamma){
  if ((u1+gamma)<=pl){
    return (u1-u)*gamma;
  }
  else if (((u1+gamma)>pl)&&((-1+u1+gamma)<=pl1)){
    if((u+gamma)<=pl){
      return gamma*(pl-gamma-u)+(gamma+pl-u1)*(u1-pl+gamma)/2;
    }
    else{
      return (u1-u)*(pl-u+pl-u1)/2;
    }
  }
  else{
    if((u+gamma)<=pl){
      return gamma*(pl-gamma-u)+(1-pl+pl1)*(2*gamma+pl-pl1-1)/2+(pl-pl1+gamma-1)*(-1+u1+gamma-pl1);
    }
    else if (((u+gamma)>pl)&&((-1+u+gamma)<=pl1)){
      return (2*pl-pl1-u+gamma-1)*(pl1-gamma-u+1)/2+(pl-pl1+gamma-1)*(-1+u1-pl1+gamma);
    }
    else{
      return (pl-pl1+gamma-1)*(u1-u);
    }
  }
}

//[[Rcpp::export]]
double jlmore (const double &pi1, const double &pi, const double &pl1, const double &pl, const double &gamma){
  if ((-1+pi+gamma)<=pl1){
    return 0;
  }
  else if ((pl1<=(-1+pi+gamma))&&((-1+pi+gamma)<pl)){
    if ((-1+pi1+gamma)<=pl1){
      return pow((-1+pi+gamma-pl1),2)/2;
    }
    else{
      return pow((-1+pi+gamma-pl1),2)/2-pow((-1+pi1+gamma-pl1),2)/2;
    }
  }
  else{
    if ((-1+pi1+gamma)<=pl1){
      return pow((pl-pl1),2)/2+(pl-pl1)*(-1+pi+gamma-pl);
    }
    else if ((pl1<=(-1+pi1+gamma))&&((-1+pi1+gamma)<pl)){
      return pow((pl-pl1),2)/2-pow((-1+pi1+gamma-pl1),2)/2+(pl-pl1)*(-1+pi+gamma-pl);
    }
    else{
      return (pl-pl1)*(pi-pi1);
    }
  }
}

//[[Rcpp::export]]
double tranint(const double &u, const double &u1, const double &pl1, const double &pl, const double &gamma){
  if (u<=pl1){
    if (u1<=pl1){
      return jlless(u, u1, pl1, pl, gamma);
    }
    else if ((pl1<u1)&&(u1<=pl)){
      return jlless(u, pl1, pl1, pl, gamma)+jlwithin(pl1, u1, pl1, pl, gamma);
    }
    else{
      return jlless(u, pl1, pl1, pl, gamma)+jlwithin(pl1, pl, pl1, pl, gamma)+jlmore(pl, u1, pl1, pl, gamma);
    }
  }
  else if ((pl1<u)&&(u<=pl)){
    if ((pl1<u1)&&(u1<=pl)){
      return jlwithin(u, u1, pl1, pl, gamma);
    }
    else{
      return jlwithin(u, pl, pl1, pl, gamma)+jlmore(pl, u1, pl1, pl, gamma);
    }
  }
  else{
    return jlmore(u, u1, pl1, pl, gamma);
  }
}

// This function computes the transition probability p(x_1 | x_0) 
// for the over-relaxed sampler. In this context, `alpha` represents 
// the sign preceding w_0, and `gamma` corresponds to the β parameter 
// in the paper. The over-relaxation update is given by:
//     w_1 = (alpha * w_0 + gamma * w_tilde) % 1

//[[Rcpp::export]]
double trans_prob(const int &i, const int &l, const vec &cumd, const int &alpha, const double &gamma){
  vec a(1, fill::zeros);
  vec c;
  c = join_cols(a, cumd);
  if(gamma>=0){
    if(alpha==1){
      if(i<l){
        return jlless(c((i-1)), c(i), c((l-1)), c(l), gamma)/(gamma*(c(i)-c((i-1))));
      }
      else if(i==l){
        return jlwithin(c((l-1)), c(l),c((l-1)), c(l), gamma)/(gamma*(c(i)-c((i-1))));
      }
      else{
        return jlmore(c((i-1)), c(i), c((l-1)), c(l), gamma)/(gamma*(c(i)-c((i-1))));
      }
    }
    else{
      return tranint(1-c(i), 1-c((i-1)), c((l-1)), c(l), gamma)/(gamma*(c(i)-c((i-1))));
    }
  }
  else{
    if (alpha==1){
      if((c(i)+gamma) <= 0){
        return tranint(1+c((i-1))+gamma, 1+c(i)+gamma, c((l-1)), c(l), -gamma)/(-gamma*(c(i)-c((i-1))));
      }
      else if ((c((i-1))+gamma)>=0){
        return tranint(c((i-1))+gamma, c(i)+gamma, c((l-1)), c(l), -gamma)/(-gamma*(c(i)-c((i-1))));
      }
      else{
        return (tranint(c((i-1))+gamma+1, 1, c((l-1)), c(l), -gamma)+tranint(0, c(i)+gamma, c((l-1)), c(l), -gamma))/(-gamma*(c(i)-c((i-1))));
      }
    }
    else{
      if ((1-c((i-1))+gamma)<=0){
        return tranint(2-c(i)+gamma, 2-c((i-1))+gamma, c((l-1)), c(l), -gamma)/(-gamma*(c(i)-c((i-1))));
      }
      else if ((1-c(i)+gamma)>0){
        return tranint(1-c(i)+gamma, 1-c((i-1))+gamma, c((l-1)), c(l), -gamma)/(-gamma*(c(i)-c((i-1))));
      }
      else{
        return (tranint(2-c(i)+gamma,1, c((l-1)), c(l), -gamma)+tranint(0,  1-c((i-1))+gamma, c((l-1)), c(l), -gamma))/(-gamma*(c(i)-c((i-1))));
      }
    }
  }
}

