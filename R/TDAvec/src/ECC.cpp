#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

NumericVector computeVAB(NumericMatrix D, int homDim, NumericVector scaleSeq);

// [[Rcpp::export]]
NumericVector computeECC(NumericMatrix D, int maxhomDim, NumericVector scaleSeq){
  NumericMatrix ecc(scaleSeq.size()-1,maxhomDim+1);
  for (int d=0;d<=maxhomDim;++d){
    ecc(_,d) = pow(-1,d)*computeVAB(D,d,scaleSeq);
  }
  return rowSums(ecc);
}

