#include <Rcpp.h>
using namespace Rcpp;

//' Multiply a Number by Two
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}


